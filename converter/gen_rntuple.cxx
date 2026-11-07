/* Based on Jakob Blomer's ATLAS conversion example from iotools:
 * https://github.com/jblomer/iotools/blob/master/gen_atlas.cxx
 * Adapted here for split Go4 detector trees using ROOT's RNTupleImporter.
 * This we keep out of any MONAD dependency. Is just straight-up ROOT+C++ */

/* Expose the detector leaves of a split Go4 tree to ROOT's standard RNTuple importer. */
#include <ROOT/RNTupleImporter.hxx>
#include <ROOT/RField.hxx>

#include <TBranch.h>
#include <TBranchElement.h>
#include <TDataType.h>
#include <TClass.h>
#include <TFile.h>
#include <TKey.h>
#include <TLeaf.h>
#include <TLeafC.h>
#include <TObjArray.h>
#include <TROOT.h>
#include <TStreamerInfo.h>
#include <TStreamerElement.h>
#include <TSystem.h>
#include <TTree.h>

#include <charconv>
#include <cstdio>
#include <cstring>
#include <filesystem>
#include <getopt.h>
#include <limits>
#include <memory>
#include <optional>
#include <set>
#include <stdexcept>
#include <string>
#include <vector>

namespace fs = std::filesystem;

struct LeafType {
	char code;
	std::size_t size;
};

static std::optional<LeafType> GetLeafType(const std::string& type) {
	const TDataType* dataType = gROOT->GetType(type.c_str());
	if(!dataType) return std::nullopt;
	switch(dataType->GetType()) {
		case kBool_t:    return LeafType{'O', sizeof(Bool_t)   };
		case kChar_t:    return LeafType{'B', sizeof(Char_t)   };
		case kUChar_t:   return LeafType{'b', sizeof(UChar_t)  };
		case kShort_t:   return LeafType{'S', sizeof(Short_t)  };
		case kUShort_t:  return LeafType{'s', sizeof(UShort_t) };
		case kInt_t:     return LeafType{'I', sizeof(Int_t)    };
		case kUInt_t:    return LeafType{'i', sizeof(UInt_t)   };
		case kLong_t:    return LeafType{'G', sizeof(Long_t)   };
		case kULong_t:   return LeafType{'g', sizeof(ULong_t)  };
		case kLong64_t:  return LeafType{'L', sizeof(Long64_t) };
		case kULong64_t: return LeafType{'l', sizeof(ULong64_t)};
		case kFloat_t:   return LeafType{'F', sizeof(Float_t)  };
		case kDouble_t:  return LeafType{'D', sizeof(Double_t) };
		default: return std::nullopt;
	}
}

// This tree contains only branch descriptions and buffers, never intermediate baskets on disk.
// ROOT's importer supplies the destination buffers; GetEntry fills them from the source event.
class DetectorTree : public TTree {
	struct Column {
		TBranch* source = nullptr;
		TBranch* destination = nullptr;
		std::vector<unsigned char> initial;
		std::unique_ptr<void*> object;
		std::shared_ptr<void> emptyObject;
	};

	TTree& source;
	std::vector<Column> columns;
	std::size_t skipped = 0;
	bool sourceReady = false;

	void Skip(TBranch& branch, const std::string& type, const std::string& reason) {
		++skipped;
		std::fprintf(stderr, "Skipping %s [%s]: %s\n", branch.GetName(), type.c_str(), reason.c_str());
	}

	static void* SourceAddress(TBranch& branch) {
		if(const auto element = dynamic_cast<TBranchElement*>(&branch)) {
			if(element->GetID() < 0) return element->GetObject();
			if(!element->GetAddress() || !element->GetInfo())
				throw std::runtime_error("Cannot locate source data: " + std::string{branch.GetName()});
			return element->GetAddress() + element->GetInfo()->GetElementOffset(element->GetID());
		}
		return static_cast<TLeaf*>(branch.GetListOfLeaves()->At(0))->GetValuePointer();
	}

	void AddObject(TBranch& branch, TLeaf& leaf, const std::string& name) {
		const std::string type = leaf.GetTypeName();
		auto cl = TClass::GetClass(type.c_str(), true, true);
		if(!cl || !cl->HasDictionary()) {
			Skip(branch, type, "no dictionary available");
			return;
		}
		if(const auto element = dynamic_cast<TBranchElement*>(&branch); element && element->GetID() >= 0) {
			const auto member = static_cast<TStreamerElement*>(element->GetInfo()->GetElements()->At(element->GetID()));
			if(!member || member->IsaPointer()) {
				Skip(branch, type, "raw object pointers are not supported");
				return;
			}
		}
		const auto count = leaf.GetLenStatic();
		const std::string fieldType = count > 1
			? "std::array<" + type + ',' + std::to_string(count) + '>' : type;
		try {
			ROOT::RFieldBase::Create("probe", fieldType).Unwrap();
		} catch(const std::exception& error) {
			const std::string reason = error.what();
			Skip(branch, type, reason.substr(0, reason.find('\n')));
			return;
		}
		if(!sourceReady && source.GetEntries() != 0) {
			if(source.GetEntry(0) <= 0) throw std::runtime_error("Cannot read first source entry.");
			sourceReady = true;
		}
		Column column;
		column.source = &branch;
		if(sourceReady) {
			column.object = std::make_unique<void*>(SourceAddress(branch));
		} else {
			// Empty trees still need a valid object pointer while the importer prepares its schema.
			void* object = count > 1 ? cl->NewArray(count) : cl->New();
			if(!object) {
				Skip(branch, type, "cannot construct an object for the empty input");
				return;
			}
			column.emptyObject = std::shared_ptr<void>(object, [cl, count](void* p) {
				if(count > 1) cl->DeleteArray(p); else cl->Destructor(p);
			});
			column.object = std::make_unique<void*>(object);
		}
		column.destination = Branch(name.c_str(), cl->GetName(), column.object.get(), 32000, 0);
		if(!column.destination) throw std::runtime_error("Cannot create object branch: " + name);
		static_cast<TLeaf*>(column.destination->GetListOfLeaves()->First())->SetLen(count);
		columns.push_back(std::move(column));
	}

	void AddBranches(TObjArray& branches) {
		for(auto object : branches) {
			auto branch = static_cast<TBranch*>(object);
			if(const auto element = dynamic_cast<TBranchElement*>(branch)) {
				if(element->GetStreamerType() == TStreamerInfo::kBase) {
					Skip(*branch, "BASE", "base-class fields are omitted");
					continue;
				}
			}
			if(branch->GetListOfBranches()->GetEntriesFast() != 0) {
				AddBranches(*branch->GetListOfBranches());
				continue;
			}
			if(branch->GetNleaves() != 1) {
				Skip(*branch, "leaf list", "expected a single terminal leaf");
				continue;
			}
			auto leaf = static_cast<TLeaf*>(branch->GetListOfLeaves()->At(0));
			if(dynamic_cast<TLeafC*>(leaf)) {
				Skip(*branch, leaf->GetTypeName(), "C-string leaves are not supported by the adapter");
				continue;
			}
			if(leaf->GetLeafCount()) {
				Skip(*branch, leaf->GetTypeName(), "variable-length leaf arrays are not supported");
				continue;
			}
			std::string name = branch->GetName();
			if(const auto bracket = name.find('['); bracket != std::string::npos) name.erase(bracket);
			const auto type = GetLeafType(leaf->GetTypeName());
			if(!type) {
				AddObject(*branch, *leaf, name);
				continue;
			}
			std::string dimensions;
			const std::string title = leaf->GetTitle();
			if(const auto bracket = title.find('['); bracket != std::string::npos)
				dimensions = title.substr(bracket);
			const std::string declaration = name + dimensions + '/' + type->code;
			columns.emplace_back();
			auto& column = columns.back();
			column.source = branch;
			column.initial.resize(static_cast<std::size_t>(leaf->GetLenStatic()) * type->size);
			column.destination = Branch(name.c_str(), column.initial.data(), declaration.c_str());
			if(!column.destination) throw std::runtime_error("Cannot create detector branch: " + name);
		}
	}

public:
	explicit DetectorTree(TTree& input)
		: TTree(input.GetName(), input.GetTitle(), 99, nullptr), source(input)
	{
		for(auto object : *input.GetListOfBranches()) static_cast<TBranch*>(object)->SetAutoDelete(false);
		AddBranches(*input.GetListOfBranches());
		if(columns.empty()) throw std::runtime_error("No supported fields found.");
		SetEntries(input.GetEntries());
		std::fprintf(stderr, "Selected %zu fields; skipped %zu branches.\n", columns.size(), skipped);
	}

	using TTree::SetBranchAddress;
	Int_t SetBranchAddress(const char* name, void* address, TBranch** branch = nullptr) override {
		for(const auto& column : columns) {
			if(column.object && std::string{column.destination->GetName()} == name) {
				// The importer writes directly from the source object, with no byte-copy of C++ containers.
				*static_cast<void**>(address) = *column.object;
				if(branch) *branch = column.destination;
				return TTree::kMatch;
			}
		}
		return TTree::SetBranchAddress(name, address, branch);
	}

	Int_t SetBranchAddress(const char* name, void* address, TClass*, EDataType, bool) override {
		return SetBranchAddress(name, address);
	}

	Int_t SetBranchAddress(const char* name, void* address, TBranch** branch, TClass*, EDataType, bool) override {
		return SetBranchAddress(name, address, branch);
	}

	Int_t GetEntry(Long64_t entry, Int_t getall = 0) override {
		const auto bytes = source.GetEntry(entry, getall);
		if(bytes <= 0) throw std::runtime_error("Cannot read source entry " + std::to_string(entry));
		for(const auto& column : columns) {
			const void* data = SourceAddress(*column.source);
			if(column.object) {
				if(data != *column.object) throw std::runtime_error("Source object moved: " + std::string{column.source->GetName()});
			} else {
				if(!data || !column.destination->GetAddress())
					throw std::runtime_error("Null buffer for detector branch: " + std::string{column.source->GetName()});
				std::memcpy(column.destination->GetAddress(), data, column.initial.size());
			}
		}
		return bytes;
	}
};

static int Compression(const std::string& name) {
	if(name == "none") return   0;
	if(name == "zlib") return 101;
	if(name == "lz4")  return 404;
	if(name == "lzma") return 207;
	if(name == "zstd") return 505;
	throw std::runtime_error("Unknown compression: " + name);
}

static void FindTrees(TDirectory& directory, const std::string& prefix, std::vector<std::string>& names) {
	std::set<std::string> visited;
	for(auto object : *directory.GetListOfKeys()) {
		auto key = static_cast<TKey*>(object);
		if(!visited.insert(key->GetName()).second) continue; // Multiple key cycles describe the same object.
		auto cl = TClass::GetClass(key->GetClassName());
		if(!cl) continue;
		const std::string path = prefix + key->GetName();
		if(cl->InheritsFrom(TTree::Class())) names.push_back(path);
		else if(cl->InheritsFrom(TDirectory::Class())) {
			auto child = directory.GetDirectory(key->GetName());
			if(!child) throw std::runtime_error("Cannot read directory: " + path);
			FindTrees(*child, path + '/', names);
		}
	}
}

static void Usage(const char* program) {
	std::printf("Usage: %s -i input.root [-o output-directory] [-c none|zlib|lz4|lzma|zstd] [-t tree-name] [-n entries] [-l library.so]\n\n"
		"Defaults: output-directory=<input-directory>/rn, compression=zstd (level 5).\n"
		"Without -t, the input must contain exactly one TTree (including subdirectories).\n"
		"Writes <input-filename>; an existing output file is never overwritten.\n"
		"Base-class fields are omitted. Unsupported fields are skipped with a reason.\n"
		"Use repeated -l options to load dictionaries. Dots in names become underscores.\n\n"
		"Compression | Level | Decompression speed           | Typical file size\n"
		"------------|-------|-------------------------------|--------\n"
		"none        | -     | No decompression cost         | Largest\n"
		"lz4         | 4     | Fastest compressed option     | Larger\n"
		"zstd        | 5     | Very fast                     | Small\n"
		"zlib        | 1     | Slower than LZ4/ZSTD          | Medium\n"
		"lzma        | 7     | Usually slowest               | Often smallest\n\n"
		"These are general tendencies; results depend on the data and storage.\n"
		"Start with zstd for balance, or lz4 for decompression speed.\n"
		"Compressed files can read faster than none by reducing disk/network I/O.\n"
		"If -c is repeated, the last setting wins.\n", program);
}

int main(int argc, char** argv) {
	try {
		fs::path input, outputDirectory;
		std::string treeName;
		int compression = Compression("zstd");
		std::optional<std::uint64_t> maxEntries;
		std::vector<std::string> libraries;
		int option;
		while((option = getopt(argc, argv, "hi:o:c:t:n:l:")) != -1) {
			switch(option) {
				case 'h': Usage(argv[0]); return 0;
				case 'i': input = optarg; break;
				case 'o': outputDirectory = optarg; break;
				case 'c': compression = Compression(optarg); break;
				case 't': treeName = optarg; break;
				case 'l': libraries.emplace_back(optarg); break;
				case 'n': {
					std::uint64_t value;
					const auto end = optarg + std::strlen(optarg);
					const auto result = std::from_chars(optarg, end, value);
					if(result.ec != std::errc{} || result.ptr != end
						|| value > static_cast<std::uint64_t>(std::numeric_limits<Long64_t>::max()))
						throw std::runtime_error("Invalid entry limit: " + std::string{optarg});
					maxEntries = value;
					break;
				}
				default: Usage(argv[0]); return 2;
			}
		}
		if(input.empty() || optind != argc) {
			Usage(argv[0]);
			return 2;
		}
		if(outputDirectory.empty()) outputDirectory = input.parent_path() / "rn";
		const fs::path output = outputDirectory / input.filename().string();
		if(fs::weakly_canonical(input) == fs::weakly_canonical(output)
			|| (fs::exists(output) && fs::equivalent(input, output)))
			throw std::runtime_error("Input and output must be different files.");
		if(fs::exists(output)) throw std::runtime_error("Output already exists: " + output.string());
		for(const auto& library : libraries)
			if(gSystem->Load(library.c_str()) < 0)
				throw std::runtime_error("Cannot load library: " + library);
		std::unique_ptr<TFile> file{TFile::Open(input.c_str(), "READ")};
		if(!file || file->IsZombie()) throw std::runtime_error("Cannot open input file: " + input.string());
		if(treeName.empty()) {
			std::vector<std::string> names;
			FindTrees(*file, "", names);
			if(names.empty()) throw std::runtime_error("No TTree found in input file.");
			if(names.size() > 1) {
				std::string message = "Multiple TTrees found; select one with -t:";
				for(const auto& name : names) message += " " + name;
				throw std::runtime_error(message);
			}
			treeName = names.front();
		}
		auto tree = file->Get<TTree>(treeName.c_str());
		if(!tree) throw std::runtime_error("Cannot find TTree: " + treeName);
		DetectorTree detectorTree{*tree};
		fs::create_directories(outputDirectory);
		auto importer = ROOT::Experimental::RNTupleImporter::Create(&detectorTree, output.string());
		importer->SetConvertDotsInBranchNames(true);
		auto options = importer->GetWriteOptions();
		options.SetCompression(compression);
		importer->SetWriteOptions(options);
		if(maxEntries) importer->SetMaxEntries(*maxEntries);
		importer->Import();
		return 0;
	} catch(const std::exception& error) {
		std::fprintf(stderr, "Conversion failed: %s\n", error.what());
		return 1;
	}
}
