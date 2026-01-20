#include "ROOT/RNTupleModel.hxx"
#include "ROOT/RNTupleReader.hxx"
#include "ROOT/RNTupleWriter.hxx"
#include "ROOT/RCanvas.hxx"
#include "ROOT/RNTupleDS.hxx"
#include "ROOT/RDataFrame.hxx"
using namespace ROOT;
using namespace ROOT::Experimental;
using nlohmann::json;

template<typename T>
void PrintVec(const std::vector<T>& );
template<typename T>
T AvgVec(const std::vector<T>& );

void dump_json_custom(const json& j, std::ostream& out, int indent, int level = 0);

static const char* labels[] = {
	"21", "22", "23", "24", "31", "41", "42"
};

void tpc_alignment(std::string prefix = "", std::vector<std::string> fileNames = {}, int i_tpc=0) {
	if(i_tpc < 0 or i_tpc >= RNFRSCal::N_VALID_TPC)
		throw std::runtime_error("Bad i_tpc arg");

	ROOT::EnableImplicitMT();

	using Measurement = RNTPCCal::Measurement;

	TH1I *h1_x[2], *h1_y[4];
	for(int i=0; i<2; ++i)
		h1_x[i] = new TH1I(Form("h1_TPC%d_x%d", i_tpc, i), Form("h1_TPC%d_x%d", i_tpc, i), 550, -100, 100);

	for(int i=0; i<4; ++i)
		h1_y[i] = new TH1I(Form("h1_TPC%d_y%d", i_tpc, i), Form("h1_TPC%d_y%d", i_tpc, i), 380, -60, 60);

	std::vector<double> xs[2], ys[4];
	
	std::vector<TH1I*> ext_h_x[2];
	std::vector<TH1I*> ext_h_y[4];

	for(const auto& fileName: fileNames) { 
		auto model = RNTupleModel::Create();
		auto frs = model->MakeField<RNFRSCal>("FRS"); // shared_ptr.
		auto ntuple = RNTupleReader::Open(std::move(model), "h103", prefix+fileName+".root");

		for(auto entryId : *ntuple) {
			ntuple->LoadEntry(entryId);
			auto& tpc = frs->tpc.at(i_tpc);
			
			const std::array<std::vector<Measurement>, 2>& tpc_hits = tpc.hits;
			for(int d : {0,1}) {
				const std::vector<Measurement>& m = tpc_hits[d];
				if(m.size() != 1)
					continue;
				
				if(!std::isnan(m[0].x))
					h1_x[d]->Fill(m[0].x);

				for(int a : {0,1}) {
					if(m[0].mask >> a)
						h1_y[2*d + a]->Fill(m[0].y);
				}
			}
		}

		/* Fit Gauss around those values. */
		for(int i=0; i<2; ++i) {
			auto* h = h1_x[i];
			int maxBin = h->GetMaximumBin();
			double fitMin = h->GetBinCenter(maxBin - 20);
			double fitMax = h->GetBinCenter(maxBin + 20);
			h->Fit("gaus", "Q", "", fitMin, fitMax);
			TF1* fitF = h->GetFunction("gaus");
	
			xs[i].push_back(fitF->GetParameter(1));
			
			ext_h_x[i].push_back( (TH1I*)h->Clone() );

			h->GetListOfFunctions()->Delete();
			h->Reset("ICESM");
		}

		for(int i=0; i<4; ++i) {
			auto* h = h1_y[i];
			int maxBin = h->GetMaximumBin();
			double fitMin = h->GetBinCenter(maxBin - 30);
			double fitMax = h->GetBinCenter(maxBin + 30);
			h->Fit("gaus", "Q", "", fitMin, fitMax);
			TF1* fitF = h->GetFunction("gaus");

			ys[i].push_back(fitF->GetParameter(1));

			ext_h_y[i].push_back( (TH1I*)h->Clone() );
			
			h->GetListOfFunctions()->Delete();
			h->Reset("ICESM");
		}
	}

	printf("Corresponding values: \n");
	for(int i=0; i<2; ++i) {
		printf("x[%d]: ---\n", i);
		PrintVec(xs[i]);
	}
	printf("\n");
	for(int i=0; i<4; ++i) {
		printf("y[%d]: ---\n", i);
		PrintVec(ys[i]);
	}

	// Draw a canvas for each TOF calibration file.
	const int N = (int)ext_h_x[0].size();
	for(int i=0; i<N; ++i) {
		TCanvas* cX = new TCanvas(Form("cTOF%d_x", i), Form("cTOF%d_x", i),1200, 800);
		cX->Divide(1,2);
		for(int xi = 0; xi < 2; ++xi) {
			cX->cd(xi+1);
			ext_h_x[xi].at(i)->Draw();
		}

		TCanvas* cY = new TCanvas(Form("cTOF%d_y", i), Form("cTOF%d_y", i), 1200, 800);
		cY->Divide(2,2);
		for(int yi = 0; yi < 4; ++yi) {
			cY->cd(yi+1);
			ext_h_y[yi].at(i)->Draw();
		}
	}

	/* Also output the suggested JSON file: */
	std::string oneFile = prefix + fileNames[0] + ".root";
	std::unique_ptr<TFile> f = std::make_unique<TFile>(oneFile.c_str(), "READ");
	
	std::string setup = *f->Get<std::string>("FRS_setup_file");
	
	std::ifstream input(setup);
	json jfull = json::parse(input);
	json& j = jfull.at("TPC").at( labels[i_tpc] );
	for(int i=0; i<2; ++i) {
		double& off = j["x_offset"][i].get_ref<double&>();
		off -= AvgVec( xs[i] );
	}
	for(int i=0; i<4; ++i) {
		double& off = j["y_offset"][i].get_ref<double&>();
		off -= AvgVec( ys[i] );
	}

	printf("\n~~ Suggested changes : \n");
	std::cout << "\"" << labels[i_tpc] << "\": ";

	dump_json_custom(j, std::cout, 4);
	std::cout << std::endl;
}

template<typename T>
void PrintVec(const std::vector<T>& v) {
	std::cout << "[";
	for(auto it = v.begin(); it < v.end()-1; ++it) 
		std::cout << *it << ", ";
	if(v.size() > 0)
		std::cout<< v.back();
	std::cout << "]\n";
}

template<typename T>
T AvgVec(const std::vector<T>& v) {
	return static_cast<T> (
		std::accumulate(v.begin(), v.end(), static_cast<T>(0)) / v.size()
	);
}

void dump_json_custom(const json& j, std::ostream& out, int indent, int level) {
	auto pad = [&]() {
		for(int i = 0; i < level * indent; ++i) out.put(' ');
	};

	if(j.is_object()) {
		out << "{\n";
        bool first = true;
        for(auto it = j.begin(); it != j.end(); ++it) {
            if(!first) out << ",\n";
            first = false;

            pad(); 
			out << std::string(indent, ' ');
            
			out << json(it.key()).dump() << ": ";
            dump_json_custom(it.value(), out, indent, level + 1);
        }
        out << "\n"; pad(); out << "}";
	}
    else if(j.is_array()) {
        out << j.dump(); // inline entire numeric array tree (including nested arrays)
    }
	else {
		out << j.dump();
	}
}
