#pragma once

#include <filesystem>
#include <sstream>
#include <stdexcept>

#include "TROOT.h"
#include "TCanvas.h"
#include "TInterpreter.h"

using A2 = std::array<double, 2>;
using A3 = std::array<double, 3>;

enum class DoSave { yes, no };

namespace canvas {
	enum struct Extension { png, jpeg, pdf, C, root };

	inline void save_all (
		Extension extension = Extension::png, 
		std::vector<std::string_view> extra_tag = {}
	) { 
		const char* ext;
#define HANDLE_CASE(e) case(Extension::e): { ext = #e; break; }

		switch(extension) {
			HANDLE_CASE(png)
			HANDLE_CASE(jpeg)
			HANDLE_CASE(pdf)
			HANDLE_CASE(C)
			HANDLE_CASE(root)
		}

		const char* macro_name = gInterpreter->GetCurrentMacroName(); 
		std::string stem = std::filesystem::path(macro_name).stem().string();

		std::filesystem::path p = "autosave";
		p = p / stem / ext;
		for(auto& tag: extra_tag)
			p /= tag;

		try {
			std::filesystem::create_directories(p); 
		} catch(const std::filesystem::filesystem_error& e) {
			fprintf(stderr, "canvas::save_all : Error in creating directories: \'%s\', err: %s\n", p.c_str(), e.what());  
			return;
		}

		std::vector<TCanvas*> cs;
		for(TObject* k_ : *gROOT->GetListOfCanvases()) {
			if(TCanvas* c = dynamic_cast<TCanvas*>(k_))
				cs.push_back(c);
		}
		if(cs.empty()) return;

		for(auto* c : cs) {
			auto name = std::string(c->GetName()); 
			auto outfile = p / (name + "." + ext);
			// Force rendering
			c->Modified();
			c->Update();
			c->SaveAs( outfile.c_str() ); 
		}
		if(cs.size() == 1 or (extension != Extension::png and extension != Extension::jpeg)) return;

		// In case of two or more canvases saved, also collect them into a .pdf
		std::string outpdf = (p / "all.pdf").string();

		cs.front()->Print(Form("%s(", outpdf.c_str()));
		for(size_t i=1; i < cs.size() - 1; ++i) {
			cs[i]->Print( outpdf.c_str() );
		} 
		cs.back()->Print(Form("%s)", outpdf.c_str()));
	}
};

struct MeanStddev {
	double mean, stddev;
	inline operator std::pair<double,double>() noexcept { return { mean, stddev }; }

	friend std::ostream& operator<<(std::ostream&, const MeanStddev& ) noexcept;

	inline std::string string() const noexcept { 
		std::stringstream ss;
		ss << *this;
		return ss.str();
	}

	inline std::string lstring() const noexcept {
		std::stringstream ss;
		ss << mean << " #pm " << stddev;
		return ss.str();
	}
};

inline std::ostream& operator<<(std::ostream& os, const MeanStddev& rhs) noexcept {
	os << rhs.mean << " ± " << rhs.stddev;
	return os;
}

template<typename T, typename = void>
struct is_range : std::false_type {};

template<typename T>
struct is_range<T, std::void_t <
	decltype( std::cbegin(std::declval<T>()) ),
	decltype( std::cend(std::declval<T>()) )
>> : std::true_type {};

template <typename Range,
	typename = std::enable_if_t<is_range<Range>::value>
> MeanStddev mean_stddev(const Range& r) {
	using std::end;
	auto first = std::cbegin(r);
	auto last = std::cend(r);
	const auto n = std::distance(first, last);
	if(n <= 0) return { NAN, NAN };

	double sum = 0.0;
	for(auto it = first; it != last; ++it) {
		sum += *it;
	}
	double mean = sum / static_cast<double>(n);

	double sq_sum = 0.0;
	for(auto it = first; it != last; ++it) {
		double d = *it - mean;
		sq_sum += d * d;
	}
	double stddev = std::sqrt(sq_sum / static_cast<double>(n-1));
	return {mean, stddev};
}
