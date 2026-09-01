#pragma once

#include <cstring>
#include <sstream>
#include <stdexcept>
#include <string_view>
#include <type_traits>
#include <variant>
#include <filesystem>

#include "RtypesCore.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TColor.h"
#include "TPad.h"
#include "TLatex.h"
#include "TGClient.h"
#include "TLine.h"
#include "TCanvas.h"

#include "GaussFitMax.hxx"

#if !defined(MND_INCLUDE_SPAN_IS_DEFINED)
#define MND_INCLUDE_SPAN_IS_DEFINED
#if __cplusplus >= 202000L 
#	include <span>
	namespace mnd {
		template<typename T>
		using span = std::span<T>;
	}
#elif __has_include("boost/beast/core/span.hpp")
#	include "boost/beast/core/span.hpp"
	namespace mnd {
		template<typename T>
		using span = boost::beast::span<T>;
	}
#else
#	error "Neither C++20 given, nor boost span library found. Cannot proceed."
#endif
#endif // MND_INCLUDE_SPAN_IS_DEFINED

struct ORGB { uint32_t v; };

namespace mnd::detail {

inline const char* skip_whitespace(const char* str) noexcept {
	while(*str and isspace(*str)) {
		++str;
	}
	return str;
}

inline std::string trim(std::string_view s) noexcept {
	std::stringstream ss{};
	auto is_ws = [](const char c) -> bool { return std::isspace(c); };
	while(!s.empty() && !is_ws(s.front())) ss << s;
	return ss.str();
}

inline std::string strip_square_brackets(std::string_view s) noexcept {
	std::string out;
	out.reserve(s.size());

	bool in_brackets = false;
	for(char ch : s) {
		if(!in_brackets) {
			if(ch == '[') {
				in_brackets = true;
			} else {
				out.push_back(ch);
			}
		} 
		else {
			if(ch == ']') in_brackets = false;
		}
	}
	return out;
}

/* Collapse all consecutive sequences of non-alphanumberic char to a single '_'.
 * Also trims leading/trailing whitespace by default. */
inline std::string nonalnum_to_underscore(std::string_view s) noexcept {
	std::string out;
	out.reserve(s.size());

	auto is_non_alnum = [](const char c) -> bool { return !std::isalnum(c); };

	bool prev_was_ws = true; // treat leading ws as "already ws" -> no leading underscore
	for(char ch : s) {
		if(is_non_alnum(ch)) {
			if(!prev_was_ws) out.push_back('_');
			prev_was_ws = true;
		} else {
			out.push_back( (char)std::tolower(ch) );
			prev_was_ws = false;
		}
	}

	// remove trailing underscores
	while(!out.empty() && out.back() == '_') out.pop_back();
	
	return out;
}

} // namespace mnd::detail

namespace mnd::col {

struct RGBA {
	double r = 0.0;
	double g = 0.0;
	double b = 0.0;
	double a = 1.0;

	RGBA() = default;
	RGBA(double r_, double g_, double b_, double a_ = 1.0) :
		r(r_), g(g_), b(b_), a(a_) {}

	/* ROOT palette/index color, e.g. kRed+1. */
	RGBA(Color_t );
	/* No RBGA(uint32_t ) ctor, as that would be ambiguous. */

	/* Preserve the old ORGB packed representation: 0xTTRRGGBB, where
	 * TT is transparency (00 = opaque, ff = fully transparent). */
	static RGBA from_packed(uint32_t value) noexcept;

	Int_t GetColorCode() const;
	void ApplyFill(TH1* h) const;
};
namespace literals {
	RGBA operator""_c(unsigned long long int );
}

static inline Int_t hex_to_col(uint32_t val) {
	Float_t b = (val & 0xff) / 255.0; val >>= 8;
	Float_t g = (val & 0xff) / 255.0; val >>= 8;
	Float_t r = (val & 0xff) / 255.0; val >>= 8;
	return TColor::GetColor(r,g,b);
}

inline Int_t Col(uint32_t i) {
	static uint32_t cols[] = {
		0xC41E3A, 0xA330C9, 0xFF7C0A, 0x33937F,
		0xAAD372, 0x3FC7EB, 0x00FF98, 0xF48CBA,
		0xFFF468, 0x0070DD, 0x8788EE, 0xC69B6D
	};
	constexpr static size_t Ncols = sizeof cols / sizeof *cols;
	
	return hex_to_col( cols[i % Ncols] );
}

} // namespace mnd::col

#define MND_FWD_DRAW(inner) \
	template<typename... Ts> \
	auto Draw(Ts&&... args) { \
		h.Draw( std::forward<Ts>(args)... ); \
		if(gPad) gPad->SetGrid(); \
	}

#define MND_FWD_FCN(fcn, inner) \
	template<typename... Ts> \
	auto fcn(Ts&&... args) { return inner.fcn( std::forward<Ts>(args)... ); }

struct TH1P {
	using inner_type = TH1D;
	
	inner_type h;

	template<typename... Ts>
	TH1P(const char* label, mnd::col::RGBA col, Ts&&... args) {
		const char* semicolon = strchr(label, ';');
		if(semicolon)
			throw std::invalid_argument(Form("Label input to 'TH1P' must not contain ';' semicolon delimiter! Received: \'%s\'", label));

		std::string xlabel, hname_extra, title_extra;
		
		// Possible `(( ))` block
		const char* bracket_open = strstr(label, "((");
		if(bracket_open !=  nullptr) {
			const char* bracket_close = strstr(bracket_open, "))");
			if(bracket_close == nullptr)
				throw std::invalid_argument(Form("Label input to 'TH1P' contains bracket open \"((\" "
					"delimiter but not the closing one \"))\" after the opening bracket. Received: \'%s\'", label));
			if(std::distance(bracket_open, bracket_close) <= 2)
				throw std::invalid_argument(Form("Label input to 'TH1P' has empty block inside \"((\" and \"))\" brackets? Received: \'%s\'", label));
			hname_extra = std::string(bracket_open+2, bracket_close);
			label = bracket_close+2;
			label = mnd::detail::skip_whitespace(label);
		}

		// Possible `@` separator 
		const char* at = strchr(label, '@');
		if(at != nullptr) {
			xlabel = std::string(label, at);
			title_extra = std::string(at+1);
		} else {
			xlabel = std::string(label);
		}
		
		std::string xs = mnd::detail::strip_square_brackets( xlabel );
		std::string hname = hname_extra + "_"
			+ mnd::detail::nonalnum_to_underscore(xs);
		
		h = inner_type(hname.c_str(), "", std::forward<Ts>(args)...); 
		std::stringstream title;
		title << xs;
		if(title_extra.length() > 0)
			title << " (" << title_extra << ')';
		auto title_materialied = title.str();
		h.SetTitle(Form("%s;%s;%s", title_materialied.c_str(), xlabel.c_str(), "Count") );
		h.GetYaxis()->SetTitleOffset(1.0);

		col.ApplyFill(&h);
		h.SetLineColor(kBlack);
		h.SetLineWidth(2);
	}

	/* Forward only Fill and Draw methods. Don't care about others. */
	MND_FWD_DRAW(h);
	MND_FWD_FCN(Fill, h);

	/* Draw while also fitting a small gauss-chan 🥺 👉👈 around the peak value. */
	inline auto DrawAndFit (
		double side_ratio = GAUSS_FIT_SIDE_RATIO_DEFAULT,
		mnd::col::RGBA col = kRed,
		Width_t line_width = 2, // type: short
		uint32_t niter = 2,
		Verbosity v = Verbosity::SILENT
	) -> decltype( GaussFitMax(std::declval<TH1D*>()) ) {
		auto fitresult   = GaussFitMax(&h, side_ratio, niter, v);
		const auto& res = fitresult.first;
		const double A     = res[0];
		const double mu    = res[1];
		const double sigma = res[2];

		const double m_  = h.GetXaxis()->GetBinCenter( h.GetMaximumBin() );
		const double s_  = h.GetStdDev();
		const double xlo = m_ - side_ratio * s_;
		const double xhi = m_ + side_ratio * s_;

		auto* f = new TF1 (
			Form("extern_gaus[%.1f]:%s", side_ratio, h.GetName()),
			"[0]*exp(-0.5*((x-[1])/[2])^2)",
			xlo, xhi
		);
		
		f->SetParameters(A, mu, sigma);
		f->SetParNames("A", "#mu", "#sigma");
		f->SetLineColor( col.GetColorCode() );
		f->SetLineWidth( line_width );

		this->Draw();
		f->Draw("same");

		return fitresult;
	}

	inline __attribute__((always_inline))
	bool IsInside(double x) const noexcept {
		return (
			x >= h.GetXaxis()->GetXmin() &&
			x <  h.GetXaxis()->GetXmax()
		);
	}
	
	inline auto FillInside(double x) {
		if(IsInside(x)) return h.Fill(x);
		else return -1;
	}
	inline auto FillInside(double x, double w) {
		if(IsInside(x)) return h.Fill(x, w);
		else return -1;
	}

	/* Implicit ref cvt */
	inline operator inner_type&()             noexcept { return h; }
	inline operator const inner_type&() const noexcept { return h; }

	/* Implicit pointer cvt */
	inline operator inner_type*()             noexcept { return &h; }
	inline operator const inner_type*() const noexcept { return &h; }

	inline inner_type* operator->() noexcept { return &h; }
	inline inner_type& operator*()  noexcept { return h; }
	inline const inner_type* operator->() const noexcept { return &h; }
	inline const inner_type& operator*()  const noexcept { return h; }
};

/* Wrapper to pretty-decorate standard CERN ROOT histograms
 * based on label and its rules. 
 * Ctor API example:
 * >> TH2P h("SCI21X [mm]:SCI22X [mm]", ... )
 * will create a interally a hist:
 * TH2D("sci21x.v.sci22x", "SCI21X vs. SCI22X", ... )
 * ... and set the titles appropriately. */
struct TH2P {
	using inner_type = TH2D;
	inner_type h;
	
	template<typename... Ts>
	TH2P(const char* label, Ts&&... args) {
		const char* colon = strchr(label, ':');
		const char* semicolon = strchr(label, ';');
		if(!colon)
			throw std::invalid_argument(Form("Label input to 'TH2P' must contain ':' delimiter! Received: \'%s\'", label));
		if(semicolon)
			throw std::invalid_argument(Form("Label input to 'TH2P' must not contain ';' semicolon delimiter! Received: \'%s\'", label));

		std::string xlabel, ylabel, hname_extra, title_extra;
		
		// Possible `(( ))` block
		const char* bracket_open = strstr(label, "((");
		if(bracket_open !=  nullptr) {
			const char* bracket_close = strstr(bracket_open, "))");
			if(bracket_close == nullptr)
				throw std::invalid_argument(Form("Label input to 'TH2P' contains bracket open \"((\" "
					"delimiter but not the closing one \"))\" after the opening bracket. Received: \'%s\'", label));
			if(std::distance(bracket_open, bracket_close) <= 2)
				throw std::invalid_argument(Form("Label input to 'TH2P' has empty block inside \"((\" and \"))\" brackets? Received: \'%s\'", label));
			hname_extra = std::string(bracket_open+2, bracket_close);
			label = bracket_close+2;
			label = mnd::detail::skip_whitespace(label);
		}

		ylabel = std::string(label, colon);

		// Possible `@` separator 
		const char* at = strchr(label, '@');
		if(at != nullptr) {
			xlabel = std::string(colon+1, at);
			title_extra = std::string(at+1);
		} else {
			xlabel = std::string(colon+1);
		}
		std::string ys = mnd::detail::strip_square_brackets( ylabel );
		std::string xs = mnd::detail::strip_square_brackets( xlabel );
		
		std::string hname = hname_extra + "_" 
			+ mnd::detail::nonalnum_to_underscore(ys) + ".v." + mnd::detail::nonalnum_to_underscore(xs);
		
		h = inner_type(hname.c_str(), "", std::forward<Ts>(args)...); 
		std::stringstream title;
		title << ys << " vs. " << xs;
		if(title_extra.length() > 0)
			title << " (" << title_extra << ')';
		auto title_materialied = title.str();
		h.SetTitle(Form("%s;%s;%s", title_materialied.c_str(), xlabel.c_str(), ylabel.c_str()) );
		h.GetYaxis()->SetTitleOffset(1.0);
	}

	/* Forward only Fill and Draw methods. Don't care about others. */
	MND_FWD_DRAW(h);
	MND_FWD_FCN(Fill, h);

	inline __attribute__((always_inline))
	bool IsInside(double x, double y) const noexcept {
		return (
			x >= h.GetXaxis()->GetXmin() &&
			x <  h.GetXaxis()->GetXmax() &&
			y >= h.GetYaxis()->GetXmin() &&
			y <  h.GetYaxis()->GetXmax()
		);
	}
	auto FillInside(double x, double y) {
		if(IsInside(x,y)) return h.Fill(x,y);
		else return -1;
	}
	auto FillInside(double x, double y, double w) {
		if(IsInside(x,y)) return h.Fill(x,y, w);
		else return -1;
	}

	/* Implicit ref cvt */
	operator inner_type&()             noexcept { return h; }
    operator const inner_type&() const noexcept { return h; }

	/* Implicit pointer cvt */
	operator inner_type*()             noexcept { return &h; }
    operator const inner_type*() const noexcept { return &h; }

	inner_type* operator->() noexcept { return &h; }
	inner_type& operator*()  noexcept { return h; }
	const inner_type* operator->() const noexcept { return &h; }
	const inner_type& operator*()  const noexcept { return h; }
};

namespace mnd::type_traits {

template<typename T, typename = void>
struct is_text_sequence : std::false_type {};

template<typename T>
struct is_text_sequence<T,
		std::void_t<
			typename std::decay_t<T>::value_type
		>
	> : std::bool_constant<
			std::is_convertible_v<
				typename std::decay_t<T>::value_type const&,
				std::string_view
			>
		> 
	{};

template<typename T>
inline constexpr bool is_text_sequence_v = is_text_sequence<T>::value;

} // namespace mnd::type_traits

/* This instance is thread-unsafe! But anyway if some psycho like me uses latex'ing
 * in multithreaded/async, they should start touching grass. Honest opinion.
 * Also - very badly written API. I apologise.. */
struct PLatex {
	static constexpr double MAX_H = 0.9;
	static constexpr double MIN_H = 0.1;
	static constexpr double START_X = 0.1;
	uint32_t N;
	uint32_t current_text_index;
	TLatex latex;
    std::string __temporary {};

	/* Ts is one the: std::string/std::string_view/const char* or a range of it. */
	template<typename... Ts>
	PLatex(double textsize, Ts&&... texts) : N((...+ TextCount(texts)))
	{
		assert(textsize > 0 && "First argument to PLatex ctor must be textsize > 0.");

		latex.SetNDC();
		latex.SetTextSize(textsize);
		current_text_index = 0;
		(..., DrawOne( std::forward<Ts>(texts) ));
		gPad->Modified();
		gPad->Update();
	}
	template<
		typename T,
	    typename std::enable_if_t<std::is_convertible_v<T, std::string_view>>* = nullptr
	> void DrawOne(T&& text) {
		double height = MAX_H - current_text_index * (MAX_H - MIN_H) / (N-1);
		if(!std::isfinite(height)) height = (MAX_H + MIN_H) / 2;

		latex.DrawLatex(START_X, height, to_const_char( std::forward<T>(text) ));
		++current_text_index;
	}
	/* Special overload for a sequence of texts. */
	template<typename T>
	void DrawOne(mnd::span<T const> text_sequence) {
		static_assert(std::is_convertible_v<T const&, std::string_view>);
		for(const auto& text : text_sequence) {
			DrawOne(std::string_view{text});
		}
	}

	/* Small routine to convert random-C++ text back into non-dangling `const char*` range.
	 * Also doing things such as replacing different color palettes and ansi codes
	 * to ROOT-valid styles. */
	template<typename T> 
	const char* to_const_char(T&& text) {
		/* Tricky, because std::string_view isn't null-terminated by default, and it cannot mutate
		 * the underlying object before taking a local copy first. However, if we materialize a proper std::string
		 * locally here, its `.c_str()` reference dangles when this call returns, and the caller needs a proper reference. */ 
		__temporary = std::string( std::forward<T>(text) );
		
		ParseProperTLatex();
		return __temporary.c_str();
	}

	inline void ParseProperTLatex() {
		static std::pair<const char*, const char*> rules[] {
			{ "\e[0m", "}" }, // normal
			{"\e[0;30m", Form("#color[%d]{", kBlack )},
			{"\e[0;31m", Form("#color[%d]{", kRed    + 2)},
			{"\e[0;32m", Form("#color[%d]{", kGreen  + 2)},
			{"\e[0;33m", Form("#color[%d]{", kYellow + 2)},
			{"\e[0;34m", Form("#color[%d]{", kBlue   + 2)},
			{"\e[0;35m", Form("#color[%d]{", kMagenta+ 2)},
			{"\e[0;36m", Form("#color[%d]{", kCyan   + 2)}
		};
		for(const auto& [raw, rep] : rules)
			PLatex::replace_all(__temporary, raw, rep);
	}
	inline static void replace_all(std::string& s, const char* sub, const char* rep) {
		if(!sub || !*sub) return;
		size_t pos = 0;
		
		while((pos = s.find(sub, pos)) != std::string::npos) {
			s.replace(pos, strlen(sub), rep);
			pos += strlen(rep);
		}
	}

private:
	template <typename T>
	static uint32_t TextCount(const T& x) {
		if constexpr(mnd::type_traits::is_text_sequence_v<T>) {
			return static_cast<uint32_t>(x.size());
		} else {
			return 1;
		}
	}
	
	/* Phew, why are humans still using that god-awful `const char*` API..? */
};

namespace mnd::hist {

inline double lo_y(TH2D* h2, double r=0) {
	if(r > 1 || r <= 0)
		throw std::runtime_error("lo_y::second arg must be <0, 1]");
	return  (1-r)/2 * h2->GetYaxis()->GetXmax() + (1+r)/2 * h2->GetYaxis()->GetXmin();
}
inline double hi_y(TH2D* h2, double r=0) {
	if(r > 1 || r <= 0)
		throw std::runtime_error("hi_y::second arg must be <0, 1]");
	return  (1+r)/2 * h2->GetYaxis()->GetXmax() + (1-r)/2 * h2->GetYaxis()->GetXmin();
}
inline double lo_x(TH2D* h2, double r=0) {
	if(r > 1 || r <= 0)
		throw std::runtime_error("lo_x::second arg must be <0, 1]");
	return  (1-r)/2 * h2->GetXaxis()->GetXmax() + (1+r)/2 * h2->GetXaxis()->GetXmin();
}
inline double hi_x(TH2D* h2, double r=0) {
	if(r > 1 || r <= 0)
		throw std::runtime_error("hi_x::second arg must be <0, 1]");
	return  (1+r)/2 * h2->GetXaxis()->GetXmax() + (1-r)/2 * h2->GetXaxis()->GetXmin();
}

inline double lo_y(TH2P* h2, double r=0) { return lo_y(&h2->h, r); };
inline double hi_y(TH2P* h2, double r=0) { return hi_y(&h2->h, r); };
inline double lo_x(TH2P* h2, double r=0) { return lo_x(&h2->h, r); };
inline double hi_x(TH2P* h2, double r=0) { return hi_x(&h2->h, r); };

[[nodiscard ]] 
inline TLine* vline(TH2D* h2, double x, double r = 0) {
	TLine* line = new TLine( x, lo_y(h2,r), x, hi_y(h2,r) );
	line->SetLineColor(kBlack);
	line->SetLineStyle(3);
	line->SetLineWidth(6);
	return line;
}
[[nodiscard ]] inline TLine* vline(TH2P* h2, double x, double r = 0) {
	return vline(&h2->h, x, r);
}

} // namespace mnd::hist

class TGraph;
class TGraphErrors;

namespace mnd::plot {

template<typename T>
using Maybe = std::optional<T>;
using mnd::col::RGBA;

enum class HistMode {
	stairs,
	errors,
	points,
};
enum class Hist2DMode {
	mesh,
	image,
	contour,
};
enum class GraphMode {
	plot,
	scatter,
	scatterplot,
};

struct HistStyle {
	Maybe<std::string> label_;
	Maybe<double>      line_width_;
	Maybe<std::string> line_style_;
	Maybe<std::string> marker_;
	Maybe<RGBA>        line_color_;
	Maybe<RGBA>        fill_color_;
	HistMode mode_ = HistMode::stairs;

	bool do_fill = false;

	auto label(std::string ) && -> HistStyle;
	auto line_width(double ) && -> HistStyle;
	auto line_style(std::string ) && -> HistStyle;
	auto marker(std::string ) && -> HistStyle;
	auto fill(bool = true) && -> HistStyle;
	auto facecolor(RGBA ) && -> HistStyle;
	auto edgecolor(RGBA ) && -> HistStyle;

	auto stairs() && -> HistStyle;
	auto errors() && -> HistStyle;
	auto points() && -> HistStyle;
};

struct Hist2DStyle {
	Maybe<std::string> label_;
	bool colorbar_ = true;
	Hist2DMode mode_ = Hist2DMode::mesh;

	auto label(std::string ) && -> Hist2DStyle;
	auto colorbar(bool = true) && -> Hist2DStyle;

	auto mesh() && -> Hist2DStyle;
	auto image() && -> Hist2DStyle;
	auto contour() && -> Hist2DStyle;
};

struct GraphStyle {
	Maybe<std::string> label_;
	Maybe<double> line_width_;
	Maybe<std::string> line_style_;
	Maybe<std::string> marker_;
	Maybe<double> marker_size_;
	Maybe<RGBA> color_;
	GraphMode mode_ = GraphMode::plot;

	auto label(std::string ) && -> GraphStyle;
	auto line_width(double ) && -> GraphStyle;
	auto line_style(std::string ) && -> GraphStyle;
	auto marker(std::string ) && -> GraphStyle;
	auto marker_size(double ) && -> GraphStyle;

	auto plot() && -> GraphStyle;
	auto scatter() && -> GraphStyle;
	auto scatterplot() && -> GraphStyle;
};

struct Hist1DData {
	std::vector<double> edges;
	std::vector<double> values;
	std::vector<double> errors;

	HistStyle style;
};

struct Hist2DData {
	std::vector<double> x_edges;
	std::vector<double> y_edges;

	// row-major:
	// values[iy * nx + ix]
	std::vector<double> values;

	std::size_t nx = 0;
	std::size_t ny = 0;

	Hist2DStyle style;
};

struct GraphData {
	std::vector<double> x;
	std::vector<double> y;

	std::vector<double> xerr;
	std::vector<double> yerr;

	GraphStyle style;
};

using Drawable = std::variant<
    Hist1DData,
    Hist2DData,
    GraphData
>;

struct Figure {
	Figure() = default;

	auto plot(const TH1& , HistStyle = {}) && -> Figure;
	auto plot(const TH2& , Hist2DStyle = {}) && -> Figure;
	auto plot(const TGraph& , GraphStyle = {}) && -> Figure;
	auto plot(const TGraphErrors& , GraphStyle = {}) && -> Figure;

	auto xlabel(std::string ) && -> Figure;
	auto xlabel(const TObject* ) && -> Figure;
	auto ylabel(std::string ) && -> Figure;
	auto ylabel(const TObject* ) && -> Figure;
	auto title(std::string ) && -> Figure;
	auto title(const TObject*) && -> Figure;

	auto xlim(double , double ) && -> Figure;
	auto ylim(double , double ) && -> Figure;

	auto logx(bool = true) && -> Figure;
	auto logy(bool = true) && -> Figure;
	auto grid(bool = true) && -> Figure;
	auto enable_right_top_spline(bool = true) && -> Figure;
	auto legend(bool = true) && -> Figure;

	auto save(const std::filesystem::path& ) const -> void;

private:
	std::vector<Drawable> objects_;

	std::string xlabel_;
	std::string ylabel_;
	std::string title_;

	Maybe<std::pair<double, double>> xlim_;
	Maybe<std::pair<double, double>> ylim_;
	bool hide_right_and_top_spline_ = true;

	bool logx_ = false;
	bool logy_ = false;
	bool grid_ = false;
	bool legend_ = false;
};

} // namespace mnd::plot

namespace mnd::python {

/* Poke the embedded Python interpreter. */
void poke(bool verbose = false);

}
