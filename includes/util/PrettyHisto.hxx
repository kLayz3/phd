#pragma once

#include <cstring>
#include <sstream>
#include "TH1D.h"
#include "TH2D.h"
#include "TColor.h"
#include <stdexcept>
#include <string_view>
#include <type_traits>
#include "TPad.h"
#include "TLatex.h"
#include "TGClient.h"
#include "TLine.h"
#include "TCanvas.h"

struct ORGB { uint32_t v; };

namespace ph_detail {

static inline const char* skip_whitespace(const char* str) {
	while(*str and isspace(*str)) {
		++str;
	}
	return str;
}

static inline std::string trim(std::string_view s) {
	std::stringstream ss{};
	auto is_ws = [](const char c) -> bool { return std::isspace(c); };
	while(!s.empty() && !is_ws(s.front())) ss << s;
	return ss.str();
}

static inline std::string strip_square_brackets(std::string_view s) {
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
static inline std::string nonalnum_to_underscore(std::string_view s) {
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

struct __ARGB {
	double a = 1.0, r, g, b;
	Color_t id;

	__ARGB() = default;
	__ARGB(ORGB val) {
		b = (val.v & 0xff) / 255.0; val.v >>= 8;
		g = (val.v & 0xff) / 255.0; val.v >>= 8;
		r = (val.v & 0xff) / 255.0; val.v >>= 8;
		a = std::max(1.0 - (val.v & 0xff) / 255.0, 0.0);
		mode_ = Mode::argb;
	}
	__ARGB(Color_t rc) {
		id = rc;
		mode_ = Mode::root_col;
	}

	void ApplyFill(TH1* h) {
		Int_t idx;
		switch(mode_) {
			case Mode::empty: return;
			case Mode::argb: {
				idx = TColor::GetColor((Float_t)r, (Float_t)g, (Float_t)b);
				break;
			}
			case Mode::root_col: {
				idx = id;
			}
		}
		/* Sometimes setting alpha doesn't work. The transparency is implicitly enabled through
		 * $ROOTSYS/etc/system.rootrc
		 * OpenGL.CanvasPreferGL = 1
		 * , if it is 0, then alpha will always be total (1.0). */
		h->SetFillColorAlpha(idx, a);
	}

	enum class Mode { empty, argb, root_col } mode_ = Mode::empty;
};

}
#define FWD_DRAW(inner) \
	template<typename... Ts> \
	auto Draw(Ts&&... args) { \
		h.Draw( std::forward<Ts>(args)... ); \
		if(gPad) gPad->SetGrid(); \
	}

#define FWD_FCN(fcn, inner) \
	template<typename... Ts> \
	auto fcn(Ts&&... args) { return inner.fcn( std::forward<Ts>(args)... ); }

struct TH1P {
	using inner_type = TH1D;
	
	inner_type h;
		
	template<typename... Ts>
	TH1P(const char* label, ph_detail::__ARGB col, Ts&&... args) {
		const char* semicolon = strchr(label, ';');
		if(semicolon)
			throw std::invalid_argument("Label input to 'TH1P' must not contain ';' semicolon delimiter!");

		std::string xlabel, hname_extra, title_extra;
		
		// Possible `(( ))` block
		const char* bracket_open = strstr(label, "((");
		if(bracket_open !=  nullptr) {
			const char* bracket_close = strstr(bracket_open, "))");
			if(bracket_close == nullptr)
				throw std::invalid_argument("Label input to 'TH1P' contains bracket open \"((\" "
					"delimiter but not the closing one \"))\" after the opening bracket.");
			if(std::distance(bracket_open, bracket_close) <= 2)
				throw std::invalid_argument("Label input to 'TH1P' has empty block inside \"((\" and \"))\" brackets?");
			hname_extra = std::string(bracket_open+2, bracket_close);
			label = bracket_close+2;
			label = ph_detail::skip_whitespace(label);
		}

		// Possible `@` separator 
		const char* at = strchr(label, '@');
		if(at != nullptr) {
			xlabel = std::string(label, at);
			title_extra = std::string(at+1);
		} else {
			xlabel = std::string(label);
		}
		
		std::string xs = ph_detail::strip_square_brackets( xlabel );
		std::string hname = hname_extra + "_"
			+ ph_detail::nonalnum_to_underscore(xs);
		
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
	FWD_DRAW(h);
	FWD_FCN(Fill, h);
	
	inline bool IsInside(const double x) const {
		return (
			x >= h.GetXaxis()->GetXmin() &&
			x <  h.GetXaxis()->GetXmax()
		);
	}

	/* Implicit ref cvt */
	operator inner_type&()             noexcept { return h; }
	operator const inner_type&() const noexcept { return h; }

	/* Implicit pointer cvt */
	operator inner_type*()             noexcept { return &h; }
	operator const inner_type*() const noexcept { return &h; }

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
			throw std::invalid_argument("Label input to 'TH2P' must contain ':' delimiter!");
		if(semicolon)
			throw std::invalid_argument("Label input to 'TH2P' must not contain ';' semicolon delimiter!");

		std::string xlabel, ylabel, hname_extra, title_extra;
		
		// Possible `(( ))` block
		const char* bracket_open = strstr(label, "((");
		if(bracket_open !=  nullptr) {
			const char* bracket_close = strstr(bracket_open, "))");
			if(bracket_close == nullptr)
				throw std::invalid_argument("Label input to 'TH1P' contains bracket open \"((\" "
					"delimiter but not the closing one \"))\" after the opening bracket.");
			if(std::distance(bracket_open, bracket_close) <= 2)
				throw std::invalid_argument("Label input to 'TH1P' has empty block inside \"((\" and \"))\" brackets?");
			hname_extra = std::string(bracket_open+2, bracket_close);
			label = bracket_close+2;
			label = ph_detail::skip_whitespace(label);
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
		std::string ys = ph_detail::strip_square_brackets( ylabel );
		std::string xs = ph_detail::strip_square_brackets( xlabel );
		
		std::string hname = hname_extra + "_" 
			+ ph_detail::nonalnum_to_underscore(ys) + ".v." + ph_detail::nonalnum_to_underscore(xs);
		
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
	FWD_DRAW(h);
	FWD_FCN(Fill, h);

	inline bool IsInside(const double x, const double y) const {
		return (
			x >= h.GetXaxis()->GetXmin() &&
			x <  h.GetXaxis()->GetXmax() &&
			y >= h.GetYaxis()->GetXmin() &&
			y <  h.GetYaxis()->GetXmax()
		);
	}

	/* Implicit ref cvt */
	operator inner_type&()             noexcept { return h; }
    operator const inner_type&() const noexcept { return h; }

	/* Implicit pointer cvt */
	operator inner_type*()             noexcept { return &h; }
    operator const inner_type*() const noexcept { return &h; }

	inline inner_type* operator->() noexcept { return &h; }
	inline inner_type& operator*()  noexcept { return h; }
	inline const inner_type* operator->() const noexcept { return &h; }
	inline const inner_type& operator*()  const noexcept { return h; }
};

/* This instance is thread-unsafe! But anyway if some psycho like me uses latex'ing
 * in multithreaded/async, they should start touching grass. Honest opinion. */
struct PLatex {
	static constexpr double MAX_H = 0.9;
	static constexpr double MIN_H = 0.1;
	static constexpr double START_X = 0.1;
	uint32_t N;
	uint32_t current_text_index;
	TLatex latex;
    std::string __temporary {};

	/* Ts is one the: std::string/std::string_view/const char* */
	template<typename... Ts>
	PLatex(double textsize, Ts&&... texts) : N( static_cast<uint32_t>(sizeof...(Ts))) 
	{
		static_assert(sizeof...(Ts) > 0, "Cannot create PLatex from empty text box!");
		assert(textsize > 0 && "First argument to PLatex ctor must be textsize > 0.");

		latex.SetNDC();
		latex.SetTextSize(textsize);
		current_text_index = 0;
		(..., DrawOne( std::forward<Ts>(texts) ));
		gPad->Modified();
		gPad->Update();
	}
	template<typename T>
	void DrawOne(T&& text) {
		double height = MAX_H - current_text_index * (MAX_H - MIN_H) / (N-1);
		if(!std::isfinite(height)) height = (MAX_H + MIN_H) / 2;

		latex.DrawLatex(START_X, height, to_const_char( std::forward<T>(text) ));
		++current_text_index;
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
		for(const auto [raw, rep] : rules)
			PLatex::replace_all(__temporary, raw, rep);
	}
	static void replace_all(std::string& s, const char* sub, const char* rep) {
		if(!sub || !*sub) return;
		size_t pos = 0;
		
		while((pos = s.find(sub, pos)) != std::string::npos) {
			s.replace(pos, strlen(sub), rep);
			pos += strlen(rep);
		}
	}

	/* Phew, why are humans still using that god-awful `const char*` API..? */
};

inline void MaximizeCanvas(TCanvas* c) {
	c->SetWindowSize (
		gClient->GetDisplayWidth(),
		gClient->GetDisplayHeight()
	);
}
/* TODO */
struct PCanvas {
	using inner_type = TCanvas;
	inner_type c;
	
	template<typename... Ts>
	PCanvas(const char* label, Ts&&... args) {
		const char* semicolon = strchr(label, ';');
		if(!semicolon)
			throw std::invalid_argument("Label input to 'PCanvas' must contain ';' semicolon delimiter");
		
		std::string obj_name = std::string(label, semicolon);
		std::string title    = std::string(semicolon + 1);	
	}

	/* Implicit ref cvt */
	operator inner_type&()             noexcept { return c; }
	operator const inner_type&() const noexcept { return c; }

	/* Implicit pointer cvt */
	operator inner_type*()             noexcept { return &c; }
	operator const inner_type*() const noexcept { return &c; }

	inline inner_type* operator->() noexcept { return &c; }
	inline inner_type& operator*()  noexcept { return c; }
	inline const inner_type* operator->() const noexcept { return &c; }
	inline const inner_type& operator*()  const noexcept { return c; }

};

namespace hist {

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

}
