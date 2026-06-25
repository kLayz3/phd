#pragma once
#include "../monad/monad.hxx"

#include "TFile.h"
#include "TObject.h"
#include <iostream>
#include <variant>
#include <filesystem>
#include <sstream>
#include <cmath>
#include <stdexcept>

#include "TROOT.h"
#include "TCanvas.h"
#include "TInterpreter.h"

namespace _detail {
inline TFile* file_ptr(TFile* f) noexcept {
    return f;
}
inline TFile* file_ptr(std::unique_ptr<TFile> const& f) noexcept {
    return f.get();
}
}
/* Handle can be either unique ptr, or standard pointer. 
 * `var` must be a raw pointer! */
template<typename F, typename P>
void get_obj(F&& fhandle, P& var, const char* label) {
	static_assert(std::is_pointer_v<P>, "get_obj(): var must be a raw pointer");

	using T = std::remove_pointer_t<P>;
	
	TFile* f = _detail::file_ptr(fhandle);

	if constexpr (std::is_base_of_v<TObject, T>) {
		var = dynamic_cast<T*>(f->Get(label));
	} else {
		var = f->Get<T>(label);
	}

	if (!var) {
		std::cerr << "get_obj(): cannot extract object: '" << label << "'\n";
		std::abort();
	}
}

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
#define HANDLE_CASE_SAVE_ALL(e) case(Extension::e): { ext = #e; break; }

		switch(extension) {
			HANDLE_CASE_SAVE_ALL(png)
			HANDLE_CASE_SAVE_ALL(jpeg)
			HANDLE_CASE_SAVE_ALL(pdf)
			HANDLE_CASE_SAVE_ALL(C)
			HANDLE_CASE_SAVE_ALL(root)
		}
#undef HANDLE_CASE_SAVE_ALL

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

namespace mnd {

/* Nicer API to name different inputs. */

template<typename T, typename Tag = void, bool CanInherit = std::is_class_v<T>>
struct InputWrapper;

template<typename T, typename Tag>
struct InputWrapper<T,Tag,true> : T {
	using T::T;
	using underlying_type = T;
	
	T& get() noexcept { return *this; }
	const T& get() const noexcept { return *this; }
};

/* Arrays are aggregates and don't provide any ctors by default. How nice. */
template<typename U, std::size_t N, typename Tag>
struct InputWrapper<std::array<U,N>, Tag, true> : std::array<U,N> {
	using underlying_type = std::array<U,N>;

	InputWrapper() = default;
	template<typename... Args,
		 typename = typename std::enable_if_t<sizeof...(Args) == N>
	>
    InputWrapper(Args&&... args)
        : underlying_type{ { static_cast<U>(std::forward<Args>(args))... } }
    {}

	underlying_type& get() noexcept { return *this; }
	const underlying_type& get() const noexcept { return *this; }
};

template<typename T, typename Tag>
struct InputWrapper<T, Tag, false> {
	using underlying_type = T;
	T value;

	InputWrapper() = default;

	InputWrapper(T v) : value(v) {}

	operator T&() noexcept { return value; }
	operator const T&() const noexcept { return value; }

	T& get() noexcept { return value; }
	const T& get() const noexcept { return value; }
};

} // namespace mnd

/* Parse a file first thru the GCC preprocessor, and then
 * try to parse the output as a sequence of lines. 
 * Is not thread safe! */
std::vector<std::string> ParseFile(const std::string& );
extern std::vector<std::string> ParseFile(const std::string& );
