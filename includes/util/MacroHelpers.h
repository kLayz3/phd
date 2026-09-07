#pragma once

#include <iostream>
#include <string_view>
#include <variant>
#include <filesystem>
#include <sstream>
#include <cmath>
#include <stdexcept>

#include "TROOT.h"
#include "TCanvas.h"
#include "TInterpreter.h"
#include "TSystem.h"
#include "TImage.h"

#include "TFile.h"

#include "../monad/monad.hxx"
#include "../magic_enum/magic_enum.hpp"
#include "CLI.h"

namespace _detail {
inline TFile* file_ptr(TFile* f) noexcept {
    return f;
}
inline TFile* file_ptr(std::unique_ptr<TFile> const& f) noexcept {
    return f.get();
}
}
/* Handle can be either unique ptr, or standard pointer.
 * `var` must be a raw pointer! You own the object now, never the ROOT. */
template<typename F, typename P>
void get_obj(F&& fhandle, P& var, const char* label) {
	static_assert(std::is_pointer_v<P>, "get_obj(): var must be a raw pointer");

	using T = std::remove_pointer_t<P>;
	
	TFile* f = _detail::file_ptr(fhandle);

	if constexpr(std::is_base_of_v<TObject, T>) {
		var = dynamic_cast<T*>(f->Get(label));
		if constexpr(std::is_base_of_v<TH1, T>) {
			var->SetDirectory(nullptr);
		}
	} else {
		var = f->Get<T>(label);
	}

	if(!var) ERROR("get_obj(): cannot extract object: '%s'\n", label);
}

using A2 = std::array<double, 2>;
using A3 = std::array<double, 3>;
template<typename T, size_t M, size_t N>
using Arr2 = std::array<std::array<T,N>, M>;

enum class DoSave { yes, no };

/* File names are often of the form: `main_0XXX_0YYY.root`, as such
 * metadata'ing multiple files can be concatenated e.g.:
 * => main_0123_0144.root
 * +  main_0145_0174.root
 * +  main_0175_0178.root
 * ----------------------
 * =  main_0123_0178.root
 * */

namespace mnd::fs {

constexpr std::string_view FILE_PREFIX    = "main";
constexpr std::string_view FILE_EXTENSION = ".root";

/* If `main_0023_0144.root` => "0023"sv */
std::string_view file_start_number(const std::string& );

/* If `main_0023_0144.root` => "0023"sv */
std::string_view file_end_number(const std::string& );

/* If `main_0023_0144.root` => { "0023"sv, "0144"sv } */
std::pair <
    std::string_view,
    std::string_view
> file_number_bounds(const std::string& );

/* If `{ main_0023_0144.root, main_0145_0174.root }` => { "0023"sv, "0174"sv }.
 * It doesn't internally sort the sequence. Assumes sequence comes already sorted. */
std::pair <
    std::string_view,
    std::string_view
> file_number_bounds(const std::vector<std::string>& );

/* If `{ main_0023_0144.root, main_0145_0174.root }` => "main_0023_0174". Doesn't have the extension.
 * It doesn't internally sort the sequence. Assumes sequence comes already sorted. */
std::string file_names_concatenated(const std::vector<std::string>& );


} // namespace mnd::fs

namespace canvas {

/* Get all TObject-derived elements from a Canvas. */
void DumpPrimitives(TVirtualPad* , int = 0);
/* Extract all histogram objects from a Canvas. */
std::vector<TH1*> GetHistograms(TVirtualPad* );

enum struct Extension { png, jpeg, pdf, C, root, nil };

/* ROOT implements GetCurrentMacroName() in pre 6.38 as a simple forward to:
 * return fCurExecutingMacros.back();
 * Which, if called in a standalone program, simply segfaults on the spot. And there's no
 * public API to make a check. Thanks for this API friends. Take care. */
enum WhereAmI { Macro, Exe };

template<enum WhereAmI loc>
void save_all (
	Extension extension = Extension::png,
	std::vector<std::string_view> extra_tag = {}
) {
	const char* ext = "";
#define HANDLE_CASE_SAVE_ALL(e) case(Extension::e): { ext = #e; break; }

	switch(extension) {
		HANDLE_CASE_SAVE_ALL(png)
		HANDLE_CASE_SAVE_ALL(jpeg)
		HANDLE_CASE_SAVE_ALL(pdf)
		HANDLE_CASE_SAVE_ALL(C)
		HANDLE_CASE_SAVE_ALL(root)
		case(Extension::nil): return;
	}
#undef HANDLE_CASE_SAVE_ALL

	std::string stem{};
	if constexpr(loc == Macro) {
		const char* macro_name = gInterpreter->GetCurrentMacroName();
		stem = std::filesystem::path(macro_name).stem().string();
	}
	else { // standalone executable, or no interpreted macro currently active
		stem = mnd::fs::current_executable_path().stem().string();
	}

	std::filesystem::path p = "autosave";
	p = p / stem / ext;
	for(auto& tag: extra_tag)
		p /= tag;

	try {
		std::filesystem::create_directories(p);
	} catch(const std::filesystem::filesystem_error& e) {
		WARN("canvas::save_all : Error in creating directories: \'%s\', err: %s\n", p.c_str(), e.what());
		return;
	}

	std::vector<TCanvas*> cs;
	for(TObject* k_ : *gROOT->GetListOfCanvases()) {
		if(TCanvas* c = dynamic_cast<TCanvas*>(k_))
			cs.push_back(c);
	}
	if(cs.empty()) return;
	WARN("Will stash %zu canvases as *.%s in path: \'%s\'\n", cs.size(), ext, p.c_str());

	for(auto* c : cs) {
		auto name = std::string(c->GetName());
		auto outfile = p / (name + "." + ext);
		// Force rendering
		c->cd();
		c->Modified();
		c->Update();
		gSystem->ProcessEvents();
		if(extension == Extension::png) {
			auto img = std::unique_ptr<TImage>{TImage::Create()};
			if(!img) {
				WARN("canvas::save_all: TImage::Create() failed for Canvas: \'%s\' (title: \'%s\')\n",
					c->GetName(), c->GetTitle());
				continue;
			}
			img->FromPad(c);
			img->WriteImage(outfile.c_str());
		} else {
			c->Print(outfile.c_str());
		}
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

template<enum WhereAmI loc>
void save_all (
	std::vector<Extension> extension,
	std::vector<std::string_view> extra_tag = {}
) {
	for(auto const e : extension)
		save_all<loc>(e, extra_tag);
}
/* All vectors passed by value, but it's small so overhead is whatever,
 * and this fnc won't ever get called in a loop. */

}; // namespace canvas
extern template void canvas::save_all<canvas::Macro>(canvas::Extension , std::vector<std::string_view> );
extern template void canvas::save_all<canvas::Exe  >(canvas::Extension , std::vector<std::string_view> );
extern template void canvas::save_all<canvas::Macro>(std::vector<canvas::Extension> , std::vector<std::string_view> );
extern template void canvas::save_all<canvas::Exe  >(std::vector<canvas::Extension> , std::vector<std::string_view> );

inline std::ostream& operator<<(std::ostream& os, canvas::Extension e) {
	return os << magic_enum::enum_name(e);
}
inline std::ostream& operator<<(std::ostream& os, DoSave e) {
	return os << magic_enum::enum_name(e);
}

namespace mnd {

/* Nicer API to allow strong typing.
 * This is basically a zero-cost abstraction that allows really
 * pretty API's to directly name the positional arguments. */
template<typename T, typename Tag = void>
struct InputWrapper {
	using value_type = T;
	T value;

	InputWrapper() = default;
	InputWrapper(T v) : value(std::move(v)) {}

	operator T&() noexcept { return value; }
	operator const T&() const noexcept { return value; }

    T&       get() &       noexcept { return value; }
    T const& get() const & noexcept { return value; }
    T&&      get() &&      noexcept { return std::move(value); }
};

template<typename T, typename = void>
struct is_istreamable : std::false_type {};
template<typename T>
struct is_istreamable<T,
    std::void_t <
        decltype(std::declval<std::istream&>() >> std::declval<T&>())
    >
> : std::true_type {};

template<typename T, typename = void>
struct is_ostreamable : std::false_type {};
template<typename T>
struct is_ostreamable<T,
    std::void_t<
        decltype(std::declval<std::ostream&>() << std::declval<T const&>())
    >
> : std::true_type {};

/* Why is Rust amazing? Well, simple:
 * enum Option<T> {
 *   Some(T),
 *   None
 * }
 * A true algebraic sum type! Nullability isn't tied to
 * a self-defined `nil` subset within `T` itself. */

template<typename T>
class Option {
	struct NoTag {};

public:
	using value_type = T; // needed for CLI11

	struct Yes { T value; };

	static constexpr NoTag No{};

	Option(Yes y) : data(std::move(y)) {}
	Option(NoTag) : data(std::monostate{}) {};
	
	bool is_some() const noexcept { return data.index() == 0; }
	bool is_none() const noexcept { return data.index() == 1; }

	/* May panic (throw). Unlike rust, returns back a reference when called on lvalue. */
	T const& unwrap() const& { return std::get<0>(data).value; }
	T&       unwrap() &      { return std::get<0>(data).value; }
	T        unwrap() &&     { return std::get<0>(data).value; }

	decltype(auto) get() const noexcept { return (data); }
	decltype(auto) get() noexcept { return (data); }

protected:
	std::variant<Yes, std::monostate> data;
};

} // namespace mnd

/* Parse a file first thru the GCC preprocessor, and then
 * try to parse the output as a sequence of lines.
 * Is not thread safe! */
std::vector<std::string> ParseFile(const std::string& );
extern std::vector<std::string> ParseFile(const std::string& );

std::string ParseFileToString(const std::string& );
extern std::string ParseFileToString(const std::string& );

/* For the Option<T> wrapper, also expose a CLI tool template specialization
 * to parse it properly, otherwise boilerplate reeks through the code. */
template <
	typename T
> CLI::Option* add_logged_option (
	CLI::App& app,
	const std::string& name,
	mnd::Option<T>& variable,
	const std::string& description
) {
	auto state = std::make_shared<mnd::cli::detail::State>();
	auto* opt = app.add_option_function<T>(
			name,
			[&variable, name, state](const T& match) {
				if(state->current_is_authoritative) { // Respect my authoritah.
					variable = typename mnd::Option<T>::Yes{ .value = match };
					state->authoritative_seen = true;
					WARN("Parsed %sauthoritative%s sum-type option ", BOLD, KNRM);
					std::cerr << KBH_YEL << name << KNRM << " as "
						<< KBH_CYN << match << KNRM << '\n';
				} else if(!state->authoritative_seen) {
					variable = typename mnd::Option<T>::Yes{ .value = match };
					WARN("Parsed sum-type option ");
					std::cerr << KBH_YEL << name << KNRM << " as "
						<< KBH_CYN << match << KNRM << '\n';
				}
				state->current_is_authoritative = false;
			},
			description
		)
		->transform([state](std::string input) -> std::string {
			if(input == "@") return "{}";
			if(!input.empty() && input.front() == mnd::cli::detail::auth_sym) {
				state->current_is_authoritative = true;
				input.erase(input.begin());	
			}
			return input;
        }, "@ means an explicitly empty (but in the value variant)")
		->expected(0,-1)
		->trigger_on_parse()
		->default_str("none");

	return opt;
}

/* Custom char buffer streaming operations for the phantom wrapper types, if the underlying type
 * implements them. If underlying type's definitions are not found at this point, then this
 * template is sfinae'd out. E.g. vector|array overload is in `json_struct_def.hh`, and won't be
 * automatically detected here, if that header is included *after* this one.
 *
 * Non-templated specialized overloads can still be defined and compiler will like them more. Obviously. */
template<typename T, typename Tag,
    typename = std::enable_if_t<mnd::is_istreamable<T>::value>
> std::istream& operator>>(std::istream& in, mnd::InputWrapper<T, Tag>& value) {
    return in >> value.get();
}

template<typename T, typename Tag,
    typename = std::enable_if_t<mnd::is_ostreamable<T>::value>
> std::ostream& operator<<(std::ostream& out, mnd::InputWrapper<T, Tag> const& value) {
    return out << value.get();
}

namespace mnd {

/* Invoke a function `func` over a range of objects, over nthreads.
 * `Range` here binds here to any type anything that is indexable such as array/vector/span.
 * Function invocation can carry mutable (outside) state, and each thread
 * gets a copy of the functor. */
template<typename Range, typename F>
void parallel_process(
	Range& objects,
	size_t nthreads,
	F&& func
) {
	/* Idea of this call is the following - outside objects will get captured
	 * by value and each thread gets its own copy, invokes the functor over its given objects,
	 * and then merges its copy's objects back to the original one sitting in the main thread.
	 *
	 * How we achieve this, is that the objects that are explicitly mutated in the lambda must have
	 * structure similar to TH1P/TH2P, where each cloned object (via T(const T&) copy-ctor) 'remembers' its parent
	 * as a simple pointer, and during destruction gives its acquired contents back. */
	if(objects.empty() || nthreads == 0)
		return;

	using Fn  = std::decay_t<F>;
	using Ref = decltype(objects[size_t{}]);

	static_assert(
		std::is_copy_constructible_v<Fn>,
		"parallel_process requires a copyable callable"
	);

	static_assert(
		std::is_invocable_v<Fn&, size_t, Ref> ||
		std::is_invocable_v<Fn&, Ref>,
		"func must accept either (size_t, T&) or (T&)"
	);

	nthreads = std::min(nthreads, objects.size());

	/* Stable copies that live in this stack frame.
	 * In particular, these functors will NOT be owned by the
	 * std::thread/jthread callable objects. Inside, the threads just
	 * touch the underlying raw ptr. */
	std::vector<std::unique_ptr<Fn>> worker_funcs;
	worker_funcs.reserve(nthreads);

    /* Materialize here *exactly* one callable object first.
     * Each worker below gets a COPY of this object. This original
	 * functors is to be pinned to this stack frame and does not move. */
    worker_funcs.emplace_back(
		std::make_unique<Fn>( std::forward<F>(func) )
	);

	for(size_t i = 1; i < nthreads; ++i)
		worker_funcs.emplace_back(
			std::make_unique<Fn>( *worker_funcs.front() )
		);

	std::vector<jthread> threads;
	threads.reserve(nthreads);

	std::atomic<size_t> next{0};
	for(size_t t = 0; t < nthreads; ++t) {
		threads.emplace_back([&, t]() {
			auto& func = *worker_funcs[t];
			while(true) {
				const auto i = next.fetch_add(1, std::memory_order_relaxed);

				if(i >= objects.size())
					return;
				
				if constexpr(std::is_invocable_v<Fn&, size_t, Ref>) {
					std::invoke(func, i, objects[i]);
				} else {
					std::invoke(func, objects[i]);
				}
			}
		});
	}
	
	threads.clear();

	/* Worker functors that will be destroyed are serialized back to the original copy.
	 * Thus TH1P::~TH1P() merges worker histograms into
	 * *worker_funcs.front() one at a time. No data remains in the clones!
	 * Calling vector<T>::clear() does not guarantee sequential dtor calls from back to front
	 * element. */
	while(worker_funcs.size() > 1)
		worker_funcs.pop_back();
	
	// *worker_funcs[0] gets destructed here and returns its contents back to the main thread.
}

} // namespace mnd
