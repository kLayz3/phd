#pragma once

#include <iostream>
#include <istream>
#include <optional>
#include <ostream>
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
#include "cli/CLI11.hpp"

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

/* Make a sequence of paths based on an iterable container.
 * E.g.: path_sequence{"foo", "bar", "baz.txt"} -> "foo/bar/baz.txt" */
template<typename Range>
std::filesystem::path path_sequence(Range const& parts) {
	std::filesystem::path p;

	for(auto const& part : parts)
		p /= part;

	return p;
}

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
namespace type_traits {

template<typename T, typename S = std::istream>
using is_istreamable = CLI::detail::is_istreamable<T, S>;
template<typename T, typename S = std::ostream>
using is_ostreamable = CLI::detail::is_ostreamable<T, S>;

} // namespace type_traits

template<typename T>
std::string to_string(const T& val) {
	static_assert(type_traits::is_istreamable<T>::value,
		"Type T must be input-streamable to stringstream. AKA: there must be at least "
		"istream& operator>>(..) overload (works also for std::stringstream).");
	
	std::stringstream ss{};
	ss << val;
	return ss.str();
}

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

/* Rust fanboy? Well, simple:
 * enum Option<T> {
 *   Some(T),
 *   None
 * }
 * A true algebraic sum type! Nullability isn't tied to
 * a self-defined `nil` subset within `T` itself.
 *
 * NOTE: this type shouldn't be used as direct replacement of std::optional. */

template<typename T>
struct Some {
	T value;
};
template<typename T>
Some(T) -> Some<T>;

template<typename T>
class Option {
public:
	using value_type = T; // needed for CLI11

	constexpr Option()      : data(None) {};
	constexpr Option(std::nullopt_t) : data(None) {};
	
	template<typename U>
	constexpr Option(Some<U> some) : data( T{ std::move(some.value) } ) {}

	constexpr bool is_some() const noexcept { return data.has_value(); }
	constexpr bool is_none() const noexcept { return !is_some(); }

	/* May panic (throw). Unlike rust, returns back a reference when called on lvalue. */
	constexpr T const& unwrap() const& { return data.value(); }
	constexpr T&       unwrap() &      { return data.value(); }
	constexpr T&&      unwrap() &&     { return std::move(data.value()); }

	constexpr decltype(auto) get() const noexcept { return (data); }
	constexpr decltype(auto) get() noexcept { return (data); }

	/* Normally in STL, the functor type `F` is constrained by different concepts. */
	template<typename F>
	constexpr auto and_then(F&& f) & {
		if(is_some())
			return std::invoke(std::forward<F>(f), data.value());
		else
			return mnd::remove_cvref_t<std::invoke_result_t<F, T&>>{};
	}
	template<typename F>
	constexpr auto and_then(F&& f) const& {
		if(is_some())
			return std::invoke(std::forward<F>(f), data.value());
		else
			return mnd::remove_cvref_t<std::invoke_result_t<F, T const&>>{};
	}
	template<typename F>
	constexpr auto and_then(F&& f) && {
		if(is_some())
			return std::invoke(std::forward<F>(f), std::move(data.value()));
		else
			return mnd::remove_cvref_t<std::invoke_result_t<F, T>>{};
	}
	template<typename F>
	constexpr auto and_then(F&& f) const&& {
		if(is_some())
			return std::invoke(std::forward<F>(f), std::move(data.value()));
		else
			return mnd::remove_cvref_t<std::invoke_result_t<F, T const>>{};
	}
	template<typename F>
	constexpr Option or_else( F&& f ) const& {
		return *this ? *this : std::forward<F>(f)();
	};
	template<typename F>
	constexpr Option or_else( F&& f ) && {
		return *this ? std::move(*this) : std::forward<F>(f)();
	};

	/* Reset the state back to the `No` variant. */
	constexpr void reset() noexcept { data.reset(); }
	
	template<typename U>
	Option& operator=(U&& rhs) {
		data = std::forward<U>(rhs);
		return *this;
	}

protected:
	std::optional<T> data;
};

/* Predicate if the value is inside a range spanned by last 2 elements of some array.
 * Note, variant state is *assumed* to be valued here, and isn't checked! */
template<typename T, typename U, std::size_t N>
bool IsInside(const T& value, const Option<std::array<U,N>>& bounds) {
	static_assert(N >= 2, "Array size must be >= 2");
	static_assert (
		std::is_convertible_v<decltype(std::declval<const U&>() <= std::declval<const T&>()), bool> &&
		std::is_convertible_v<decltype(std::declval<const T&>() <  std::declval<const U&>()), bool>,
		"T and U must support comparison operators U <= T and T < U."
	);
	const auto& bval = bounds.unwrap();
	return bval[N-2] <= value and value < bval[N-1];
}

/* This option is nullable by value type not by NAN boundary. */
template<typename T, std::size_t N>
bool IsValid(const Option<std::array<T,N>>& bounds) {
	static_assert(N >= 2, "Array size must be >= 2");
	return bounds.is_some();
}

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
namespace mnd::cli::detail {

inline constexpr char empty_sym = '@';
inline constexpr char reset_sym = '~';
}

template <
	typename T
> CLI::Option* add_logged_option (
	CLI::App& app,
	const std::string& name,
	mnd::Option<T>& variable,
	const std::string& description
) {
	static_assert(
		std::is_default_constructible_v<T>,
		"mnd::Option<T> CLI parsing requires default-constructible T"
	);
	/* God I love undocumented API. So basically the tokenising begins *before*
	 * transformers take place. E.g. `--flag=!@` cannot be parsed for Option<array> ...
	 * Just had to dissect the library like usual. HINT: people please write your docs. */
	 
	/* Inside CLI11.hpp:
	 *   using results_t = std::vector<std::string>;
	 *   using callback_t = std::function<bool(const results_t &)>;
	 */

	auto state = std::make_shared<mnd::cli::detail::State>();
	
	CLI::callback_t callback =
	[&variable, name, state](const CLI::results_t& raw) -> bool {
		if(raw.empty())
			return false;

		/* Copy because we're going to strip the authoritative prefix
		 * from the first vector's element before delegating to CLI11. */
		auto input = raw;
		auto& first = input.front();
	
		/* Problem is that for containers, passing a single flag e.g. '~' will result 
		 * in its results vector being padded by empty strings... */
		auto is_magic = [&input, &first](char symbol) {
			return first.size() == 1 &&
			first.front() == symbol &&
			std::all_of(
				std::next(input.begin()),
				input.end(),
				[](const std::string& s) { return s.empty(); }
			);
		};

		const bool authoritative =
			!first.empty() &&
			first.front() == mnd::cli::detail::auth_sym;

		if(authoritative)
			first.erase(first.begin());

		/* A lone '!' isn't a valid value. */
		if(first.empty())
			return false;

		const bool reset = is_magic(mnd::cli::detail::reset_sym);
		const bool empty = is_magic(mnd::cli::detail::empty_sym);

		T match{}; // local placeholder for value

		/* '@' and '~' both yield T{},
		 * otherwise use CLI11's normal conversion for T. */
		if(!reset && !empty) {
			if(!CLI::detail::lexical_conversion<T,T>(input, match))
				return false;
			/* On successful conversion, keep going. */
		}

		/* '!x' always wins
		 * ordinary x is ignored after any !x
		 * a later !x can replace an earlier !x */
		if(!authoritative && state->authoritative_seen)
			return true;

		if(authoritative)
			state->authoritative_seen = true;

		variable = std::move(match);

		if(reset) {
			variable.reset();

			WARN("Parsed %ssum-type reset%s option ",
				authoritative ?  BOLD "authoritative ": "",
				authoritative ? KNRM : ""
			);

			std::cerr << KBH_YEL << name << KNRM << " as "
				<< MND_RGB_COL(204,102,0) << "none" << KNRM << '\n';

			return true;
		}

		if(authoritative) {
			WARN("Parsed %sauthoritative%s sum-type option ", BOLD, KNRM);
		} else {
			WARN("Parsed sum-type option ");
		}

		std::cerr << KBH_YEL << name << KNRM << " as "
			<< (empty ? MND_RGB_COL(255,102,255) "defaulted value: " : "")
			<< KBH_CYN << variable.unwrap() << KNRM << '\n';

		return true;
	};
	
	return app.add_option(
		name,
		std::move(callback),
		description
		 + mnd::msg("\n%s\'%c\'%s flag requests an explicitly empty object (but in the value-given variant), "
			        "\n%s\'%c\'%s requests a reset back to the no-value-given variant.",
		            MND_RGB_COL(255,102,255), mnd::cli::detail::empty_sym, KNRM,
		            MND_RGB_COL(204,102,0  ), mnd::cli::detail::reset_sym, KNRM)
	)
	->type_name(CLI::detail::type_name<T>())
	->type_size(
		1, CLI::detail::type_count<T>::value
	)
	->expected(
		CLI::detail::expected_count<T>::value
	)
	->trigger_on_parse()
	->default_str("none");
}

/* Custom char buffer streaming operations for the phantom wrapper types, if the underlying type
 * implements them. If underlying type's definitions are not found at this point, then this
 * template is sfinae'd out. E.g. vector|array overload is in `json_struct_def.hh`, and won't be
 * automatically detected here, if that header is included *after* this one.
 *
 * Non-templated specialized overloads can still be defined and compiler will like them more. Obviously. */
template<typename T, typename Tag,
    typename = std::enable_if_t<mnd::type_traits::is_istreamable<T>::value>
> std::istream& operator>>(std::istream& in, mnd::InputWrapper<T, Tag>& value) {
    return in >> value.get();
}

template<typename T, typename Tag,
    typename = std::enable_if_t<mnd::type_traits::is_ostreamable<T>::value>
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
	/* Idea of this call is the following: outside objects will get captured
	 * by value and each thread gets its own copy, invokes the functor over the state of its given objects,
	 * and then merges its copy's mutated objects back to the original one which is sitting in the main thread.
	 *
	 * How we achieve this, is that the objects that are explicitly mutated in the lambda must have
	 * structure similar to TH1P/TH2P, where each cloned object (via T(const T&) copy-ctor) 'remembers' its parent
	 * as a simple pointer, and during destruction gives its acquired contents back to the parent. */
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
	 * Thus THXP::~THXP() merges worker histograms into
	 * *worker_funcs.front() one at a time. No data remains in the clones!
	 * Calling vector<T>::clear() does not guarantee sequential dtor calls from back to front
	 * element. */
	while(worker_funcs.size() > 1)
		worker_funcs.pop_back();
	
	// *worker_funcs[0] gets destructed here and returns its contents back to the main thread.
}

} // namespace mnd
