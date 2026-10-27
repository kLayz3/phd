#include "PythonCLI.h"

#include "util/PrettyHisto.h"
#include "util/TextUtils.h"
#include "util/FromChars.h"

#include <iomanip>
#include <iterator>
#include <optional>
#include <stdexcept>

#include "monad/accessories/formatting.hxx"

using namespace mnd::plot;

namespace {

using View = std::string_view;

// Plot grammar rejects empty fields and ignores surrounding whitespace.
std::vector<std::string> split_fields(View text) {
	auto parts = mnd::split(std::string{text}, ',');
	for(const auto& part : parts)
		if(part.empty()) MND_THROW("empty item (or trailing equals '=' symbol? )");

	/* Parts now must have either:
	 * [1]: key=value syntax, for either name= or figure=
	 * [2]: key:value syntax, for specific underblocks. */

	return parts;
}

} // namespace


// Uniform adapters handle the different arities and argument types.
void PyInfo::add_figure_cfg(View kv_cfg) {
	const auto colon = kv_cfg.find(':');
	if(colon == View::npos)
		MND_THROW("parse_figure_cfg(): the string '%.*s' isn't of the form key:value ?",
			(int)kv_cfg.size(), kv_cfg.data());

	const View key   = mnd::trim(kv_cfg.substr(0, colon));
	const View value = kv_cfg.substr(colon + 1);

	if(key == "xlabel") {
		fig_style.xlabel_ = std::string{value};
	} else if(key == "ylabel") {
		fig_style.ylabel_ = std::string{value};
	} else if(key == "title") {
		fig_style.title_ = std::string{value};
	} else if(key == "save_dpi") {
		double v = mnd::stof<double>(value).unwrap(mnd::msg(
			"save_dpi flag value: '%.*s' not parsable to double", (int)value.size(), value.data())
		);
		fig_style.save_dpi_ = Some{v};
	} else if(key == "figsize") {
		auto fs = mnd::split_view(value, ',');
		if(fs.size() != 2)
			MND_THROW("Found %.*s:%.*s  pair, but rhs cannot be destructured into integer,integer .",
				(int)key.size(),key.data(), (int)value.size(), value.data());

		u32 d0 = mnd::stou(fs[0]).unwrap( mnd::msg(
			"parse_figure_cfg(): '%.*s:%.*s' parsing cannot parse 1st arg: '%.*s' as integer.",
			 (int)key.size(),key.data(), (int)value.size(), value.data(), (int)fs[0].size(), fs[0].data())
		);
		u32 d1 = mnd::stou(fs[1]).unwrap( mnd::msg(
			"parse_figure_cfg(): '%.*s:%.*s' parsing cannot parse 2nd arg: '%.*s' as integer.",
			(int)key.size(),key.data(), (int)value.size(), value.data(), (int)fs[1].size(), fs[1].data())
		);

		fig_style.figsize_ = Some{ std::pair<u32,u32>{d0,d1} };
	} else if(key == "xlim" || key == "ylim") {
		auto fs = mnd::split_view(value, ',');
		if(fs.size() != 2)
			MND_THROW("Found %.*s:%.*s  pair, but rhs cannot be destructured into double,double .",
				(int)key.size(),key.data(), (int)value.size(), value.data()
			);

		double d0 = mnd::stof<double>(fs[0]).unwrap( mnd::msg(
			"parse_figure_cfg(): '%.*s:%.*s' parsing cannot parse 1st arg: '%.*s' as double.",
				(int)key.size(),key.data(), (int)value.size(), value.data(), (int)fs[0].size(), fs[0].data())
			);
		double d1 = mnd::stof<double>(fs[1]).unwrap( mnd::msg(
			"parse_figure_cfg(): '%.*s:%.*s' parsing cannot parse 2nd arg: '%.*s' as double.",
				(int)key.size(),key.data(), (int)value.size(), value.data(), (int)fs[1].size(), fs[1].data())
			);

		if(key == "xlim")
			fig_style.xlim_ = Some{ std::pair<double,double>{d0,d1} };
		else 
			fig_style.ylim_ = Some{ std::pair<double,double>{d0,d1} };
	}
	else if(key == "logx" || key == "logy" || key == "logz"
		|| key == "enable_right_top_spline" || key == "legend" || key == "grid") {
		bool b = false;
		if(!mnd::parse(value, b))
			MND_THROW("parse_figure_cfg(): '%.*s:%.*s', rhs not parsable to a boolean.",
				(int)key.size(),key.data(), (int)value.size(), value.data());

		if(key == "logx")      fig_style.logx_ = b;
		else if(key == "logy") fig_style.logy_ = b;
		else if(key == "logz") fig_style.logz_ = b;
		else if(key == "enable_right_top_spline") fig_style.hide_right_and_top_spline_ = !b;
		else if(key == "legend") fig_style.legend_ = b;
		else if(key == "grid") fig_style.grid_ = b;
		else MND_THROW("should have been unreachable?");
	}
	else {
		MND_THROW("unknown figure key? Parsed 'key:value' as ...  '%.*s':'%.*s'",
			(int)key.size(),key.data(), (int)value.size(), value.data());
	}
}

void PyInfo::ObjectInfo::add_object_cfg(View kv_cfg) {
	const auto colon = kv_cfg.find(':');
	if(colon == View::npos)
		MND_THROW("add_object_cfg(): the string '%.*s' isn't of the form key:value ?",
			(int)kv_cfg.size(), kv_cfg.data());

	const View key   = mnd::trim(kv_cfg.substr(0, colon));
	const View value = kv_cfg.substr(colon + 1);
	
	options.emplace_back(key, value);
}

PyInfo PyInfo::parse(View text) {
	PyInfo info;

	auto comma_sep_fields = split_fields(text);
	bool is_in_figure_parsing = false;
	
	for(View subtext: comma_sep_fields) {
		/* Subtext has either:
		 * [1]: key=value syntax, for either name= or figure=
		 * [2]: key:value syntax, for specific underblocks. */
		
		if( has_prefix(subtext, "name=") ) {
			subtext.remove_prefix(5);
			
			if(subtext.empty()) MND_THROW("name= empty rhs not allowed.");
			info.objects.emplace_back();
			info.objects.back().name = subtext;

			is_in_figure_parsing = false;
		}
		else if( has_prefix(subtext, "figure=") ) {
			subtext.remove_prefix(7);
			if(subtext.empty()) MND_THROW("figure= empty rhs not allowed.");

			is_in_figure_parsing = true;

			info.add_figure_cfg(subtext);
		}
		else if(is_in_figure_parsing) {
			info.add_figure_cfg(subtext);
		}
		else {
			if(info.objects.size() == 0)
				MND_THROW("didn't pass name=value initially to specify an entry object.");

			info.objects.back().add_object_cfg(subtext);
		}
	}
	return info;
}

std::istream& mnd::plot::operator>>(std::istream& in, PyInfo& info) {
	std::istream::sentry sentry{in};
	if(!sentry) return in;

	const std::string text{
		std::istreambuf_iterator<char>{in}, std::istreambuf_iterator<char>{}
	};
	try {
		auto parsed = PyInfo::parse(text);
		info = std::move(parsed);
	} catch(const std::runtime_error&) {
		in.setstate(std::ios::failbit);
	} catch(const std::bad_optional_access&) {
		in.setstate(std::ios::failbit);
	}
	return in;
}

std::ostream& mnd::plot::operator<<(std::ostream& out, const PyInfo::ObjectInfo& o) {
	using mnd::fmt::operator<<;
	return out << "name: " << o.name << "; opts: " << o.options;
}

std::ostream& mnd::plot::operator<<(std::ostream& out, const PyInfo& info) {
	using mnd::fmt::operator<<;

	out << "PyInfo{ objects: ";
	out << info.objects;

	const auto& f = info.fig_style;
	out << ", figure: {xlabel:" << std::quoted(f.xlabel_)
		<< ", ylabel:" << std::quoted(f.ylabel_)
		<< ", title:" << std::quoted(f.title_);

	const auto pair = [&](const char* key, const auto& value) {
		if(value.is_some())
			out << ", " << key << ":(" << value->first << ',' << value->second << ')';
	};
	pair("figsize", f.figsize_);
	pair("xlim", f.xlim_);
	pair("ylim", f.ylim_);
	if(f.save_dpi_.is_some()) out << ", save_dpi:" << f.save_dpi_.unwrap();

	out << ", logx:" << (f.logx_ ? "true" : "false")
		<< ", logy:" << (f.logy_ ? "true" : "false")
		<< ", logz:" << (f.logz_ ? "true" : "false")
		<< ", grid:" << (f.grid_ ? "true" : "false")
		<< ", legend:" << (f.legend_ ? "true" : "false")
		<< ", enable_right_top_spline:"
		<< (f.hide_right_and_top_spline_ ? "false" : "true")
		<< "}}";
	return out;
}

const char* helper_ = R"(
This small DSL serves as a placeholder to manually select plotting of different histograms from ROOT also into Python.
Basically API is something like: the users will pass in:

--python=name=histname

Will default-plot the histogram with the API from PrettyHisto.h|cxx

If histname corresponds to something dynamically castable to, say, TH2, then users can also select options such as:

--python=\
   name=histname,\
     mode:contour,\
     colorbar:false,\
     label:"Foo Bar TH2",\
   name=graphname,\
     mode:scatterplot,\
     line_style:'--',\
     line_width:2,\
     label:'Foo Bar Graph',\
   figure=\
     xlabel:"X label: FooBar",\
     ylabel:"Y label: FooBar",\
     title:"Does it work?",\
     figsize:'10,12',\
     save_dpi:200,\
     logy:true

Where methods (same name) from Style classes come as key:pair after name=histname tag,
and methods for Figure come after the figure= tag.
Key:value items not given are selected as defaulted, same way as the PrettyHisto API does it.

Whats also acceptable is something like:
--python=\
   name=histname,\
     mode:contour,\
     colorbar:false,\
     label:"Foo Bar TH2",\
   name=graphname,\
   figure=\
     xlabel:"X label FooBar",\
     ylabel:"Y label FooBar",\
     title:"Does it work?",\
     save_dpi:200

Which means that the TH2 hist has some user-passed style selection, and graph just inherits
the ROOT's style. Figure also inherits the non-given fields from default. E.g., this also works:
--python=\
   name=histname,\
     mode:contour,\
     colorbar:false,\
     label:"Foo Bar TH2",\
   name=graphname

--python=\
   name=histname,\
   name=graphname

Plots both the hist and the graph onto the same canvas.

)";
