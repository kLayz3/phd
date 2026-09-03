/* pybind11 stuff must be first to be included. */
#include <pybind11/numpy.h>
#include <pybind11/embed.h>
#include <pybind11/pytypes.h>
#include <stdexcept>
#include <csignal>
#include <string>

#include "PrettyHisto.h"
#include "RtypesCore.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TROOT.h"

/* First histogramming API thta we left undefined.. */

#define APPEND_TO_TITLE_IMPL(THXP) \
void THXP::AppendToTitle(std::string_view v) { \
	auto current_title = std::string{ h.GetTitle() }; \
	if(current_title.empty()) { \
		h.SetTitle( std::string(v).c_str() ); \
		return; \
	} \
	/* Current title can end with a comment: \
	 * (text) \
	 * Means strip out the last bracket, paste in the text and add \
	 * back the bracket */ \
	\
	switch(current_title.back()) { \
		case ')': { \
			current_title.pop_back(); \
			current_title += std::string{' '} + std::string{v} + std::string{')'}; \
			break; \
		} \
		default: { \
			current_title += std::string{' '} + std::string{v}; \
		} \
	} \
	h.SetTitle( current_title.c_str() ); \
}

APPEND_TO_TITLE_IMPL(TH1P)
APPEND_TO_TITLE_IMPL(TH2P)

namespace py = pybind11;
using namespace pybind11::literals;

namespace {

/* explicit scoped_interpreter(
 *   bool init_signal_handlers = true,
 *   int argc = 0,
 *   const char* const* argv = nullptr,
 *   bool add_program_dir_to_path = true
 * );
 * But even with this API, somehow signal handler gets remapped to a custom handler within
 * CPython... This was a fucking nightmare to debug.
 * Usually we're a normal C++/ROOT process and Python isn't
 * initialized yet. But don't try to initialize it twice if
 * somebody else already embedded CPython. */
struct PythonRuntime {
	std::optional<py::scoped_interpreter> owned;

	PythonRuntime() {
		if(!Py_IsInitialized())
			owned.emplace(false);

		py::gil_scoped_acquire gil;

		auto matplotlib = py::module_::import("matplotlib");

		/* No GUI required. Must happen before importing pyplot. */
		matplotlib.attr("use")("Agg");
		
		py::module_::import("matplotlib.pyplot");

		PyOS_setsig(SIGINT, SIG_DFL); // hehehehehehe losing my mind
	}
};

void ensure_python() {
	static PythonRuntime runtime;
}

} // namespace {anonymous}

void mnd::python::poke(bool verbose) {
	ensure_python();

	py::gil_scoped_acquire gil;

	try {
		auto sys = py::module_::import("sys");
		auto np  = py::module_::import("numpy");
		auto mpl = py::module_::import("matplotlib");
		mpl.attr("use")("Agg");
		auto plt = py::module_::import("matplotlib.pyplot");

		PyOS_setsig(SIGINT, SIG_DFL);

		// Force a minimal matplotlib object
		auto result = plt.attr("subplots")().cast<py::tuple>();
		py::object fig = result[0];

		plt.attr("close")(fig);

		if(verbose) {
			std::cout
				<< "[PYTHON]: OK, executable='"
				<< py::str(sys.attr("executable")).cast<std::string>()
				<< "', numpy='"
				<< py::str(np.attr("__version__")).cast<std::string>()
				<< "', matplotlib='"
				<< py::str(mpl.attr("__version__")).cast<std::string>()
				<< "'\n";
		}
	} catch(py::error_already_set const& e) {
		throw std::runtime_error(
			std::string{"Python/matplotlib environment not usable:\n"} + e.what()
		);
	}
}

mnd::col::RGBA::RGBA(Color_t root_color) {
	auto* c = gROOT->GetColor(root_color);
	if(!c) return;

	r = c->GetRed();
	g = c->GetGreen();
	b = c->GetBlue();
	a = c->GetAlpha();
}

mnd::col::RGBA mnd::col::RGBA::from_packed(uint32_t value) noexcept {
	const double b = (value & 0xffu) / 255.0; value >>= 8;
	const double g = (value & 0xffu) / 255.0; value >>= 8;
	const double r = (value & 0xffu) / 255.0; value >>= 8;
	const double a = std::max(1.0 - (value & 0xffu) / 255.0, 0.0);
	return {r, g, b, a};
}

Int_t mnd::col::RGBA::GetColorCode() const {
	return TColor::GetColor(
		static_cast<Float_t>(r),
		static_cast<Float_t>(g),
		static_cast<Float_t>(b)
	);
}

void mnd::col::RGBA::ApplyFill(TH1* h) const {
	Int_t idx = this->GetColorCode();
	
	/* Sometimes setting alpha doesn't work. The transparency is implicitly enabled through
	 * $ROOTSYS/etc/system.rootrc
	 * OpenGL.CanvasPreferGL = 1
	 * , if it is 0, then alpha will always be total (1.0). */
	h->SetFillColorAlpha(idx, a);
}

mnd::col::RGBA mnd::col::literals::operator""_c(unsigned long long int value) {
	if(value > std::numeric_limits<uint32_t>::max())
		throw std::out_of_range("packed color literal does not fit in 32 bits");
	return mnd::col::RGBA::from_packed(static_cast<uint32_t>(value));
}

using namespace mnd::plot;

static mnd::col::RGBA root_color(Color_t id) {
	auto* c = gROOT->GetColor(id);
	if(!c) return {0., 0., 0., 1.};

	return {
		c->GetRed(),
		c->GetGreen(),
		c->GetBlue(),
		c->GetAlpha()
	};
}

static std::string root_line_style(Style_t s) {
	switch(s) {
		case 1: return "-";
		case 2: return "--";
		case 3: return ":";
		case 4: return "-.";
		default: return "-";
	}
}

static std::string root_marker_style(Style_t s) {
	s = TAttMarker::GetMarkerStyleBase(s);

	switch(s) {
		/* Basic markers. */
		case 1:       return ".";   // dot
		case 2:       return "+";   // plus
		case 3:       return "*";   // star / asterisk
		case 4:       return "o";   // open circle
		case 5:       return "x";   // multiply
		case 6 ... 19:return ".";   // fixed/scalable dots

		/* Filled basic shapes. */
		case 20:      return "o";   // full circle
		case 21:      return "s";   // full square
		case 22:      return "^";   // full triangle up
		case 23:      return "v";   // full triangle down

		/* Open basic shapes.
		 *
		 * Shape is the same in matplotlib; whether it is open/filled
		 * should be handled separately using markerfacecolor.
		 */
		case 24:      return "o";   // open circle
		case 25:      return "s";   // open square
		case 26:      return "^";   // open triangle up
		case 27:      return "D";   // open diamond
		case 28:      return "+";   // open cross
		case 29:      return "*";   // full star
		case 30:      return "*";   // open star
		case 31:      return "*";   // asterisk
		case 32:      return "v";   // open triangle down
		case 33:      return "D";   // full diamond
		case 34:      return "P";   // full cross / filled plus

		/*
		 * ROOT's more exotic compound markers have no exact
		 * matplotlib shorthand equivalent. Pick the closest
		 * recognizable geometry.
		 */
		case 35:      return "X";   // open diamond cross
		case 36:      return "s";   // open diagonal square
		case 37:      return "^";   // open three triangles
		case 38:      return "8";   // octagon + cross
		case 39:      return "^";   // full three triangles
		case 40:      return "X";   // open four triangles X
		case 41:      return "X";   // full four triangles X
		case 42:      return "D";   // open double diamond
		case 43:      return "D";   // full double diamond
		case 44:      return "P";   // open four triangles +
		case 45:      return "P";   // full four triangles +
		case 46:      return "x";   // open cross X
		case 47:      return "X";   // full cross X
		case 48:      return "X";   // four squares X
		case 49:      return "P";   // four squares +

		default:      return ".";
	}
}

[[ maybe_unused ]]
static bool root_marker_is_open(Style_t s) {
	s = TAttMarker::GetMarkerStyleBase(s);

	switch(s) {
		case 4:
		case 24 ... 28:
		case 30:
		case 32:
		case 35 ... 38:
		case 40:
		case 42:
		case 44:
		case 46:
			return true;

		default:
			return false;
	}
}

static GraphStyle resolve_graph_style(
	const TGraph& g,
	GraphStyle style
) {
	style.label_ =
		style.label_.value_or(g.GetTitle());

	style.line_width_ =
		style.line_width_.value_or(g.GetLineWidth());

	style.line_style_ =
		style.line_style_.value_or(
			root_line_style(g.GetLineStyle())
		);

	style.marker_ =
		style.marker_.value_or(
			root_marker_style(g.GetMarkerStyle())
		);
	
	style.marker_size_ =
		style.marker_size_.value_or(g.GetMarkerSize());

	style.color_ =
		style.color_.value_or(
			root_color(g.GetLineColor())
		);
	
	return style;
}

static auto rgba(const RGBA& c) -> py::tuple {
    return py::make_tuple(c.r, c.g, c.b, c.a);
}
static auto numpy_view(const std::vector<double>& v) -> py::array_t<double> {
    return py::array_t<double>(
        static_cast<py::ssize_t>(v.size()),
        v.data()
    );
}

[[ maybe_unused ]]
static auto numpy_view_2d(const Hist2DData& h) -> py::array_t<double> {
	return py::array_t<double>(
		{
			static_cast<py::ssize_t>(h.ny),
			static_cast<py::ssize_t>(h.nx)
		},
		{
			static_cast<py::ssize_t>(h.nx * sizeof(double)),
			static_cast<py::ssize_t>(sizeof(double))
		},
		h.values.data()
	);
}
auto bin_centers(const Hist1DData& h) -> std::vector<double> {
	std::vector<double> centers;
	centers.reserve(h.values.size());

	for(std::size_t i = 0; i < h.values.size(); ++i)
		centers.push_back(
			0.5 * (h.edges[i] + h.edges[i + 1])
		);

	return centers;
}

/* =============================================================== */

auto HistStyle::label(std::string s) && -> HistStyle {
	label_ = std::move(s);
	return std::move(*this);
}
auto HistStyle::line_width(double w) && -> HistStyle {
	line_width_ = w;
	return std::move(*this);
}
auto HistStyle::line_style(std::string s) && -> HistStyle {
	line_style_ = std::move(s);
	return std::move(*this);
}
auto HistStyle::marker(std::string s) && -> HistStyle {
	marker_ = std::move(s);
	return std::move(*this);
}
auto HistStyle::stairs() && -> HistStyle {
	mode_ = HistMode::stairs;
	return std::move(*this);
}
auto HistStyle::errors() && -> HistStyle {
	mode_ = HistMode::errors;
	return std::move(*this);
}
auto HistStyle::points() && -> HistStyle {
	mode_ = HistMode::points;
	return std::move(*this);
}
auto HistStyle::fill(bool b) && -> HistStyle {
	do_fill = b;
	return std::move(*this);
}
auto HistStyle::facecolor(mnd::col::RGBA col) && -> HistStyle {
	fill_color_ = std::move(col);
	return std::move(*this);
}
auto HistStyle::edgecolor(mnd::col::RGBA col) && -> HistStyle {
	line_color_ = std::move(col);
	return std::move(*this);
}

auto Hist2DStyle::label(std::string s) && -> Hist2DStyle {
	label_ = std::move(s);
	return std::move(*this);
}
auto Hist2DStyle::colorbar(bool enabled) && -> Hist2DStyle {
	colorbar_ = enabled;
	return std::move(*this);
}
auto Hist2DStyle::mesh() && -> Hist2DStyle {
	mode_ = Hist2DMode::mesh;
	return std::move(*this);
}
auto Hist2DStyle::image() && -> Hist2DStyle {
	mode_ = Hist2DMode::image;
	return std::move(*this);
}
auto Hist2DStyle::contour() && -> Hist2DStyle {
	mode_ = Hist2DMode::contour;
	return std::move(*this);
}

/* =============================================================== */

auto GraphStyle::label(std::string s) && -> GraphStyle {
	label_ = std::move(s);
	return std::move(*this);
}
auto GraphStyle::line_width(double w) && -> GraphStyle {
	line_width_ = w;
	return std::move(*this);
}
auto GraphStyle::line_style(std::string s) && -> GraphStyle {
	line_style_ = std::move(s);
	return std::move(*this);
}
auto GraphStyle::marker(std::string s) && -> GraphStyle {
	marker_ = std::move(s);
	return std::move(*this);
}
auto GraphStyle::marker_size(double w) && -> GraphStyle {
	marker_size_ = w;
	return std::move(*this);
}
auto GraphStyle::plot() && -> GraphStyle {
	mode_ = GraphMode::plot;
	return std::move(*this);
}
auto GraphStyle::scatter() && -> GraphStyle {
	mode_ = GraphMode::scatter;
	return std::move(*this);
}
auto GraphStyle::scatterplot() && -> GraphStyle {
	mode_ = GraphMode::scatterplot;
	return std::move(*this);
}

/* =============================================================== */

auto Figure::plot(const TH1& h, HistStyle style) && -> Figure {
	Hist1DData out;

	const int n = h.GetNbinsX();

	out.edges.reserve(n + 1);
	out.values.reserve(n);
	out.errors.reserve(n);

	for(int i = 1; i <= n; ++i) {
		out.edges.push_back(
			h.GetXaxis()->GetBinLowEdge(i)
		);

		out.values.push_back(
			h.GetBinContent(i)
		);

		out.errors.push_back(
			h.GetBinError(i)
		);
	}

	out.edges.push_back(
		h.GetXaxis()->GetBinUpEdge(n)
	);

	out.style.label_ =
		style.label_.value_or(h.GetTitle());

	out.style.line_width_ =
		style.line_width_.value_or(h.GetLineWidth());

	out.style.line_style_ =
		style.line_style_.value_or(
			root_line_style(h.GetLineStyle())
		);

	out.style.marker_ =
		style.marker_.value_or(
			root_marker_style(h.GetMarkerStyle())
		);

	out.style.line_color_ =
		style.line_color_.value_or(
			root_color(h.GetLineColor())
		);
	
	/* Passed-in style color takes priority */
	if(style.do_fill) {
		out.style.do_fill = style.do_fill;
		out.style.fill_color_ = style.fill_color_; // can be nullopt

		const bool has_fill = (h.GetFillStyle() != 0 && h.GetFillColor() != 0);

		/* Fallback case, in case facecolor attribute not supplied to the style. */
		if(has_fill and !out.style.fill_color_.has_value())
			out.style.fill_color_ = root_color(h.GetFillColor());
		
		/* In this case, just reject the fill request. */
		if(!out.style.fill_color_.has_value()) {
			fprintf(stderr,
				"mnd::plot::Figure::plot(): requested `fill()` for (ROOT) hist: \'%s\' "
				" but facecolor not provided and also histogram isn't filled.\n", h.GetTitle());
			out.style.do_fill = false;
		}
	}

	out.style.mode_ = style.mode_;

	objects_.emplace_back( std::move(out) );
	return std::move(*this);
}

auto Figure::plot(const TH2& h, Hist2DStyle style) && -> Figure {
	Hist2DData out;

	const int nx = h.GetNbinsX();
	const int ny = h.GetNbinsY();

	out.nx = nx;
	out.ny = ny;

	out.x_edges.reserve(nx + 1);
	out.y_edges.reserve(ny + 1);
	out.values.reserve(nx * ny);

	for(int ix = 1; ix <= nx; ++ix)
		out.x_edges.push_back(
			h.GetXaxis()->GetBinLowEdge(ix)
		);

	out.x_edges.push_back(
		h.GetXaxis()->GetBinUpEdge(nx)
	);

	for(int iy = 1; iy <= ny; ++iy)
		out.y_edges.push_back(
			h.GetYaxis()->GetBinLowEdge(iy)
		);

	out.y_edges.push_back(
		h.GetYaxis()->GetBinUpEdge(ny)
	);

	for(int iy = 1; iy <= ny; ++iy)
		for(int ix = 1; ix <= nx; ++ix)
			out.values.push_back(
				h.GetBinContent(ix, iy)
			);

	style.label_ =
		style.label_.value_or(h.GetTitle());

	out.style = std::move(style);

	objects_.emplace_back(std::move(out));
	return std::move(*this);
}

auto Figure::plot(const TGraph& g, GraphStyle style) && -> Figure {
	GraphData out;

	const int n = g.GetN();

	out.x.assign(g.GetX(), g.GetX() + n);
	out.y.assign(g.GetY(), g.GetY() + n);

	out.style = resolve_graph_style(g, style);

	objects_.emplace_back(std::move(out));
	return std::move(*this);
}

auto Figure::plot (
	const TGraphErrors& g,
	GraphStyle style
) && -> Figure {
	GraphData out;

	const int n = g.GetN();

	out.x.assign(g.GetX(), g.GetX() + n);
	out.y.assign(g.GetY(), g.GetY() + n);

	out.xerr.assign(g.GetEX(), g.GetEX() + n);
	out.yerr.assign(g.GetEY(), g.GetEY() + n);

	out.style = resolve_graph_style(g, style);

	objects_.emplace_back(std::move(out));
	return std::move(*this);
}

auto Figure::xlabel(std::string s) && -> Figure {
	xlabel_ = std::move(s);
	return std::move(*this);
}
auto Figure::xlabel(const TObject* obj) && -> Figure {
	if(auto const* g = dynamic_cast<TGraph const*>(obj)) {
		xlabel_ = g->GetXaxis()->GetTitle();
	}
	else if(auto const* g = dynamic_cast<TH1 const*>(obj)) {
		xlabel_ = g->GetXaxis()->GetTitle();
	}
	else {
		fprintf(stderr, "Figure::xlabel(const TObject*): passed [%s] obj, but expects something castable to TGraph/TH1\n",
			obj->ClassName());
	}
	return std::move(*this);
}
auto Figure::ylabel(std::string s) && -> Figure {
	ylabel_ = std::move(s);
	return std::move(*this);
}
auto Figure::ylabel(const TObject* obj) && -> Figure {
	if(auto const* g = dynamic_cast<TGraph const*>(obj)) {
		ylabel_ = g->GetYaxis()->GetTitle();
	}
	else if(auto const* g = dynamic_cast<TH1 const*>(obj)) {
		ylabel_ = g->GetYaxis()->GetTitle();
	}
	else {
		fprintf(stderr, "Figure::ylabel(const TObject*): passed [%s] obj, but expects something castable to TGraph/TH1\n",
			obj->ClassName());
	}
	return std::move(*this);
}

auto Figure::title(std::string s) && -> Figure {
	title_ = std::move(s);
	return std::move(*this);
}
auto Figure::title(const TObject* obj) && -> Figure {
	if(obj)
		title_ = obj->GetTitle();
	return std::move(*this);
}
auto Figure::legend(bool enabled) && -> Figure {
	legend_ = enabled;
	return std::move(*this);
}
auto Figure::logy(bool enabled) && -> Figure {
	logy_ = enabled;
	return std::move(*this);
}
auto Figure::xlim(double lo, double hi) && -> Figure {
	xlim_ = {lo, hi};
	return std::move(*this);
}
auto Figure::ylim(double lo, double hi) && -> Figure {
	ylim_ = {lo, hi};
	return std::move(*this);
}
auto Figure::logx(bool v) && -> Figure {
	logx_ = v;
	return std::move(*this);
}
auto Figure::grid(bool v) && -> Figure {
	grid_ = v;
	return std::move(*this);
}
auto Figure::enable_right_top_spline(bool v) && -> Figure {
	hide_right_and_top_spline_ = !v;
	return std::move(*this);
}

/* Rendering of the Python objects. */
static auto make_stairs_kwargs(const HistStyle& style) -> py::dict {
	py::dict kwargs;

	if(style.label_ && !style.label_->empty())
		kwargs["label"] = *style.label_;

	if(style.line_width_)
		kwargs["linewidth"] = *style.line_width_;

	if(style.line_style_)
		kwargs["linestyle"] = *style.line_style_;

	if(style.line_color_) {
		const auto& c = *style.line_color_;

		kwargs["edgecolor"] = py::make_tuple(
			c.r, c.g, c.b, c.a
		);
	}
	if(style.do_fill) {
		const auto& c = *style.fill_color_; // shouldn't throw!

		kwargs["facecolor"] = py::make_tuple(
			c.r, c.g, c.b, c.a
		);
		kwargs["fill"] = py::bool_(true);
	}

    return kwargs;
}

static auto make_line_kwargs(const HistStyle& style) -> py::dict {
	py::dict kwargs;

	if(style.label_ && !style.label_->empty())
		kwargs["label"] = *style.label_;

	if(style.line_width_)
		kwargs["linewidth"] = *style.line_width_;

	if(style.line_style_)
		kwargs["linestyle"] = *style.line_style_;

	if(style.line_color_) {
		const auto& c = *style.line_color_;
		kwargs["color"] = rgba(c);
	}

	return kwargs;
}

static void render_stairs(
	py::object& ax,
	const Hist1DData& h
) {
	auto values = numpy_view(h.values);
	auto edges  = numpy_view(h.edges);

	auto kwargs = make_stairs_kwargs(h.style);

	ax.attr("stairs")(
		values,
		edges,
		**kwargs
	);
}

static void render_errors(
	py::object& ax,
	const Hist1DData& h
) {
	auto centers_data = bin_centers(h);

	auto centers = numpy_view(centers_data);
	auto values  = numpy_view(h.values);
	auto errors  = numpy_view(h.errors);

	auto kwargs = make_line_kwargs(h.style);
	kwargs["marker"] = h.style.marker_.value_or("o");

	ax.attr("errorbar")(
		centers,
		values,
		py::arg("yerr") = errors,
		**kwargs
	);
}

static void render_points(
	py::object& ax,
	const Hist1DData& h
) {
	auto centers_data = bin_centers(h);

	auto centers = numpy_view(centers_data);
	auto values  = numpy_view(h.values);

	auto kwargs = make_line_kwargs(h.style);

	kwargs["linestyle"] = "None";
	kwargs["marker"] = h.style.marker_.value_or("o");

	ax.attr("plot")(
		centers,
		values,
		**kwargs
	);
}

static void render_hist1d(
	py::object& ax,
	const Hist1DData& h
) {
	assert(h.edges.size()  == h.values.size() + 1);
	assert(h.errors.size() == h.values.size());

	switch(h.style.mode_) {
		case HistMode::stairs:
			return render_stairs(ax, h);

		case HistMode::errors:
			return render_errors(ax, h);

		case HistMode::points:
			return render_points(ax, h);
	}
}

void Figure::save(const std::filesystem::path& path) const {
	/* First: create the output directories. */
	std::filesystem::create_directories( path.parent_path() );
	ensure_python();

	/* Important if save() is ever called from a thread other than
	 * the one which initialized Python. */
	py::gil_scoped_acquire gil;

	auto plt = py::module_::import("matplotlib.pyplot");

	auto result =
		plt.attr("subplots")().cast<py::tuple>();

	py::object fig = result[0];
	py::object ax  = result[1];

	try {
		/* For now we intentionally support only TH1-derived data. */
		for(const auto& drawable : objects_) {
			if(const auto* hist = std::get_if<Hist1DData>(&drawable)) {
				render_hist1d(ax, *hist);
			}
			else {
				throw std::runtime_error(
					"mnd::plot::Figure::save(): "
					"TH2/TGraph rendering is not implemented yet."
				);
			}
		}

		if(!xlabel_.empty())
			ax.attr("set_xlabel")(xlabel_);

		if(!ylabel_.empty())
			ax.attr("set_ylabel")(ylabel_);

		if(!title_.empty())
			ax.attr("set_title")(title_);

		if(xlim_) {
			ax.attr("set_xlim")(
				xlim_->first,
				xlim_->second
			);
		}

		if(ylim_) {
			ax.attr("set_ylim")(
				ylim_->first,
				ylim_->second
			);
		}

		if(logx_)
			ax.attr("set_xscale")("log");

		if(logy_)
			ax.attr("set_yscale")("log");

		if(grid_) {
			ax.attr("grid")(true);
			ax.attr("set_axisbelow")(true);
			ax.attr("grid")(
				"color"_a = "gray",
				"linestyle"_a = "dashed"
			);
		}

		if(hide_right_and_top_spline_) {
			ax.attr("spines")["right"].attr("set_visible")(false);
			ax.attr("spines")["top"].attr("set_visible")(false);
		}

		if(legend_)
			ax.attr("legend")();

		fig.attr("tight_layout")();

		/*
		 * matplotlib infers PDF/PNG/SVG/etc. from the extension.
		 */
		fig.attr("savefig")(path.string());

		plt.attr("close")(fig);
	}
	catch(...) {
		/* Don't leak matplotlib figures if a Python call throws. */
		plt.attr("close")(fig);
		
		throw;
	}

	fprintf(stdout, "[PYTHON]: Figure::save(): output file \'%s\' saved as: \'%s\'\n",
		path.filename().c_str(), path.c_str());
}


