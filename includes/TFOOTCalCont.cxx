#include "TFOOTCalCont.h"
#include "util/JSONParser.h"
#include "nlohmann/json.hpp"
#include "util/json_struct_def.hh"
#include "util/FFT.h"
#include "filesystem"

using json = nlohmann::json;
namespace fs = std::filesystem;

std::mutex FOOTParam::mtx_{};

const FMultiPoly* FOOTAsicGainParam::GetPoly(int Z) const noexcept {
	for(const auto& p : multi_poly) { if(p.Z == Z) return &p; }
	return nullptr;
}
FMultiPoly* FOOTAsicGainParam::GetPoly(int Z) noexcept {
	for(auto& p : multi_poly) { if(p.Z == Z) return &p; }
	return nullptr;
}

double FOOTAsicGainParam::GetReferentMeasurement(int Z, double x) const { // slow
	const auto* p = GetPoly(Z);
	if(!p) return NAN;
	return poly::Eval(x, p->pol);
}

bool FOOTAsicGainParam::IsSaneZ() const {
	return std::is_sorted( multi_poly.begin(), multi_poly.end(),
		[](const auto& lhs, const auto& rhs) { return lhs.Z < rhs.Z; }
	);	
}

std::vector<int> FOOTAsicGainParam::GetZ() const {
	std::vector<int> Zs{}; Zs.reserve(multi_poly.size());
	std::transform(multi_poly.begin(), multi_poly.end(), std::back_inserter(Zs), [](const auto& x) { return x.Z; });
	return Zs;
}

bool FOOTAsicGainParam::IsSane(double x, std::vector<double>& evals) const noexcept {
	evals.clear();
	for(const auto& mp : multi_poly) {
		evals.push_back(poly::Eval(x, mp.pol));
	}
	return evals.size() > 0 and std::is_sorted(evals.begin(), evals.end());
}

/* Return an iterator to the first referent measurment line (along E; for set x) that has ref >= e */
auto FOOTAsicGainParam::GetMeasurementBound(double x, double e) const -> typename Vec::const_iterator {
	return std::lower_bound( multi_poly.begin(), multi_poly.end(), std::pair{x,e},
		[](const FMultiPoly& p, const std::pair<double, double>& q) noexcept {
			const auto [x, e] = q;
			return poly::Eval(x, p.pol) < e;
		}
	);
}

std::vector<int> FOOTGainParam::GetNominalZ() const noexcept {
	std::vector<int> Zs{}; Zs.reserve(nominal_value.size());
	std::transform(nominal_value.begin(), nominal_value.end(), std::back_inserter(Zs), [](const auto& x) { return x.Z; });
	return Zs;
}
bool FOOTGainParam::IsSane() const noexcept {
	std::vector<double> evals(20);
	const auto nominalZ = GetNominalZ();

	for(int i=0; i<TFOOTMapCont::N_ASIC; ++i) {
		const auto& asic = fit[i];
		if(!asic.IsSaneZ()) return false;

		const auto Zs = asic.GetZ();
		if(Zs != nominalZ) return false;

		for(int n=0; n < TFOOTMapCont::N_STRIPS_PER_ASIC; ++n) {
			double x = TFOOTMapCont::N_STRIPS_PER_ASIC*i + n + 0.5; // centre of the strip.
			if(!asic.IsSane(x, evals)) return false;
		}
	}
	/* Check also if the referent Z measurements are sorted, and all are unique. */	
	return (
		std::is_sorted(nominal_value.begin(), nominal_value.end()) and
		std::adjacent_find(nominal_value.begin(), nominal_value.end(),
			[](const auto& lo, const auto& hi) { return lo.Z == hi.Z; }) == nominal_value.end()
	);
}

double FOOTGainParam::GetNominalValue(i32 Z) const noexcept {
	auto it = std::find_if( nominal_value.begin(), nominal_value.end(), [Z](const auto& x) { return x.Z == Z; });
	if(it == nominal_value.end()) return NAN;
	return it->val;
}

double FOOTGainParam::CorrectionFactor(double x, double e) const {
	const auto& asic = GetASIC(x);
	double gain;

	auto ref_hi = asic.GetMeasurementBound(x, e);
	if(ref_hi == asic.multi_poly.begin()) { // `e` smaller than the lower bound.
		const auto& [Z, poly] = *ref_hi;
		double avg_val = GetNominalValue(Z);
		double val_x = poly::Eval(x, poly);
		gain = avg_val / val_x;
	} else if(ref_hi == asic.multi_poly.end()) { // `e` bigger than the upper bound.
		const auto& [Z, poly] = asic.multi_poly.back();
		double avg_val = GetNominalValue(Z);
		double val_x = poly::Eval(x, poly);
		gain = avg_val / val_x;
	}
	else { // Iterator valid and not pointed at the range head.
		const auto& [Z_lo, poly_lo] = *(ref_hi - 1);
		const auto& [Z_hi, poly_hi] = *ref_hi;

		double nom_lo = GetNominalValue(Z_lo);
		double nom_hi = GetNominalValue(Z_hi);
		double evaluated_lo = poly::Eval(x, poly_lo);
		double evaluated_hi = poly::Eval(x, poly_hi);
		double gain_lo = nom_lo / evaluated_lo;
		double gain_hi = nom_hi / evaluated_hi;
		auto line = mnd::geom::GetLine2D( {evaluated_lo, gain_lo}, {evaluated_hi, gain_hi} );
		gain = poly::Eval(e, line.array());
	}
	return gain;
}

TH2D* FOOTGainParam::GetHisto (
	int nbins_x, int nbins_y,
	int lo_y, int hi_y
) const {
	TH2D* h = new TH2D("_hFOOTGainParam", "FOOT Gain Parameter",
		nbins_x, 0, TFOOTMapCont::N_STRIPS, nbins_y, lo_y, hi_y );

	for(int ix = 1; ix <= h->GetNbinsX(); ++ix) {
		double x = h->GetXaxis()->GetBinCenter(ix);

		for(int iy = 1; iy <= h->GetNbinsY(); ++iy) {
			double y = h->GetYaxis()->GetBinCenter(iy);

			double value = this->CorrectionFactor(x, y);
			h->SetBinContent(ix, iy, value);
		}
	}
	h->GetXaxis()->SetTitle("Strip #");
	h->GetYaxis()->SetTitle("Cluster ADC");

	return h; /* Draw via: `h->Draw("SURF1")`  :-)  */
}
std::pair<TGraph*, TGraph*> FOOTGainParam::GetGraph (
	const double x, int nbins_y,
	int lo_y, int hi_y
) const {
	TGraph* g = new TGraph(nbins_y);
	const auto& asic = GetASIC(x);
	for(int i=0; i<nbins_y; ++i) {
		double y = lo_y + (i+0.5) * (hi_y - lo_y) / nbins_y; // centre
		double value = this->CorrectionFactor(x, y);
		g->SetPoint(i, y, value);
	}
	g->SetLineWidth(4);
	g->SetLineColor(kRed + 2);
	const std::vector<i32> Zs = this->GetNominalZ();
	TGraph* gpts = new TGraph{};
	for(i32 Z : Zs) {
		double value = asic.GetReferentMeasurement(Z, x);
		double gain = this->CorrectionFactor(x, value);
		gpts->AddPoint(value, gain);
	}

	gpts->SetMarkerStyle(20);
	gpts->SetMarkerSize(1.35);
	return { gpts, g };
	/* Draw via:
	 * auto [pts, graph] = GetGraph();
	 * graph->Draw("AL"); pts->Draw("P SAME"), ; */
}

std::pair<TGraph*, TLine*> FOOTGainParam::GetRefZGraph(
	int Z, int Npts
) const {
	const double nm = this->GetNominalValue(Z);
	if(!std::isfinite(nm))
		throw std::invalid_argument(Form("FOOTGainParam::GetRefZGraph, Passed Z=%d, which isn't found in the config.", Z));

	double xlo = 0, xhi = 640;
	TGraph* g = new TGraph(Npts);
	TLine *lref = new TLine(xlo, nm, xhi, nm);
	for(int i=0; i<Npts; ++i) {
		double x = xlo + (i+0.5)*(xhi - xlo) / Npts;
		const auto& asic = GetASIC(x);
		double val = asic.GetReferentMeasurement(Z, x);
		g->SetPoint(i, x, val);
	}
	constexpr double lwidth_ = 4.0;
	g->SetLineWidth(lwidth_);
	g->SetLineColor(kRed + 2);
	lref->SetLineStyle(2);
	lref->SetLineColor(kMagenta + 2);
	lref->SetLineWidth(lwidth_);

	return {g, lref};
}
[[ nodiscard ]] std::vector<std::tuple<i32, TGraph*, TLine*>>
FOOTGainParam::GetAllRefZGraph(int Npts) const {
	std::vector<std::tuple<i32, TGraph*, TLine*>> r;
	const auto nominalZ = GetNominalZ();
	for(const i32 Z : nominalZ) {
		auto [g,l] = GetRefZGraph(Z, Npts);
		r.emplace_back(Z,g,l);
	}
	return r;
}


double FOOTDeltaFFT::Evaluate(const double x) const {
	return FFTW::Evaluate(this->c, this->n, x, -0.5, 0.5);
}

TGraph* FOOTDeltaFFT::GetGraph(int Npts) const {
	constexpr double x_lo = -0.5;
	constexpr double x_hi =  0.5;

	TGraph* g = new TGraph(Npts);
	for(int i=0; i<Npts; ++i) {
		double x0 = x_lo + (i+0.5) * (x_hi - x_lo) / Npts;
		double y0 = this->Evaluate(x0);
		g->SetPoint(i, x0, y0);
	}
	g->SetLineColor(kRed);
	g->SetLineWidth(3);
	return g;
}

double FOOTDeltaParam::CorrectionBasic(double delta) const noexcept {
	double x = std::min( std::abs(delta), s );
	return 1 - (1-f) / s * x;
}
/* Total factor. Corrected cluster energy is then: e' = e / CorrectionFactor(δ) */
double FOOTDeltaParam::CorrectionFactor(double delta) const noexcept {

	double corr1 = this->CorrectionBasic(delta); // Interval [f, 1]
	double corr2 = this->f2.Evaluate(delta); // ~approx. [-0.3, 0.3]

	return corr1 * (1 + corr2);
}

/* Small helper function to plot what we're actually matching upon. */
[[ nodiscard ]] TGraph* FOOTDeltaParam::GetGraph(int Npts, double x_lo, double x_hi) const {
	TGraph* g = new TGraph(Npts);
	for(int i=0; i<Npts; ++i) {
		double x0 = x_lo + (i+0.5) * (x_hi - x_lo) / Npts;
		double y0 = this->CorrectionFactor(x0);
		g->SetPoint(i, x0, y0);
	}
	g->SetLineColor(kRed);
	g->SetLineWidth(3);
	return g;
}

Orientation FOOTParam::GetOrientation() const noexcept {
	MND_ASSERT(orientation.size() > 0 && "FOOTParam::GetOrientation() but orientation is not specified. Not allowed!");
	if(orientation[0] == 'x' or (orientation.length() > 1 and orientation[1] == 'x')) return Orientation::X;
	if(orientation[0] == 'y' or (orientation.length() > 1 and orientation[1] == 'y')) return Orientation::Y;
	return Orientation::UNKNOWN;
};
double FOOTParam::R() const noexcept {
	MND_ASSERT(orientation.size() > 0 && "FOOTParam::R() but orientation is not specified. Not allowed!");
	if(orientation[0] == '-') return -1.0;
	return 1.0;
}
double FOOTParam::Q(double E) const noexcept {
	if( MND_UNLIKELY(!fit_.load_acq()) )
		Init();
	return std::pow(fit_->f*E, fit_->c);
}

double FOOTParam::E(const RNFOOTCluster& clust) const noexcept {
	double ce = clust.fCE;
	double const d  = clust.Delta();
	double cx = clust.fCX;

	ce *= this->gain.CorrectionFactor(cx, ce);
	ce /= this->de.CorrectionFactor(d);

	return ce;
};

double FOOTParam::Q(const RNFOOTCluster& clust) const noexcept {
	const double e = this->E(clust);
	return this->Q(e);
}

double FOOTParam::BarePosition(const RNFOOTCluster& clust) const noexcept {
	return R() * (clust.fCX - DETECTOR_MIDPOINT) * STRIP_TO_MM;
}
double FOOTParam::X0(const RNFOOTCluster& clust) const noexcept {
	return BarePosition(clust) + delta_p;
}

double FOOTParam::BareToStripnum(double m) const noexcept {
	return DETECTOR_MIDPOINT + R() * m * MM_TO_STRIP;
}
double FOOTParam::X0ToStripnum(double x0) const noexcept {
	double b = x0 - delta_p;
	return BareToStripnum(b);
}

void FOOTParam::Init() const {
	auto lock = std::lock_guard{mtx_};
	
	// Another thread may have initialized it while waiting
	if(fit_.load() == true)
		return;
	std::vector<double> x, y;
	for(auto [Q,E] : gain.nominal_value) {
		x.push_back( std::log(Q) );
		y.push_back( std::log(E) );
	}
	auto r = PolyFit<1>(x,y);
	fit_->f = std::exp(-r[0]);
	fit_->c = 1.0 / r[1];
	fit_.publish(true);
}

RNFOOTCluster::RNFOOTCluster(double x, double e, u32 m, ClusterType t) :
	fCX(x), fCE(e), fCM(m), fCT(t) {}

/* ------------------------------------------------------- */

RNFOOTCal::RNFOOTCal() {
	fCl.reserve(INIT_CAPACITY);
	fRaw.reserve(N_STRIPS);
}
void RNFOOTCal::Clean() noexcept {
	fCl.clear();
	fRaw.clear();
    w_ev_type = EventType::Undetermined;
}
std::vector<double> RNFOOTCal::E() const noexcept {
	std::vector<double> res;
	res.reserve(fCl.size());
	for(auto const& c : fCl) res.push_back(c.fCE);
	return res;
}
std::vector<double> RNFOOTCal::X() const noexcept {
	std::vector<double> res;
	res.reserve(fCl.size());
	for(auto const& c : fCl) res.push_back(c.fCX);
	return res;
}

/* ------------------------------------------------------- */

static std::string full_setup_file_name {};

/* Initialized with 2 keys: "ID" and "Setup" */
void TFOOTCalCont::Init(TDictInfo info) {
	constexpr const char* id_key = "ID";
	auto it = info.find(id_key);
	if(it == info.end())
		ERROR("\"%s\" key not found in the info hashmap.", id_key);
	try {
		FOOT_N = std::stoi(it->second);
	} catch(const std::exception& e) {
		ERROR("For key \"%s\", value found: " EMPH(%s) " but unparsable to integer. Err: %s", id_key, it->second.c_str(), e.what());
	}

	it = info.find("Setup");	
	if(it == info.end())
		ERROR("TFOOTCalCont::Init(): \"Setup\" key not found for info (%s).\n" 
			, mnd::type_name<TDictInfo>().c_str());
	const std::string& file_name = it->second;
	json j = ParseJSON(file_name);
	
	auto j_it = j.find(Form("FOOT%d", FOOT_N));
	if(j_it == j.end()) 
		ERROR("Trying to set up FOOT[%d] cal, but setup file doesn't contain the key \"FOOT%d\".", FOOT_N, FOOT_N);
	
	json& jf = j_it.value();
	UNROLL_JSON_PARAM(par, jf, 9);

	par.de10_index_ = FOOT_N;
	if(par.N == -1)
		ERROR("Parsed the setup file fine, but the \"N\" table entry for FOOT%d not found. "
			"It is mandatory to label the FOOT's!", FOOT_N);
	/* Note, if two FOOT's in the setup have identical `.N` field, it will throw a
	 * wildest RNTuple error 'cause of columns of identical name... */
	this->SetName(Form("FOOT%d", par.N));

	if(should_register_box_) {
		if(!j.contains("box"))
			ERROR("Container meant to also register the FOOT box, but inside %s file, \"box\" key not found?",
				file_name.c_str());
		
		UNROLL_JSON_PARAM(bpar, j["box"], 7);
		
		/* `file_name` could've been passed relative. Fetch the full path and bake it in. */
		try {
			fs::path pfull = fs::canonical(file_name);
			full_setup_file_name = pfull.c_str();
		} catch(std::filesystem::filesystem_error const& e) {
			ERROR("Fetching full path from: \'%s\' failed? Reason: %s", file_name.c_str(), e.what());
		}
	}
}

void TFOOTCalCont::Setup() {
	if(strlen(GetName()) == 0) ERROR("Setup called before Init() ? Or name field set to '' ?");
	h1_mult        = RegisterObject<TH1I>("h1_mult", Form("FOOT (%2d: %d) multiplicity", FOOT_N, par.N), 200,0,50);
	h1_dE          = RegisterObject<TH1I>("h1_dE", Form("FOOT(%2d: %d) energy", FOOT_N, par.N), 5000,0,1000);
	h1_X           = RegisterObject<TH1I>("h1_x", Form("FOOT(%2d: %d) clust position", FOOT_N, par.N), 5*N_STRIPS, 0, N_STRIPS);
	h1_cl_type     = RegisterObject<TH1I>("h1_cl_type", Form("(%2d:%d) clust type", FOOT_N, par.N), 7, 0, 7);
	h1_dE_m3_light = RegisterObject<TH2I>("h1_dE_m3_light", Form("(%2d:%d) energy with multiplicity 3 (bare) and below heavy-ion cutoff",
		FOOT_N, par.N), N_STRIPS,0,N_STRIPS,400,0,200);
	h1_dE_m3_heavy = RegisterObject<TH2I>("h1_dE_m3_heavy", Form("(%2d:%d) energy with multiplicity 3 (bare) and above heavy-ion cutoff",
		FOOT_N, par.N), N_STRIPS,0,N_STRIPS,2000,0,2000);
	h2_mult_e = RegisterObject<TH2I>("h2_mult_e", Form("FOOT(%2d:%d) multiplicity and energy", FOOT_N, par.N),
		2000,0,2000, 10, 1, 10);
	h2_mult_e->GetXaxis()->SetTitle("Individual energy per strip hit");
	h2_mult_e->GetYaxis()->SetTitle("Multiplicity");

	setup = RegisterObject<FOOTParam>("setup", this->par /* copy ctor */);
	
	/* In the cal step, box object not used. It's just written down for helper scripts to process stuff. */
	if(should_register_box_) {
		box = RegisterObject<FOOTBoxParam>("box", this->bpar);
		setupName = RegisterObject<std::string>("setup_file", mnd::noop_fn<std::string>(), full_setup_file_name); 
	}
}

ClassImp(FOOTClusterFit);
ClassImp(RNFOOTCluster);
ClassImp(RNFOOTCal);
ClassImp(FMultiPoly);
ClassImp(FOOTAsicGainParam);
ClassImp(FOOTReferentADCMeasurement);
ClassImp(FOOTGainParam);
ClassImp(FOOTDeltaFFT);
ClassImp(FOOTDeltaParam);
ClassImp(FOOTParam);
ClassImp(ExpertTarget);
ClassImp(FOOTBoxParam);
