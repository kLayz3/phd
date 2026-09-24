#include "ExperimentAssumptions.hh"
#include "TFRSCalCont.h"
#include "TFRSHitCont.h"
#include "TFRSHitProc.h"
#include "util/PolyFitter.h"
#include <cmath>

#include "util/Geometry.h"
#include "util/MPhysics.h"
#include "util/RunsheetParser.h"
#include "util/json_struct_def.hh"

thread_local mnd::geom::Line3D g_upstream_track {};

/* This ctor only gets called once, at the creation. Later on, the clones
 * call the implicit copy-ctor. */
TFRSHitProc::TFRSHitProc (
	TFRSHitCont& out,
	const TFRSCalCont& in,
	int s2_bt_mask,
	std::string_view runsheet_name
) : TFRSHitProc::Base(out, in),
	s2_bt_tracking_mask(s2_bt_mask)
{
	x.reserve(8); y.reserve(8);
	zx.reserve(8); zy.reserve(8);
	
	extern double g_expert_target_z; // extern'ed from: `includes/TFOOTHitCont.cxx`
	z0 = g_expert_target_z;

	out.h2_target_xy->SetTitle(Form("%s@z0=%.1f", out.h2_target_xy->GetTitle(), z0));

	/* Next, load the ToF parameter objects from the cal's setup json object. */
	if(in.setup.empty())
		ERROR("TFRSHitProc::TFRSHitProc(..): needing to init the ToF parameter, but cal's setup JSON object is empty? "
			"Did you call TFRSCalCont::Init(..) ?");
	
	constexpr auto tof_param_json_name = FRSToFParam::get_name<0>();
	if(in.setup.find(tof_param_json_name) == in.setup.end()) {
		ERROR("'%s' parameter not found in the input JSON. Cannot continue (no ToF).\n",
			tof_param_json_name);
	} else {
		UNROLL_JSON_PARAM( (*out.sTof), in.setup, 0);
		WARN("'%s' JSON successfully parsed as: %s\n", tof_param_json_name, mnd::streamable(*out.sTof).c_str());
	}

	const std::string_view fname = mnd::g_input_file.name();
	WARN("TFRSHitProc::TFRSHitProc(..) found current input file from MONAD: %s%.*s%s\n",
		KBH_MAG, (int)fname.size(), fname.data(), KNRM);

	/* Self clarification: `in` is the input container, owning the underlying `_vc` vector
	 * pinned on the heap via shared_ptr API. AKA: the raw pointers such as in.sci_param will just get copied around
	 * during the ping-ponging of constructing the full TAnalysisProcess<..> and TAnalysisPool<..> but
	 * still refer to the correct `TOnce<T>` inside the _vc vector.
	 *
	 * It's only important to refer to the CAL step's converter, not output current HIT step's. */
#define INIT_DE_TO_Q(SCI_LABEL) \
    { \
        const SCIParam& s = in.sci_param->at( RNFRSCal::SCI##SCI_LABEL##_I ); \
        u32 r = s.SetConverter(fname); \
        if(r > 1) { \
            WARN("SCI%s parameter: found >1 (out of %zu) \"de_to_q\" matches for current input file: '%.*s'. Is OK, will take last match.\n", \
                 #SCI_LABEL, s.de_to_q.size(), (int)fname.size(), fname.data()); \
        } else if(r == 0) { \
            WARN("SCI%s parameter: found =0 (out of %zu) \"de_to_q\" matches for current input file: '%.*s'. Is OK, charge values will be NAN.\n", \
                 #SCI_LABEL, s.de_to_q.size(), (int)fname.size(), fname.data()); \
        } \
		else { \
			WARN("SCI%s parameter: " BOLD "successfully attached the charge converter.\n" KNRM, #SCI_LABEL); \
		} \
    }
    INIT_DE_TO_Q(21);
    INIT_DE_TO_Q(22);
    INIT_DE_TO_Q(31);

	mnd::fs::load_runsheet(runsheet_name);
	const mnd::RunsheetState runsheet_row = mnd::QueryRunsheet(fname);

	if(!runsheet_row.primary.valid() || (runsheet_row.primary != mnd::assume::primary))
		WARN("Primary beam identified as %s%s%s but is mismatched with assumption %s%s%s. Is fine, just to be noted.\n",
			KBH_RED, mnd::streamable(runsheet_row.primary).c_str(), KNRM,
			KBH_MAG, mnd::streamable(mnd::assume::primary).c_str(), KRNM);

	if(!runsheet_row.secondary.valid())
		ERROR("Runsheet parsed fine, but secondary beam: %s is marked invalid?",
			runsheet_row.secondary.chem_to_string().c_str());
	
	WARN("Runsheet path loaded and identified as %s%.*s%s, with corresponding row being:\n  %s\n",
		BOLD, (int)runsheet_name.size(), runsheet_name.data(), KNRM, mnd::streamable(runsheet_row).c_str());

	BeamInfo* binfo = out.binfo;
	if(!binfo) ERROR("Forgot to call TFRSHitCont::Setup() ?");
	binfo->A0 = runsheet_row.secondary.A;
	binfo->Z0 = runsheet_row.secondary.Z;
	binfo->R0 = runsheet_row.brho.s1_s2;
	binfo->R1 = runsheet_row.brho.s2_s3;
	binfo->R2 = runsheet_row.brho.s3_s4;
	WARN("Beam impinging on S2 target identified as %s%s%s\n", KBH_MAG,
	  phy::Nucleus{
		.A = binfo->A0,
		.Z = binfo->Z0
	  }.chem_to_string().c_str(), KNRM
	);
}

void TFRSHitProc::ProcessEntry() noexcept {
	RNFRSHit& out = (this->out).inner();
	out.Clean();

	out.cal = std::get<0>(this->in).inner(); // RNFRSCal& operator=(RNFRSCal& )
	ProcessS2BT();
    ProcessS2AT();
    ProcessS3();
}

/* Before target we don't have full Q vs. A/Q measurement, since there's no 
 * ToF S1-S2. A we just blindly take from the setup in this case. */
void TFRSHitProc::ProcessS2BT() noexcept {
	x.clear(); y.clear();
	zx.clear(); zy.clear();

	const TFRSCalCont& cal = std::get<0>( this->in );
	const RNFRSCal& in = cal.inner();

	i64 code{0};
	if(s2_bt_tracking_mask & RNFRSHit::S2_BT_TRACKING_INCLUDE_SCI21_MASK) {
		static const double sci21_z = cal.sci_param->at(0).z0;
		const std::vector<RNSciCal::Measurement>& hits = in.sci[0].hits;
		if(hits.size() == 1) { // otherwise can't resolve.
			x.push_back( hits[0].x );
			zx.push_back( sci21_z );
			code |= RNFRSHit::S2_BT_TRACKING_INCLUDE_SCI21_MASK;
		}
	}

#define TRY_INCLUDE_TPC_INTO_BT_TRACKING(LABEL, INDEX) \
	if(s2_bt_tracking_mask & RNFRSHit::S2_BT_TRACKING_INCLUDE_TPC##LABEL##_MASK) { \
		const double ztpc = cal.tpc_param->at(INDEX).z0; \
		const RNTPCCal& tpc = in.tpc[INDEX]; \
		const double xtpc = tpc.X0(); \
		const double ytpc = tpc.Y0(); \
		if( std::isfinite(xtpc) and std::isfinite(ytpc) ) { \
			x.push_back( xtpc ); \
			y.push_back( ytpc ); \
			zx.push_back( ztpc ), zy.push_back( ztpc ); \
			\
			code |= RNFRSHit::S2_BT_TRACKING_INCLUDE_TPC##LABEL##_MASK; \
		} \
	} \
	MND_EMPTY_MACRO(INDEX)
	
	TRY_INCLUDE_TPC_INTO_BT_TRACKING(21, 0)
	TRY_INCLUDE_TPC_INTO_BT_TRACKING(22, 1)
	TRY_INCLUDE_TPC_INTO_BT_TRACKING(23, 2)
	
	double x0{NAN}, ax{NAN}, y0{NAN}, ay{NAN};
	double& xT = out.inner().xT;
	double& yT = out.inner().yT;

	if(x.size() >= 2) {
		PolyFit<1>(zx, x, this->fit_result);
		x0 = fit_result[0];
		ax = fit_result[1];
		xT = poly::Eval(z0, fit_result);
	}
	if(y.size() >= 2) {
		PolyFit<1>(zx, y, this->fit_result);
		y0 = fit_result[0];
		ay = fit_result[1];
		yT = poly::Eval(z0, fit_result);
	}
	
	out.h2_target_xy->Fill(xT, yT);

    const SCIParam& sci21_p   = cal.sci_param->operator[]( RNFRSCal::SCI21_I );
    const RNSciCal& sci21_data = in.sci[ RNFRSCal::SCI21_I ];

	double Q = sci21_p.Q( sci21_data );
	
	/* Since S2 'Q' measurement can vary, plugging in a constant 'A'
	 * number into beta calculation would be wrong.
	 * C++11: 'If control enters the declaration concurrently while the variable is being initialized,
	 *         the concurrent execution shall wait for completion of the initialization'
	 * https://timsong-cpp.github.io/cppwp/n3337/stmt.dcl?utm_source=chatgpt.com , Paragraph [4] */
	static double const s2_beta_assume = phy::Beta(
		out.binfo->A0,
		out.binfo->Z0,
		phy::Brho_t{ out.binfo->R0 }
	);
	
	RNFRSHit::Id& bt = out.inner().s2_bt;
	/* bt.A ==> should be measured, not assumed! */
	bt.Q    = Q;
	bt.x0   = x0;
	bt.y0   = y0;
	bt.ax   = ax;
	bt.ay   = ay;
	bt.beta = s2_beta_assume;
	bt.code = code;

	g_upstream_track = RNTrackToLine3D(bt); // could be null.
}

void TFRSHitProc::ProcessS2AT() noexcept {
    const TFRSCalCont& cal = std::get<0>( this->in );
	const RNFRSCal& in = cal.inner();
	
    const RNSciCal& sci21_data = in.sci[ RNFRSCal::SCI21_I ];

    const SCIParam& sci22_p   = cal.sci_param->operator[]( RNFRSCal::SCI22_I );
    const RNSciCal& sci22_data = in.sci[ RNFRSCal::SCI22_I ];

	RNFRSHit::Id& at = out.inner().s2_at;
    double Q = sci22_p.Q( sci22_data );
	
	out.h2_s2_q->Fill( out.inner().s2_bt.Q, Q );

	at.Q    = Q;
	at.code = (sci21_data.hits.size() << 32)
		| sci22_data.hits.size();

	/* Other fields are unmeasurable or not measurable to wanted precision, like
	 * position and angles. */
}

void TFRSHitProc::ProcessS3() noexcept {
    const TFRSCalCont& cal = std::get<0>( this->in );
	const RNFRSCal& in = cal.inner();
	
    const SCIParam& sci31_p   = cal.sci_param->operator[]( RNFRSCal::SCI31_I );
    const RNSciCal& sci31_data = in.sci[ RNFRSCal::SCI31_I ];

    const RNSciCal& sci22_data = in.sci[ RNFRSCal::SCI22_I ];

	double Q = sci31_p.Q( sci31_data );

	const double beta = tofp_s3_s2
		? tofp_s3_s2->Beta(in)
		: NAN;

	const double gamma = phy::Gamma(beta);

	double AoQ = phy::AoQ(
		phy::Brho_t{ out.binfo->R1 },
		phy::BetaGamma_t{ beta*gamma }
	);

	out.h2_s3_id->Fill(AoQ, Q);
	out.h2_s3_q->Fill( out.inner().s2_bt.Q, Q );
	out.h1_s3_beta->Fill(beta);

	/* Local references. */
	RNFRSHit::Id& s3 = out.inner().s3;

	s3.Q    = Q;
	s3.A    = AoQ * Q;
	s3.beta = beta;
	s3.code = (sci22_data.hits.size() << 32)
		| sci31_data.hits.size();
}

void TFRSHitProc::FinalInit() {
	/* Is const-qualified and once loaded, won't change. */
	tofp_s3_s2 = out.sTof
		? out.sTof->Get(
			RNFRSCal::SCI22_I,
			RNFRSCal::SCI31_I
		)
		: nullptr;

}
