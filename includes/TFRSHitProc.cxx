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
double TFRSHitProc::z0 = NAN;

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
		for(const auto& tofp : out.sTof->ToF) {
			if(!tofp.IsValid())
				ERROR("Parsed in ToF combination marked as invalid:\nIt is: %s\n",
					mnd::streamable(tofp).c_str());
		}
		WARN("'%s' JSON successfully parsed as: %s\n", tof_param_json_name, mnd::streamable(*out.sTof).c_str());
	}

	const std::string_view fname = mnd::g_input_file.name();
	WARN("TFRSHitProc::TFRSHitProc(..) found current input file from MONAD: %s%.*s%s\n",
		KBH_MAG, (int)fname.size(), fname.data(), KNRM);

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
	binfo->GetFrom(runsheet_row);
	
	WARN("Beam impinging on S2 target identified as %s%s%s\n", KBH_MAG,
		binfo->SecondaryIon()->chem_to_string().c_str(), KNRM
	);
}

void TFRSHitProc::ProcessEntry() noexcept {
	RNFRSHit& out = (this->out).inner();
	out.Clean();

	out.cal = std::get<0>(this->in).inner(); // RNFRSCal& operator=(RNFRSCal& )
	
	ProcessToF();
	ProcessS2BT();
    ProcessS2AT();
    ProcessS3();
}

void TFRSHitProc::ProcessToF() noexcept {
	const TFRSCalCont& cal = std::get<0>( this->in );
	const RNFRSCal& in = cal.inner();

    const SCIParam& sci21_p   = cal.sci_param->operator[]( RNFRSCal::SCI21_I );
    const SCIParam& sci22_p   = cal.sci_param->operator[]( RNFRSCal::SCI22_I );

	const double beta_21_22 = tofp_s21_s22
		? tofp_s21_s22->Beta(in)
		: NAN;

	const double beta_22_31 = tofp_s22_s31
		? tofp_s22_s31->Beta(in)
		: NAN;
	
	/* C++11: 'If control enters the declaration concurrently while the variable is being initialized,
	 *         the concurrent execution shall wait for completion of the initialization'
	 * https://timsong-cpp.github.io/cppwp/n3337/stmt.dcl , §4 */
	static const double lambda = (
		(sci21_p.z0 > 0 && sci22_p.z0 > 0 && z0 > 0 &&
		 mnd::is_strictly_ascending(sci21_p.z0, z0, sci22_p.z0))
		? (z0 - sci21_p.z0)/(sci22_p.z0 - sci21_p.z0)
		: NAN
	);
	/* SCI21 - SCI22 ToF:
	 * λ/β₁ = 1/β - (1-λ)/β₂ */
	const double lambda_over_b1 = 1.0 / beta_21_22 - (1-lambda)/beta_22_31;
	const double beta_21 = lambda/lambda_over_b1;

	out->s2_bt.beta = beta_21;
	out->s2_at.beta = beta_22_31;
	out->s3.beta    = beta_22_31;
}

/* Before target we don't have full Q vs. A/Q measurement, since there's no 
 * ToF S1-S2. A we just blindly take from the setup in this case. */
void TFRSHitProc::ProcessS2BT() noexcept {
	x.clear(); y.clear();
	zx.clear(); zy.clear();

	const TFRSCalCont& cal = std::get<0>( this->in );
	const RNFRSCal& in = cal.inner();
	RNFRSHit::Id& bt = out.inner().s2_bt;

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

	const double beta = bt.beta; /*from ProcessToF()*/
	const double E0 = sci21_p.E0( sci21_data );
	const double Q0 = sci21_p.Q(E0, NAN);
	const double Qc = sci21_p.Q(E0, beta);
	this->s2_q0[0] = Q0;
	this->s2_qc[0] = Qc;
	
	const double AoQ = phy::AoQ(
		phy::Brho_t{ out.binfo->BrhoAt(RNFRSCal::SCI22_I) },
		phy::Beta_t{ beta }
	);

	bt.AoQ = AoQ;
	bt.Q   = Qc;
	bt.x0  = x0;
	bt.y0  = y0;
	bt.ax  = ax;
	bt.ay  = ay;
	// bt.beta assigned in `ProcessTof()`.
	bt.code = code;

	g_upstream_track = RNTrackToLine3D(bt); // could be null.
}

void TFRSHitProc::ProcessS2AT() noexcept {
    const TFRSCalCont& cal = std::get<0>( this->in );
	const RNFRSCal& in = cal.inner();
	RNFRSHit::Id& at = out.inner().s2_at;

	const RNSciCal& sci21_data = in.sci[ RNFRSCal::SCI21_I ];
	
    const SCIParam& sci22_p   = cal.sci_param->operator[]( RNFRSCal::SCI22_I );
    const RNSciCal& sci22_data = in.sci[ RNFRSCal::SCI22_I ];
	
	const double beta = at.beta; /*from ProcessToF()*/
	const double E0 = sci22_p.E0( sci22_data );
	const double Q0 = sci22_p.Q(E0, NAN);
	const double Qc = sci22_p.Q(E0, beta);
	this->s2_q0[1] = Q0;
	this->s2_qc[1] = Qc;

	out.h2_s2_q ->Fill( s2_q0[0], s2_q0[1] );
	out.h2_s2_qc->Fill( s2_qc[0], s2_qc[1] );

	at.Q    = Qc;
	// at.beta assigned in `ProcessTof()`
	at.code |= (sci21_data.hits.size() << 32);
	at.code |=  sci22_data.hits.size();

	/* Other fields are unmeasurable or not measurable to wanted precision, like
	 * position and angles. */
}

void TFRSHitProc::ProcessS3() noexcept {
    const TFRSCalCont& cal = std::get<0>( this->in );
	const RNFRSCal& in = cal.inner();
	RNFRSHit::Id& s3 = out.inner().s3;

	const RNSciCal& sci22_data = in.sci[ RNFRSCal::SCI22_I ];
	
    const SCIParam& sci31_p   = cal.sci_param->operator[]( RNFRSCal::SCI31_I );
    const RNSciCal& sci31_data = in.sci[ RNFRSCal::SCI31_I ];

	const double beta = out->s3.beta;
	const double E0 = sci31_p.E0( sci31_data );
	const double Q0 = sci31_p.Q(E0, NAN);
	const double Qc = sci31_p.Q(E0, beta);

	const double AoQ = phy::AoQ(
		phy::Brho_t{ out.binfo->BrhoAt(RNFRSCal::SCI31_I) },
		phy::Beta_t{ beta }
	);

	out.h2_s3_id->Fill(AoQ, Qc);
	out.h2_s3_q ->Fill( s2_q0[1], Q0 );
	out.h2_s3_qc->Fill( s2_qc[1], Qc );
	out.h1_s3_beta->Fill(beta);

	s3.AoQ  = AoQ;
	s3.Q    = Qc;
	// s3.beta assigned in `ProcessTof()`
	s3.code |= (sci22_data.hits.size() << 32);
	s3.code |=  sci31_data.hits.size();
}

void TFRSHitProc::FinalInit() {
	/* Is const-qualified and once loaded, won't change. */
	tofp_s22_s31 = out.sTof
		? out.sTof->Get(
			RNFRSCal::SCI22_I,
			RNFRSCal::SCI31_I
		)
		: nullptr;
	
	tofp_s21_s22 = out.sTof
		? out.sTof->Get(
			RNFRSCal::SCI21_I,
			RNFRSCal::SCI22_I
		)
		: nullptr;
}
