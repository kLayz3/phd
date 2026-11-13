#include "TFOOTHitCont.h"

#include "monad/monad.hxx"
#include "util/JSONParser.h"
#include "util/Geometry.h"
#include "util/FloatHacks.hxx"

#include "TParameter.h"

using json = nlohmann::json;

static nlohmann::json setup {};
static FOOTBoxParam _box {};

double g_expert_target_z {};

std::ostream& operator<<(std::ostream& os, const FOOTHit& rhs) noexcept {
	os << mnd::msg("(%s%.2f%s,%d,%s%5.1f%s)",
		KBH_CYN, rhs.Q.q, KNRM,
		rhs.Q.fCM,
		KBH_RED, rhs.m, KNRM);
	return os;
}

bool RNFOOTPair::HasDataX() const noexcept {
    return mnd::get_bit<_DATA_PRESENT_BITINDEX>(zx);
}
bool RNFOOTPair::HasDataY() const noexcept {
    return mnd::get_bit<_DATA_PRESENT_BITINDEX>(zy);
}
bool RNFOOTPair::HasData() const noexcept {
	return HasDataX() && HasDataY();
}

u32 RNFOOTHit::DataPresentCount() const noexcept {
	u32 sum = 0;
	mnd::static_for_each(pair, [&](const RNFOOTPair& fp) {
		sum += static_cast<u32>(fp.HasData());
	});
    return sum;
}
bool RNFOOTHit::HasData() const noexcept {
    return DataPresentCount() == pair.size();
}

TFOOTHitCont::TFOOTHitCont() : TContainer("FOOT") {}

RNFOOTTrack::RNFOOTTrack (
	const std::array<double, 2>& xline,
	const std::array<double, 2>& yline,
	double Q_,
	double score_,
	std::size_t n_
#ifdef MND_FOOTTRACK_DEBUG
		,
		const std::array<double, N_PAIRS>& u_,
		const std::array<double, N_PAIRS>& v_,
		const std::array<double, N_PAIRS>& zx_,
		const std::array<double, N_PAIRS>& zy_,
		const std::array<double, N_PAIRS>& qx_,
		const std::array<double, N_PAIRS>& qy_,
		const std::array<double, N_PAIRS>& e0_x_,
		const std::array<double, N_PAIRS>& delta_x_,
		const std::array<u32,    N_PAIRS>& c0_x_,
		const std::array<double, N_PAIRS>& e0_y_,
		const std::array<double, N_PAIRS>& delta_y_,
		const std::array<u32,    N_PAIRS>& c0_y_
#endif
) : x0(xline[0]), ax(xline[1]), y0(yline[0]), ay(yline[1]),
	Q(Q_), score(score_), n(n_)
#ifdef MND_FOOTTRACK_DEBUG
		,
		_u(u_),
		_v(v_),
		_zx(zx_),
		_zy(zy_),
		_qx(qx_),
		_qy(qy_),
		_e0_x(e0_x_),
		_delta_x(delta_x_),
		_c0_x(c0_x_),
		_e0_y(e0_y_),
		_delta_y(delta_y_),
		_c0_y(c0_y_)
#endif
{}

/* Ok, hacks. Usually two tracks compare equal if their geometry is the same.
 * We can trust it's almost impossible for two distinct tracks to have exactly the same geometry.
 * Then this minimizes to comparison of 4 doubles. As long as we don't fuck around with the member order! */
bool operator==(const RNFOOTTrack& l, const RNFOOTTrack& r) noexcept {
#ifdef MND_HAS_INTEL_INTRIN
	const __m256d a = _mm256_loadu_pd(&l.x0);
	const __m256d b = _mm256_loadu_pd(&r.x0);
	const __m256d cmp = _mm256_cmp_pd(a, b, _CMP_EQ_OQ); // ordered comparison

	return _mm256_movemask_pd(cmp) == 0xf;
#else
	return l.x0 == r.x0 &&
	       l.y0 == r.y0 &&
	       l.ax == r.ax &&
	       l.ay == r.ay;
#endif
}

constexpr static u32 NP = RNFOOTHit::N_PAIRS;

std::unique_ptr<TrackingGeometry<NP>>
tg_from_stg(const STrackingGeometry& stg) {
	auto result = std::make_unique<TrackingGeometry<NP>>();

	for(u32 i = 0; i < NP; ++i) {
		const auto& saved = stg.layer[i];
		auto& layer = result->layer[i];
		for(Eigen::Index row = 0; row < 2; ++row) {
			for(Eigen::Index col = 0; col < 4; ++col)
				layer.A(row, col) = saved.A[4 * row + col];
			layer.dxy(row) = saved.dxy[row];
			layer.z(row) = saved.z[row];
		}
	}

	mnd::static_for<2, NP + 1>([&](auto K_) {
		constexpr u32 K = decltype(K_)::value;
		using DMatrix = DesignMatrix<NP, K>;
		auto& cached = result->fit.template get<K>();

		for(u32 index = 0; index < DMatrix::n_combinations; ++index) {
			typename DMatrix::Matrix matrix;
			typename DMatrix::Vector offsets;
			const u32 mask = mnd::combi::combo_lookup_table<NP, K>[index];
			Eigen::Index row = 0;
			for(u32 i = 0; i < NP; ++i) {
				if(mask & (u32{1} << i)) {
					matrix.template middleRows<2>(row) = result->layer[i].A;
					offsets.template segment<2>(row) = result->layer[i].dxy;
					row += 2;
				}
			}

			auto& parameters = cached.lookup_table[index];
			parameters.first.compute(matrix);
			parameters.second = offsets;
			if(parameters.first.rank() != 4)
				ERROR("tg_from_stg(): rank() !=4 for mask %#08x.", mask);
		}
	});

	return result;
}

void TFOOTHitCont::Init(TDictInfo info) {
auto it = info.find("Setup");
	if(it == info.end())
		ERROR("Setup key not found for info (%s).\n", mnd::type_name<TDictInfo>().c_str());
	const std::string& file_name = it->second;
	
	setup = ParseJSON(file_name);
	setup["file_name"] = file_name;
	
	if(!setup.contains("box")) ERROR("\'box\' key not found in JSON file: %s\n", file_name.c_str());

	UNROLL_JSON_PARAM(_box, setup["box"], 7);

	/* Single global to be shipped over to TFRSHitProc */
	g_expert_target_z = _box.GetTargetZ();
}

using FOOTParams = std::array<FOOTParam, 2>;
using A3 = std::array<double,3>;
template<> void Add(FOOTParams&, const FOOTParams&) {}
template<> void Add(A3&, const A3& ) {}

void TFOOTHitCont::Setup() {
	h1_qtrack = RegisterObject<TH1I>("h1_qtrack", "Charge (Q) of recognized tracks",
		100, 0, 8);
	h1_track_nsampled = RegisterObject<TH1I>("h1_qtrack_nsampled", "Points sampled (N) in the track",
		12, -0.5, 5.5);
	diff_heavy_frag_vs_upstream = RegisterObject<TH1I>("diff_heavy_frag_vs_upstream",
		"Distance heaviest fragment track to TPC track;d [mm]", 500, 0.0, 5);
	
	box = RegisterObject<FOOTBoxParam>("box", _box);
	setupName = RegisterObject<std::string>("distinctup_file", mnd::noop_fn<std::string>(), setup["file_name"].get_ref<const std::string&>());
	geometry = RegisterObject<STrackingGeometry>("geometry", mnd::noop_fn<STrackingGeometry>(), STrackingGeometry{});

	cost_coeff = RegisterObject<A3>("cost_coeff", {}); // filled by Proc ctor
	max_cost = RegisterObject<TParameter<double>>("max_cost", mnd::noop_fn<TParameter<double>>(), NAN); // filled by Proc ctor
	max_cost_f = RegisterObject<TParameter<double>>("max_cost_f", mnd::noop_fn<TParameter<double>>(), NAN); // filled by Proc ctor

	/* TFOOTHitCont cannot know ahead of time which FOOT will get put into which position,
	 * for the tracking. It's the Proc which decides that.
	 * Therefore, we allocate the `foot_param` objects here, but leave them defaulted.
	 * The TFOOTHitProc ctor will fill the rest. */

	for(u32 i=0; i<N_PAIRS; ++i) {
		foot_param[i] = RegisterObject<FOOTParams>(Form("%u_setup", i), {});
	#ifdef MND_FOOTTRACK_DEBUG
		h1_diff_q[i] = RegisterObject<TH1I>(Form("h1_diff_q_%u", i),
			Form("Innovation difference dq [e] after layer %u", i), 400,-2,2);
		h1_diff_r[i] = RegisterObject<TH1I>(Form("h1_diff_r_%u", i),
			Form("Innovation difference |dr| [mm] after layer %u", i), 1000,0,2);
		h1_diff_t[i] = RegisterObject<TH1I>(Form("h1_diff_t_%u", i),
			Form("Projected track's matching to upstream at target [mm] after layer %u", i), 500,0,100);

		h2_diff_q[i] = RegisterObject<TH2I>(Form("h2_diff_q_%u", i),
			Form("Innovation difference de [q] after layer %u vs. dQ in layer %u", i, i), 400,-2,2, 400,-2,2);
		h2_diff_r[i] = RegisterObject<TH2I>(Form("h2_diff_r_%u", i),
			Form("Innovation difference dy [mm] vs. dx [mm] after layer %u", i), 500,-2,2, 500,-2,2);
		h2_diff_t[i] = RegisterObject<TH2I>(Form("h2_diff_t_%u", i),
			Form("Innovation in projected track's matching dx vs. dy [mm] after layer %u", i), 500,-20, 20, 500,-20,20);
	#endif
	}
}

mnd::geom::Line3D RNTrackToLine3D(const RNFOOTTrack& t) noexcept {
	return *t;
}

/* Instantiate the massive template. */
template mnd::geom::VertexingResult<RNFOOTTrack>
mnd::geom::FindVertexingTracksMut(std::vector<RNFOOTTrack> &, double );

ClassImp(FOOTQ);
ClassImp(FOOTHit);
ClassImp(RNFOOTPair);
ClassImp(RNFOOTTrack);
ClassImp(RNFOOTHit);
ClassImp(RNFOOTHit::Vertex);

ClassImp(SFOOTLayerGeometry);
ClassImp(STrackingGeometry);
