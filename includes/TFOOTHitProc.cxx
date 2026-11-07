#include "TFOOTHitProc.h"
#include "TFOOTCalCont.h"
#include "TFOOTHitCont.h"

#include <algorithm>
#include <cmath>
#include <tuple>
#include <type_traits>
#include <utility>

#include "TParameter.h"

#include "util/DirectedAGraph.hxx"
#include "util/FTrack.h"
#include "util/Geometry.h"
#include "util/HitMatrix.hxx"
#include "util/FloatHacks.hxx"
#include "util/Combinatorics.hxx"

static constexpr u32 N_PAIRS = TFOOTHitCont::N_PAIRS;

extern thread_local mnd::geom::Line3D g_upstream_track; // extern'ed from `includes/TFRSHitProc.cxx`
using mnd::fmt::operator<<; // to shut up the whiney-ass LSP :X

/* These following values are assigned only once; at the TFOOTHitProc ctor.
 * Later access is read-only, as such safe to be non thread_local. */
static std::array<Eigen::Vector2d, N_PAIRS> refl; // +-1 based on the on each of the `orientation` params
static TrackingGeometry<N_PAIRS> geometry;

#define GEN_ARG_INSTANCE_FOOT(z, n, data) \
	const TFOOTCalCont& in_##n
#define GEN_ARG_NAME_FOOT(z, n, data) in_##n

/* This part below is a bit stupid, since anyway once the compiler sees the
 * field inside TFOOTHitProc, it can implicitly instantiate the templated type.
 * However, if we leave the *full* instantiation there, we will be
 * running into the Cling's good old Eigen non-digestability.
 * Hence, we extern it there, and explicitly instantiate here. Once and for all.
 * Cute little fact: when the compiler first sees `TFOOTHitProc::hm`, it will instantiate just enough
 * of the template to answer: alignof, alignas, sizeof, and all the ABI stuff
 * (vptr/vtable / dispatches). But it does not need to implicitly define all the
 * member functions/inlined statics/dtor/ctors. That part is extern'ed.
 * These actually get their definition and home right here. :-) */
template struct HitMatrix<RNFOOTPair>;
template struct Track<N_PAIRS, RNFOOTPair>;

using FHitMatrix = TFOOTHitProc::FHitMatrix;
using FTrackOnline = TFOOTHitProc::FTrackOnline;

namespace mnd {
template<typename Tuple, typename BinaryOp, std::size_t... Is>
void _for_pair_in_tuple_impl(Tuple&& t, BinaryOp&& f, std::index_sequence<Is...>) {
	(..., f (std::get<2*Is>(std::forward<Tuple>(t)),
			 std::get<2*Is+1>(std::forward<Tuple>(t))) );
}
template<typename Tuple, typename BinaryOp>
void for_pair_in_tuple(Tuple&& t, BinaryOp&& f) {
	constexpr std::size_t N = std::tuple_size_v<std::decay_t<Tuple>>;
	static_assert(N % 2 == 0, "Tuple size must be an even number.");

	_for_pair_in_tuple_impl(std::forward<Tuple>(t), std::forward<BinaryOp>(f),
		std::make_index_sequence<N/2>{});
}
};

/* Some global analysis states, which are set-up during TFOOTHitProc ctor, and then
 * everafter only read from. */
double TrackCost::Cr = TrackCost::DEFAULT_COST_R; // cost per mm^2 difference
double TrackCost::Cq = TrackCost::DEFAULT_COST_Q; // cost per charge^2 difference
double TrackCost::Ct = TrackCost::DEFAULT_COST_T; // cost per mm^2 of upstream @target difference
double TrackCost::max_cost = TrackCost::DEFAULT_MAX_CANDIDATE_COST; // max amount that a candidate point can accumulate
double TrackCost::max_cost_final_track = TrackCost::DEFAULT_MAX_FINAL_COST; // max amount that a final track can accumulate

/* Bind references for easier access... */
double& Cr = TrackCost::Cr;
double& Cq = TrackCost::Cq;
double& Ct = TrackCost::Ct;
double& max_cost = TrackCost::max_cost;
double& max_cost_final_track = TrackCost::max_cost_final_track;

Verbosity TFOOTHitProc::v = Verbosity::SILENT;
bool TFOOTHitProc::requires_valid_upstream_track;

std::ostream& operator<<(std::ostream& os, const TrackCost& rhs) {
	return os << "(kr: " << rhs.kr_ << ','
	          << " kq: " << rhs.kq_ << ','
	          << " kt: " << rhs.kt_ << ')';
}

/* Maximum cost to allow a point to branch the Kalman cannot be constant,
 * as it depends on the layer/points numbered:
 * Layer [1] => only `kq` cost
 * Layer [2] => `kq` + `kt` cost
 * Layer [3] => `kq` + `kt` + `kr` cost
 * Layer [4] => `kq` + `kt` + `kr` (3-point-fit) cost. */
bool operator>(const TrackCost& c, double max_cost) {
	// Use upstream availability even before the target cost has been evaluated.
	if(c.test_track_size >= 1 && !g_upstream_track.HasValue())
		max_cost *= 0.666667; // reduce by 2/3
	switch(c.test_track_size) {
		case 0: return c.sum() > max_cost * 0.5; // Only kq measurement possible (without track contribution)
		case 1: return c.sum() > max_cost * 1.4; // kq             + kt
		case 2: return c.sum() > max_cost * 1.2; // kq + kr        + kt(3 pts)
		case 3: return c.sum() > max_cost * 1.0; // kq + kr(3 pts) + kt(4 pts)
		default: mnd::unreachable();
	}
}
inline bool operator<(const TrackCost& c, double max_cost) noexcept { return !(c > max_cost); }

void TFOOTHitProc::ValidateParameters() {
	u32 i = 0;

	std::vector<
		std::map<Orientation, double> // orientation: z
	> test_vec(N_PAIRS);
	const FOOTBoxParam* box = out.box;
	if(!box) ERROR("`box` pointer left as null? Did you call TFOOTHitCont::Init(..) before constructing the processor?\n");
	
	mnd::for_each_in_tuple(this->in, [this, &i, &test_vec, box](const TFOOTCalCont& cfoot) {
		const auto* s = cfoot.setup;
		if(!s) ERROR("FOOT%u setup is nullptr?", i);
		int n = s->N;
		if(n >= N_FOOT_DETECTORS)
			ERROR("Found FOOT index: %d, and is out of range [0,%d> ?", n, N_FOOT_DETECTORS);

		double z = box->GetFOOTZ(n, s);
		Orientation o = s->GetOrientation();
		if(o == Orientation::UNKNOWN) ERROR("FOOT%u, orientation unknown.\n", i);
		
		const u32 ipair = i/2;
		if(ipair >= N_PAIRS)
			ERROR("Checking for input setup validity, encountered %i >= %i ?", ipair, N_PAIRS);
		auto& map = test_vec[ipair];
		map.insert({o,z});

		++i;
	});
	MND_ASSERT(i == 2*N_PAIRS && "Paranoia V1");
	
	/* Check that the parameter file has sane input. */
	if(!std::all_of(test_vec.begin(), test_vec.end(),
		[](const auto& layer) {
			return layer.count(Orientation::X) == 1
				and layer.count(Orientation::Y) == 1;
		})
	) ERROR("Validity check failed. In FOOT setup file it's marked that some pair "
		"does not contain both X- and Y- measuring detectors ('orientation' key).\n");

	if(std::adjacent_find(test_vec.begin(), test_vec.end(),
		[](const auto& previous, const auto& next) {
          return !( std::max(previous.at(Orientation::X), previous.at(Orientation::Y))
            < std::min(next.at(Orientation::X), next.at(Orientation::Y)));
		}
	) != test_vec.end())
		ERROR("Validity check failed. In FOOT setup file it's marked that "
			"the two coupled xy's aren't closer to each-other than to the following pair.\n");

	MND_ASSERT(i == 2*N_PAIRS && "Paranoia V1");
}

static bool _called_set_cv_mat = false;

void TFOOTHitProc::AssignParameters() {
	if(_called_set_cv_mat)
		WARN("Calling a heavy-load init function more than once. Is not necessary ..!\n");

	ValidateParameters();
	const FOOTBoxParam* box = out.box;
	
	using std::cos; using std::sin;

	/* Assign the centre-point z-values and the offsets for each pair. */
	u32 i = 0;
	mnd::for_pair_in_tuple(this->in, [this, &i, box](
		const TFOOTCalCont& f1, // X or Y- measuring FOOT
		const TFOOTCalCont& f2  // Y or X- measuring FOOT
		) {
			const Orientation o1 = f1.setup->GetOrientation();
			const Orientation o2 = f2.setup->GetOrientation();
			if(o1 == o2) ERROR("Pair: %u, orientations marked as identical?", i);
			
			const FOOTParam& px = (o1 == Orientation::X) ? *f1.setup : *f2.setup;
			const FOOTParam& py = (o1 == Orientation::Y) ? *f1.setup : *f2.setup;

			this->SetPairConversionMatrices(i, px, py);

			this->out.foot_param[i]->at(0) = px;
			this->out.foot_param[i]->at(1) = py;
			++i;
		}
	);

	MND_ASSERT(i == N_PAIRS && "Paranoia V2");
	SetFitConversionMatrices();

	*out.geometry = stg_from_tg(geometry);
	_called_set_cv_mat = true;
}

/* Conversion matrix from a track vector normal repr: ν = (x0, ax, y0, ay)ᵀ into FOOT pair's measurement (u,v).
 * Where u|v = R(cx - 319.5) , for X|Y- measuring FOOT.
 * Each detector is rotated by an angle θx, and placed at different positions: zu and zv.
 * A = ( cosθx ,  zu*cosθx , sinθx , zu*sinθx
 *      -sinθy , -zv*sinθy , cosθy , zv*cosθy ) .
 * Then:  (u,v)ᵀ = Aν - dxy . */
void TFOOTHitProc::SetPairConversionMatrices(
	int ipair, const FOOTParam& px, const FOOTParam& py
) {

	using std::cos; using std::sin;
	const FOOTBoxParam* box = out.box;

	const double tx = px.delta_a * M_PI / 180.0;
	const double ty = py.delta_a * M_PI / 180.0;

	const double zu = box->GetFOOTZRel(&px) + TARGET_Z;
	const double zv = box->GetFOOTZRel(&py) + TARGET_Z;
	
	geometry.layer[ipair].A << cos(tx),  zu*cos(tx), sin(tx), zu*sin(tx),
	                          -sin(ty), -zv*sin(ty), cos(ty), zv*cos(ty);

	geometry.layer[ipair].dxy << px.delta_p,
	                             py.delta_p;
	
	geometry.layer[ipair].z << zu,
	                           zv;

	refl[ipair] << px.R(),
	               py.R();
}

STrackingGeometry TFOOTHitProc::stg_from_tg(const TrackingGeometry<N_PAIRS>& tg) {
	STrackingGeometry result;
	
	for(u32 i = 0; i < N_PAIRS; ++i) {
		const auto& layer = tg.layer[i];
		auto& saved = result.layer[i];
		
		for(Eigen::Index row = 0; row < 2; ++row) {
			for(Eigen::Index col = 0; col < 4; ++col)
				saved.A[4 * row + col] = layer.A(row, col);
			saved.dxy[row] = layer.dxy(row);
			saved.z[row] = layer.z(row);
		}
	}
	return result;
}

void TFOOTHitProc::SetFitConversionMatrices() {
	constexpr static u32 N = N_PAIRS;
	auto& g_fits_lookup_table = geometry.fit;
	static_assert(g_fits_lookup_table.size() == N, "Paranoia");

	mnd::static_for<2, N_PAIRS+1>( [&](const auto K_) {
		constexpr u32 K = decltype(K_)::value;
		constexpr u32 n_combos = (u32)mnd::combi::n_choose_k(N,K);
		
		using DMatrix = DesignMatrix<N,K>;
		DMatrix& dm = g_fits_lookup_table.get<K>();

		static_assert(n_combos == DMatrix::n_combinations);
		
		for(u32 icombo = 0; icombo < n_combos; ++icombo) {
			typename DMatrix::Matrix underlying_matrix {};
			typename DMatrix::Vector offset_vector {};
			
			const u32 mask = mnd::combi::combo_lookup_table<N, K>[icombo];
			Eigen::Index row = 0;
			for(u32 i=0; i<N; ++i) {
				if(mask & (1U << i)) {
					underlying_matrix.template middleRows<2>(row) =
						geometry.layer[i].A;
					offset_vector.template middleRows<2>(row) =
						geometry.layer[i].dxy;
					row += 2;
				}
			}
			if(row != 2*K)
				ERROR("N=%u, K=%u, row=%u => check row == 2*K failed, e.g. the "
					"underlying_matrix for QR is entirely wrong.\n", N, K, (u32)row);

			auto& qr = dm.lookup_table[icombo].first;
			qr.compute(underlying_matrix);
			
			dm.lookup_table[icombo].second = offset_vector;

			if(qr.rank() != 4)
				ERROR("Rank-deficient geometry for mask %#08x. QR decomposition failed. "
					"Are the parameters fine or degenerate for this case.. ?", mask);
		}
	});
}

TFOOTHitProc::TFOOTHitProc (
	TFOOTHitCont& out,
	BOOST_PP_ENUM(N_FOOT_DETECTORS, GEN_ARG_INSTANCE_FOOT, ~),
	double max_cost_,
	double max_cost_final_track_,
	const std::array<double, 3>& cost,
	bool req_,
	Verbosity v_
) : TFOOTHitProc::Base (
		out,
		BOOST_PP_ENUM(N_FOOT_DETECTORS, GEN_ARG_NAME_FOOT, ~)
	)
{
	TFOOTHitProc::v = v_;
	TFOOTHitProc::requires_valid_upstream_track = req_;
	TrackCost::max_cost = max_cost_;
	TrackCost::max_cost_final_track = max_cost_final_track_;

#ifndef MND_FOOTTRACK_DEBUG
	if(v > 0)
		WARN("Asking for a verbose analysis, but `MND_FOOTTRACK_DEBUG` option not compiled in.. Recompile it please.");
#endif

	/* Assign the Kalman cost coefficients. */
	*out.cost_coeff = {
		TrackCost::DEFAULT_COST_R,
		TrackCost::DEFAULT_COST_Q,
		TrackCost::DEFAULT_COST_T
	};
	out.max_cost->SetVal(max_cost);
	out.max_cost_f->SetVal(max_cost_final_track);

	for(int i=0; i<3; ++i) {
		double val = cost[i];
		if( !std::isfinite(val) or val < 0 ) continue;
		switch(i) {
			case(0): ::Cr = val; break;
			case(1): ::Cq = val; break;
			case(2): ::Ct = val; break;
		}
		out.cost_coeff->at(i) = val;
	}
	
	WARN("Cost coefficients [cr, cq, ct] = [%.2f, %.2f, %.2f]\n", Cr, Cq, Ct);
	WARN("Max candidate allowed cost: %s%.2f%s , max complete track allowed cost %s%.2f%s\n",
		BOLD, max_cost, KNRM, BOLD, max_cost_final_track, KRNM);
	WARN("Cluster size==1 cutoff point is %.2f units of charge. Anything below that considered "
		"pure noise and ignored\n", CLUSTER_SIZE_ONE_Q_CUTOFF);

	AssignParameters();

	if(v > 2) {
		WARN("FOOT pairs placed at: "); for(u32 n=0; n<N_PAIRS; ++n) {
			std::cerr << ", " << geometry.layer[n].z << std::endl;
		}
		WARN("EXPERT target placed at: %.2f\n", TFOOTHitProc::TARGET_Z);
	}
} // TFOOTHitProc::TFOOTHitProc( ... )

void TFOOTHitProc::ProcessEntry() noexcept {
	out.Clean();
	if(requires_valid_upstream_track && !g_upstream_track.HasValue()) return;

#ifdef MND_FOOTTRACK_DEBUG
	if(v > 1) {
		static u64 ev_num = 0;
		fprintf(stderr, "\n%s>>> Entry[%lu] <<<%s\n", KBH_GRN, ++ev_num, KNRM);
	}
#endif

	PreProcess();
	
	int ipair = 0;
	mnd::for_pair_in_tuple(this->in, [this, &ipair](const auto& f1, const auto& f2) {
		/* Pepega code, but can't do much ... */
		using PairRef = std::pair<const TFOOTCalCont&, const TFOOTCalCont&>;
		
		Orientation e1 = f1.setup->GetOrientation();
		switch(e1) {
			case Orientation::X:
				this->ProcessPair(PairRef{f1,f2}, ipair);
				break;
			case Orientation::Y:
				this->ProcessPair(PairRef{f2,f1}, ipair);
				break;
			default:
				mnd::unreachable();
		}
		++ipair;
	});

	mnd::static_for<0, N_PAIRS>([this](auto I) {
		constexpr size_t i = decltype(I)::value;
		this->hm[i].InitEvent( this->out.inner().pair[i] );
	});

	ConstructDAG();
	AnalyseDAG();
	PostProcess();
}

/* First std::pair member is `x`, second is `y` */
void TFOOTHitProc::ProcessPair (
	const std::pair<const TFOOTCalCont&, const TFOOTCalCont&>& f, i32 ipair
) noexcept {
	const FOOTParam &px = *f.first.setup, &py = *f.second.setup;
	const RNFOOTCal &fx = f.first.inner(), &fy = f.second.inner();

	RNFOOTPair& output = out.inner().pair[ipair];

	for(const RNFOOTCluster& hit : fx.fCl) {
		auto [cx, ce, mult, ctype] = hit;

		float q = static_cast<float>( px.Q(hit) );
		
		const double u = refl[ipair].x() * (cx - FOOTParam::DETECTOR_MIDPOINT) * FOOTParam::STRIP_TO_MM;
		// Cluster size 1 fucks with everything above Z >~ 1,
		// so only care about it if its sitting at low energies..
		if(mult > 1 or q < CLUSTER_SIZE_ONE_Q_CUTOFF)
			output.x.emplace_back(q, mult
			#ifdef MND_FOOTTRACK_DEBUG
				, ctype, hit.Delta(), ce
			#endif
				, u);
	}

	for(const RNFOOTCluster& hit : fy.fCl) {
		auto [cy, ce, mult, ctype] = hit;

		float q = static_cast<float>( py.Q(hit) );

		const double v = refl[ipair].y() * (cy - FOOTParam::DETECTOR_MIDPOINT) * FOOTParam::STRIP_TO_MM;
		if(mult > 1 or q < CLUSTER_SIZE_ONE_Q_CUTOFF)
			output.y.emplace_back(q, mult
			#ifdef MND_FOOTTRACK_DEBUG
				, ctype, hit.Delta(), ce
			#endif
				, v);
	}

	/* Sort these vectors, in ascending values of average charge (Q) */
	thread_local const auto comparator = [](const auto& lhs, const auto& rhs) { return lhs.Q.q < rhs.Q.q; };
	std::sort(output.x.begin(), output.x.end(), comparator);
	std::sort(output.y.begin(), output.y.end(), comparator);
	
	output.zx = geometry.layer[ipair].z.x();
	output.zy = geometry.layer[ipair].z.y();

	constexpr u32 B = RNFOOTPair::_DATA_PRESENT_BITINDEX;
	fx.IsEmptyEvent() ? mnd::set_bit<B,0>(output.zx) : mnd::set_bit<B,1>(output.zx);
	fy.IsEmptyEvent() ? mnd::set_bit<B,0>(output.zy) : mnd::set_bit<B,1>(output.zy);
}

constexpr auto X = FHitMatrix::X;
constexpr auto Y = FHitMatrix::Y;
using Entry = FHitMatrix::Entry;

double TFOOTHitProc::kr(const FLayerInfo::PredictionP& predicted, const FHitMatrix::Entry& candidate) const noexcept {
	const Eigen::Vector2d& measured = candidate.pos; // (u,v);

	if(! predicted.array().isFinite().all() ) // power of expression templates! :)
		return TrackCost::NIL_VALUE;
	return Cr * ((predicted - measured).squaredNorm());
}

/* kQ = Cq || Qij - Qn ||^2
 * ==== Cq( [ mean(Qij) - mean(Qtrack) ]^2 + variance_ij ) 
 * Can *never* return a nil value. */
double TFOOTHitProc::kq(const FLayerInfo::PredictionQ& predicted, const FHitMatrix::Entry& candidate) const noexcept {

	// Cost starts at the difference for the layers x,y measurements
	double cost = candidate.q.var();
	if(std::isfinite(predicted)) { // Track is non-null
		double diff = predicted - candidate.q.mean();
		cost += diff*diff;
	}
	cost *= Cq;
	return cost;
}

double TFOOTHitProc::kt(FTrackOnline& ft, const FHitMatrix::Entry& candidate, u32 layer) const noexcept {
	/* If event has no upstream track or
	 * if currently 0 measurements are available in the online track, can't do anything. */
	if( !upstream_hit_loc.eigen_view().allFinite() || ft.N() == 0)
		return TrackCost::NIL_VALUE;

	FTrack current_track = ft.get(geometry);
	
	ft.Add(candidate, layer);
	const auto target_point = ft.get(geometry).l.Eval(TARGET_Z); // This next call cannot throw, so pop_back() always gets executed.

#ifdef MND_FOOTTRACK_DEBUG
	if(v > 3) {
		fprintf(stderr, "TFOOTHitProc::kt: Updated test track: ");
		std::cerr << ft << " :: value at target: " << target_point << std::endl;
	}
#endif
	ft.pop_back( std::move(current_track) );

	if(!target_point.eigen_view().allFinite())
		return TrackCost::NIL_VALUE;

	return Ct * target_point.Distance2(upstream_hit_loc);
}

/* Once identified a valid track, remove the corresponding col/row from the
 * individual hit matrices. */
void TFOOTHitProc::PoisonEntriesFromHMs(const DAG::DAGPath& path) noexcept {
	static_assert(DAG::depth == std::tuple_size_v<decltype(path.node)>, "Paranoia V3");
	mnd::static_for<0, DAG::depth>([this, &path](auto i_) {
		constexpr size_t i = decltype(i_)::value;
		const DAG::Index& index = path.node[i];
		if(!index) return; // operator bool; node is nil

		hm[i].poison(index[0], index[1]);
	});
}

using CostIndex =  TFOOTHitProc::CostIndex;
static constexpr auto sort_by_cost = [](const CostIndex& lhs, const CostIndex& rhs) {
	return lhs.first < rhs.first;
};
static_assert(std::tuple_size<TFOOTHitProc::CandidatesBuffer>::value >= 2*TFOOTHitProc::MAX_CANDIDATES);

void TFOOTHitProc::ConstructDAG() noexcept {
	/* Idea is explained in the PhD writeup.
	 * If you don't have it, ask Klayze. */

	/* To keep the algorithm invariant between layers, tracks can also be nullable (e.g. containing 0 or 1 point). */	
	dag.Initialize();

	for(u32 n = 0; n < N_PAIRS; ++n) { /* Don't unroll this. Hurts L1I locality. */
		const FHitMatrix& h = hm[n];

		const size_t nx = h.GetN<X>();
		const size_t ny = h.GetN<Y>();
		if(nx == 0 or ny == 0) continue;

		std::vector<DAG::DAGPath>& current_paths = dag.path;
		std::vector<DAG::DAGPath>  new_paths;

		/* If every single combination, for every single existing path is good,
		 * could have at most 1 path branching out to nx*ny paths extra. */
		size_t graph_branch_factor = std::min(nx*ny, (size_t)MAX_CANDIDATES);
		new_paths.reserve( graph_branch_factor * current_paths.size() ); // could be larger.

	#ifdef MND_FOOTTRACK_DEBUG
		if(v > 0)
			WARN("[#%u] ConstructDAG(): new_paths to be reserved at size: %zu (size currently:%zu)\n", n, new_paths.capacity(), new_paths.size());
	#endif
		/* Each path already draws (an optional) preliminary track.
		 * Try to match some of the current layer's hitmatrix elements against that
		 * specific preliminary test track. */
		for(const DAG::DAGPath& path : current_paths) {
			
			// Fetch the preliminary track that the path describes.
			FTrackOnline tau = this->GetPrelimTrackFromPath(path);
			const FTrack& ftrack = tau.get(geometry);
			
			// Prediction (u,v) at the following layer
			FLayerInfo::Prediction prediction = geometry.layer[n].prediction(ftrack);
			const auto& prediction_q = prediction.first;
			const auto& prediction_p = prediction.second;
			u32 ncounted = 0;
			
			// FHitMatrix::Cached is column-major (Eigen::Matrix convention).
			for(size_t j=0; j<ny; ++j) {
				for(size_t i=0; i<nx; ++i) {
					const mnd::hm::Data& e = h(i,j);
			
			        TrackCost cost(tau.N());
					cost.set<TrackCost::KQ>( kq(prediction_q, e) );
					cost.set<TrackCost::KR>( kr(prediction_p, e) );
					if(cost > max_cost) continue;
					/* Even though the other cost(s) aren't yet evaluated, `kr` is the stictest.
					 * If that one already fails here, no need evaluate further costs. */

					cost.set<TrackCost::KT>( kt(tau, e, n) );
					if(cost > max_cost) continue;

					// If the flow survives til this point means that
					// the candidate is decent. Add it to the `path_specific_candidates_buf` buffer.
					path_specific_candidates_buf[ncounted++] = { cost.sum(), DAG::Index(i,j) };
					
					// If the buffer is full, sort it and hold on to only the first `MAX_CANDIDATES`
					if( ncounted == (u32)path_specific_candidates_buf.size() ) {
						std::sort(
							path_specific_candidates_buf.begin(),
							path_specific_candidates_buf.begin() + ncounted,
							sort_by_cost
						);
						ncounted = MAX_CANDIDATES;
					}
				}
			}
			
			/* Sort the candidates buffer. */
			std::sort(
				path_specific_candidates_buf.begin(),
				path_specific_candidates_buf.begin() + ncounted,
				sort_by_cost
			);

			/* Take at most the best `MAX_CANDIDATES` from `path_specific_candidates_buf`
             * as branching points. */
			ncounted = std::min( ncounted, MAX_CANDIDATES );
			for(u32 i=0; i<ncounted; ++i) {
				new_paths.emplace_back(path); // copy-ctor.
				new_paths.back().node[n] = path_specific_candidates_buf[i].second;
			}
			/* Add a null node only if current candidate path has no pre-existing "layer holes" */
			if(path.Rank() >= static_cast<int>(n))
				new_paths.emplace_back(path);

		} // loop over existing paths
		
		dag.path = std::move( new_paths );
	} // for(u32 n = 0; n < N_PAIRS; ++n)
} // void TFOOTHitProc::ConstructDAG()

void TFOOTHitProc::AnalyseDAG() noexcept {
	static_assert(N_PAIRS - 1 >= 3, "Maximal paranoia");

	/* Allow only one layer to be missing. Must have 3 or more
	 * valid measurements and be below total allowed cost. */
	dag.KeepIf( [this](const DAG::DAGPath& p) -> bool {
		return (
			(p.Rank() >= (int)N_PAIRS - 1) and
			this->GetPrelimTrackFromPath(p).GetScore(geometry) < max_cost_final_track
		);
	});	

	#ifdef MND_FOOTTRACK_DEBUG
		if(v > 0) {
			WARN(KBH_CYN "~~~> Valid %zu paths thru the tree.\n" KNRM, dag.path.size());
		}
	#endif
	/* Sort the sequence of paths through the graph according to their score.
	 * Each paths's score anyway has to be evaluated to be robust. There's no way around this.
	 * Heavy lifting call. */
	std::sort (
		dag.path.begin(),
		dag.path.end(),
		[this](const DAG::DAGPath& pl, const DAG::DAGPath& pr) -> bool {
			/* OK this requires explaination. Basically we 'weigh' according to reduced chi2 from FTrack::GetScore(),
			 * but a path making a 4-point track will also very very likely make a 3x 3-point tracks. Now the score 
			 * R3 = 2*R4, but assuming the worst measurement, in like 2/3 of the cases, a 3-point track overtakes the 
			 * main 4-point track and poisons the hitmatrix entry, thus forever pruning it. */
			if(pl.Rank() != pr.Rank())
				return pl.Rank() > pr.Rank();
			
			const FTrackOnline tl = this->GetPrelimTrackFromPath(pl);
			const FTrackOnline tr = this->GetPrelimTrackFromPath(pr);
			return tl.GetScore(geometry) < tr.GetScore(geometry);
		}
	);

	/* Go over the score table, pluck the entries in order, and poison the occupied hitmatrix row/cols.
	 * If an entry's matched indices aren't poisoned, they are encoded as valid tracks. */
	for(const DAG::DAGPath& path : dag.path) {
		bool is_valid = true;
		for(u32 i=0; i < N_PAIRS; ++i) {
			const DAG::Index& index = path.node[i];
			if(!index) continue; // is fine.
			
			FHitMatrix& hm = this->hm[i];
			if( hm.is_poisoned(index[0] , index[1]) ) {
				/* In this case, a better candidate yoinked this element. Nothing I can do. 🤷 */
				is_valid = false;
			}
		}
		if(!is_valid) continue;

		/* Export the track. */
		FTrackOnline tau = this->GetPrelimTrackFromPath(path);
		const size_t N = tau.N();
		const double score = tau.GetScore(geometry);
		
		const FTrack& t = tau.get(geometry); // Evaluate the actual track fit.

	#ifdef MND_FOOTTRACK_DEBUG
		auto us      = mnd::make_filled_array<double, N_PAIRS>(NAN);
		auto vs      = mnd::make_filled_array<double, N_PAIRS>(NAN);
		auto zxs     = mnd::make_filled_array<double, N_PAIRS>(NAN);
		auto zys     = mnd::make_filled_array<double, N_PAIRS>(NAN);
		auto qm      = mnd::make_filled_array<double, N_PAIRS>(NAN);
		auto sqm     = mnd::make_filled_array<double, N_PAIRS>(NAN);
		auto e0_x    = mnd::make_filled_array<double, N_PAIRS>(NAN);
		auto delta_x = mnd::make_filled_array<double, N_PAIRS>(NAN);
		auto c0_x    = mnd::make_filled_array<u32,    N_PAIRS>(-1);
		auto e0_y    = mnd::make_filled_array<double, N_PAIRS>(NAN);
		auto delta_y = mnd::make_filled_array<double, N_PAIRS>(NAN);
		auto c0_y    = mnd::make_filled_array<u32,    N_PAIRS>(-1);

		for(size_t layer=0, valid=0; layer < N_PAIRS; ++layer) {
			if( path.node[layer] ) {
				us[layer] = tau.m[valid].x();
				vs[layer] = tau.m[valid].y();
				zxs[layer] = geometry.layer[layer].z.x();
				zys[layer] = geometry.layer[layer].z.y();
				qm[layer] = tau.q[valid].mean();
				sqm[layer] = tau.q[valid].s();
				
				// Is asserted to return a value based on the `if` block, but CBA changing the API.
				mnd::Maybe<FHitMatrix::RawHitPairRef> p = GetHitFromPath(layer, path);
				if(!p) ERROR("huh?");
				FHitMatrix::RawHitPairRef& val = p.value();
				e0_x[layer]    = val.x->Q.e0;
				delta_x[layer] = val.x->Q.delta;
				c0_x[layer] = static_cast<u32>(refl[layer].x() * val.x->m * FOOTParam::MM_TO_STRIP + FOOTParam::DETECTOR_MIDPOINT);
				e0_y[layer]    = val.y->Q.e0;
				delta_y[layer] = val.y->Q.delta;
				c0_y[layer] = static_cast<u32>(refl[layer].y() * val.y->m * FOOTParam::MM_TO_STRIP + FOOTParam::DETECTOR_MIDPOINT);

				++valid;
			}
		}
	#endif // MND_FOOTTRACK_DEBUG

		out.inner().track.emplace_back (
			t.l.xarray(), t.l.yarray(), t.q.mean(), score, N
		#ifdef MND_FOOTTRACK_DEBUG
				,
				us, vs, zxs, zys,
				qm, sqm,
				e0_x, delta_x, c0_x,
				e0_y, delta_y, c0_y
		#endif
		);
		
		/* Poison the exported hitmatrix elements along their rows and columns.
		 * This way, no other track can reuse either the same `x` or the same `y` entry. */
		PoisonEntriesFromHMs(path);
	
	} // for(const DAG::DAGPath& path : dag.path)

#ifdef MND_FOOTTRACK_DEBUG
	/* At the end, replay the best scoring path, if it contains N_PAIRS measurements. */
	if(!dag.path.empty() and dag.path.front().Rank() == N_PAIRS)
		ReplayCostCalculation(dag.path.front());
#endif

} // void TFOOTHitProc::AnalyseDAG()

/* High frequency call, decoding a DAGPath sequence into an online track. */
FTrackOnline TFOOTHitProc::GetPrelimTrackFromPath(const DAG::DAGPath& p) const noexcept {
	static_assert(DAG::depth == N_PAIRS, "Just in case. Must pass");
	FTrackOnline t{};
	
	mnd::static_for<0, N_PAIRS>([this, &p, &t](auto layer_) {
		constexpr u32 layer = static_cast<u32>( decltype(layer_)::value ); // layer index
		const DAG::Index& i = p.node[layer];
		const FHitMatrix& h = hm[layer];
		if(i) { // explicit operator bool()
			t.Add( h(i[0], i[1]) , layer );
		}
	});

	return t;
}

#ifdef MND_FOOTTRACK_DEBUG
mnd::Maybe<FHitMatrix::RawHitPairRef> TFOOTHitProc::GetHitFromPath(size_t layer, const DAG::DAGPath& p) const {
	const DAG::Index& i = p.node.at(layer);
	const FHitMatrix& h = hm[layer];

	if(i) { // explicit operator bool()
		u16 xi = i[0];
		u16 yi = i[1];
		return h.GetRaw(xi, yi);
	} else {
		return {};
	}
}
#endif

void TFOOTHitProc::PreProcess() noexcept {
	/* Upstream track is in FRS coordinates, and must be represented in FOOT array coordinates. */
	g_upstream_track %= ( out.box->GetTargetZ() + TFOOTHitProc::TARGET_Z );

	/* Can be nil. */
	upstream_hit_loc = g_upstream_track.Eval( TFOOTHitProc::TARGET_Z );
};

void TFOOTHitProc::PostProcess() noexcept {
	/* For recognised tracks, try to find their vertex, together with the upstream track. */

	std::sort( /* Sort in descending charge (.Q) attribute. */
		out.inner().track.begin(),
		out.inner().track.end()
	);
	
	this->lines.clear();
	for(const auto& t : out.inner().track) {
		lines.emplace_back( RNTrackToLine3D(t) );
	}
	lines.push_back( g_upstream_track ); // copy-ctor

	mnd::geom::Point3D vertex = FindVertex( mnd::as_span(lines) );
	out.inner().vertex = RNFOOTHit::Vertex{ vertex };

	if(lines.size() >= 2) {	
		out.diff_heavy_frag_vs_upstream->Fill ( /* Hardly converged, fucked up due to rounding. */
			lines.front().DistanceTo( g_upstream_track )
		);
	}
};

#ifdef MND_FOOTTRACK_DEBUG
void TFOOTHitProc::ReplayCostCalculation(const DAG::DAGPath& path) noexcept {
	DAG::DAGPath p{};
	for(u32 n = 0; n < N_PAIRS; ++n) {
		const DAG::Index& i = path.node[n];
		if(!i) continue; // possible that it is null?
		
		const mnd::hm::Data& e = hm[n]( i[0], i[1] );
		FTrackOnline tau = this->GetPrelimTrackFromPath(p);
		const FTrack& ftrack = tau.get(geometry);
	
		// Prediction (u,v) at the following layer
		FLayerInfo::Prediction prediction = geometry.layer[n].prediction(ftrack);
		const auto& prediction_q = prediction.first;
		const auto& prediction_p = prediction.second;

		/* Calculate cost of the current path `p`, if I add the node `index` from i'th layer. */
		const double kq = this->kq(prediction_q, e);
		const double kr = this->kr(prediction_p, e);
		const double kt = this->kt(tau, e, n);

		if(std::isfinite(kq)) out.h1_diff_q[n]->Fill( kq / Cq );
		if(std::isfinite(kr)) out.h1_diff_r[n]->Fill( kr / Cr );
		if(std::isfinite(kt)) out.h1_diff_t[n]->Fill( kt / Ct );

		// Add the node to the current 'fake' online track.
		p.node[n] = path.node[n];
	}
}
#endif // MND_FOOTTRACK_DEBUG
