#include "TFOOTHitProc.h"

#include "TFOOTCalCont.h"
#include "TFOOTHitCont.h"
#include "TFOOTMapCont.h"
#include "util/DirectedAGraph.hxx"
#include "util/Geometry.h"
#include <algorithm>
#include <cmath>
#include <tuple>
#include <utility>
#include "TParameter.h"
#include "util/HitMatrix.hxx"

extern thread_local mnd::geom::Line3D g_upstream_track; // extern'ed from `includes/TFRSHitProc.cxx`

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
template struct Track<TFOOTHitCont::N_PAIRS, RNFOOTPair>;

using mnd::Maybe;
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
	// Sometimes `kt` measurement missing because the upstream has no track. Adjust the `max_cost` in this case.
	if(c.test_track_size >= 1 and !c.is_set<TrackCost::KT>())
		max_cost *= 0.666667; // reduce by 2/3
	switch(c.test_track_size) {
		case 0: return c.sum() > max_cost * 0.5; // Only kq measurement possible (without track contribution)
		case 1: return c.sum() > max_cost * 1.4; // kq             + kt       
		case 2: return c.sum() > max_cost * 1.2; // kq + kr        + kt(3 pts)
		case 3: return c.sum() > max_cost * 1.0; // kq + kr(3 pts) + kt(4 pts)
		default: __builtin_unreachable();
	}
}
inline bool operator<(const TrackCost& c, double max_cost) noexcept { return !(c > max_cost); } 

/* Read the param file and create the matrix `A`:
 *  1/cos(θx-θy) * (  cos(θx)  sin(θx) )
 *                 ( -sin(θy)  cos(θy) ) 
 *  and the offset vector `d_xy`. */
void TFOOTHitProc::SetConversionMatrices(int ipair, const FOOTParam& px, const FOOTParam& py ) {
	using std::cos; using std::sin;
	double tx = px.delta_a * M_PI / 180.0, ty = py.delta_a * M_PI / 180.0;	
	hm[ipair].A << cos(ty), -sin(tx),
	               sin(ty),  cos(tx);
	hm[ipair].A *= 1.0/cos(tx-ty);

	hm[ipair].dxy << px.delta_p, // already in [mm] scale, don't need to convert.
	                 py.delta_p;
	
	refl[ipair] << px.R(), 
	               py.R();
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

	u32 i = 0;

	std::vector< 
		std::map<Orientation, double> // orientation vs. z
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
		if(ipair > N_PAIRS) 
			ERROR("Checking for input setup validity, encountered %i>%i ?", ipair, N_PAIRS);	
		auto& map = test_vec[ipair];
		map.insert({o,z});

		++i;
	});
	MND_ASSERT(i == 2*N_PAIRS && "Paranoia V1");

	/* Check that the parameter file has sane input. */
	bool is_fine = std::is_sorted(test_vec.begin(), test_vec.end(),
		[](const std::map<Orientation,double>& lhs, const std::map<Orientation,double>& rhs) {
			return (
				lhs.find(Orientation::X) != lhs.end() and // x orientation entry found.
				lhs.find(Orientation::Y) != lhs.end() and // y orientation entry found.
				rhs.find(Orientation::X) != rhs.end() and // x orientation entry found.
				rhs.find(Orientation::Y) != rhs.end() and // y orientation entry found.
				lhs.at(Orientation::X) < rhs.at(Orientation::X) and // x's come in a sequence
				lhs.at(Orientation::Y) < rhs.at(Orientation::Y) and // y's come in a sequence
				( std::abs(lhs.at(Orientation::X) - lhs.at(Orientation::Y)) // y from following won't get mixed
				  < std::abs(lhs.at(Orientation::X) - rhs.at(Orientation::Y))) and
				( std::abs(lhs.at(Orientation::X) - lhs.at(Orientation::Y)) // x from following won't get mixed
				  < std::abs(rhs.at(Orientation::X) - lhs.at(Orientation::Y)))
			);
		}
	);
	if(!is_fine)
		ERROR("Validity check failed. In FOOT setup file it's either not marked:\n"
			"(1): sequence xy xy xy (with `-` optional modifier)\n"
			"(2): the two coupled xy's aren't closer to eachother than the following pair.\n");

	/* Assign the average z-values of each of the pairs. */
	i = 0;
	mnd::for_pair_in_tuple(this->in, [this, &i, box](const TFOOTCalCont& f1, const TFOOTCalCont& f2) {
			this->pair_z[i] = ( 
				box->GetFOOTZRel(f1.setup->N) + f1.setup->dz + 
				box->GetFOOTZRel(f2.setup->N) + f2.setup->dz
			) / 2.0 + TFOOTHitProc::TARGET_Z;
			Orientation o1 = (f1.setup->orientation[1] == 'x') ? Orientation::X : Orientation::Y;
			if(o1 == Orientation::X) { 
				this->SetConversionMatrices(i, *f1.setup, *f2.setup);

				this->out.foot_param[i]->at(0) = *f1.setup;
				this->out.foot_param[i]->at(1) = *f2.setup;
			} else { 
				this->SetConversionMatrices(i, *f2.setup, *f1.setup);

				this->out.foot_param[i]->at(0) = *f2.setup;
				this->out.foot_param[i]->at(1) = *f1.setup;
			}
			++i;
		}
	);
	MND_ASSERT(i == N_PAIRS && "Paranoia V2");
	const ExpertTarget& target = box->target;
	
	/* Assign S2 Be target parameters (EXPERT target). */
	this->target_xy = mnd::geom::Rectangle2D (
		{target.dx, target.dy}, target.width_x, target.width_y
	);

	if(v > 2) {
		WARN("FOOT pairs placed at: "); std::cerr << this->pair_z << std::endl;
		WARN("EXPERT target placed at: %.2f\n", TFOOTHitProc::TARGET_Z);
		WARN("EXPERT target dimensions: "); std::cerr << target_xy << std::endl;
	}
}

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
				__builtin_unreachable();
		}
		++ipair;
	});

	mnd::static_for<0, N_PAIRS>([this](auto I) {
		constexpr size_t i = decltype(I)::value;
		this->hm[i].InitEvent( this->out.inner().pair[i] );
	});

#ifdef MND_FOOTTRACK_DEBUG
	ConstructObviousTracks();
#else
	ConstructDAG();
	AnalyseDAG();
#endif
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
		
		double xprime = refl[ipair].x() * (cx - FOOTParam::DETECTOR_MIDPOINT) * FOOTParam::STRIP_TO_MM; 
		// Cluster size 1 fucks with everything above Z >~ 1,
		// so only care about it if its sitting at low energies.. 
		if(mult > 1 or q < CLUSTER_SIZE_ONE_Q_CUTOFF)
			output.x.emplace_back(q, mult 
#ifdef MND_FOOTTRACK_DEBUG
				, ctype, ce, hit.Delta() 
#endif
				, xprime);
	}

	for(const RNFOOTCluster& hit : fy.fCl) {
		auto [cy, ce, mult, ctype] = hit;

		float q = static_cast<float>( py.Q(hit) );

		double yprime = refl[ipair].y() * (cy - FOOTParam::DETECTOR_MIDPOINT) * FOOTParam::STRIP_TO_MM; 
		if(mult > 1 or q < CLUSTER_SIZE_ONE_Q_CUTOFF)
			output.y.emplace_back(q, mult 
#ifdef MND_FOOTTRACK_DEBUG
				, ctype, ce, hit.Delta() 
#endif
				, yprime);
	}

	/* Sort these vectors, in ascending values of average charge (Q) */
	thread_local const auto comparator = [](const auto& lhs, const auto& rhs) { return lhs.Q.q < rhs.Q.q; };
	std::sort(output.x.begin(), output.x.end(), comparator);
	std::sort(output.y.begin(), output.y.end(), comparator);
	
	output.z = pair_z[ipair];
}

constexpr auto X = FHitMatrix::X; 
constexpr auto Y = FHitMatrix::Y;
using Entry = FHitMatrix::Entry;

double TFOOTHitProc::kr(const FTrackOnline& ft, const FHitMatrix::Entry& candidate, u32 k) const noexcept {
	const Eigen::Vector2d& measured = candidate.v;
	mnd::geom::Point2D extrapolated = ft.extrapolate_to( pair_z[k] );

	if(! extrapolated.eigen_view().array().isFinite().all() ) // power of expression templates! :)
		return TrackCost::NIL_VALUE;
	else
		return Cr * ((extrapolated.eigen_view() - measured).squaredNorm());
}

/* kQ = Cq || Qij - Qn ||^2 
 * ==== Cq( [ mean(Qij) - mean(Qtrack) ]^2 + variance_ij )  
 * Can *never* return a nil value. */
double TFOOTHitProc::kq(const FTrackOnline& ft, const FHitMatrix::Entry& candidate, u32 k) const noexcept {
	(void)k;

	double mean_track_q = ft.q.mean();
	double candidate_var = candidate.q.var();

	double cost = candidate_var;
	if(std::isfinite(mean_track_q)) { // Track is non-null
		double diff = mean_track_q - candidate.q.mean();
		cost += diff*diff;
	}
	cost *= Cq;
	return cost;
}

double TFOOTHitProc::kt(const FTrackOnline& ft, const FHitMatrix::Entry& candidate, u32 k) const noexcept {
	FTrackOnline& mft = const_cast<FTrackOnline&>(ft);
	mft.Add(candidate, pair_z[k]);
	const FTrack& tt = mft.get();

	double cost_t = TrackCost::NIL_VALUE;

	/* In case the current online track has only 1 point, then we can't squeeze out a value. */
	if( tt.l.HasValue() ) {
		mnd::geom::Point2D extrapolate_to_target = mft.extrapolate_to( TARGET_Z );

#ifdef MND_FOOTTRACK_DEBUG
		if(v > 3) {
			fprintf(stderr, "TFOOTHitProc::kt: Test track: ");
			std::cerr << mft << " :: value at target: " << extrapolate_to_target << std::endl;
		}
#endif
		cost_t = Ct * extrapolate_to_target.Distance2( upstream_hit_loc );
		
	}

	mft.pop_back();
	return cost_t;
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

/* Some heavy particle tracks can be obviously taken out of the full glory
 * CKF algorithm. E.g. if each FOOT in-order measures Z=4,5,6 is fine. 
 * NOTE: the only obvious particle we can find at this point, is the largest Q one. 
 * Isn't very efficient, so its called only in debug mode. */
void TFOOTHitProc::ConstructObviousTracks() noexcept {
	static_assert(DAG::depth == N_PAIRS, "Paranoia V4");

	/* Do the same CKF but with reduced phase space: meaning that 
	 * we only take last row/column as candidates.
	 * Taking only the last entry (largest Qx,Qy) is a bit lazy, since the hit could go through
	 * noisy/dead strip, and the corresponding entry won't be the largest Z. But anyway, this should
	 * get caught by the full CKF. */

	/* Only one path is viable here, no branching possible. Always take the best candidate... */
	DAG::DAGPath path{};
	
	for(u32 n = 0; n < N_PAIRS; ++n) {
		const FHitMatrix& h = hm[n];
		const size_t nx = h.GetN<X>();
		const size_t ny = h.GetN<Y>();
		const size_t n_last_row_col = nx+ny-1;

		if(nx == 0 or ny == 0) continue; 

		// Fetch the preliminary track that the path describes.
		FTrackOnline tau = this->GetPrelimTrackFromPath(path);
		TrackCost cost(n);

		double cost_min_current = INFINITY;
		DAG::Index best_i;

#ifdef MND_FOOTTRACK_DEBUG
		if(v > 1) {
			WARN(KBH_BLU "Track entering layer: #%u: " KNRM, n);
			std::cerr << tau << std::endl;
		}
		if(v > 2) std::cerr << h << std::endl;
#endif
		                      /* nx + ny - 1 */
		for(size_t k = 0; k < n_last_row_col; ++k) {
			size_t i, j;
			if(k < ny) {
				i = nx - 1; // last row
				j = k;
			} else {
				j = ny - 1; // last col
				i = nx+ny - k - 2;
			}
			// FHitMatrix::Cached is column-major (Eigen convention).
			const mnd::hm::Data& e = h(i,j);
			if(e.q.mean() < 2.5) continue; 
			
			cost.reset();
			const double kq = this->kq(tau, e, n);
			const double kr = this->kr(tau, e, n);
			const double kt = this->kt(tau, e, n);

			cost.set<TrackCost::KQ>(kq);
			cost.set<TrackCost::KR>(kr);
			cost.set<TrackCost::KT>(kt);

#ifdef MND_FOOTTRACK_DEBUG
			if(v > 3) std::cerr << DAG::Index{i,j} << " : " << e << " :: cost: " << cost << std::endl;
#endif
			if(cost.sum() < cost_min_current) {
				cost_min_current = cost.sum();
				best_i = {i,j};
			}

			out.h1_diff_q[n]->Fill( kq / Cq ); 
			out.h1_diff_r[n]->Fill( kr / Cr );	 
			out.h1_diff_t[n]->Fill( kt / Ct );
		}
		// After processing the layer, check if total cost acquired 
		// for this specific {i,j} entry of layer `n` is within the limits of `max_cost`.
		if(best_i and cost < max_cost) {
			path.node[n] = best_i;

#ifdef MND_FOOTTRACK_DEBUG
			if(v > 2) std::cerr << "Found best index: " << best_i << std::endl; 
#endif
		}
	}

	FTrackOnline tau = this->GetPrelimTrackFromPath(path);
	const size_t N = tau.N();
	
	// Demand all layers
	if(N != N_PAIRS or tau.GetScore() > max_cost_final_track) return;

	/* This is only a single track, and should be matched as well in the main Kalman.
	 * So don't kick its entries out of the hit matrices. */	
	// PoisonEntriesFromHMs(path);
#ifdef MND_FOOTTRACK_DEBUG
	ReplayCostCalculation( path );
#endif

	const double score = tau.GetScore();
	const FTrack& t = tau.get(); // evaluate the actual track fit.
	
#ifdef MND_FOOTTRACK_DEBUG
	std::array<double, N_PAIRS> qm, sqm; 
	auto e0_x    = mnd::make_filled_array<double, N_PAIRS>(NAN);
	auto e0_y    = mnd::make_filled_array<double, N_PAIRS>(NAN);
	auto delta_x = mnd::make_filled_array<double, N_PAIRS>(NAN);
	auto delta_y = mnd::make_filled_array<double, N_PAIRS>(NAN);
	for(size_t layer=0; layer<N; ++layer) {
		qm[layer] = tau.q[layer].mean();
		sqm[layer] = tau.q[layer].s();
		Maybe<FHitMatrix::RawHitPairRef> p = GetHitFromPath(layer, path);
		if(p) {
			FHitMatrix::RawHitPairRef& val = p.value();
			e0_x[layer]    = val[0]->Q.e0;
			delta_x[layer] = val[0]->Q.delta;
			e0_y[layer]    = val[1]->Q.e0;
			delta_y[layer] = val[1]->Q.delta;
		}
	}
#endif

	out.inner().heavy_fragment = RNFOOTTrack (
		t.l.xarray(), t.l.yarray(), t.q.mean(), score, N
#ifdef MND_FOOTTRACK_DEBUG
		,
		tau.xs,
		tau.ys,
		tau.zs,
		qm,
		sqm,
		e0_x, delta_x, e0_y, delta_y
#endif
	);
 
#ifdef MND_FOOTTRACK_DEBUG
	if(v > 2) {
		fprintf(stderr, BKH_GRN KBH_RED KBLINK "🐄🐄🐄 found a matching heavy fragment track:" KNRM BOLD);
		std::cerr << ' ' << tau << std::endl;
	}
#endif
}

using CostIndex =  std::pair<double, TFOOTHitProc::DAG::Index>;
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
			TrackCost cost(tau.N());
			u32 ncounted = 0;
			
			// FHitMatrix::Cached is column-major (Eigen::Matrix convention).
			for(size_t j=0; j<ny; ++j) {
				for(size_t i=0; i<nx; ++i) {
					const mnd::hm::Data& e = h(i,j);
			
					cost.set<TrackCost::KQ>( kq(tau, e, n) );
					cost.set<TrackCost::KR>( kr(tau, e, n) );
					if(cost > max_cost) continue;
					/* Even though other 2 costs aren't yet evaluated, `kr` is the stictest.
					 * If that one already fails here, no need evaluate it. */

					const double kt  = this->kt(tau, e, n);
					cost.set<TrackCost::KT>(kt);
					if(cost > max_cost) continue;

					// If the flow survives til this point means that
					// the candidate is good. Add it to the `path_specific_candidates_buf` buffer.
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
					new_paths.emplace_back(path); // copy-ctor.
					new_paths.back().node[n] = DAG::Index(i,j);
				}
			}
			
			/* Sort the candidates buffer. */
			std::sort(
				path_specific_candidates_buf.begin(),	
				path_specific_candidates_buf.begin() + ncounted,
				sort_by_cost
			);

			/* Take at most the best `MAX_CANDIDATES` as branching point. */
			ncounted = std::min( ncounted, MAX_CANDIDATES );
			for(u32 i=0; i<ncounted; ++i) {
				new_paths.emplace_back(path);
				new_paths.back().node[n] = path_specific_candidates_buf[i].second;
			}
			/* Add a null node only if current candidate path has no pre-existing "layer holes" */ 
			if(path.Rank() >= static_cast<int>(n))
				new_paths.emplace_back(path);

		} // loop over existing paths
		
		dag.path = std::move( new_paths );
	} // for(u32 n = 0; n < N_PAIRS; ++n)
}

void TFOOTHitProc::AnalyseDAG() noexcept {
	static_assert(N_PAIRS - 1 >= 3, "Maximal paranoia");

	/* Allow only one layer to be missing. Must have 3 or more 
	 * valid measurements and be below total allowed cost. */
	dag.KeepIf( [this](const DAG::DAGPath& p) -> bool {
		return (
			(p.Rank() > (int)N_PAIRS - 1) and
			this->GetPrelimTrackFromPath(p).GetScore() < max_cost_final_track 
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
			const FTrackOnline tl = this->GetPrelimTrackFromPath(pl);
			const FTrackOnline tr = this->GetPrelimTrackFromPath(pr);
			return tl.GetScore() < tr.GetScore();
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
		const double score = tau.GetScore();
		
		const FTrack& t = tau.get(); // Evaluate the actual track fit.

#ifdef MND_FOOTTRACK_DEBUG
		auto xs      = mnd::make_filled_array<double, N_PAIRS>(NAN); 
		auto ys      = mnd::make_filled_array<double, N_PAIRS>(NAN); 
		auto zs      = mnd::make_filled_array<double, N_PAIRS>(NAN); 
		auto qm      = mnd::make_filled_array<double, N_PAIRS>(NAN); 
		auto sqm     = mnd::make_filled_array<double, N_PAIRS>(NAN); 
		auto e0_x    = mnd::make_filled_array<double, N_PAIRS>(NAN);
		auto e0_y    = mnd::make_filled_array<double, N_PAIRS>(NAN);
		auto delta_x = mnd::make_filled_array<double, N_PAIRS>(NAN);
		auto delta_y = mnd::make_filled_array<double, N_PAIRS>(NAN);

		for(size_t layer = 0, valid=0; layer < N; ++layer) {
			if( path.node[layer] ) {
				xs[layer] = tau.xs[valid];
				ys[layer] = tau.ys[valid];
				zs[layer] = tau.zs[valid];
				qm[layer] = tau.q[valid].mean();
				sqm[layer] = tau.q[valid].s();
				
				// Is asserted to return a value based on the `if` block, but CBA changing the API.
				Maybe<FHitMatrix::RawHitPairRef> p = GetHitFromPath(layer, path);
				if(!p) ERROR("huh?");
				FHitMatrix::RawHitPairRef& val = p.value();
				e0_x[layer]    = val[0]->Q.e0;
				delta_x[layer] = val[0]->Q.delta;
				e0_y[layer]    = val[1]->Q.e0;
				delta_y[layer] = val[1]->Q.delta;

				++valid;
			}
		}
#endif

		out.inner().track.emplace_back (
			t.l.xarray(), t.l.yarray(), t.q.mean(), score, N
#ifdef MND_FOOTTRACK_DEBUG
				, xs, ys, zs, qm, sqm, e0_x, delta_x, e0_y, delta_y
#endif
		);
		
		/* Poison the exported hitmatrix elements along their rows and columns.
		 * This way, no other track can reuse either the same `x` or the same `y` entry. */
		PoisonEntriesFromHMs(path);
	}
}

/* High frequency call, decoding a DAGPath sequence into an online track. */
FTrackOnline TFOOTHitProc::GetPrelimTrackFromPath(const DAG::DAGPath& p) const noexcept {
	static_assert(DAG::depth == N_PAIRS, "Just in case. Must pass");
	FTrackOnline t{};
	
	mnd::static_for<0, N_PAIRS>([this, &p, &t](auto layer_) {
		constexpr size_t i_ = decltype(layer_)::value; // layer index
		const DAG::Index& i = p.node[i_];
		const FHitMatrix& h = hm[i_];
		if(i) { // explicit operator bool()
			t.Add( h(i[0], i[1]) , pair_z[i_] );
		}
	});

	return t;
}

#ifdef MND_FOOTTRACK_DEBUG
Maybe<FHitMatrix::RawHitPairRef> TFOOTHitProc::GetHitFromPath(size_t layer, const DAG::DAGPath& p) const {
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
		out.diff_heavy_frag_vs_upstream->Fill (	/* Hardly converged, fucked up due to rounding. */
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
		const FTrackOnline tau = this->GetPrelimTrackFromPath(p);

		/* Calculate cost of the current path `p`, if I add the node `index` from i'th layer. */
		const double kq = this->kr(tau, e, n);
		const double kr = this->kr(tau, e, n);
		const double kt = this->kt(tau, e, n);
		out.h1_acc_q[n]->Fill( kq / Cq ); 
		out.h1_acc_r[n]->Fill( kr / Cr );	 
		out.h1_acc_t[n]->Fill( kt / Ct );

		p.node[n] = path.node[n];
	}
}
#endif
