#include "TFOOTHitProc.h"
#include "TFOOTCalCont.h"
#include "TFOOTHitCont.h"
#include "TFOOTMapCont.h"
#include "util/DirectedAGraph.hxx"
#include "util/PolyFitter.h"
#include <algorithm>
#include <cmath>

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
 * member functions/statics/dtor/ctors. These actually get their definition and home right here. :-) */
template struct HitMatrix<RNFOOTPair>;
template struct Track<TFOOTHitCont::N_PAIRS + 1, RNFOOTPair>;

using FHitMatrix = TFOOTHitProc::FHitMatrix;
using FTrack = TFOOTHitProc::FTrack; 

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

Verbosity TFOOTHitProc::v = Verbosity::SILENT;


std::ostream& operator<<(std::ostream& os, const TrackCost& rhs) {
	return os << "(kr: " << rhs.kr_ << ','
	          << " kq: " << rhs.kq_ << ','
	          << " kp: " << rhs.kp_ << ','
	          << " kt: " << rhs.kt_ << ')';
}

/* Read the param file and create the matrix `A`:
 *  1/cos(θx-θy) * ( cos(θx)  -sin(θx) )
 *                 ( sin(θy)   cos(θy) ) 
 *  and the offset vector `dxy`. */
void TFOOTHitProc::SetConversionMatrices(int ipair, const FOOTParam& px, const FOOTParam& py ) {
	double tx = px.delta_a * M_PI / 180.0, ty = py.delta_a * M_PI / 180.0;	
	hm[ipair].A << std::cos(ty), sin(tx),
	              -std::sin(ty), cos(tx);
	hm[ipair].A *= 1.0/std::cos(tx-ty);

	A_inv[ipair] << std::cos(tx), -sin(tx),
	                std::sin(ty), cos(ty);

	hm[ipair].dxy << px.delta_p, // already in [mm] scale, don't need to convert.
	                 py.delta_p;
	
	refl[ipair] << ((px.orientation == "x") ? 1.0 : -1.0),
		           ((py.orientation == "y") ? 1.0 : -1.0);
}

void TFOOTHitProc::e_to_z_t::Init(const FOOTParam& p) {
	const auto& values = p.gain.nominal_value;
	std::vector<double> x, y;
	for(auto [Z,E] : values) {
		x.push_back( std::log(Z) );
		y.push_back( std::log(E) );
	}
	auto r = PolyFit<1>(x,y);
	if(TFOOTHitProc::v > 0) {
		WARN("FOOT%d energy dependence: E(Z) = %.1f * Z^%.2f\n", p.N, std::exp(r[0]), r[1]);
	}
	this->f = std::exp(-r[0]);
	this->c = 1.0 / r[1];
};

TFOOTHitProc::TFOOTHitProc (
	TFOOTHitCont& out, 
	BOOST_PP_ENUM(N_FOOT_DETECTORS, GEN_ARG_INSTANCE_FOOT, ~),
	double q_tolerance_,
	double max_cost_,
	Verbosity v_
) : TFOOTHitProc::Base (
		out,
		BOOST_PP_ENUM(N_FOOT_DETECTORS, GEN_ARG_NAME_FOOT, ~)
	), q_tolerance(q_tolerance_), max_cost(max_cost_)
{
	TFOOTHitProc::v = v_;
	u32 i = 0;

	/* Do some verification + E-to-Q converter init */
	enum class Orientation { X, Y };
	std::vector< 
		std::map<Orientation, double> // orientation: z` 
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
		Orientation o = (s->orientation == "x" || s->orientation == "-x") ? Orientation::X : Orientation::Y;

		const u32 ipair = i/2;
		if(ipair > N_PAIRS) 
			ERROR("Checking for input setup validity, encountered %i>%i ?", ipair, N_PAIRS);	
		auto& map = test_vec[ipair];
		map.insert({o,z});

		this->e_to_z[i].Init(*s);
		++i;	
	});

	/* Check that the parameter file has sane input. */
	bool is_fine = std::is_sorted(test_vec.begin(), test_vec.end(),
		[](const auto& lhs, const auto& rhs) {
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
	int ipair = 0;
	const FOOTBoxParam* b = out.box;
	mnd::for_pair_in_tuple(this->in, [this, &ipair, b](const auto& f1, const auto& f2) {
			this->pair_z[ipair] = ( 
				b->GetFOOTZ(f1.setup->N) + f1.setup->dz + 
				b->GetFOOTZ(f2.setup->N) + f2.setup->dz
			) / 2.0;
			Orientation o1 = (f1.setup->orientation == "x" || f1.setup->orientation == "-x") ? Orientation::X : Orientation::Y;
			if(o1 == Orientation::X) { 
				this->SetConversionMatrices(ipair, *f1.setup, *f2.setup);
			} else { 
				this->SetConversionMatrices(ipair, *f2.setup, *f1.setup);
			}
			ipair++;
		}
	);
	const ExpertTarget& target = box->target;
	
	/* Assign S2 Be target parameters (EXPERT target). */
	this->target_z = 0.0; // by convention
	this->target_xy = mnd::geom::Rectangle2D (
		{target.dx, target.dy}, target.width_x, target.width_y
	);
}

void TFOOTHitProc::ProcessEntry() noexcept {
	out.Clean();
	
	if(v > 1) {
		static u64 ev_num = 0;
		fprintf(stderr, "\n%s>>> Entry[%lu] <<<%s\n", KBH_GRN, ++ev_num, KNRM);
	}
	int ipair = 0;
	mnd::for_pair_in_tuple(this->in, [this, &ipair](const auto& f1, const auto& f2) {
		const std::pair<const TFOOTCalCont&, const TFOOTCalCont&> 
			pair_xy = (f1.setup->orientation == "x" || f1.setup->orientation == "-x") ? std::pair{f1,f2} : std::pair{f2,f1};

		this->ProcessPair(pair_xy, ipair);
		++ipair;
	});

	for(u32 i=0; i<N_PAIRS; ++i)
		hm[i].InitEvent( out->pair[i] );

	ConstructObviousTracks();
	//ConstructDAG();
}

void TFOOTHitProc::ProcessPair (
	const std::pair<const TFOOTCalCont&, const TFOOTCalCont&>& f, i32 ipair
) noexcept {
	const FOOTParam *px = f.first.setup, *py = f.second.setup;
	const int nx = px->N, ny = py->N;
	const RNFOOTCal &fx = f.first.inner(), &fy = f.second.inner();

	RNFOOTPair& output = out.inner().pair[ipair];

	for(const RNFOOTCluster& hitx : fx.fCl) {
		double delta_x = hitx.Delta();
		auto [x, Ex, mult, ctype] = hitx;
		Ex *= px->gain.CorrectionFactor(x, Ex);
		Ex /= px->de.CorrectionFactor(delta_x);

		double Zx = e_to_z[nx](Ex);

		// Cluster size 1 fucks with everything above Z >~ 1,
		// so only care about it if its sitting at low energies.. 
		if(mult > 1 or Zx < CLUSTER_SIZE_ONE_Q_CUTOFF)
			output.x.emplace_back(Zx, mult, ctype, x);
	}
	for(const RNFOOTCluster& hity : fy.fCl) {
		double delta_y = hity.Delta();
		auto [y, Ey, mult, ctype] = hity;
		Ey *= py->gain.CorrectionFactor(y, Ey);
		Ey /= py->de.CorrectionFactor(delta_y);

		double Zy = e_to_z[ny](Ey);
		if(mult > 1 or Zy < CLUSTER_SIZE_ONE_Q_CUTOFF)
			output.y.emplace_back(Zy, mult, ctype, y);
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

double TFOOTHitProc::kr(const FTrack& ft, const FHitMatrix::Entry& candidate, u32 k) const {
	Eigen::Vector2d extrapolated = ft.extrapolate_to( pair_z[k] );
	const Eigen::Vector2d& measured = candidate.v;
	return Cr * ((extrapolated - measured).norm());
	return 0;
}

double TFOOTHitProc::kq(const FTrack& ft, const FHitMatrix::Entry& candidate, u32 k) const {
	(void)k;

	double mean_track_q = ft.q.mean();
	double candidate_q_mean = candidate.q.mean();

	double diff = mean_track_q - candidate_q_mean;
	return Cq * diff*diff;
}

/* Bundle these two cost functions together since both need to calculate the track update. */
std::pair<double,double> TFOOTHitProc::kt_kp(const FTrack& ft, const FHitMatrix::Entry& candidate, u32 k) const {
	FTrack& mft = const_cast<FTrack&>(ft);
	mft.Add(candidate, pair_z[k]);

	double cost_p = 0;
	double cost_t = 0;

	// Target sits nominally at z=0.
	Eigen::Vector2d extrapolate_to_target = mft.extrapolate_to( target_z ).array();
	if( !target_xy.IsInside( extrapolate_to_target.x(), extrapolate_to_target.y() ) )
		cost_t = Ct;

	// Check that track goes through next layer `k+1`. If we are in last layer, it's a no-op.
	if(k < static_cast<u32>(pair_z.size() - 1)) {
		const Eigen::Vector2d extrapolate_to_next_layer = mft.extrapolate_to( pair_z[k+1] );
		
		Eigen::Array2d pair_coords = MM_TO_STRIP * refl[k+1].cwiseProduct( 
				A_inv[k+1] * (extrapolate_to_next_layer - hm[k+1].dxy)
			).array() + DETECTOR_MIDPOINT;
		
		double cx = pair_coords.x();
		double cy = pair_coords.y();
		
		if(cx < 0 || cx > N_STRIPS || cy < 0 || cy > N_STRIPS)
			cost_p = Cp;
	} 

	mft.pop_back();
	return { cost_t, cost_p };
}

void TFOOTHitProc::ConstructDAG() noexcept {
	/* Idea is explained in the PhD writeup. 
	 * If you don't have it, ask Klayze. */

	/* To keep the algorithm invariant between layers, tracks can also be nullable
	 * (for now). And we add an extra `-1`th layer, which can only ever have an empty hitmatrix. */
	
	dag.Initialize();
	TrackCost cost{};

	for(u32 n = 0; n < N_PAIRS; ++n) {
		const FHitMatrix& h = hm[n];

		const size_t nx = h.GetN<X>();
		const size_t ny = h.GetN<Y>();
		if(nx == 0 or ny == 0) continue;

		std::vector<DAG::DAGPath>& current_paths = dag.path;
		std::vector<DAG::DAGPath>  new_paths( 2*current_paths.size() ); // could be larger.
		
		/* Each path already draws (an optional) preliminary track. 
		 * Try to match some of the current layer's hitmatrix elements against that */
		for(const DAG::DAGPath& path : current_paths) {

			// Fetch the preliminary track that the path describes.
			FTrack tau = this->GetPrelimTrackFromPath(path);
			const FTrack::Status status = tau.GetStatus();
			
			// FHitMatrix::Cached is column-major (Eigen convention).
			for(size_t j=0; j<ny; ++j) {
				for(size_t i=0; i<nx; ++i) {
					const mnd::hm::Data& e = h(i,j);
					
					cost = TrackCost{};
					
					switch(status) {
						case FTrack::Status::WellDefined: {
							cost.set<TrackCost::KR>( kr(tau, e, n) );
							if(cost.sum() > max_cost) continue;
							[[ fallthrough ]];
						}
						case FTrack::Status::SinglePoint: {
							cost.set<TrackCost::KQ>( kq(tau, e, n) );
							if(cost.sum() > max_cost) continue;
							auto [kt, kp] = kt_kp(tau, e, n);
							cost.set<TrackCost::KP>(kp);
							cost.set<TrackCost::KT>(kt);
							if(cost.sum() > max_cost) continue;
							[[ fallthrough ]];
						}
						case FTrack::Status::Bare: {
							double qv = e.q.var();
							cost.set<TrackCost::KQ>( Cq * qv*qv );
							if(cost.sum() > max_cost) continue;
						}
					}
					
					// If the flow survives til this point means that
					// the candidate is stellar. Add it to the list of paths.
					new_paths.emplace_back(path); // copy-ctor.
					new_paths.back().node[n] = DAG::Index(i,j);
				}
			}

			// Add a null node.
			new_paths.emplace_back(path);
		}
		
		dag.path = std::move( new_paths );
	}
}

FTrack TFOOTHitProc::GetPrelimTrackFromPath(const DAG::DAGPath& p) const {
	FTrack t{};
	static_assert(DAG::depth == N_PAIRS, "Just in case. Must pass");
	for(size_t i_ = 0; i_ < N_PAIRS; ++i_) {
		const DAG::Index& i = p.node[i_];
		const FHitMatrix& h = hm[i_];
		if(i) {
			t.Add( h( i[0], i[1] ) , pair_z[i_] );
		}
	}
	return t;
}

/* Some heavy particle tracks can be obviously taken out of the full glory
 * CKF algorithm. E.g. if each FOOT in-order measures Z=4,5,6 is fine. 
 * NOTE: the only obvious particle we can find at this point, is the largest Q one. */
void TFOOTHitProc::ConstructObviousTracks() noexcept {
	static_assert(DAG::depth == N_PAIRS, "Hello there");

	/* Do the same CKF but with reduced phase space: meaning that 
	 * we only take last row/column as candidates.
	 * Taking only the last entry (largest Qx,Qy) is a bit lazy, since the hit could go through
	 * noisy/dead strip, and the corresponding entry won't be the largest Z. */

	/* Only one path is viable here, no branching possible. Always take the best candidate... */
	DAG::DAGPath path{};

	TrackCost cost {};

	for(u32 n = 0; n < N_PAIRS; ++n) {
		const FHitMatrix& h = hm[n];
		if(v > 1)
			std::cerr << "PFOOT" << n << ":\n" << h << std::endl;
		const size_t nx = h.GetN<X>();
		const size_t ny = h.GetN<Y>();

		/* If even a single layer shows no hits, the whole call returns, basically. */
		if(nx == 0 or ny == 0) break; 

		// Fetch the preliminary track that the path describes.
		FTrack tau = this->GetPrelimTrackFromPath(path);
		const FTrack::Status status = tau.GetStatus();

		double cost_min_current = INFINITY;
		DAG::Index best_i;

#define EXPAND_COST_CALC \
	{ \
		const mnd::hm::Data& e = h(i,j); \
		cost = TrackCost{}; \
		switch(status) { \
			case FTrack::Status::WellDefined: { \
				cost.set<TrackCost::KR>( kr(tau, e, n) ); \
				/* if(cost.sum() > max_cost) continue; */\
				[[ fallthrough ]]; \
			} \
			case FTrack::Status::SinglePoint: { \
				cost.set<TrackCost::KQ>( kq(tau, e, n) ); \
				/* if(cost.sum() > max_cost) continue; */\
				auto [kt, kp] = kt_kp(tau, e, n); \
				cost.set<TrackCost::KP>(kp); \
				cost.set<TrackCost::KT>(kt); \
				/* if(cost.sum() > max_cost) continue; */\
				[[ fallthrough ]]; \
			} \
			case FTrack::Status::Bare: { \
				double qv = e.q.var(); \
				cost.set<TrackCost::KQ>( Cq * qv*qv ); \
				/* if(cost.sum() > max_cost) continue; */\
			} \
		} \
		if(v > 2) { \
			std::cerr << DAG::Index{i,j} << " : " << cost << std::endl; \
		} \
		if(cost.sum() < cost_min_current) { \
			cost_min_current = cost.sum(); \
			best_i = {i,j}; \
		} \
	}
		// FHitMatrix::Cached is column-major (Eigen convention).
		size_t j = ny-1, i = nx-1;
		for(size_t i=0; i<nx; ++i) // shadows `i` outside
			EXPAND_COST_CALC
		for(size_t j=0; j<ny-1; ++j) // shadows `j` outside
			EXPAND_COST_CALC

		if(best_i) { // operator bool(); checks if it isn't null
			path.node[n] = best_i;
			if(v > 2)
				std::cerr << "Found best index: " << best_i << std::endl; 
		}
	}

	FTrack tau = this->GetPrelimTrackFromPath(path);
	const size_t N = tau.N();
	
	// Demand all the layers
	if(N != N_PAIRS) return;

	//double score = tau.GetScore();
	//const auto& t = tau.get();
	//out.inner().track.emplace_back (
	//	t.l.a.array(), t.l.b.array(), t.q.mean(), score, N
	//);

	//if(v > 2) {
	//	std::cerr << "Found a matching track: " << tau << std::endl;
	//}

}
