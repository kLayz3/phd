#include "TFOOTHitProc.h"

#include "TFOOTCalCont.h"
#include "TFOOTHitCont.h"
#include "TFOOTMapCont.h"
#include "util/DirectedAGraph.hxx"
#include "util/Geometry.h"
#include <algorithm>
#include <cmath>

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

Verbosity TFOOTHitProc::v = Verbosity::SILENT;

std::ostream& operator<<(std::ostream& os, const TrackCost& rhs) {
	return os << "(kr: " << rhs.kr_ << ','
	          << " kq: " << rhs.kq_ << ','
	          << " kp: " << rhs.kp_ << ','
	          << " kt: " << rhs.kt_ << ')';
}

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

	A_inv[ipair] << cos(tx), sin(tx),
	               -sin(ty), cos(ty);

	hm[ipair].dxy << px.delta_p, // already in [mm] scale, don't need to convert.
	                 py.delta_p;
	
	refl[ipair] << px.R(), 
	               py.R();
}

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
		Orientation o = (s->orientation[1] == 'x') ? Orientation::X : Orientation::Y;
		
		const u32 ipair = i/2;
		if(ipair > N_PAIRS) 
			ERROR("Checking for input setup validity, encountered %i>%i ?", ipair, N_PAIRS);	
		auto& map = test_vec[ipair];
		map.insert({o,z});

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
	mnd::for_pair_in_tuple(this->in, [this, &ipair, box](const TFOOTCalCont& f1, const TFOOTCalCont& f2) {
			this->pair_z[ipair] = ( 
				box->GetFOOTZRel(f1.setup->N) + f1.setup->dz + 
				box->GetFOOTZRel(f2.setup->N) + f2.setup->dz
			) / 2.0 + TFOOTHitProc::TARGET_Z;
			Orientation o1 = (f1.setup->orientation[1] == 'x') ? Orientation::X : Orientation::Y;
			if(o1 == Orientation::X) { 
				this->SetConversionMatrices(ipair, *f1.setup, *f2.setup);

				this->out.foot_param[ipair]->at(0) = *f1.setup;
				this->out.foot_param[ipair]->at(1) = *f2.setup;
			} else { 
				this->SetConversionMatrices(ipair, *f2.setup, *f1.setup);

				this->out.foot_param[ipair]->at(0) = *f2.setup;
				this->out.foot_param[ipair]->at(1) = *f1.setup;
			}
			++ipair;
		}
	);
	MND_ASSERT(ipair == N_PAIRS && "Paranoia");
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
		auto [cx, _, mult, ctype] = hit;

		double q = px.Q(hit);
		
		double xprime = refl[ipair].x() * (cx - FOOTParam::DETECTOR_MIDPOINT) * FOOTParam::STRIP_TO_MM; 
		// Cluster size 1 fucks with everything above Z >~ 1,
		// so only care about it if its sitting at low energies.. 
		if(mult > 1 or q < CLUSTER_SIZE_ONE_Q_CUTOFF)
			output.x.emplace_back(q, mult, ctype, xprime);
	}

	for(const RNFOOTCluster& hit : fy.fCl) {
		auto [cy, _, mult, ctype] = hit;

		double q = py.Q(hit);

		double yprime = refl[ipair].y() * (cy - FOOTParam::DETECTOR_MIDPOINT) * FOOTParam::STRIP_TO_MM; 
		if(mult > 1 or q < CLUSTER_SIZE_ONE_Q_CUTOFF)
			output.y.emplace_back(q, mult, ctype, yprime);
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

double TFOOTHitProc::kr(const FTrackOnline& ft, const FHitMatrix::Entry& candidate, u32 k) const {
	const Eigen::Vector2d& measured = candidate.v;

	Eigen::Vector2d extrapolated = ft.extrapolate_to( pair_z[k] );
	if(! extrapolated.array().isFinite().all() ) // power of expression templates :)
		return 0.0;
	else
		return Cr * ((extrapolated - measured).norm());
}

/* kQ = Cq || Qij - Qn ||^2 
 * ==== Cq( (mean(Qij) - mean(Qtrack))^2 + variance_ij )  */
double TFOOTHitProc::kq(const FTrackOnline& ft, const FHitMatrix::Entry& candidate, u32 k) const {
	(void)k;

	double mean_track_q = ft.q.mean();
	double candidate_var = candidate.q.var();

	double cost = Cq * candidate_var;
	if(std::isfinite(mean_track_q)) { // Track is non-null
		double diff = mean_track_q - candidate.q.mean();
		cost += Cq*diff*diff;
	} 
	return cost;
}

/* Bundle these two cost fncs together since both need to calculate the track update. */
std::pair<double,double> TFOOTHitProc::kt_kp(const FTrackOnline& ft, const FHitMatrix::Entry& candidate, u32 k) const {
	FTrackOnline& mft = const_cast<FTrackOnline&>(ft);
	mft.Add(candidate, pair_z[k]);
	const auto& tt = mft.get();

	double cost_p = 0;
	double cost_t = 0;

	if( tt.l.HasValue() ) {
	
		// Target sits nominally at z=0.
		Eigen::Vector2d extrapolate_to_target = mft.extrapolate_to( TARGET_Z );
		if(v > 3) {
			fprintf(stderr, "TFOOTHitProc::kt_kp: Test track: ");
			std::cerr << mft << " :: value at target: (" << extrapolate_to_target.transpose() << ')' << std::endl;
		}
		if( !target_xy.IsInside( extrapolate_to_target.x(), extrapolate_to_target.y() ) )
			cost_t = Ct;

		// Check that track goes through next layer `k+1`. If we are in last layer, it's a no-op.
		if(k < static_cast<u32>(pair_z.size() - 1)) {
			Eigen::Vector2d extrapolate_to_next_layer = mft.extrapolate_to( pair_z[k+1] );
		
			// Convert the [mm x mm] measurement from next layer, the hitmatrix entry,
			// back into strip units.
			Eigen::Array2d pair_coords = FOOTParam::MM_TO_STRIP * refl[k+1].cwiseProduct( 
					A_inv[k+1] * (extrapolate_to_next_layer) - hm[k+1].dxy
				).array() + FOOTParam::DETECTOR_MIDPOINT;
			
			double cx = pair_coords.x();
			double cy = pair_coords.y();
			
			if(cx < 0 || cx > FOOTParam::N_STRIPS || cy < 0 || cy > FOOTParam::N_STRIPS)
				cost_p = Cp;
		}
	}

	mft.pop_back();
	return { cost_t, cost_p };
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
		const size_t nx = h.GetN<X>();
		const size_t ny = h.GetN<Y>();
		const size_t n_last_row_col = nx+ny-1;

		/* If even a single layer shows no hits, the whole call abruptly returns. */
		if(nx == 0 or ny == 0) break; 

		// Fetch the preliminary track that the path describes.
		FTrackOnline tau = this->GetPrelimTrackFromPath(path);

		double cost_min_current = INFINITY;
		DAG::Index best_i;

		if(v > 1) {
			WARN(KBH_BLU "Track entering layer: #%u: " KNRM, n);
			std::cerr << tau << std::endl;
		}
		if(v > 2) {
			std::cerr << h << std::endl;
		}
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

			cost = {};
			
			cost.set<TrackCost::KQ>( kq(tau, e, n) );
			//if(cost.sum() > max_cost) continue;

			cost.set<TrackCost::KR>( kr(tau, e, n) );
			//if(cost.sum() > max_cost) continue;

			auto [kt,kp] = kt_kp(tau, e, n);
			cost.set<TrackCost::KP>(kp);
			cost.set<TrackCost::KT>(kt);
			//if(cost.sum() > max_cost) continue;

			if(v > 3) {
				std::cerr << DAG::Index{i,j} << " : " << e << " :: cost: " << cost << std::endl;
			}
			if(cost.sum() < cost_min_current) {
				cost_min_current = cost.sum();
				best_i = {i,j};
			}
		}

		if(best_i) { // operator bool(); checks if the index object is non-null
			path.node[n] = best_i;
			if(v > 2)
				std::cerr << "Found best index: " << best_i << std::endl; 
		}
	}

	FTrackOnline tau = this->GetPrelimTrackFromPath(path);
	const size_t N = tau.N();
	
	// Demand all the layers
	if(N != N_PAIRS) return;

	double score = tau.GetScore();
	const auto& t = tau.get(); // evaluate the actual track fit.
	
#ifdef MND_FOOTTRACK_DEBUG
	std::array<double, N_PAIRS> qm, sqm;
	for(size_t i=0; i<N; ++i) {
		qm[i] = tau.q[i].mean();
		sqm[i] = tau.q[i].s();
	}
#endif

	out.inner().track.emplace_back (
		t.l.xarray(), t.l.yarray(), t.q.mean(), score, N
#ifdef MND_FOOTTRACK_DEBUG
		,
		tau.xs,
		tau.ys,
		tau.zs,
		qm,
		sqm
#endif
	);
	out.h1_qtrack->Fill(t.q.mean());
	out.h1_track_nsampled->Fill(N);

	if(v > 2) {
		fprintf(stderr, BKH_GRN KBH_RED KBLINK "🐄🐄🐄 found a matching track:" KNRM BOLD);
		std::cerr << ' ' << tau << std::endl;
	}
}

void TFOOTHitProc::ConstructDAG() noexcept {
#if 0
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
#endif
}

FTrackOnline TFOOTHitProc::GetPrelimTrackFromPath(const DAG::DAGPath& p) const {
	FTrackOnline t{};
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

void TFOOTHitProc::PostProcess() noexcept {
	using namespace mnd::geom;

	/* For recognised tracks, try to find their vertex, together with the upstream track.
	 * Upstream track, however is in FRS coordinates, and must be represented in FOOT coordinates. */
	g_upstream_track %= ( out.box->GetTargetZ() + TFOOTHitProc::TARGET_Z );

	std::sort( /* Sort in descending charge (.Q) attribute. */ 
		out.inner().track.begin(), 
		out.inner().track.end() 
	);
	
	lines.clear();
	for(const auto& t : out.inner().track) {
		lines.emplace_back( RNTrackToLine3D(t) );
	}
	lines.push_back( g_upstream_track ); // copy-ctor

	mnd::geom::Point3D vertex = FindVertex( mnd::as_span(lines) );
	out.inner().vertex = RNFOOTHit::Vertex{ vertex };

	if(lines.size() >= 2) {
		out.diff_heavy_frag_vs_upstream->Fill (	
			lines.front().DistanceTo( g_upstream_track )
		);
	}
};

