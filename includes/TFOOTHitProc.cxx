#include "TFOOTHitProc.h"
#include "TFOOTCalCont.h"
#include "TFOOTHitCont.h"
#include "TFOOTMapCont.h"
#include "util/PolyFitter.h"
#include <algorithm>
#include <cmath>

#define GEN_ARG_INSTANCE_FOOT(z, n, data) \
	const TFOOTCalCont& in_##n
#define GEN_ARG_NAME_FOOT(z, n, data) in_##n

/* This part is a bit stupid, since anyway the field inside TFOOTHitProc will instantiate
 * the template, but somehow Cling misses this (???) and (re)instantiates it, running into the 
 * good old Eigen non-digestability. Hence, we *explicitly* instantiate it here. Once and for all. */
template struct HitMatrix<RNFOOTPair>;
template struct Track<RNFOOTPair>;

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

/* Read the param file and create the matrix `A`:
 *  1/cos(θx-θy) * ( cos(θx)  -sin(θx) )
 *                 ( sin(θy)   cos(θy) ) 
 *  and the offset vector `dxy`. */
static void SetConversionMatrices(FHitMatrix& hm, const FOOTParam& px, const FOOTParam& py ) {
	double tx = px.delta_a * M_PI / 180.0, ty = py.delta_a * M_PI / 180.0;	
	hm.A << std::cos(ty), sin(tx),
	    -std::sin(ty), cos(ty);
	hm.A *= 1.0/std::cos(tx-ty);

	double dx = px.delta_p, dy = py.delta_p; // already in [mm] scale, don't need to convert.
	hm.dxy << dx,
	       dy;
}

void TFOOTHitProc::e_to_z_t::Init(const FOOTParam& p) {
	const auto& values = p.gain.nominal_value;
	std::vector<double> x, y;
	for(auto [Z,E] : values) {
		x.push_back( std::log(Z) );
		x.push_back( std::log(E) );
	}
	auto r = PolyFit<1>(x,y);
	this->f = std::exp(-r[0]);
	this->c = 1.0 / r[1];
};

TFOOTHitProc::TFOOTHitProc (
	TFOOTHitCont& out, 
	BOOST_PP_ENUM(N_FOOT_DETECTORS, GEN_ARG_INSTANCE_FOOT, ~),
	double e_diff_tolerance
) : TFOOTHitProc::Base (
		out,
		BOOST_PP_ENUM(N_FOOT_DETECTORS, GEN_ARG_NAME_FOOT, ~)
	)
{	
	u32 i = 0;
	/* Do some verification + E-to-Q converter init */
	enum class Orientation { X, Y };
	std::vector< 
		std::map<Orientation, double>
		> test_vec(N_PAIRS);
	const FOOTBoxParam* box = out.box;
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
		}
	);
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

	/* Assign the z-values of the pairs. */
	int ipair = 0;
	const FOOTBoxParam* b = out.box;
	mnd::for_pair_in_tuple(this->in, [this, &ipair, b](const auto& f1, const auto& f2) {
			this->pair_z[ipair] = ( 
				b->GetFOOTZ(f1.setup->N) + f1.setup->dz + 
				b->GetFOOTZ(f2.setup->N) + f2.setup->dz
			) / 2.0;
			SetConversionMatrices(this->hm[ipair], *f1.setup, *f2.setup);
			ipair++;
		}
	);

	pair_z.back() = 0.0; // Dummy z- coordinate.
}

void TFOOTHitProc::ProcessEntry() noexcept {
	out.Clean();

	int ipair = 0;
	mnd::for_pair_in_tuple(this->in, [this, &ipair](const auto& f1, const auto& f2) {
		auto pair_xy = (f1.setup->orientation == "x" || f1.setup->orientation == "-x") ? std::pair{f1,f2} : std::pair{f2,f1};
		this->ProcessPair(pair_xy, ipair);
		++ipair;
	});

	ProcessTracks();
}

void TFOOTHitProc::ProcessPair (
	const std::pair<const TFOOTCalCont&, const TFOOTCalCont&>& f, i32 ipair
) noexcept {
	const FOOTParam *px = f.first.setup, *py = f.second.setup;
	const int nx = px->N, ny = py->N;
	const RNFOOTCal& fx = f.first.inner(), fy = f.second.inner();

	RNFOOTPair& output = out.inner().pair[ipair];
	for(const RNFOOTCluster& hitx : fx.fCl) {
		double delta_x = hitx.Delta();
		auto [x, Ex, mx, _] = hitx;
		Ex *= px->gain.CorrectionFactor(x, Ex);
		Ex /= px->de.CorrectionFactor(delta_x);

		double Zx = e_to_z[nx](Ex);
		output.x.emplace_back(Zx, x);
	}
	for(const RNFOOTCluster& hity : fy.fCl) {
		double delta_y = hity.Delta();
		auto [y, Ey, my, _] = hity;
		Ey *= py->gain.CorrectionFactor(y, Ey);
		Ey /= py->de.CorrectionFactor(delta_y);

		double Zy = e_to_z[ny](Ey);
		output.y.emplace_back(Zy, y);
	}
	/* Sort these vectors, in ascending values of charge (Q) */
	thread_local auto comparator = [](const auto& lhs, const auto& rhs) { return lhs.Q < rhs.Q; };
	std::sort(output.x.begin(), output.x.end(), comparator);
	std::sort(output.y.begin(), output.y.end(), comparator);
	
	output.z = pair_z[ipair];
}

constexpr auto X = FHitMatrix::X; 
constexpr auto Y = FHitMatrix::Y;
using Entry = FHitMatrix::Entry;


double TFOOTHitProc::kr(const FTrack& ft, const FHitMatrix::Entry& candidate, int k) {
	[[maybe_unused]] const auto& track = ft.get();
	FHitMatrix::Entry::xy_type xy_track = ft.extrapolate_to( pair_z[k] );
	return Cr * (xy_track - candidate.v).norm();
}

double TFOOTHitProc::kQ(const FTrack& ft, const FHitMatrix::Entry& candidate, int k) {
	(void)k;

	double mean_track_q = ft.q.mean();
	double candidate_q_mean = candidate.q.mean();
	double candidate_q_var = candidate.q.var();

	double diff = mean_track_q - candidate_q_mean;
	return CQ * (diff*diff + candidate_q_var);
}

/* Bundle these two cost functions together since both need to calculate the track update. */
double TFOOTHitProc::kt_kp(const FTrack& ft, const FHitMatrix::Entry& candidate, int k) {
	FTrack& mft = const_cast<FTrack&>(ft);
	mft.Add(candidate, pair_z[k]);
	mft.pop_back();
	return 0;
}

double TFOOTHitProc::k(const FTrack& ft, const FHitMatrix::Entry& candidate, int k) {
	return (
		kr(ft, candidate, k) +
		kQ(ft, candidate, k) +
		kt_kp(ft, candidate, k)
	);
}
void TFOOTHitProc::ProcessTracks() noexcept {
#if 0
	for(u32 i=0; i<N_PAIRS; ++i) {
		hm[i].InitEvent( out->pair[i] );
	}

	/* Idea is explained in the PhD writeup. 
	 * If you don't have it, ask Klayze. */
	auto& first = hm[0];

	[[ maybe_unused ]] size_t nx, ny;

	/* Start from highest Z values.
	 * From the initial pair's hit matrix try to match an entry with something.. */
	while((nx = first.GetN<X>()) > 0 and (ny = first.GetN<Y>()) > 0) {
		for(int j=ny-1; j>=0; --j) {
			const Entry& candidate_hit0 = first(nx-1, j);
			
			/* Take this candidate and look for candidates in the following layer. */
		}
		first.pop_back<X>();
	}
#endif
}
