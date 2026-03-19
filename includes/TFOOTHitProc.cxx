#include "TFOOTHitProc.h"
#include "TFOOTCalCont.h"
#include "TFOOTHitCont.h"
#include "TFOOTMapCont.h"
#include "TH2I.h"

#define GEN_ARG_INSTANCE_FOOT(z, n, data) \
	const TFOOTCalCont& in_##n
#define GEN_ARG_NAME_FOOT(z, n, data) in_##n

static constexpr int FOOTHIT_BUFFER_MAX_SIZE = 50;

TFOOTHitProc::TFOOTHitProc (
	TFOOTHitCont& out, 
	BOOST_PP_ENUM(N_FOOT_DETECTORS, GEN_ARG_INSTANCE_FOOT, ~),
	double e_diff_tolerance
) : TFOOTHitProc::Base (
		out,
		BOOST_PP_ENUM(N_FOOT_DETECTORS, GEN_ARG_NAME_FOOT, ~)
	),
	e_diff_tolerance(e_diff_tolerance)
{	
	for(auto& footbuf : buf) {
		footbuf.xs.reserve(FOOTHIT_BUFFER_MAX_SIZE);
		footbuf.ys.reserve(FOOTHIT_BUFFER_MAX_SIZE);
	}
}

static double delta_factor_basic(double delta, double s, double M) {
	double x = std::min( std::abs(delta), s );
	return M - (M-1) / s * x;
}

void TFOOTHitProc::ProcessEntry() noexcept {
	out.Clean();
	
	for(auto& footbuf : buf)
		footbuf.Clear();

	mnd::for_each_in_tuple(this->in, [this](const auto& cfoot) {
		this->ProcessSingle(cfoot);
	});
}

void TFOOTHitProc::ProcessSingle(const TFOOTCalCont& cfoot) noexcept {
	FOOTParam* p = cfoot.setup;
	FOOTBoxParam* b = out.box;
	const i32 n = p->N;
	const double z = b->GetFOOTZ(n, p); /* Second argument is for `.dz` field that p could expose. */
	const std::string& orientation = p->orientation;
	
	if(n >= N_FOOT_DETECTORS)
		ERROR("Found FOOT index: %d, and is out of range [0,%d> ?", n, N_FOOT_DETECTORS);
	
	int pn = n / 2;
	HitsBuffer& pair = this->buf[pn];
	bool is_x = (orientation == "x") || (orientation == "-x");
	std::vector<HitCandidate>& v = (is_x) ? pair.xs : pair.ys;
	if(is_x) { pair.nx = n; pair.zx = z; }
	else     { pair.ny = n; pair.zy = z; }
	
	const RNFOOTCal& foot = cfoot.inner();

	for(const RNFOOTCluster& hit : foot.fCl) {
		double delta = hit.Delta();
		double E = hit.fCE / delta_factor_basic(delta, 0.32, 2.0);
		v.emplace_back(E, hit.fCX);
		out.h_single_all[n]->Fill(E);
	}
	
	if( pair.IsValid())  {
		for(const auto& hit_x : pair.xs) {
			for(const auto& hit_y : pair.ys) {
				out.h_corr_all[pn]->Fill(hit_x.e, hit_y.e);
			}
		}

		ConstructHits(pair, pn);
	}
}

void TFOOTHitProc::ConstructHits(HitsBuffer& pair, int n) noexcept {
	auto& xs = pair.xs; auto& ys = pair.ys;
	using TV = decltype(pair.xs)::value_type;
	static auto lambda =  [](const TV& l, const TV& r) { return l.e < r.e; };

	std::sort(xs.begin(), xs.end(), lambda);
	std::sort(ys.begin(), ys.end(), lambda);
	

	while(xs.size() > 0 and ys.size() > 0) {
		const auto& hit_x = xs.back();
		const auto& hit_y = ys.back();
		
		if( IsCompatible(hit_x.e, hit_y.e) ) {
			out.inner().hits.emplace_back (
				(hit_x.e + hit_y.e) / 2, hit_x.pos, hit_y.pos, pair.Z()
			);

			out.h_corr_gud[n]->Fill(hit_x.e, hit_y.e);
			out.h_single_gud[pair.nx]->Fill(hit_x.e);
			out.h_single_gud[pair.ny]->Fill(hit_y.e);

			xs.pop_back(); ys.pop_back();
		}
		/* If not compatible, throw away the higher energy hit.
		 * References don't dangle as `pop_back` calls are final in the loop. */
		else if(hit_x.e > hit_y.e)
			xs.pop_back();
		else
			ys.pop_back();
	}
}

bool TFOOTHitProc::IsCompatible(const double ex, const double ey) noexcept {
	auto [min, max] = std::minmax(ex, ey);
	if(min > (1 - e_diff_tolerance) * max) 
		return true;
	
	return false;
}
