#pragma once

#include "monad/monad.hxx"
#include "TFOOTMapCont.h"
#include "TFOOTCalCont.h"
#include <type_traits>

struct TFOOTMapCont;

struct TFOOTCalProc : TProcessor <
	TFOOTCalCont
	(TFOOTMapCont)
> {
	using Base = TProcessor<TFOOTCalCont(TFOOTMapCont)>;
	static constexpr int N_STRIPS    = TFOOTCalCont::N_STRIPS; /* 640 */
	static constexpr int MAX_CL_SIZE = 40; /* Maximal allowed cluster size. */
	static constexpr int MASSIVE_CLUSTER_CUTOFF = 20; /* After which multiplicity a cluster is called 'massive'. */
    static constexpr f64 HEAVY_ISOLATED_CLUSTER_CUTOFF = 4.0; /* Central strip ADC [in units of c_thr] w/o neightbours to be called 'heavy isolated'. */

	static constexpr f64 BAD_STRIP_THRESHOLD = INFINITY; /* Bad strips labelled in map step get assigned this thresholds. */

	static_assert(MAX_CL_SIZE > MASSIVE_CLUSTER_CUTOFF);

	using ClusterType = RNFOOTCluster::ClusterType;
	using ClusterFullType = std::pair<ClusterType, ClusterType>;

	friend bool operator==(const ClusterFullType&, ClusterType );
	friend bool operator!=(const ClusterFullType&, ClusterType );

	enum LookAhead { kNEG, kPOS, };
	
	/* Single hit belonging to a cluster. */
	struct TClustHit {
		double e; /* ADC value. */
		i32 x;    /* Position.  */
		enum Flag : u32 {
			PosEdge    = 1<<0,
			NegEdge    = 1<<1,
			Breakpoint = 1<<2,
			Seed       = 1<<3,
			BelowThr   = 1<<4,
		};
		u32 flag; // default: all false
	};
	static_assert(sizeof(TClustHit) == 16);

	double X_CENTRE_THR;
	double X_NEIGHB_THR;

	TFOOTCalProc(TFOOTCalCont& out, TFOOTMapCont& in);
	TFOOTCalProc() = default;

	void ProcessEntry() noexcept;

	template<LookAhead> bool _IsAddibleToCluster(const int );

	void MakeACluster(int& );

	Double_t* _e;
	Double_t e[N_STRIPS]{0};
	
	MND_FORCE_INLINE Double_t c_thr(int i) const noexcept { return _hidden_c_thr[i+1]; }
	MND_FORCE_INLINE Double_t n_thr(int i) const noexcept { return _hidden_n_thr[i+1]; }

private:
	std::array<Double_t, N_STRIPS+2> _hidden_c_thr
		= mnd::make_filled_array<Double_t, N_STRIPS+2>( INFINITY );
	std::array<Double_t, N_STRIPS+2> _hidden_n_thr
		= mnd::make_filled_array<Double_t, N_STRIPS+2>( INFINITY );

public:
	/* Temporarily storing the information
	 * for a specific collected cluster. */
	TClustHit _buf[MAX_CL_SIZE]{};
	using TClustHitIterator = std::remove_reference<decltype( std::begin(_buf) )>::type;

	u32 _cl_cnt = 0; /* Points to one-ahead valid index in `_buf`. Is also size of `_buf`. */
	ClusterFullType _ct { ClusterType::kUNKNOWN, ClusterType::kUNKNOWN };
	ClusterType GetClusterType();

	void CharacterizeCurrentCluster() noexcept;

	inline void _ClearClust() noexcept { _cl_cnt = 0; _ct = {ClusterType::kGOOD, ClusterType::kGOOD}; }
	inline void _AddHit(double e, int x, u32 flag=0) noexcept {
		if(_cl_cnt >= MAX_CL_SIZE) {
			WARN("(%s => %s) cluster size exceeded: " EMPH(max = %d) ". Ignoring the hit (strip:ADC): " EMPH(%d: %7.2f; N-Thr: %5.2f (X=%.1f)) ". Consider increasing the buffer size.\n", 
				std::get<0>(in).GetName(), out.GetName(), MAX_CL_SIZE,
				x, e, n_thr(x), X_NEIGHB_THR);
			return;
		}
		_buf[_cl_cnt].e = e;
		_buf[_cl_cnt].x = x;
		_buf[_cl_cnt].flag = flag;
		++_cl_cnt;
	}

	void PrintBuff();
	
	/* Save the state of where the possible cluster fragmentation took place. */
	TClustHit _frag{};
};
