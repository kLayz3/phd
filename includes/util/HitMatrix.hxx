#pragma once

#include "../monad/monad.hxx"

#include <cmath>
#include <type_traits>
#include "../Eigen.h"


#ifndef MND_HITMATRIX_DO_BOUNDS_CHECKmessed //#	define MND_HITMATRIX_DO_BOUNDS_CHECK
#endif

template<typename FOOTPair> struct HitMatrix;

namespace mnd { namespace hm {

struct Q {
	float qx, qy;

	inline double mean() const noexcept { return 0.5*(qx+qy); }
	inline double var() const noexcept {
		double d = qx-qy;
		return 0.5*d*d;
	}
	inline double s() const noexcept { return std::sqrt( var() ); } 

	inline bool HasValue() const noexcept { return std::isfinite(qx); }
	inline explicit operator bool() const noexcept { return HasValue(); }

	friend std::ostream& operator<<(std::ostream& os, const Q& rhs) {
		return os << rhs.mean() << "±" << rhs.s();
	}
};
static_assert(sizeof(Q) == 2*sizeof(float));

/* Data cannot be packed as just (v,q) into std::optional, as also
 * the third state is required. So the type is extended a'la optional,
 * but with the third state present. Furthermore, std::optional<Data> would be padded to 64-byte wide. */
struct Data {
	using q_type = Q;
	using xy_type = Eigen::Vector2d;
	
	enum State : u32 {
		UNEVALUATED, // ditto
		READY,       // if the current entry is fully evaluated
		POISONED     // if the current entry got plucked into a well defined track. 
	};

	xy_type v;   // 16B
	q_type q;    //  8B
	State state; //  4B

	inline void reset() noexcept { state = UNEVALUATED; }
	inline bool has_value() const noexcept { return state != UNEVALUATED; }
	inline bool is_poisoned() const noexcept { return state == POISONED; }
		
	inline double X() const noexcept { return v(0); }
	inline double Y() const noexcept { return v(1); }
};
static_assert(std::is_aggregate_v<Data> && sizeof(Data) == 32);

inline std::ostream& operator<<(std::ostream& os, const Data& d) {
	return os << '(' << d.v.transpose() << "; " << d.q << ')';
}

struct Cached {
	template<typename> friend struct ::HitMatrix; 
	using Storage = Eigen::Matrix<
		Data,
		Eigen::Dynamic,
		Eigen::Dynamic
	>;

	decltype(auto) operator()(size_t i, size_t j) const noexcept { return ( cache(i,j) ); }
	decltype(auto) operator()(size_t i, size_t j)       noexcept { return ( cache(i,j) ); }

	inline void resize(size_t nx, size_t ny) __attribute__((always_inline)) {
		cache.resize(static_cast<Eigen::Index>(nx),
		             static_cast<Eigen::Index>(ny));
	}
	inline void clear() __attribute__((always_inline)) {
		for(auto j = 0; j < cache.cols(); ++j)
			for(auto i = 0; i < cache.rows(); ++i) {
				cache(i,j).reset(); 
			}
	}
	inline void resize_and_clear(size_t nx, size_t ny) { 
		this->resize(nx, ny);
		this->clear();
	}
	
private:
	Storage cache;
};
} /* namespace hm */ 

#define DECL_TYPE_TRAIT_HAS_FIELD(FIELD) \
	template<typename, typename = void> \
	struct has_##FIELD##_field : std::false_type {}; \
	\
	template<typename T> \
	struct has_##FIELD##_field<T, std::void_t<decltype(std::declval<T&>().FIELD)>> \
		: std::true_type {}; \

DECL_TYPE_TRAIT_HAS_FIELD(x)
DECL_TYPE_TRAIT_HAS_FIELD(y)
DECL_TYPE_TRAIT_HAS_FIELD(z)
DECL_TYPE_TRAIT_HAS_FIELD(Q)
DECL_TYPE_TRAIT_HAS_FIELD(q)
DECL_TYPE_TRAIT_HAS_FIELD(m)

} /* namespace mnd */ 

template<typename FOOTPair>
struct HitMatrix {
private:
	static_assert(mnd::has_x_field<FOOTPair>::value,
		"Underlying type `T` must have `x` field.");
	static_assert(mnd::has_x_field<FOOTPair>::value,
		"Underlying type `T` must have `y` field.");
	static_assert(mnd::has_x_field<FOOTPair>::value,
		"Underlying type `T` must have `z` field.");
	using hit_type_x = typename mnd::is_indexable_range <
		std::decay_t<decltype( std::declval<FOOTPair&>().x )>
	>::underlying_type;
	using hit_type_y = typename mnd::is_indexable_range <
		std::remove_reference_t<std::remove_cv_t<
			decltype( std::declval<FOOTPair&>().y )
		>>
	>::underlying_type;

	static_assert(std::is_same<hit_type_x, hit_type_y>::value,
		"FOOTPair::x and FOOTPair::y must contain the same hit type");
	static_assert(mnd::has_m_field<hit_type_x>::value,
		"hit_type must have field .m");
	static_assert(mnd::has_Q_field<hit_type_x>::value,
		"hit_type must have field .Q");
	
	using q_type_x = std::decay_t<decltype( std::declval<hit_type_x>().Q )>;
	static_assert(mnd::has_q_field<q_type_x>::value,
		"Q_type must have field .q");

public:
	using hit_type = hit_type_x;
	using q_type = q_type_x;
	using Entry = mnd::hm::Data;
	using Cached  = mnd::hm::Cached;
	using ActiveIndex = uint32_t;
	static constexpr size_t X = 0;
	static constexpr size_t Y = 1;
	
	HitMatrix() = default;
	HitMatrix(FOOTPair const& rhs): p(&rhs) {}

	/* Evaluate c(i,j) element. */
	Entry const& operator()(ActiveIndex i, ActiveIndex j) const noexcept {
#if defined(MND_HITMATRIX_DO_BOUNDS_CHECK) 
		if(i > static_cast<ActiveIndex>(GetN<X>()))
			ERROR("HitMatrix::operator(i,j): Requested x-index %u out of bounds (%zu)", static_cast<u32>(i), GetN<X>());
		if(j > static_cast<ActiveIndex>(GetN<Y>()))
			ERROR("HitMatrix::operator(i,j): Requested y-index %u out of bounds (%zu)", static_cast<u32>(j), GetN<Y>());
#endif

		auto& slot = cache(i,j);
		if(! slot.has_value() ) {
			const hit_type& hx = p->x[i];
			const hit_type& hy = p->y[j];
			slot = Entry {
				A * ( Eigen::Vector2d(hx.m, hy.m) + dxy ),
				{hx.Q.q, hy.Q.q}, 
				Entry::READY
			}; 
		}
		return slot;
	}
	
	/* Called at the start of tracking, to initialize the entry. */
	void InitEvent(const FOOTPair& cont ) noexcept {
		p = &cont;
		size_t const nx = p->x.size(), ny = p->y.size();
		cache.resize_and_clear(nx, ny);
	}

	bool is_poisoned(ActiveIndex i, ActiveIndex j) const noexcept {
#if defined(MND_HITMATRIX_DO_BOUNDS_CHECK) 
		if(i > static_cast<ActiveIndex>(GetN<X>()) or j > static_cast<ActiveIndex>(GetN<Y>())) {
			ERROR("HitMatrix::is_poisoned(): Requested %s-index %d/%d out of bounds (%zu/%zu)", 
					"X/Y", i, j, GetN<X>(), GetN<Y>() );
		}
#endif
		return cache(i,j).is_poisoned();
	}

	/* Flag the entire row `i` and column `j` as 'poisoned'. It won't get removed from the index redirection table, but will
	 * rather just carry the byteflag. */
	void poison(ActiveIndex i, ActiveIndex j) noexcept {
		const size_t nx = GetN<X>();
		const size_t ny = GetN<Y>();

#if defined(MND_HITMATRIX_DO_BOUNDS_CHECK) 
		if(i > static_cast<ActiveIndex>(nx) or j > static_cast<ActiveIndex>(ny)) {
			ERROR("HitMatrix::poison(): Requested %s-index %d/%d out of bounds (%zu/%zu)", 
				"X/Y", i, j, nx, ny);
		}
#endif
		auto col = cache.cache.col(j); // Eigen::ColXpr
		for(size_t i_ = 0; i_ < nx; ++i_) 
			col(i_).state = Entry::POISONED;

		auto row = cache.cache.row(i); // Eigen::ColXpr
		for(size_t j_ = 0; j_< ny; ++j_) 
			row(j_).state = Entry::POISONED;

		// Exact (i,j) entry addressed twice, but its fine.	
	}
	
	/* Eager function, force the computation of the entire matrix. */
	const Cached& EvalAll() const noexcept {
		const size_t nx = GetN<X>();
		const size_t ny = GetN<Y>();
		cache.resize(nx, ny);

		for(size_t j=0; j<ny; ++j) {
			const hit_type& hy = p->y[j];
			for(size_t i=0; i<nx; ++i) {
				if(cache(i,j).has_value()) continue;
				
				const hit_type& hx = p->x[i];

				cache(i,j) = Entry {
					A * (Eigen::Vector2d(hx.m, hy.m) + dxy),
					{ hx.Q.q, hy.Q.q }, 
					Entry::READY
				};
			}
		}
		return cache;
	}

	/* Dimension of the full matrix state along the axes. Invariant relative to the prior `pop` invocations. */
	template<size_t L>
	inline size_t GetN() const noexcept {
		if constexpr(L == X)
			return p->x.size(); 
		else if constexpr(L == Y)
			return p->y.size();
		else
			static_assert(L == X || L == Y, "Axis parameter supplied must be 'X' or 'Y' (or 0,1) respectively.");
	}   

	Eigen::Matrix2d A; // Rotation matrix, each detector might be rotated by an angle θx or θy
	Eigen::Vector2d dxy; // Translation vector, each detector might be translated by dx/dy

	friend std::ostream& operator<<(std::ostream& os, const HitMatrix& hm) {
		constexpr size_t X = HitMatrix<FOOTPair>::X;
		constexpr size_t Y = HitMatrix<FOOTPair>::Y;

		const size_t nx =  hm.template GetN<X>();
		const size_t ny =  hm.template GetN<Y>();

		os << KBH_YEL 
		   << " --- == |Y| ==> " << KNRM;
		for(size_t j=0; j < ny; ++j) {
			os << hm.p->y[j] << hm.template GetPoisonLabel<Y>(j);
			if(nx > 0)
			if(j != hm.template GetN<Y>()-1)
				os << ", ";
		}	
		os << std::endl << KBH_YEL 
		   << " --- vv |X| vvv";
		for(size_t i=0; i < nx; ++i) {
			os << '\n' << hm.p->x[i] << hm.template GetPoisonLabel<X>(i);
		}
		os << KBH_YEL << " -----END---- " << KNRM;
		return os;
	}

private: 
	const FOOTPair* p; /* Reassigned, event-by-event. */
	mutable Cached cache;
	
	template<size_t L>
	inline const char* GetPoisonLabel(ActiveIndex i) const {
		if      constexpr(L == X) {
			return (GetN<!L>() > 0)? (cache(i,0).is_poisoned()? "🐍": (cache(i,0).has_value()? "✅": "❓")): "";
		}
		else if constexpr(L == Y) {
			return (GetN<!L>() > 0)? (cache(0,i).is_poisoned()? "🐍": (cache(0,i).has_value()? "✅": "❓")): "";
		}
		else static_assert(L == X || L == Y, "Axis parameter supplied must be 'X' or 'Y' (or 0,1) respectively.");	
	}
}; // struct HitMatrix

struct RNFOOTPair;
struct FOOTHit;
//extern template struct HitMatrix<RNFOOTPair>;

/* No implementation file. The type is instantiated in the TFOOTHitProc.cxx. */
