#pragma once

#include <tuple>
#include <cmath>
#include <optional>
#include <type_traits>
#include "../Eigen/Dense"
#include "../monad/monad.hxx"

namespace mnd { namespace hm {

struct Q {
	double qx, qy;

	inline double mean() const noexcept { return 0.5*(qx+qy)/2; }
	inline double var() const noexcept {
		double d = qx-qy;
		return 0.5*d*d;
	}
	inline double s() const noexcept { return std::sqrt( var() ); } 
};

struct Data {
	using q_type = Q;
	using xy_type = Eigen::Vector2d;

	q_type q; 
	xy_type v; 
	inline double X() const noexcept { return v(0); }
	inline double Y() const noexcept { return v(1); }
};

struct Cached {
	using Storage = Eigen::Matrix<
		std::optional<Data>,
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
			for(auto i = 0; i < cache.rows(); ++i)
				cache(i,j).reset();
	}
	inline void resize_and_clear(size_t nx, size_t ny) { 
		this->resize(nx, ny);
		this->clear();
	}

private:
	Storage cache;
};
} /* namespace hm */ 

template<typename, typename = void>
struct is_indexable_range : std::false_type {};
template<typename T>
struct is_indexable_range <
	T,
	std::void_t <
		typename T::value_type,
		decltype(std::begin(std::declval<T&>())),
		decltype(std::end(std::declval<T&>())),
		decltype(std::declval<T&>().size()),
		decltype(std::declval<T&>()[std::declval<std::size_t>()])
	>
> : std::true_type { 
	using underlying_type = typename T::value_type;
};

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

public:
	using hit_type = hit_type_x;
	using Entry = mnd::hm::Data;
	using Cached  = mnd::hm::Cached;
	using ActiveIndex = uint32_t;

	
	HitMatrix() = default;
	HitMatrix(FOOTPair const& rhs): p(&rhs) {}
	static constexpr size_t X = 0;
	static constexpr size_t Y = 1;

	/* Evaluate c_ij element. Note there's no bounds checking! */
	Entry const& operator()(ActiveIndex i, ActiveIndex j) const noexcept {
		const auto real_i = static_cast<Eigen::Index>(active_x[i]);
		const auto real_j = static_cast<Eigen::Index>(active_y[j]);

		auto& slot = cache(real_i, real_j);
		if(! slot.has_value() ) {
			const hit_type& hx = p->x[real_i];
			const hit_type& hy = p->y[real_j];
			slot = Entry {
				{hx.Q, hy.Q}, 
				A * Eigen::Vector2d(hx.m, hy.m) + dxy
			}; 
		}
		return *slot;
	}
	
	/* Bound checking version of operator() */
	Entry const& at(ActiveIndex i, ActiveIndex j) const {
		if(i > static_cast<ActiveIndex>(GetN<X>()))
			ERROR("HitMatrix::at(): Requested x-index %u out of bounds (%zu)", static_cast<u32>(i), GetN<X>());
		if(i > static_cast<ActiveIndex>(GetN<X>()))
			ERROR("HitMatrix::at(): Requested x-index %u out of bounds (%zu)", static_cast<u32>(i), GetN<Y>());
		return this->operator()(i,j);
	}
	
	/* Called at the start of tracking, to initialize the entry. */
	void InitEvent(const FOOTPair& cont ) {
		p = &cont;
		size_t const nx = p->x.size(), ny = p->y.size();
		cache.resize_and_clear(p->x.size(), p->y.size());

		active_x.resize(nx); active_y.resize(ny);
		std::iota(active_x.begin(), active_x.end(), 0);
		std::iota(active_y.begin(), active_y.end(), 0);
	}

	/* Remove a row/column from current the hit matrix by removing that entry from the index redirection table. 
	 * Does not check for bounds! */
	template<size_t L> void pop(ActiveIndex index) noexcept {
		if      constexpr(L == X) active_x.erase(active_x.begin() + index);
		else if constexpr(L == Y) active_y.erase(active_x.begin() + index);
		else static_assert(L == X || L == Y, "Axis parameter supplied must be 'X' or 'Y' (or 0,1) respectively.");
	}
	inline void pop(ActiveIndex i, ActiveIndex j) noexcept { pop<X>(i); pop<Y>(j); }

	/* Remove last row or column. Does not check for bounds. */
	template<size_t L> void pop_back() noexcept {
		if      constexpr(L == X) active_x.pop_back();
		else if constexpr(L == Y) active_y.pop_back();
		else static_assert(L == X || L == Y, "Axis parameter supplied must be 'X' or 'Y' (or 0,1) respectively.");
	}
	inline void pop_back() noexcept { pop_back<X>(); pop_back<Y>(); }

	/* Remove last row or column while doing the bounds check. */
	template<size_t L> void pop_back_checked() noexcept {
		if(0 == GetN<L>()) ERROR("HitMatrix::pop_back_checked<%s>(): Vector is empty.", (L == X) ? "X" : "Y");
		pop_back<L>();
	}
	inline void pop_back_checked() noexcept { pop_back_checked<X>(); pop_back_checked<Y>(); }

	/* Remove a row/column from current hitmatrix by removing that entry from the index redirection table. 
	 * while doing a bounds checked. */
	template<size_t L> void pop_checked(ActiveIndex index) {
		if(index > static_cast<ActiveIndex>(GetN<L>())) {
			ERROR("HitMatrix::pop_checked<%zu>(): Requested %s-index %d out of bounds (%zu)", 
				L, (L == X) ? "X" : "Y", index, GetN<L>() );
		}
		pop<L>(index);
	}
	inline void pop_checked(ActiveIndex i, ActiveIndex j) noexcept { pop_checked<X>(i); pop_checked<Y>(j); }

	/* Eager function, force the computation of the entire matrix. */
	const Cached& EvalAll() const {
		const size_t nx = GetNBase<X>();
		const size_t ny = GetNBase<Y>();
		cache.resize(nx, ny);

		for(size_t i=0; i<nx; ++i) {
			const hit_type& hx = p->x[i];
			for(size_t j=0; j<ny; ++j) {
				if(cache(i,j).has_value()) continue;

				const hit_type& hy = p->y[j];
				cache(i,j) = Entry {
					{hx.Q, hy.Q}, 
					A * Eigen::Vector2d(hx.m, hy.m) + dxy
				};
			}
		}
		return cache;
	}

	/* Dimension of the full matrix state along the axes. Invariant relative to the prior `pop` invocations. */
	template<size_t L>
	inline size_t GetNBase() const noexcept {
		if constexpr(L == X)
			return p->x.size(); 
		else if constexpr(L == Y)
			return p->y.size();
		else
			static_assert(L == X || L == Y, "Axis parameter supplied must be 'X' or 'Y' (or 0,1) respectively.");
	}   
	/* Dimension of the current state of the matrix along the axes. */
	template<size_t L>
	inline size_t GetN() const noexcept {
		if constexpr(L == X)
			return active_x.size(); 
		else if constexpr(L == Y)
			return active_y.size();
		else
			static_assert(L == X || L == Y, "Axis parameter supplied must be 'X' or 'Y' (or 0,1) respectively.");
	}

	Eigen::Matrix2d A;
	Eigen::Vector2d dxy;

private: 
	/* Index indirection tables. */
	std::vector<ActiveIndex> active_x;
	std::vector<ActiveIndex> active_y;

	const FOOTPair* p; /* Reassigned, event-by-event. */
	mutable Cached cache;

}; // struct HitMatrix

struct RNFOOTPair;
struct FOOTHit;
//extern template struct HitMatrix<RNFOOTPair>;

/* No implementation file. The type is instantiated in the TFOOTHitProc.h. */
