#pragma once

#include "HitMatrix.hxx"
#include "util/PolyFitter.hxx"
#include <type_traits>
#include "../Eigen/Dense"
#include "Geometry.h"

namespace mnd {
namespace track {

struct Q {
	double mu, v;
	Q(): mu(NAN), v(NAN) {}
	Q(double mu_, double v_): mu(mu_), v(v_) {}

	inline double& mean() noexcept { return mu; }
	inline double const& mean() const noexcept { return mu; }
	inline double& var() noexcept { return v; }
	inline double const& var() const noexcept { return v; }

	inline double s() const noexcept { return std::sqrt( var() ); } 
};

} // namespace track
} // namespace mnd

/* A simple function that combines two base values and returns a some kind of a combined charge measurement.
 * Keep it as a higher-order functor for now, to potentially enter some state inside the combinator later on. 
 * E.g. that different FOOTs also have some other dependence here... */
template <
	typename EntryType  = mnd::hm::Q,
	typename ResultType = mnd::track::Q,
	std::size_t Capacity = 4
> struct QCombinator {
	QCombinator& operator+=(const EntryType& e) noexcept { 
		/* it's noexcept, since `ERROR(..)` simply calls std::abort */
		if(N == Capacity) 
			ERROR("Trying to add an entry but capacity %zu already reached.\n", Capacity);
		buffer[N++] = e;
		value.reset();
		return *this;
	}

	double size() const noexcept { return N; }
	void pop_back() noexcept { // no-op for N==0
		if(N != 0) {
			--N; value.reset();
		}
	}
	void pop_back_unchecked() noexcept { --N; value.reset(); };

	QCombinator& operator--() { pop_back(); return *this; }

	const ResultType& get() const noexcept {
		if(!value) {
			double sum = 0;
			for(size_t i=0; i<N; ++i)
				sum += buffer[i].qx + buffer[i].qy;
			const double mean = sum / (2*N);
			double M2 = 0;
			for(size_t i=0; i<N; ++i) {
				double dx = buffer[i].qx - mean;
				double dy = buffer[i].qy - mean;
				M2 += dx*dx + dy*dy;
			}
			value->mean() = mean;
			value->var() = M2 / (2*N - 1);
		}
		return *value;
	}

	double mean() const noexcept { return get().mean(); }
	double var() const noexcept { return get().var(); }

private:
	size_t N = 0;
	std::array<EntryType, Capacity> buffer;
	mutable std::optional<ResultType> value {}; // lazily evaluated and cached.
};

namespace mnd {
#define DECL_TYPE_TRAIT_HAS_STATIC(VAR) \
	template<typename, typename = void> \
	struct has_static_##VAR : std::false_type {}; \
	\
	template<typename T> \
	struct has_static_##VAR<T, std::void_t<decltype(T::VAR)>> \
		: std::true_type { \
			using underlying_type = std::remove_cv_t<decltype(T::VAR)>; \
		};

DECL_TYPE_TRAIT_HAS_STATIC(N_PAIRS);

}

template<typename QResultType = mnd::track::Q>
struct TrackTau {
	mnd::geom::Line3D l;
	QResultType q;

};
template<typename T> 
TrackTau(mnd::geom::Line3D, T) -> TrackTau<T>;

extern template struct TrackTau<>;

/* Type encapsulating an 'online' track object to be handed over 
 * to the Kalman filter. */
template <
	typename FOOTPair, 
	typename QResultType = mnd::track::Q,
	std::size_t Capacity = 4
>
struct Track {
	using FHitMatrix = HitMatrix<FOOTPair>;
	using Entry = typename FHitMatrix::Entry;
	using QEntryType = typename Entry::q_type;	
	using XYEntryType = typename Entry::xy_type;	

	 inline static const auto BareTrack 
		= TrackTau { 
			mnd::geom::Line3D{{NAN,NAN}, {NAN,NAN}},
			QResultType{}
		};

	std::array<double, Capacity> xs, ys, zs;
	QCombinator<QEntryType, QResultType, Capacity> q;

	inline size_t N() const noexcept { return q.size(); }

	/* If track currently holds 0 or 1 points, fill the undeterminable fit fields with NAN. */
	const TrackTau<QResultType>& get() const noexcept { 
		if(!t) 
			evaluate_track(); 
		return *t;	
	}

	void Add(const Entry& rhs, double z) {
		const size_t n_current = N();
		zs.at(n_current) = z;
		xs[n_current] = rhs.X();
		ys[n_current] = rhs.Y();
		q += rhs.q;
		t.reset();
	}

	void pop_back() {
		if(N() == 0) return;
		q.pop_back_unchecked(); // will mutate N();
		t.reset();
	}

	XYEntryType extrapolate_to(double z) const {
		const auto& fit = this->get();
		return {
			poly::Eval(z, fit.l.a),
			poly::Eval(z, fit.l.b) 
		}; /* Eigen::Vector2d */
	}

private:
	void evaluate_track() const {
		const size_t N = this->N();
		switch(N) {
			case 0: 
				t = BareTrack;
				break;
			case 1:
				t = BareTrack;
				t->q = q.get();
				break;
			default:
				t = {	
					PolyFit<1>(zs, xs, N),
					PolyFit<1>(zs, ys, N),
					q.get()
				};
				break;
		}
	}

	// Optional in the sense that maybe it's unevaluated at a specific point during runtime.
	mutable std::optional<TrackTau<QResultType>> t; 
};

template<typename T>
Track(T&& ) -> Track<T>;
