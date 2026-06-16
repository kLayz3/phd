#pragma once

#include "../monad/monad.hxx" // only used for assert stuff

#include "PolyFitter.hxx"
#include <cmath>
#include <limits>
#include <type_traits>

#include "json_struct_def.hh"
#include "HitMatrix.hxx"
#include "Geometry.h"
#include "DirectedAGraph.hxx"

#include "../Eigen/Dense"

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
	inline bool HasValue() const noexcept { return std::isfinite(mu); }
	inline explicit operator bool() const noexcept { return HasValue(); }

	friend std::ostream& operator<<(std::ostream& os, const Q& rhs) {
		return os << rhs.mean() << "±" << rhs.s();
	}
};

} // namespace track
} // namespace mnd

/* A simple function that combines two base values and returns a some kind of a combined charge measurement.
 * Keep it as a higher-order functor for now, to potentially enter some state inside the combinator later on. 
 * E.g. that different FOOTs also have some other dependence here... */
template <
	std::size_t Capacity,
	typename EntryType  = mnd::hm::Q,
	typename ResultType = mnd::track::Q
> struct QCombinator {
	QCombinator& operator+=(const EntryType& e) noexcept { 
		/* it's noexcept, since `ERROR(..)` simply calls std::abort */
		if(N == Capacity) 
			ERROR("Trying to add an entry but capacity %zu already reached.\n", Capacity);
		buffer[N++] = e;
		value.reset();
		return *this;
	}

	const EntryType& operator[](std::size_t i) const noexcept {
#if defined(MND_HITMATRIX_DO_BOUNDS_CHECK)
		if(i >= N) 
			ERROR("Requesting index: %zu, but currently combinator (%s) has %zu entries (capacity=%zu)",
				i, _SELF_TYPE_CSTR, N, Capacity);
#endif
		return buffer[i];
	}
	
	double size() const noexcept { return N; }
	
	EntryType* pop_back() noexcept { // no-op for N==0
		if(N == 0) return nullptr;
		
		// Reset statistics result only if we pop a non-null value.
		EntryType* last = &buffer[N-1];
		if(*last) value.reset(); // mnd::hm::Q::operator bool();
		--N;
		return last;
	}

	QCombinator& operator--() { pop_back(); return *this; }

	/* Result accesses are lazily calculated. We could have null values stashed in there. 
	 * Don't take them into account when calculating the statistics. */
	const ResultType& get() const noexcept {
		if(!value) {
			int nvalid = 0;
			double sum = 0;
			for(size_t i=0; i<N; ++i) {
				if(buffer[i]) {
					sum += buffer[i].qx + buffer[i].qy;
					++nvalid;
				}
			}
			const double mean = sum / (2*nvalid);
			double M2 = 0;
			for(size_t i=0; i<N; ++i) {
				if(buffer[i]) {
					double dx = buffer[i].qx - mean;
					double dy = buffer[i].qy - mean;
					M2 += dx*dx + dy*dy;
				}
			}
			value->mean() = mean;
			value->var() = M2 / (2*nvalid - 1);
		}
		return *value;
	}

	double mean() const noexcept { return get().mean(); }
	double var() const noexcept { return get().var(); }
	friend std::ostream& operator<<(std::ostream& os, const QCombinator& Q) {
		::mnd_output_homogeneous_range_(os, Q.buffer.data(), Q.N);
		return os << " => " << Q.get();	
	}

private:
	std::size_t N = 0;
	std::array<EntryType, Capacity> buffer;
	mutable std::optional<ResultType> value {}; // lazily evaluated and cached.

}; // struct QCombinator

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
} // namespace mnd

struct FTrack {
	mnd::geom::Line3D l;
	mnd::track::Q q;
	friend std::ostream& operator<<(std::ostream& , const FTrack& ); 
};

/* Type encapsulating an 'online' track object to be handed over 
 * to the Kalman filter. Is nullable! */
template <
	std::size_t Capacity,
	typename FOOTPair
>
struct Track {
	using FHitMatrix = HitMatrix<FOOTPair>;
	using Entry = typename FHitMatrix::Entry;
	using QEntryType = typename Entry::q_type;	
	using XYEntryType = typename Entry::xy_type;	

	static constexpr double SP = 0.100 * 0.100; // [ mm^2 ]
	static constexpr double SQ = 0.7 * 0.7; // [ Q*Q ; 'Q' = one charge unit ]

	inline static const auto BareTrack 
		= FTrack { 
			mnd::geom::Line3D{NAN,NAN, NAN,NAN},
			{}
		};

	// Containers
	std::array<double, Capacity> xs, ys, zs;
	QCombinator<Capacity, QEntryType, mnd::track::Q> q;
	inline size_t N() const noexcept { return q.size(); }

	/* If track currently holds 0 or 1 points, fill the undeterminable fit fields with NAN. */
	const FTrack& get() const noexcept { 
		if(!t) 
			evaluate_track(); 
		return *t;	
	}

	void Add(const Entry& rhs, double z) noexcept {
		const size_t n_current = N();
#if defined(MND_HITMATRIX_DO_BOUNDS_CHECK)
		if(n_current >= Capacity)
			ERROR("Adding a point to '%s' object. Capacity %zu exceeded", _SELF_TYPE_CSTR, Capacity);
#endif
		zs[n_current] = z;
		xs[n_current] = rhs.X();
		ys[n_current] = rhs.Y();
		q += rhs.q;
		t.reset();
	}

	void pop_back() {
		if(N() == 0) return;
		q.pop_back(); // will mutate N();
		t.reset();
	 }

	XYEntryType extrapolate_to(double z) const noexcept {
		const FTrack& fit = this->get();
		mnd::geom::Point2D val = fit.l.Eval(z);
		return { val.x, val.y }; /* Eigen::Vector2d */
	}

	double GetScore() const noexcept {
		const size_t n_current = N();
		if(n_current <= 2)
			return std::numeric_limits<double>::infinity();
		
		double sp = 0, sq = 0;
		
		for(size_t i=0; i < n_current; ++i) {
			const XYEntryType e = this->extrapolate_to( zs[i] );
			double dx = xs[i] - e(0);
			double dy = ys[i] - e(1);
			
			sp += dx*dx + dy*dy;

			double q0 = t->q.mean();
			double dqx = q[i].qx - q0;
			double dqy = q[i].qy - q0;
			
			sq += dqx*dqx + dqy*dqy;
		}
		return (sp / Track::SP + sq / Track::SQ) / (n_current - 2);	
	}

	enum class Status { Bare, SinglePoint, WellDefined };
	
	Status GetStatus() const noexcept {
		const size_t N = this->N();
		switch(N) {
			case 0: return Status::Bare;
			case 1: return Status::SinglePoint;
			default: return Status::WellDefined;
		}
	}

	friend std::ostream& operator<<(std::ostream& os, const Track& t) {
		const size_t N = t.N();

		std::streamsize old_precision = os.precision();
		os << std::setprecision(4);
		os << "{| ";
		os << "x: "; ::mnd_output_homogeneous_range_(os, t.xs.data(), N);
		os << ", y: "; ::mnd_output_homogeneous_range_(os, t.ys.data(), N);
		os << ", z: "; ::mnd_output_homogeneous_range_(os, t.zs.data(), N);
		os << ", q: " << t.q
		   << ", R: " << t.get() << " |}";
		os.precision(old_precision);
		return os;
	}

	/* Main API: a nullable track can be constructed 
	 * from the currently samples sequence of points.  */
private:
	void evaluate_track() const noexcept {
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
					{ PolyFit<1>(zs, xs, N),
					  PolyFit<1>(zs, ys, N) },
					q.get()
				};
				break;
		}
	}

	// Result track. Calculated on demand.
	mutable std::optional<FTrack> t; 
}; // struct Track
