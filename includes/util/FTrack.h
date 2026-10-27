#pragma once

#include "monad/monad.hxx" // only used for assert+fmt stuff

#include <cmath>
#include <limits>
#include <type_traits>

#include "util/PolyFitter.h"
#include "util/HitMatrix.hxx"
#include "util/Geometry.h"
#include "util/DirectedAGraph.hxx"
#include "util/Option.hxx"

#include "../Eigen.h"

namespace mnd {
namespace track {

/* Nullable according to the mu field's special NAN. */
struct Q {
	double mu, v;

	inline double& mean() noexcept { return mu; }
	inline double const& mean() const noexcept { return mu; }
	inline double& var() noexcept { return v; }
	inline double const& var() const noexcept { return v; }

	inline double s() const noexcept { return std::sqrt( var() ); }

	friend std::ostream& operator<<(std::ostream& os, const Q& rhs) {
		return os << rhs.mean() << "±" << rhs.s();
	}
};

} // namespace track

namespace detail {

/* Hack for niche-optimized layout of optional type. */
struct TrackQ_ {
	using value_type = track::Q;
	inline constexpr static double nil = none_f64;
	inline constexpr static u64 nil_bits = none_bits_f64;

	constexpr TrackQ_() : val_{nil, 0.0} {}
	constexpr TrackQ_(None_t ) : val_{nil, 0.0} {}
	constexpr TrackQ_(value_type val) : val_(std::move(val)) {}

	constexpr bool has_value() const noexcept {
		return __builtin_bit_cast(u64, val_.mu) != nil_bits;
	}
	constexpr void reset() noexcept { val_.mu = nil; }

	MND_OPTION_ADD_UNCHECKED_ACCESS

private:
	value_type val_;
};

} // namespace detail

/* Ok: I might've overkilled it. Fuck it. Basically: by parking the `struct { double x,y; };`
 * behind an std::optional<> pads the struct to a 3x64 bits. That's why std::optional is cursed.
 * My own wrapper hacks this by allowing a custom provided underlying 'sum-type' which has to
 * only implement a few constexpr methods and expose a value_type and that's basically it.
 * A simple "glorified" NAN abuser. TBH, it's probably a negligible runtime gain, but is cool AF NGL.
 * Templates are actually so beautiful.. ／人 ⌒ ‿‿ ⌒ 人＼ */
template<>
class Option<track::Q> : public detail::OptionBase<Option<track::Q>, detail::TrackQ_> {
	using base = detail::OptionBase<Option<track::Q>, detail::TrackQ_>;
public:
	using base::base;
	using base::operator=;
};
static_assert(sizeof(Option<track::Q>) == sizeof(track::Q));

} // namespace mnd

/* A simple type that abstracts away storing some 'base' charge values and returning some kind of a combined charge measurement.
 * Keep it as a higher-order functor for now, to potentially enter some state inside the combinator later on.
 * E.g. that different FOOTs also have some other dependence here... */
template <
	size_t Capacity,
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

	const EntryType& operator[](size_t i) const noexcept {
#if defined(MND_HITMATRIX_DO_BOUNDS_CHECK)
		if(i >= N)
			ERROR("Requesting index: %zu, but currently combinator (%s) has %zu entries (capacity=%zu)",
				i, MND_SELF_TYPE_CSTR, N, Capacity);
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
		if(value.is_none()) {
			int nvalid = 0;
			double sum = 0;
			for(size_t i=0; i<N; ++i) {
				if(buffer[i]) {
					sum += buffer[i].qx + buffer[i].qy;
					++nvalid;
				}
			}
			if(nvalid == 0) {
				value = mnd::Some{ResultType{NAN, NAN}};
				return *value;
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
			value = mnd::Some {
				ResultType{mean, M2 / (2*nvalid - 1)}
			};
		}
		return *value;
	}

	double mean() const noexcept { return get().mean(); }
	double var() const noexcept { return get().var(); }
	friend std::ostream& operator<<(std::ostream& os, const QCombinator& Q) {
		return os << mnd::as_span(Q.buffer.data(), Q.N) << " => " << Q.get();
	}

private:
	size_t N = 0;
	std::array<EntryType, Capacity> buffer;
	mutable mnd::Option<ResultType> value {}; // lazily evaluated and cached.

}; // struct QCombinator

struct FTrack {
	mnd::geom::Line3D l;
	mnd::track::Q q;
	friend std::ostream& operator<<(std::ostream& , const FTrack& );
};

namespace mnd {
namespace detail {

/* Hack for niche-optimized layout of optional type. */
struct FTrack_ {
	using value_type = ::FTrack;
	inline constexpr static double nil = none_f64;
	inline constexpr static u64 nil_bits = none_bits_f64;

	constexpr FTrack_() : val_{} {
		val_.q.mu = none_f64;
	}
	constexpr FTrack_(None_t ) : FTrack_() {}
	constexpr FTrack_(value_type val) : val_(std::move(val)) {}

	constexpr bool has_value() const noexcept {
		return __builtin_bit_cast(u64, val_.q.mu) != nil_bits;
	}
	constexpr void reset() noexcept { val_.q.mu = nil; }
	
	MND_OPTION_ADD_UNCHECKED_ACCESS

private:
	value_type val_;
};
} // namespace detail

template<>
class Option<::FTrack> : public detail::OptionBase<Option<::FTrack>, detail::FTrack_> {
	using base = detail::OptionBase<Option<::FTrack>, detail::FTrack_>;
public:
	using base::base;
	using base::operator=;
};
static_assert(sizeof(Option<::FTrack>) == sizeof(::FTrack));

} // namespace mnd

/* Type encapsulating an 'online' track object to be handed over
 * to the Kalman filter. Is nullable! */
template <
	size_t Capacity,
	typename FOOTPair
> struct Track {
	using FHitMatrix = HitMatrix<FOOTPair>;
	using Entry = typename FHitMatrix::Entry;
	using QEntryType = typename Entry::q_type;	

	static constexpr double SP = 0.100 * 0.100; // [ mm^2 ]
	static constexpr double SQ = 0.5 * 0.5; // [ Q*Q ; 'Q' = one charge unit ]

	inline static constexpr auto BareTrack
		= FTrack {
			.l = mnd::geom::Line3D{NAN,NAN, NAN,NAN},
			.q = { .mu = NAN, .v = NAN }
		};

	// Containers
	std::array<double, Capacity> xs, ys, zs;
	QCombinator<Capacity, QEntryType, mnd::track::Q> q;

	inline size_t N() const noexcept { return q.size(); }

	/* If track currently holds 0 or 1 points, fill the undeterminable fit fields with ordinary qnan */
	const FTrack& get() const noexcept {
		if(t.is_none())
			evaluate_track();
		return *t;
	}

	void Add(const Entry& rhs, double z) noexcept {
		const size_t n_current = N();
#if defined(MND_HITMATRIX_DO_BOUNDS_CHECK)
		if(n_current >= Capacity)
			ERROR("Adding a point to '%s' object. Capacity %zu exceeded", MND_SELF_TYPE_CSTR, Capacity);
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

	mnd::geom::Point2D extrapolate_to(double z) const noexcept {
		const FTrack& fit = this->get();
		mnd::geom::Point2D val = fit.l.Eval(z);
		return val;
	}

	double GetScore() const noexcept {
		const size_t n_current = N();
		if(n_current <= 2)
			return std::numeric_limits<double>::infinity();
		
		double sp = 0, sq = 0;
		
		for(size_t i=0; i < n_current; ++i) {
			const mnd::geom::Point2D e = this->extrapolate_to( zs[i] );
			double dx = xs[i] - e.x;
			double dy = ys[i] - e.y;
			
			sp += dx*dx + dy*dy;

			double q0 = this->get().q.mean();
			double dqx = q[i].qx - q0;
			double dqy = q[i].qy - q0;
			
			sq += dqx*dqx + dqy*dqy;
		}
		return (sp / Track::SP + sq / Track::SQ) / (n_current - 2);
	}

	friend std::ostream& operator<<(std::ostream& os, const Track& rhs) {
		using mnd::fmt::operator<<;
		const size_t N = rhs.N();

		std::streamsize old_precision = os.precision();
		os << std::setprecision(4);
		os << "{| "
		   << "x: "   << mnd::span<const double>{rhs.xs.data(), N}
		   << ", y: " << mnd::span<const double>{rhs.ys.data(), N}
		   << ", z: " << mnd::span<const double>{rhs.zs.data(), N}
		   << ", q: " << rhs.q
		   << ", R: " << rhs.get() << " |}";
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
				t = mnd::Some{BareTrack};
				break;
			case 1:
				t = mnd::Some{BareTrack};
				(*t).q = q.get();
				break;
			default:
				t = mnd::Some<FTrack>{
					.value = {
						.l = { PolyFit<1>(zs, xs, N),
						       PolyFit<1>(zs, ys, N) },
						.q = q.get()
					}
				};
				break;
		}
	}

	// Result track. Calculated on demand.
	mutable mnd::Option<FTrack> t;

}; // struct Track
