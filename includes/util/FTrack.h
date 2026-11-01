#pragma once

/* This TU serves to encapsulate API of "online tracking" for the CKF algorithm.
 * Basically, will encapsulate everything necessary for updating, caching and going back and forth
 * while tranversing the DAGs. */

#include "Eigen/Core"
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
	using value_type = double;

	value_type mu, v;

	inline value_type& mean() noexcept { return mu; }
	inline value_type const& mean() const noexcept { return mu; }
	inline value_type& var() noexcept { return v; }
	inline value_type const& var() const noexcept { return v; }

	inline value_type s() const noexcept { return std::sqrt( var() ); }

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

	inline constexpr bool has_value() const noexcept {
		return __builtin_bit_cast(u64, val_.mu) != nil_bits;
	}
	inline constexpr void reset() noexcept { val_.mu = nil; }

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

	inline constexpr bool has_value() const noexcept {
		return __builtin_bit_cast(u64, val_.q.mu) != nil_bits;
	}
	inline constexpr void reset() noexcept { val_.q.mu = nil; }
	
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

/* Static information about a FOOT layer, namely its geometry: zx,zy, θx, θy */
struct FLayerInfo {
	using Entry = mnd::hm::Data;
	using EntryXYUnderlyingType = Entry::xy_underlying_type; // double | long double

	using PredictionQ = mnd::track::Q::value_type;
	using PredictionP = Entry::xy_type;

	/* What the (q,u,v) of the following layer should be. */
	using Prediction = std::pair<PredictionQ, PredictionP>;

	/* Conversion matrix from a track vector normal repr: ν = (x0, ax, y0, ay)ᵀ into FOOT pair's measurement (u,v).
	 * Each detector is rotated by an angle θx, and placed at different positions: zu and zv.
	 * A = ( cosθx ,  zu*cosθx , sinθx , zu*sinθx
	 *      -sinθy , -zv*sinθy , cosθy , zv*cosθy ) .
	 * Then:  (u,v)ᵀ = Aν - dxy .
	 * Since this becomes a submatrix of a bigger Nx4 design problem's A matrix, keep these rows
	 * row major, then the bigger Nx4 can just construct itself (as a row major) trivially from a sequence of references
	 * to these smaller blocks and stack them together as rows. Also transpose(A) never really gets called, so
	 * SIMD likes the RowMajor alignment more. */
	Eigen::Matrix<EntryXYUnderlyingType, 2, 4, Eigen::RowMajor> A;
	Eigen::Vector<EntryXYUnderlyingType, 2> dxy; // Translation vector, each detector might be translated by dx/dy
	Eigen::Vector<EntryXYUnderlyingType, 2> z;   // Raw z's of the two detectors in the layer.

	/* This doesn't touch or evaluate any of the underlying mutable cache entries. So it is perfectly const.
	 * Evaluate the ftrack vector normal repr: ν = (x0, ax, y0, ay)ᵀ, and charge combinator's mean from
	 * the sequence (q1,q2,..) into "prediction" for this layers' measurement (Q,(u,v)). */
	inline Prediction prediction(const FTrack& ft) const noexcept {
		const auto nl = ft.l.normal();
		ft.q.mean();
		static_assert(mnd::size<decltype(nl)>() == 4, "Normal representation must be statically sized with 4 elements.");
		Eigen::Map<const Eigen::Vector4d> mapv = Eigen::Map<const Eigen::Vector4d>{ nl.data() };

		return {
			ft.q.mean(),
			A*mapv - dxy
		};
	}
};

/* The sequence of FLayerInfo objects corresponding to subsequent layers, builds up design matrices which are
 * used in the "heaviest" part of the Kalman algorithm. Nice thing is that
 * all these matrices can be mostly calculated up front and cached up. */
template<size_t N, size_t K>
struct DesignMatrix {
	using Matrix = Eigen::Matrix<
		typename FLayerInfo::EntryXYUnderlyingType,
		2*K, 4
	>;
	using QR = Eigen::ColPivHouseholderQR<Matrix>;
	using Vector = Eigen::Vector<
		typename FLayerInfo::EntryXYUnderlyingType,
		2*K
	>;

	constexpr static size_t n_combinations = mnd::combi::n_choose_k(N,K);
	constexpr static auto index_lookup_table = /* size is 2 pow k */
		mnd::combi::combo_index_table<N,K>;

	std::array<
		std::pair<QR,Vector>, n_combinations
	> lookup_table;

	constexpr decltype(auto) get_from_mask(u32 mask) const noexcept {
		if(mask >= mnd::size(index_lookup_table))
			ERROR("DesignMatrix::get_from_mask() : invalid access: %u > %zu.", mask, mnd::size(index_lookup_table));
		u32 retval =  index_lookup_table[mask];
		if(retval == mnd::combi::invalid)
			ERROR("DesignMatrix::get_from_mask() : invalid combination 0x%08x yields index -1.", mask);
		
		return lookup_table[retval];
	}
	/* Does not check for any bounds. */
	constexpr decltype(auto) get_from_mask_unchecked(u32 mask) const noexcept {
		return lookup_table[index_lookup_table[mask]];
	}
};
template<size_t N> struct DesignMatrix<N,0> {};
template<size_t N> struct DesignMatrix<N,1> {};

/* First two tuple entries <0> and <1> are a simple unit structs. */
template<size_t N>
struct DesignMatrixCache {
	static_assert(N >= 2);
	constexpr static size_t n_layers = N;

private:
	template<size_t... Is>
	static auto tuple_type(std::index_sequence<Is...> ) -> std::tuple<DesignMatrix<N, Is>...>;

	using Storage = decltype( tuple_type(std::make_index_sequence<N+1>{}) );

    Storage subsets;

public:
	template<size_t K>
	constexpr decltype(auto) get() noexcept {
		static_assert(K >= 2 && K <= N);
		return std::get<K>(subsets);
	}

	template<size_t K>
	constexpr decltype(auto) get() const noexcept {
		static_assert(K >= 2 && K <= N);
		return std::get<K>(subsets);
	}
	constexpr size_t size() const noexcept { return N; }
};

template<size_t Nlayers>
struct TrackingGeometry {
	std::array<FLayerInfo, Nlayers> layer;
	DesignMatrixCache<Nlayers> fit;
};

/* Type encapsulating an 'online' track object to be handed over
 * to the Kalman filter. Is nullable! */
template <
	size_t Capacity,
	typename FOOTPair
> struct Track {
	using FHitMatrix     = HitMatrix<FOOTPair>;
	using FHitMatrixCont = std::array<FHitMatrix, Capacity>;

	using Entry = typename FHitMatrix::Entry;
	using entry_xy_type = typename Entry::xy_type;
	using entry_q_type  = typename Entry::q_type;
	using hm_key_type = u8;

	using Point2D = mnd::geom::Point2D;
	using Line3D = mnd::geom::Line3D;
	
	static constexpr double SP = 0.100 * 0.100; // [ mm^2 ]
	static constexpr double SQ = 0.500 * 0.500; // [ Q^2 ; [Q] = one charge unit ]

	static constexpr auto BareTrack
		= FTrack {
			.l = Line3D{NAN,NAN, NAN,NAN},
			.q = { .mu = NAN, .v = NAN }
		};

	/* ======================== FIELDS ======================== */
	// Sequence of detector position measurements { (u1,v1), (u2,v2), ..., (uN,vN) }
	std::array<entry_xy_type, Capacity> m;
	
	// Sequence of charge measurements {q1, q2, ... qN}
	QCombinator<Capacity, entry_q_type, mnd::track::Q> q;

	/* Hit-matrices of layers are saved in a lookup container, e.g. std::array<FHitMatrix, N_PAIRS>
	 * Here keep the keys (indices) of the sequence of layer indices that constituted the track.
	 * This also serves to quickly look-up the correct design matrix via a combination hash
	 * that gets involved with the track evaluation. */
	std::array<hm_key_type, Capacity> layer_i;

	/* ======================== METHODS ======================== */
	size_t N() const noexcept { return q.size(); }
	
	u32 layer_mask() const noexcept {
		const size_t n_current = N();
		u32 mask = 0x0;
		for(size_t i=0; i < n_current; ++i) {
			mask |= (u32{1} << layer_i[i]);
		}
		return mask;
	}

	/* If track currently holds 0 or 1 points, fill the undeterminable fit fields with ordinary qnan */
	const FTrack& get(const TrackingGeometry<Capacity>& geom) const noexcept {
		if(t.is_none())
			evaluate_track(geom);
		return *t;
	}

	void Add(const Entry& entry, u32 layer) noexcept {
		const size_t n_current = N();
#if defined(MND_HITMATRIX_DO_BOUNDS_CHECK)
		if(n_current >= Capacity)
			ERROR("Adding a point to '%s' object. Capacity %zu exceeded", MND_SELF_TYPE_CSTR, Capacity);
#endif
		layer_i[n_current] = static_cast<hm_key_type>(layer);
		m[n_current] = entry.pos;
		q += entry.q; // increments n_current

		t.reset();
	}

	void pop_back() {
		if(N() == 0) return;
		q.pop_back(); // will mutate N();
		t.reset();
	}
	void pop_back(FTrack&& replacement) {
		if(N() != 0)
			q.pop_back(); // will mutate N();
		t = mnd::Some{std::move(replacement)};
	}

	Point2D extrapolate_to(const TrackingGeometry<Capacity>& geom, double z) const noexcept {
		const FTrack& fit = this->get(geom);
		Point2D val = fit.l.Eval(z);
		return val;
	}

	double GetScore(const TrackingGeometry<Capacity>& geom) const noexcept {
		const size_t n_current = N();
		if(n_current <= 2)
			return std::numeric_limits<double>::infinity();
		
		double sp = 0, sq = 0;
		const FTrack& ftrack = get(geom);
		
		for(size_t i=0; i < n_current; ++i) {
			const auto predicted = /* std::pair<PredictionQ, PredictionP> */
				geom.layer[ layer_i[i]  ].prediction(ftrack);

			const double dr2 = (predicted.second - m[i]).squaredNorm();
			sp += dr2;

			double dqx = q[i].qx - predicted.first;
			double dqy = q[i].qy - predicted.first;
			
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
		   << "(u,v): " << mnd::span<const entry_xy_type>{rhs.m.data(), N}
		   << ", q: " << rhs.q
		   << ", i: " << mnd::span<const hm_key_type>{rhs.layer_i.data(), N}
			<< " |}";
		os.precision(old_precision);
		return os;
	}

private:
	/* Main API: a nullable track can be constructed
	 * from the currently samples sequence of points.
	 * This call is "heavy" since it must solve the linear design problem.
	 * The prediction matrix which converts a normal-repr track: (x0,ax,y0,ay) to measurements:
	 * A := ( cosθx_0     ,  zu*cosθx_0     , sinθx_0     , zu*sinθx_0
	 *       -sinθy_0     , -zv*sinθy_0     , cosθy_0     , zv*cosθy_0
	 *        cosθx_1     ,  zu*cosθx_1     , sinθx_1     , zu*sinθx_1
	 *       -sinθy_1     , -zv*sinθy_1     , cosθy_1     , zv*cosθy_1
	 *           ⋮                 ⋮             ⋮                ⋮
	 *        cosθx_(N-1) ,  zu*cosθx_(N-1) , sinθx_(N-1) , zu*sinθx_(N-1)
	 *       -sinθy_(N-1) , -zv*sinθy_(N-1) , cosθy_(N-1) , zv*cosθy_(N-1)
	 *      )
	 * m := ( flatten(m[0]) , flat(m[1]) , ... flat(m[N-1]) )
	 * c := ( flatten(hm[0]->dxy) , flat(hm[1]->dxy) , ... flat(hm[N-1]) )
	 *
	 * Then solve for ν = (x0,ax,y0,ay)
	 * (AᵀA)ν = Aᵀ(m + c)
	 * Due to caching, a Track<> object itself *cannot* construct its own FTrack, it needs some
	 * parameters that are then cached and look-up by a specific bitmask key from a DesignMatrixCache. */
	void evaluate_track(const TrackingGeometry<Capacity>& geom) const noexcept {
		using value_type = typename Entry::xy_underlying_type;

		const size_t N = this->N();
		switch(N) {
			case 0:
				t = mnd::Some{BareTrack};
				break;
			case 1:
				t = mnd::Some{BareTrack};
				(*t).q = q.get();
				break;
			case 2 ... Capacity: {
				const u32 mask = layer_mask();
			
				mnd::static_for<2, Capacity + 1>([&](auto K_) {
					constexpr size_t K = decltype(K_)::value;
					if(K != N) return;

					const auto& design_matrix = geom.fit.template get<K>();
					const auto& parameters = /* std::pair<QR, Eigen::Vector> */
				#ifdef MND_DEBUG_ENABLED
						design_matrix.get_from_mask(mask);
				#else
						design_matrix.get_from_mask_unchecked(mask);
				#endif
					const auto& qr = parameters.first;
					const auto& c  = parameters.second;
					
					Eigen::Vector<value_type, 2*K> rhs;
					for(size_t i = 0; i < K; ++i) {
						rhs.template segment<2>(2 * i) =
							m[i] + c.template segment<2>(2 * i);
					}

					const Eigen::Vector<value_type,4> result = qr.solve(rhs);
					
					t = mnd::Some{ FTrack{
						.l = { result(0), result(1), result(2), result(3) },
						.q = q.get()
					}};
				});
				break;
			}
			default:
				mnd::unreachable();
		}
	}

	// Result track. Calculated on demand.
	mutable mnd::Option<FTrack> t;

}; // struct Track
