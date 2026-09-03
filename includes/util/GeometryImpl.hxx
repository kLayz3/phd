#pragma once

/* Small implementation file for the massive templates defined in the header. */

#include "Geometry.h" // Guarded
#include "Combinatorics.hxx"

namespace mnd::geom {

namespace detail {

namespace traits {

/* Disentangle when ref gets returned via deref op vs. a brand new copy. */
template<typename T>
using deref_t = decltype(*std::declval<const T&>());

template<typename T, typename = void>
struct derefs_to_line3d_ref : std::false_type {};
template<typename T>
struct derefs_to_line3d_ref <
	T,
	std::void_t<deref_t<T>>
> : std::bool_constant <
	std::is_lvalue_reference_v<deref_t<T>> &&
	std::is_same_v<
		::mnd::remove_cvref_t<deref_t<T>>,
		Line3D
	>
> {};
template<typename T>
inline constexpr bool derefs_to_line3d_ref_v =
	derefs_to_line3d_ref<T>::value;

template<typename T, typename = void>
struct derefs_to_line3d_val : std::false_type {};
template<typename T>
struct derefs_to_line3d_val<
	T,
	std::void_t<deref_t<T>>
> : std::bool_constant <
	!std::is_reference_v<deref_t<T>> &&
	std::is_same_v<
		::mnd::remove_cvref_t<deref_t<T>>,
		Line3D
	>
> {};
template<typename T>
inline constexpr bool derefs_to_line3d_val_v =
	derefs_to_line3d_val<T>::value;

} // namespace traits

template<typename, size_t, typename = void>
struct Line3DCache;

template<typename T, size_t N>
struct Line3DCache<
	T,
	N,
	std::enable_if_t<traits::derefs_to_line3d_ref_v<T>>
> {
	Line3DCache() = delete;
	Line3DCache( ::mnd::span<const T> input) {
		/* This ctor is guarded up-front by `input.size() == N`
		 * We trust it. This class anyway should not be part of public API. */
		for(size_t i=0; i<N; ++i)
			_lines[i] = &*input[i];
	}

	const Line3D* operator[](u32 n) const noexcept {
		return _lines[n];
	}

private:
	std::array<const Line3D*, N> _lines;
};

template<typename T, size_t N>
struct Line3DCache<
	T,
	N,
	std::enable_if_t<traits::derefs_to_line3d_val_v<T>>
> {
	Line3DCache() = delete;
	Line3DCache( ::mnd::span<const T> input) {
		/* This ctor is guarded up-front by `input.size() == N`
		 * We trust it. This class anyway should not be part of public API. */
		for(size_t i=0; i<N; ++i)
			_lines[i] = std::move(*input[i]);
	}

	const Line3D* operator[](u32 n) const noexcept {
		return &_lines[n];
	}

private:
	std::array<Line3D, N> _lines;
};

constexpr bool SANITY_CHECK_COMBINATORICS = 1;

struct Candidate {
	double score = std::numeric_limits<double>::infinity();
	u32 n_items = 0;
	u32 bitmask = 0;
	Point3D vertex = Point3D::null;

	bool has_value() const noexcept {
		return n_items >= 2;
	}
};

/* Basically since we're limited to around N=2,..7 tracks, just bruteforcing all the combinations is usually
 * the fastest. RANSAC not really needed. Greedy regression where we remove the farthest outlier isn't really
 * stable.
 * Bruteforcing all the combinations, I want to precompute what the available combinations actually are 
 * as a compile-time lookup table.
 * Thx God we work in C++ which is a human language, and this is hackable. :-) */
template <
	typename T,
	u32 N
> static VertexingResult<T> try_solve(
	::mnd::span<const T> lines,
	double const D
) { /* It is asserted that `N == lines.size()` */
	if constexpr(SANITY_CHECK_COMBINATORICS) {
        assert(N == lines.size() && "Paranoia combinatorics (0) hehe.");
    }
	Candidate best_candidate {};

	/* We don't know a priori if `T` directly derefs to `const Line3D&` or `Line3D`.
	 * In our case, simply taking &*lines[i] and saving that pointer might dangle on the spot. */
	detail::Line3DCache<T,N> buffer{lines};
	std::array<Line3D const*, MAX_VERTEXING_MULTP> tmp;
	
	mnd::static_for<N,1>( [&](auto _K) { // N, N-1, ... 2
		constexpr std::size_t K = decltype(_K)::value;
		
		/* If the candidate is already assigned then don't descent down to smaller combinations. */
		if(best_candidate.has_value())
			return;

		constexpr auto& table = mnd::combi::combo_lookup_table<N,K>; // std::array<u32, N_choose_K>

		for(const u32 bitmask : table) { // Should I unroll this? N=7,M=3/4 gives already 35 combos!
			u32 n = 0, m = bitmask;
			while(m) {
				/* Count trailing zeroes. The first bit '1' is at exactly this position. */
				const int i = __builtin_ctz(m);
				tmp[n++] = buffer[i]; // detail::Line3DCache::operator[]

				/* Reset the lowest bit in the bitmask.
				 * https://graphics.stanford.edu/~seander/bithacks.html#CountBitsSetKernighan */
				m &= m - 1;
			}
			if constexpr(SANITY_CHECK_COMBINATORICS) {
				assert(n == (u32)K && "Paranoia combinatorics (1) hehe.");
			}

			mnd::span<Line3D const*> selected {tmp.data(), (size_t)n};
			const Point3D vertex = FindVertex(selected);
			
			double score = 0.0;
			bool reject = false;
			
			/* The first candidate which dies vetoes this specific combination.
			 * This occasion we also use to cache-in the distances calculated. */
			for(const Line3D* line : selected) {
				const double d = line->DistanceTo(vertex);
				if(!(d < D)) { // reject NAN's too!
					reject = true;
					break;
				}
				/* Score of a valid candidate combination is simply Σd(t(i), v)^2 */
				score += d*d;
			}
			
			if(!reject && score < best_candidate.score) {
				best_candidate = {
					.score = score,
					.n_items = K,
					.bitmask = bitmask,
					.vertex = vertex
				};
			}
		}
	}); // static_for
	
	VertexingResult<T> rv;

	if(best_candidate.has_value()) {
		rv.tracks.reserve( best_candidate.n_items );
		
		/* Extract according to be bitmask, this will keep the ordering! */
		u32 m = best_candidate.bitmask;
		while(m) {
			const int i = __builtin_ctz(m);
			rv.tracks.push_back( lines[i] ); // push_back(T& )
			m &= m - 1;
		}
		rv.vertex  = best_candidate.vertex;
		rv.score   = best_candidate.score;
		rv.bitmask = best_candidate.bitmask;
	}
	return rv;
}

} // namespace detail


template<typename T>
VertexingResult<T> FindVertexingTracks(::mnd::span<const T> lines, double const D) {
	switch(lines.size()) {
		case 2: return detail::try_solve<T,2>(lines, D);
		case 3: return detail::try_solve<T,3>(lines, D);
		case 4: return detail::try_solve<T,4>(lines, D);
		case 5: return detail::try_solve<T,5>(lines, D);
		case 6: return detail::try_solve<T,6>(lines, D);
		case 7: return detail::try_solve<T,7>(lines, D);
		default: return {};
	}
}

template<typename T>
VertexingResult<T> FindVertexingTracksMut(std::vector<T>& lines, double const D) {
	auto result = FindVertexingTracks(
		::mnd::span<const T>{ lines.data(), lines.size() },
		D
	);
	if(result.valid()) {
		u32 i = 0;

		/* Kick out all tracks whose indices are marked by the bitmask.
		 * mnd::Erase anyway goes element by element and simply checks it.
		 * The actual resizing of vector is done *after* the full pass of the predicate. :) */
		mnd::Erase(
			lines,
			[&](auto const& ) {
				return result.bitmask & (1ULL << i++);
			}
		);
	}
	return result;
}

} // namespace mnd::geom

struct RNFOOTTrack;

extern template
mnd::geom::VertexingResult<RNFOOTTrack>
mnd::geom::FindVertexingTracksMut(std::vector<RNFOOTTrack> &, double );

/* Phew. */
