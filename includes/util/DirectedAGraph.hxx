#pragma once

#include "../monad/monad.hxx"
#include "json_struct_def.hh"
#include <algorithm>
#include <type_traits>
#include <array>
#include <vector>
#include <ostream>

/* Small struct representing a nullable 2-index, used to memorize the "location" (i,j) of the entry used 
 * in the Kalman online + possible null location. */
template<typename T>
struct __attribute__((packed)) DAGIndex { // attribute is redundant here, but stays to convince the compiler for DAGPath packing.
	static_assert(std::is_fundamental_v<T>, "Underlying type must be fundamental.");
	static_assert(std::is_integral_v<T>, "Underlying type must be fundamental.");
	static_assert(std::is_unsigned_v<T>, "Underlying type must be unsigned.");
	
	constexpr static T NULL_VALUE = static_cast<T>(-1);

	using value_type = T;
	inline static std::array<T,2> const null = { static_cast<T>(NULL_VALUE), static_cast<T>(NULL_VALUE) };
	inline static std::array<T,2> const& Null() noexcept { return null; }

	inline T& operator[](std::size_t pos) noexcept { return i_[pos]; }
	inline T const& operator[](std::size_t pos) const noexcept { return i_[pos]; }

	DAGIndex(): i_( Null() ) {};
	explicit DAGIndex(T first, T second) : i_{first, second} {}
	
	DAGIndex(std::size_t first, std::size_t second) : i_{static_cast<T>(first), static_cast<T>(second)} {}

	explicit operator bool() const noexcept { return i_ != null; }
	explicit operator int() const noexcept  { return static_cast<int>( bool(*this) ); } // `0` or `1` 

	friend std::ostream& operator<<(std::ostream& os, const DAGIndex& i) { 
		if(!i) return os << "(∅)";
		
		os << '(';
		if(i[0] != NULL_VALUE) os << i[0];
		else os << "∅";
		os << ',';
		if(i[1] != NULL_VALUE) os << i[1];
		else os << "∅";
		os << ')';

		return os;
	}

private:
	std::array<T,2> i_;
};

/* Entire graph is packed in a local structure to preserve heap locality. */
template<typename IType, std::size_t Depth>
struct DirectedAGraph {
	constexpr static std::size_t depth = Depth;
	using Index = DAGIndex<IType>;
	
	struct /* __attribute__((packed)) */ DAGPath {
		std::array<Index, Depth> node {};
		
		friend std::ostream& operator<<(std::ostream& os, const DAGPath& p) { return os << p.node; }
		inline int Rank() const noexcept {
			return mnd::sum<mnd::Unroll::Yes, int>(node);
		}
	}; 

	DirectedAGraph() = default;

	// Called at the start, places a singular path with all nulls. */
	void Initialize() noexcept { 
		path.clear(); 
		path.emplace_back();
	}

	DirectedAGraph(std::size_t N) { path.reserve(N); };

	void TrimRankLessThan(int N) {
		/* [1] Partition the paths container, and put the "valid" paths as preceeding. */ 
		auto it = std::partition(path.begin(), path.end(),
			[N](const DAGPath& singular_path) {
				return singular_path.Rank() >= N;
			}
		);

		/* [2] Trim the container. */
		path.erase(it, path.end());
	}
	friend std::ostream& operator<<(std::ostream& os, const DirectedAGraph& g) { 
		return os << g.path;
	}
	std::vector<DAGPath> path; 
};


