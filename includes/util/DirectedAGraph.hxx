#pragma once

#include <type_traits>
#include <array>
#include <vector>
#include <ostream>

template<typename T>
struct __attribute__((packed)) DAGIndex { // attribute is redundant here, but stays to convince the compiler for DAGPath packing.
	static_assert(std::is_fundamental_v<T>, "Underlying type must be fundamental.");
	static_assert(std::is_integral_v<T>, "Underlying type must be fundamental.");
	static_assert(std::is_unsigned_v<T>, "Underlying type must be unsigned.");

	inline static std::array<T,2> const null = { static_cast<T>(-1), static_cast<T>(-1) };
	inline static std::array<T,2> const& Null() noexcept { return null; }

	inline T& operator[](std::size_t pos) noexcept { return i[pos]; }
	inline T const& operator[](std::size_t pos) const noexcept { return i[pos]; }

	DAGIndex(): i( Null() ) {};
	DAGIndex(T first, T second) : i{first, second} {}
	
	template<typename = typename std::enable_if<!std::is_same_v<T, std::size_t>>::type> // don't compile it in, if T == std::size_t
	DAGIndex(std::size_t first, std::size_t second) : i{static_cast<T>(first), static_cast<T>(second)} {}

	explicit operator bool() const noexcept { return i != null; }

	friend std::ostream& operator<<(std::ostream& os, const DAGIndex& i) { 
		return os << '(' << i[0] << ',' << i[1] << ')';
	}

private:
	std::array<T,2> i;
};

template<typename IType, std::size_t Depth>
struct DirectedAGraph {
	constexpr static std::size_t depth = Depth;
	using Index = DAGIndex<IType>;
	
	struct /* __attribute__((packed)) */ DAGPath {
		Index __node_dummy__;
		std::array<Index, Depth> node {};
	}; 

	DirectedAGraph() = default;

	// Called at the start, places a singular path with all nulls. */
	void Initialize() noexcept { 
		path.clear(); 
		path.emplace_back();
	}

	DirectedAGraph(std::size_t N) { path.reserve(N); };
	std::vector<DAGPath> path; 
};
