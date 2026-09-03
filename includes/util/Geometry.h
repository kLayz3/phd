#pragma once

#include <array>
#include <cfloat>
#include <cmath>
#include <limits>
#include <ostream>
#include <type_traits>
#include <vector>
#include "json_struct_def.hh" // std::ostream& operator<<(std::ostream&, array<T,N> const&)
#include "../monad/monad.hxx"
#include "../Eigen.h"

#define MND_FORMAT_ANGLES_IN_RADIANS

namespace mnd::geom {

struct Point2D {
	static const Point2D null;
	inline bool is_null() const noexcept { return !std::isfinite(x); }
	
	/* Returns quiet NAN if either of the points is null. */
	inline double Distance2(const Point2D& rhs) const noexcept {
		const double dx = x - rhs.x;
		const double dy = y - rhs.y;
		return dx*dx + dy*dy;
	}
	inline double Distance(const Point2D& rhs) const noexcept {
		return std::sqrt( Distance2(rhs) );
	}
	Eigen::Map<Eigen::Vector2d> eigen_view() & noexcept { return Eigen::Map<Eigen::Vector2d>{&x}; }
	Eigen::Map<const Eigen::Vector2d> eigen_view() const& noexcept { return Eigen::Map<const Eigen::Vector2d>{&x}; }

	double x, y;
};

struct Point3D {
	static const Point3D null;
	inline bool is_null() const noexcept { return !std::isfinite(x); }
	
	/* Returns quiet NAN if either of the points is null. */
	inline double Distance2(const Point3D& rhs) const noexcept {
		const double dx = x - rhs.x;
		const double dy = y - rhs.y;
		const double dz = z - rhs.z;
		return dx*dx + dy*dy + dz*dz;
	}
	inline double Distance(const Point3D& rhs) const noexcept {
		return std::sqrt( Distance2(rhs) );
	}
	inline Eigen::Map<Eigen::Vector3d> eigen_view() & noexcept { return Eigen::Map<Eigen::Vector3d>{&x}; }
	inline Eigen::Map<const Eigen::Vector3d> eigen_view() const& noexcept { return Eigen::Map<const Eigen::Vector3d>{&x}; }

	double x, y, z;
};

// Conventiently, a vector shall just be stored as a point.
using Vector2D = Point2D;
using Vector3D = Point3D;

/* Usual way to represent lines in 2D:
 * x(z) = a[1]*z + a[0]
 */
struct Line2D {
	using Arr = std::array<double, 2>;
	using value_type = Arr::value_type;

	Line2D() = default;
	Line2D(value_type f1, value_type f2) : value{f1,f2} {}
	Line2D(Arr const& arr) : value(arr) {}
	
	static const Line2D null;

	inline value_type& operator[](Arr::size_type pos) noexcept { return value[pos]; }
	inline value_type const& operator[](Arr::size_type pos) const noexcept { return value[pos]; }

	/* Predicate to see if the line is valid or null. */
	inline bool HasValue() const noexcept { return std::isfinite(value[0]); }

	/* Translate the line by an offset. */
	Line2D& operator+=(double ) noexcept; // only along x-axis direction
	Line2D& operator+=(const Vector2D& ) noexcept;

	/* Represent the line in a coordinate system made by an offset. */
	Line2D& operator%=(double ) noexcept; // only along x-axis direction
	Line2D& operator%=(const Vector2D& ) noexcept;
	
	/* Rotate the line by an an angle (amount) */
	Line2D& Rotate(double ) noexcept;

	/* Represent the line in a coordinate system made by a unary rotation,
	 * which is represented by an angle (amount). */
	Line2D& RepresentInRotated(double ) noexcept;

	/* Offset the line in a coordinate system and rotate by a unary rotation. */
	Line2D& ShiftAndRotate(const Vector2D& , double ) noexcept;

	/* Represent the line in a coordinate system by an offset and unary rotation. */
	Line2D& RepresentInShifted(const Vector2D& , double ) noexcept;

	Arr& array() noexcept { return value; }
	Arr const& array() const noexcept { return value; }

	inline double Eval(double z) const noexcept { return value[0] + value[1]*z; }

	std::array<double, 2> value {NAN,NAN}; // {a0, a1} = {offset, slope}
};

struct SphericalAngles {
	double theta, phi;
	friend std::ostream& operator<<(std::ostream& , const SphericalAngles& );
};

/* Two ways to represent lines in 3D:
 * Usual experiment:
 * x(z) = a[0] + a[1]*z
 * y(z) = b[0] + b[1]*z
 * ... and also:
 * vec{r} = vec{p} + λ * vec{v}
 */
struct Line3D {
	static constexpr double ALMOST_PARALLEL  = 1e-4;
	static constexpr double ALMOST_PARALLEL2 = ALMOST_PARALLEL*ALMOST_PARALLEL;

	Line3D() = default;
	Line3D (
		double a0_, /* offset x:z */
		double a1_, /* slope  x:z */
		double b0_, /* offset y:z */
		double b1_  /* slope  y:z */
	) noexcept : p{a0_, b0_, 0}, v{a1_, b1_, 1} {}

	Line3D(std::array<double, 2> const& a_, std::array<double, 2> const& b_) :
		p{a_[0], b_[0], 0}, v{a_[1], b_[1], 1} {}

	explicit Line3D(const Line2D& a_, const Line2D& b_ ) :
		p{a_[0], b_[0], 0}, v{a_[1], b_[1], 1} {}

	static const Line3D null;
	inline bool HasValue() const noexcept { return std::isfinite( p[0] ); }
	
	double DistanceTo(const Point3D& ) const noexcept;
	double DistanceTo(const Line3D& ) const noexcept;
	double AngleRelativeTo(const Line3D& ) const noexcept;

	/* Translate the line by an offset. */
	Line3D& operator+=(double ) noexcept; // only along z-axis direction
	Line3D& operator+=(const Vector3D& ) noexcept;

	/* Represent the line in a coordinate system made by an offset. */
	Line3D& operator%=(double ) noexcept; // only along z-axis direction
	Line3D& operator%=(const Vector3D& ) noexcept;
	
	/* Rotate the line by an unary rotation, represented by a unit vector and an angle (amount) */
	Line3D& Rotate(const Vector3D& , double ) noexcept;

	/* Represent the line in a coordinate system made by a unary rotation,
	 * which is represented by a unit vector and an angle (amount). */
	Line3D& RepresentInRotated(const Vector3D& , double ) noexcept;

	/* Offset the line in a coordinate system and rotate by a unary rotation. */
	Line3D& ShiftAndRotate(const Vector3D& , const Vector3D&, double ) noexcept;

	/* Represent the line in a coordinate system by an offset and unary rotation. */
	Line3D& RepresentInShifted(const Vector3D& , const Vector3D&, double ) noexcept;

	inline Point2D Eval(double z) const noexcept {
		if(!HasValue() || std::abs(v[2]) < 1e-20)
			return {NAN, NAN};
		return {
			p[0] + v[0]/v[2] * (z - p[2]), // x
			p[1] + v[1]/v[2] * (z - p[2])  // y
		};
	}

	[[ nodiscard ]] inline Line2D XLine() const noexcept {
		if(!HasValue() || std::abs(v[2]) < 1e-20) {
			return {NAN, NAN};
		}
		const double a1 = v[0] / v[2];
		const double a0 = p[0] - a1 * p[2];
		return {a0, a1};
	}
	[[ nodiscard ]] inline Line2D YLine() const noexcept {
		if(!HasValue() || std::abs(v[2]) < 1e-20) {
			return {NAN, NAN};
		}
		const double b1 = v[1] / v[2];
		const double b0 = p[1] - b1 * p[2];
		return {b0, b1};
	}

	SphericalAngles GetSpherical() const noexcept;
	[[ nodiscard ]] inline std::array<double,2> xarray() const noexcept { return XLine().array(); }
	[[ nodiscard ]] inline std::array<double,2> yarray() const noexcept { return YLine().array(); }

	std::array<double, 3> p {NAN,NAN,NAN}; // {a0, b0, 0} , nominally
	std::array<double, 3> v {NAN,NAN,NAN}; // {a1, b1, 1} , nominally

	/* Is a simple no-op. Needed to vertexing API later. */
	inline Line3D&       operator*()       noexcept { return *this; }
	inline Line3D const& operator*() const noexcept { return *this; }

}; // Line3D
bool operator==(const Line3D& , const Line3D& ) noexcept;

/* Can't share the same symbol, as it will confuse the overload resolver:
 * `GetLine(<brace-enclosed initializer list>, <brace-enclosed initializer list>)` is ambiguous */
Line2D GetLine2D(const Point2D& , const Point2D& ) noexcept;
Line3D GetLine3D(const Point3D& , const Point3D& ) noexcept;

struct Rectangle2D {
	/* Spanned by bottom left and top right:
	 *   .----------------- (x1,y1)
	 *   |                     |
	 *   |                     |
	 *   |                     |
	 * (x0,y0) ----------------^ */
	Rectangle2D() = default;

	/* Take the mid point coordinates, and width in x- and y- respectively. */
	Rectangle2D(Point2D midpoint, double wx, double wy);

	Point2D p0, p1;
	Point2D Mid() const noexcept;
	bool IsInside(const Point2D& ) const noexcept;
	bool IsInside(double, double ) const noexcept;
	
	friend std::ostream& operator<<(std::ostream& os, const Rectangle2D& r) {
		return os << "{R ("
			<< r.p0.x << ',' << r.p0.y
			<< ") <> ("
			<< r.p1.x << ',' << r.p1.y
			<< ") R}";
	}
};

/* Find a point with a minimal distance to a sequence of lines (aka: a vertex). */
Point3D FindVertex(::mnd::span<const Line3D>  ) noexcept;
Point3D FindVertex(::mnd::span<const Line3D*> ) noexcept;
Point3D FindVertex(const Line3D&, const Line3D& ) noexcept;

constexpr inline double VERTEXING_MIN_DISTANCE = 1.0;
constexpr inline size_t MAX_VERTEXING_MULTP    = 7;

template<typename T = Line3D>
struct VertexingResult {
	std::vector<T> tracks {};
	Point3D vertex = Point3D::null;
	double score = std::numeric_limits<double>::infinity();
	u64 bitmask = 0x0; // 8-byte cuz anyway padded to the boundary.

	bool valid() const noexcept { return tracks.size() >= 2; }
};

/* In a series of tracks: {t0, t1,... tN}, find the largest subset {𝜏0, 𝜏1, .. 𝜏M) which forms
 * a vertex with at most `D` width. E.g. the condition: d(vertex, 𝜏(i)) < D, ∀i where i<M must hold.
 * This is generic over any type `T` that dereferences into a `Line3D` or `Line3D&`. The deref operator can
 * either return a copy or some internal reference. If it returns a reference, then
 * the returned reference's lifetime MUST match the object's lifetime. */
template<typename T = Line3D>
VertexingResult<T> FindVertexingTracks(::mnd::span<const T> , double const D = VERTEXING_MIN_DISTANCE);

/* Same as FindVertexingTracks, except we explicitly mutate the input vector. Extra requirement is that
 * the type `T` needs an explicit comparison operator.
 * We don't mutate the individual objects, just kick them out of the vector, and we also
 * don't reorder the vector! */
template<typename T = Line3D>
VertexingResult<T> FindVertexingTracksMut(std::vector<T>& , double const D = VERTEXING_MIN_DISTANCE);

} /* namespace mnd::geom */

std::ostream& operator<<(std::ostream&, const mnd::geom::Point2D& );
std::ostream& operator<<(std::ostream&, const mnd::geom::Point3D& );
std::ostream& operator<<(std::ostream&, const mnd::geom::Line2D& );
std::ostream& operator<<(std::ostream&, const mnd::geom::Line3D& );

#include "GeometryImpl.hxx"
