#pragma once

#include "../monad/monad.hxx" // only used for assert stuff
#include <array>
#include <cmath>
#include <ostream>
#include "json_struct_def.hh"

/* This is only for a temporary container,
 * lines will always be saved in Cartesian coordinates. */

#define FORMAT_ANGLES_IN_RADIANS

namespace mnd {
namespace geom {

struct Point { 
	Point() = default;
	Point(double x_, double y_) : x(x_), y(y_) {}
	Point(const std::array<double, 2>& v) : x(v[0]), y(v[1]) {}

	double x, y; 
};
struct Point3D { 
	Point3D() = default;
	Point3D(double x_, double y_, double z_) : x(x_), y(y_), z(z_) {}
	Point3D(const std::array<double, 3>& v) : x(v[0]), y(v[1]), z(v[2]) {}

	double x, y, z; 
};

struct SphericalAngles {
	double theta, phi;
	friend std::ostream& operator<<(std::ostream& os, const SphericalAngles& rhs) {
		return os << "(θ: " << rhs.theta <<  ", φ: " << rhs.phi << ')';	
	}
};

// Nullable, we don't directly expose the value as a public field,
// due to nullability (maybe idiotic design by me, but who cares) 
struct Line2D {
	using Arr = std::array<double, 2>;
	using value_type = Arr::value_type;

	Line2D(value_type f1, value_type f2) : value{f1,f2} {}
	Line2D(const Arr& arr) : value(arr) {}
	
	inline value_type& operator[](Arr::size_type pos) noexcept { return value[pos]; } 
	inline value_type const& operator[](Arr::size_type pos) const noexcept { return value[pos]; } 

	/* Predicate to see if track is valid or null. */
	inline bool HasValue() const noexcept { return std::isfinite(value[0]); }

	Arr& array() noexcept { return value; }
	Arr const& array() const noexcept { return value; }

private:
	 Arr value {NAN,NAN};
};

struct Line3D {
	Line2D a, b;
	inline SphericalAngles Spherical() const noexcept {
		return { 
			std::acos( 1.0 / std::sqrt( a[1]*a[1] + b[1]*b[1] + 1) ), // theta
			std::atan2( b[1] , a[1] ) // phi
		};
	}

	friend std::ostream& operator<<(std::ostream& os, const Line3D& rhs) {
		return os << "(Lx: " << rhs.a.array() 
		          << " , Ly: " << rhs.b.array() 
				  << " :: " << rhs.Spherical() << ')';
	}

	inline bool HasValue() const noexcept { return a.HasValue(); }
	double AngleRelativeTo(const Line3D& ref) const noexcept; 
};

struct Rectangle2D {
	/* Spanned by bottom left and top right:
	 *   .----------------- (x1,y1)
	 *   |                     | 
	 *   |                     |
	 *   |                     |
	 * (x0,y0) ----------------^ */
	Rectangle2D() = default;

	/* Take the mid point coordinates, and width in x- and y- respectively. */ 
	Rectangle2D(Point midpoint, double wx, double wy);

	Point p0, p1;
	Point Mid() const noexcept;
	bool IsInside(const Point& ) const noexcept;
	bool IsInside(double, double ) const noexcept;
	
	friend std::ostream& operator<<(std::ostream& os, const Rectangle2D& r) {
		return os << "{R (" 
			<< r.p0.x << ',' << r.p0.y 
			<< ") <> (" 
			<< r.p1.x << ',' << r.p1.y 
			<< ") R}";
	}
};
} // namespace geom
} // namespace mnd

mnd::geom::Line2D GetLine(const mnd::geom::Point& , const mnd::geom::Point& ) noexcept;
mnd::geom::Line2D GetLine(const std::array<std::array<double, 2>, 2> ) noexcept;

