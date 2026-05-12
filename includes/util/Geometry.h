#pragma once
#include <array>
#include <cmath>

/* This is only for a temporary container,
 * lines will always be saved in Cartesian coordinates. */

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
};

using Line2D = std::array<double, 2>;
struct Line3D {
	std::array<double, 2> a, b;
	inline SphericalAngles Spherical() const noexcept {
		return { 
			std::acos( 1.0 / std::sqrt( a[1]*a[1] + b[1]*b[1] + 1) ), // theta
			std::atan2( b[1] , a[1] ) // phi
		};
	}

	double AngleRelativeTo(const Line3D& ref) const noexcept; 
};

struct Rectangle2D {
	// Spanned by bottom left and top right:
	//   .----------------- (x1,y1)
	//   |                     | 
	//   |                     |
	//   |                     |
	// (x0,y0) ----------------^
	Rectangle2D() = default;
	Rectangle2D(Point midpoint, double wx, double wy);

	Point p0, p1;
	Point Mid() const noexcept;
	bool IsInside(const Point& ) const noexcept;
	bool IsInside(double, double ) const noexcept;
};

} // namespace geom
} // namespace mnd

mnd::geom::Line2D GetLine(const mnd::geom::Point& , const mnd::geom::Point& ) noexcept;
mnd::geom::Line2D GetLine(const std::array<std::array<double, 2>, 2> ) noexcept;

