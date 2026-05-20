#include "Geometry.h"
#include "../Eigen/Dense"

using namespace mnd::geom;

double Line3D::AngleRelativeTo(const Line3D& ref) const noexcept {
	using Vec3 = Eigen::Vector3d;
	Vec3 v1( this->a[1], this->b[1], 1 );
	Vec3 vref( ref.a[1], ref.b[1],   1 );

	double cosTheta = v1.normalized().dot( vref.normalized() );
	cosTheta = std::clamp(cosTheta, -1.0, 1.0);
	return std::acos(cosTheta); // radians
}

Rectangle2D::Rectangle2D(Point midpoint, double wx, double wy) :
	p0(midpoint.x - wx/2, midpoint.y - wy/2),
	p1(midpoint.x + wx/2, midpoint.y + wy/2) {}

Point Rectangle2D::Mid() const noexcept {
	return { (p0.x + p1.x)/2, (p0.y + p1.y)/2 };
}

bool Rectangle2D::IsInside(const Point& test) const noexcept {
	return 
		test.x > p0.x && test.x < p1.x && 
		test.y > p0.y && test.y < p1.y; 
}
bool Rectangle2D::IsInside(double x, double y) const noexcept {
	return 
		x > p0.x && x < p1.x && 
		y > p0.y && y < p1.y; 
}

Line2D GetLine(const mnd::geom::Point& p1, const mnd::geom::Point& p2) noexcept {
	double slope = (p2.y - p1.y) / (p2.x - p1.x);
	double offset = -slope * p1.x + p1.y;
	return { offset, slope };
}
Line2D GetLine(const std::array<std::array<double, 2>, 2> p12) noexcept {
	return GetLine( {p12[0][0],p12[0][1]}, {p12[1][0], p12[1][1]} );
}


