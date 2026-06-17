#include "Geometry.h"
#include "../Eigen/Dense"
#include <cmath>

namespace mnd { namespace geom { namespace detail {
	static_assert(std::is_standard_layout_v<Vector3D>);
	static_assert(sizeof(Vector3D) == 3*sizeof(double));
	
	static Eigen::Map<const Eigen::Vector3d> mapv(const std::array<double, 3>& a) noexcept {
		return Eigen::Map<const Eigen::Vector3d>{a.data()};
	}
	static Eigen::Map<Eigen::Vector3d> mapv(std::array<double, 3>& a) noexcept {
		return Eigen::Map<Eigen::Vector3d>{a.data()};
	}
	/* Not sure if this is legal, but works so far. */
	static Eigen::Map<const Eigen::Vector3d> mapv(const Vector3D& a) noexcept {
		return Eigen::Map<const Eigen::Vector3d>{&a.x};
	}
	[[ maybe_unused ]] static Eigen::Map<Eigen::Vector3d> mapv(Vector3D& a) noexcept {
		return Eigen::Map<Eigen::Vector3d>{&a.x};
	}
	static Eigen::Matrix3d rotation_matrix(const Vector3D& axis, double angle) noexcept {
		Eigen::Vector3d u{axis.x, axis.y, axis.z};
		return Eigen::AngleAxisd{angle, u.normalized()}.toRotationMatrix();
	}
	static void store_vec(std::array<double, 3>& dst, const Eigen::Vector3d& v) noexcept {
		dst[0] = v.x();
		dst[1] = v.y();
		dst[2] = v.z();
	}
	static Line2D from_point_direction(const Eigen::Vector2d& n, const Eigen::Vector2d& d) noexcept {
		/* `d` is asserted to be unit vector. 
		 * `n` is orthogonal to `d`. */
		constexpr double eps = 1e-20;

		if(std::abs(d.x() ) < eps) {
			return Line2D{ NAN , NAN }; // vertical +line
		}

		const double a1 = d.y() / d.x();
		const double a0 = n.y() - a1 * n.x();

		return Line2D{a0, a1};
	}
}}}

using namespace mnd::geom;

const Point2D Point2D::null { NAN, NAN };
const Point3D Point3D::null { NAN, NAN };
const Line2D Line2D::null { NAN, NAN };
const Line3D Line3D::null { NAN, NAN, NAN, NAN };

/* Translate the line by an offset vector (val  0). */
Line2D& Line2D::operator+=(double val) noexcept { 
	value[0] -= value[1] * val; 
	return *this;
}
/* Translate the line by an offset vector (rhs.x  rhs.y). */
Line2D& Line2D::operator+=(const Vector2D& rhs) noexcept {
	value[0] += rhs.y - value[1] * rhs.x;
	return *this;
}
/* Represent the line in a coordinate system made by an offset vector (val  0). */ 
Line2D& Line2D::operator%=(double val) noexcept { 
	value[0] += value[1] * val; 
	return *this;
}
/* Represent the line in a coordinate system made by an offset vector (rhs.x  rhs.y). */ 
Line2D& Line2D::operator%=(const Vector2D& rhs) noexcept {
	value[0] += -rhs.y + value[1] * rhs.x;
	return *this;
}
Line2D& Line2D::Rotate(double theta) noexcept {
	const double a0 = value[0];
    const double a1 = value[1];

	// Direction vector
    Eigen::Vector2d d{1.0, value[1]};
    const double norm2 = d.squaredNorm();
	const double norm = std::sqrt(norm2);
	d /= norm;

	// Normal vector
	Eigen::Vector2d n;
	n << -a0*a1 / norm2,
	      a0 / norm2;

	const double ct = std::cos(theta);
	const double st = std::sin(theta);
	Eigen::Matrix2d R;
	R << ct, -st,
	     st,  ct;

    n = R * n;
    d = R * d;

	return *this = detail::from_point_direction(n, d);
}
Line2D& Line2D::RepresentInRotated(double theta) noexcept {
	return Rotate(-theta);
}
Line2D& Line2D::ShiftAndRotate (
	const Vector2D& off,
	double theta
) noexcept {
    return (*this += off).Rotate(theta);
}
Line2D& Line2D::RepresentInShifted (
	const Vector2D& off,
	double theta
) noexcept {
    return (*this %= off).RepresentInRotated(theta);
}

/* ====================== 3D ====================== */

double Line3D::DistanceTo(const Point3D& pt) const noexcept {
	const auto p = detail::mapv(this->p);
	const auto v = detail::mapv(this->v);
	const auto k = detail::mapv( pt );

	const double n = v.norm();

	if(n < 1e-24) { 
		return NAN; // highly unlikely
	}
	return (k-p).cross(v).norm() / n;
}

double Line3D::DistanceTo(const Line3D& rhs) const noexcept {
	auto p1 = detail::mapv(this->p);
	auto v1 = detail::mapv(this->v);
	auto p2 = detail::mapv( rhs.p );
	auto v2 = detail::mapv( rhs.v );

	const Eigen::Vector3d v1_cross_v2 = v1.cross(v2);

	const double n = v1_cross_v2.norm();

	if(n < 1e-20) {
		return (p1-p2).cross(v1).norm() / v1.norm();
	}

	return std::abs( (p1-p2).dot(v1_cross_v2)) / n;
}

double Line3D::AngleRelativeTo(const Line3D& rhs) const noexcept {
	auto v1 = detail::mapv(this->v);
	auto v2 = detail::mapv( rhs.v );
	double c = v1.dot(v2) / (v1.norm() * v2.norm());
	c = std::clamp(c, -1.0, 1.0);
	return std::acos(c);
}

/* Translate the line by an offset vector (0  0  val). */
Line3D& Line3D::operator+=(double val) noexcept {
	p[2] += val;
	return *this;
}
/* Translate the line by an offset vector (rhs.x  rhs.y  rhs.z). */
Line3D& Line3D::operator+=(const Vector3D& rhs) noexcept {
	p[0] += rhs.x;
	p[1] += rhs.y;
	p[2] += rhs.z;
	return *this;
}
/* Represent the line in a coordinate system made by an offset vector (0  0  val). */ 
Line3D& Line3D::operator%=(double val) noexcept {
	p[2] -= val;
	return *this;
}
/* Represent the line in a coordinate system made by an offset vector (rhs.x  rhs.y  rhs.z). */ 
Line3D& Line3D::operator%=(const Point3D& rhs) noexcept {
	p[0] -= rhs.x;
	p[1] -= rhs.y;
	p[2] -= rhs.z;
	return *this;
}

/* Rotate the line by an unary rotation, represented by a unit vector and an angle (amount) */ 
Line3D& Line3D::Rotate(const Vector3D& axis, double angle) noexcept {
	const Eigen::Matrix3d R = detail::rotation_matrix(axis, angle);
	const Eigen::Vector3d p_new = R * detail::mapv(p);
	const Eigen::Vector3d v_new = R * detail::mapv(v);

	detail::store_vec(p, p_new);
	detail::store_vec(v, v_new);
	return *this;
}

/* Represent the line in a coordinate system made by a unary rotation,
 * which is represented by a unit vector and an angle (amount). */ 
Line3D& Line3D::RepresentInRotated(const Vector3D& axis, double angle) noexcept {
	const Eigen::Matrix3d R = detail::rotation_matrix(axis, angle).transpose();
	
	const Eigen::Vector3d p_new = R * detail::mapv(p);
	const Eigen::Vector3d v_new = R * detail::mapv(v);

	detail::store_vec(p, p_new);
	detail::store_vec(v, v_new);
	return *this;

}
/* Offset the line in a coordinate system and rotate by a unary rotation. */ 
Line3D& Line3D::ShiftAndRotate (
	const Vector3D& offset, 
	const Vector3D& axis,
	double angle
) noexcept {
	return (*this += offset).Rotate(axis, angle);
}

/* Represent the line in a coordinate system by an offset and unary rotation. */ 
Line3D& Line3D::RepresentInShifted (
	const Vector3D& offset, 
	const Vector3D& axis,
	double angle
) noexcept {
	return (*this %= offset).RepresentInRotated(axis, angle);
}

/* ============================================================ */

Rectangle2D::Rectangle2D(Point2D midpoint, double wx, double wy) :
	p0{midpoint.x - wx/2, midpoint.y - wy/2},
	p1{midpoint.x + wx/2, midpoint.y + wy/2} {}

Point2D Rectangle2D::Mid() const noexcept {
	return { (p0.x + p1.x)/2, (p0.y + p1.y)/2 };
}

bool Rectangle2D::IsInside(const Point2D& test) const noexcept {
	return 
		test.x > p0.x && test.x < p1.x && 
		test.y > p0.y && test.y < p1.y; 
}
bool Rectangle2D::IsInside(double x, double y) const noexcept {
	return 
		x > p0.x && x < p1.x && 
		y > p0.y && y < p1.y; 
}

Line2D GetLine2D(const mnd::geom::Point2D& p1, const mnd::geom::Point2D& p2) noexcept {
	double slope = (p2.y - p1.y) / (p2.x - p1.x);
	double offset = -slope * p1.x + p1.y;
	return { offset, slope };
}

Line3D GetLine3D(const mnd::geom::Point3D& p1, const mnd::geom::Point3D& p2) noexcept {
	const double dz = p2.z - p1.z;

	double a1 = (p2.x - p1.x) / dz;
	double a0 = p1.x - a1*p1.z;

	double b1 = (p2.y - p1.y) / dz;
	double b0 = p1.y - b1*p1.z;

	return {a0, a1, b0, b1};
}

/* Find a point with a minimal distance to a sequence of lines (aka: a vertex). 
 * If a line is described as: 
 *   vec{r} = vec{p} + l*vec{v} 
 * then the distance of a point `k` to this line is:
 * d^2 = ‖ (k-p) ⨯ v / |v| ‖^2 , or simply projection:
 * d^2 = ‖ hat{P} * (k-p ‖^2 , where `P` is a projector.
 * P = 1 - v vᵀ
 * For N lines, sum them up, take gradient relative to `k`, and equation is:
 * ( sum_i hat{P_i} ) k = sum_i ( hat{P_i} * p_i )
 * */
mnd::geom::Point3D mnd::geom::FindVertex( ::mnd::span<const Line3D> lines) noexcept {
	constexpr double eps = 1e-20;

	Eigen::Matrix3d A = Eigen::Matrix3d::Zero();
	Eigen::Vector3d b = Eigen::Vector3d::Zero();

	const Eigen::Matrix3d I = Eigen::Matrix3d::Identity();
	for(const auto& line : lines) {
		if(!line.HasValue()) {
			continue; // no trolling allowed
		}

		Eigen::Map<const Eigen::Vector3d> p( line.p.data() );
		Eigen::Map<const Eigen::Vector3d> v( line.v.data() );

		const double norm2 = v.squaredNorm();

		if(norm2 < eps) {
			continue; // super degenerate case
		}

		const Eigen::Matrix3d P = I - v * v.transpose() / norm2;

		A += P;
		b += P * p;
	}

	// https://libeigen.gitlab.io/eigen/docs-nightly/classEigen_1_1LDLT.html
	// Recommendation for symmetric (semidefinite) 3x3 linear problem: Ax = b
	// There's also LLT ... ?
	Eigen::LDLT<Eigen::Matrix3d> solver(A);

	if(solver.info() != Eigen::Success)
		return Point3D::null;

	const Eigen::Vector3d k = solver.solve(b);

	return {
		k.x(),
		k.y(),
		k.z()
	};
}

mnd::geom::Point3D mnd::geom::FindVertex(const mnd::geom::Line3D& l1, const mnd::geom::Line3D& l2) noexcept {
	const std::array<Line3D, 2> lines{l1, l2};
	return FindVertex(mnd::as_span(lines));
}
