/*
 * TriangleElement.hpp
 *
 *  Created on: Nov 18, 2012
 *      Author: petr
 */
#pragma once
#ifndef TRIANGLEELEMENT_HPP_
#define TRIANGLEELEMENT_HPP_

#undef max
#undef min

#include <cmath>
#include <limits>
#include <Eigen/Dense>
#include <Eigen/Geometry>
#include <Eigen/StdVector>
#include <vector>
#include <iostream>
#include "BoundingBox.hpp"
#include "utility.h"
#include "tribox3.hpp"

// #include <boost/archive/basic_binary_oarchive.hpp>
// #include <boost/archive/basic_binary_iarchive.hpp>


template<typename type>
struct SignedDistance {
	type            dist;
	type            angle;
	int             sign;
	Eigen::Vector3d minPoint;
};


template <typename type>
class TriangleElement {


public:

	/// Stores vertices of an element. 
	/// This implementation uses triangle element, but in future, there will be support for multiple element types.
	std::vector<Eigen::Vector3d, Eigen::aligned_allocator<Eigen::Vector3d> > vertices;
	
	/// Face normal
	Eigen::Vector3d normal;

	

	/// element ID used to identify element
	long int ID;

	TriangleElement(): ID(-1) {};

	/// Triangle element contructor, triangle normal is computed automaticaly
	/// Normal is computed standard way : n = (v1 - v0) x (v2 - v0), where x is cross product
	TriangleElement(Eigen::Vector3d& v0, Eigen::Vector3d& v1, Eigen::Vector3d& v2, long int id);

	TriangleElement(const Eigen::Vector3d& v) {
		this->vertices.reserve(3);

		this->vertices.push_back(v);
		this->vertices.push_back(v);
		this->vertices.push_back(v);

		this->ID = -1;
	}

	/// Constructor to create the triangle element from complete list of vertices with an ID
	TriangleElement(std::vector<Eigen::Vector3d, Eigen::aligned_allocator<Eigen::Vector3d> >& vertices, long int id);

	TriangleElement(const TriangleElement<type> & element) : vertices(element.vertices), normal(element.normal), ID(element.ID) {};

	/// Test whether element is inside the rectangular! region
	bool isInsideRegion(const Box<type, 3>& bb) const;

	/// Test whether element is partly inside the rectangular! region
	bool isPartlyInRegion(const Box<type, 3>& bb) const;

	/// Return centroid of element
	Eigen::Vector3d centroid();

	static SignedDistance<type> computeDistance(const TriangleElement<type>& el, const Eigen::Vector3d& point);

	static type computePerpendicularDistance(const TriangleElement<type>& el, const Eigen::Vector3d& point);

	/// Return tight Axis Aligned Bounding Box around element
	BoundingBox<type, 3, TriangleElement<double> > getBoundingBox() const;

};

template<typename type> 
TriangleElement<type>::TriangleElement(Eigen::Vector3d& v0, Eigen::Vector3d& v1, Eigen::Vector3d& v2, long int id) {

	this->vertices.reserve(3);

	this->vertices.push_back(v0);
	this->vertices.push_back(v1);
	this->vertices.push_back(v2);

	Eigen::Vector3d e0 = v1 - v0;
	Eigen::Vector3d e1 = v2 - v0;

	this->normal = e0.cross(e1);
	this->normal.normalize();

	this->ID = id;

}

template<typename type> 
TriangleElement<type>::TriangleElement(std::vector<Eigen::Vector3d, Eigen::aligned_allocator<Eigen::Vector3d> >& vertices, long int id) {

	this->vertices = vertices;

	Eigen::Vector3d e0 = vertices[1] - vertices[0];
	Eigen::Vector3d e1 = vertices[2] - vertices[0];

	this->normal = e0.cross(e1);
	this->normal.normalize();

	this->ID = id;

}

template<typename type> 
inline BoundingBox<type, 3, TriangleElement<double> > TriangleElement<type>::getBoundingBox() const {

	double minX[3], maxX[3];
	minX[0] =  std::numeric_limits<type>::max();
	minX[1] =  std::numeric_limits<type>::max();
	minX[2] =  std::numeric_limits<type>::max();

	maxX[0] = -std::numeric_limits<type>::max();
	maxX[1] = -std::numeric_limits<type>::max();
	maxX[2] = -std::numeric_limits<type>::max();


	for (auto vertex: vertices) {

		// find lower bound
		if ( vertex[0] < minX[0] ) // for x
			minX[0] = vertex[0];
		if ( vertex[1] < minX[1] ) // for y
			minX[1] = vertex[1];
		if ( vertex[2] < minX[2] ) // for z
			minX[2] = vertex[2];

		// find upper bound
		if ( vertex[0] > maxX[0] ) // for x
			maxX[0] = vertex[0];
		if ( vertex[1] > maxX[1] ) // for y
			maxX[1] = vertex[1];
		if ( vertex[2] > maxX[2] ) // for z
			maxX[2] = vertex[2];

	}

	return BoundingBox<type, 3, TriangleElement<double> >(minX, maxX, this);

}

template<typename type> 
inline bool TriangleElement<type>::isInsideRegion(const Box<type, 3>& bb) const {

	bool vertex_0 = bb.pointInside(vertices[0].array());
	bool vertex_1 = bb.pointInside(vertices[1].array());
	bool vertex_2 = bb.pointInside(vertices[2].array());
	
	return vertex_0 && vertex_1 && vertex_2;

}



/// \brief Detects if triangle element collide with cuboid region. Tests for partial collisions
///        which could pose a tricky problem: long thin triangles, one vertex collision,
///        one edge collision, edge to edge collision.
///
/// The algorithm is taken from book: (Morgan Kaufmann series in  Interactive 3d technology)
///                                   Realtime collision detection - Christer Ericson
///
template<typename type>
inline bool TriangleElement<type>::isPartlyInRegion(const Box<type, 3>& bb) const {

	float trivert[3][3] = { {(float)vertices[0][0], (float)vertices[0][1], (float)vertices[0][2]},
							{(float)vertices[1][0], (float)vertices[1][1], (float)vertices[1][2]},
							{(float)vertices[2][0], (float)vertices[2][1], (float)vertices[2][2]}};

	float center[3] = {(float)bb.center(0), (float)bb.center(1), (float)bb.center(2)};

	float extent[3] = {(float)bb.extent(0), (float)bb.extent(1), (float)bb.extent(2)};

	// utilizes code I found on the internet, check the functionality
	return triBoxOverlap(center, extent, trivert);

}

template<typename type> 
inline SignedDistance<type> TriangleElement<type>::computeDistance(const TriangleElement<type>& el, const Eigen::Vector3d &point) {
	// point is supposed to have right dimension, so do it nice and clear :)
	// this function computed the distance between triangle element and point in 3D
	// complete description of a method could be found here :
	//  	http://www.geometrictools.com/Documentation/DistancePoint3Triangle3.pdf
	// sample implementation could be found here :
	//		http://www.geometrictools.com/LibMathematics/Distance/Wm5DistPoint3Triangle3.cpp
	// and here :
	// 		http://www.mathworks.com/matlabcentral/fileexchange/22857-distance-between-a-point-and-a-triangle-in-3d/

	//double E0[3] = {el.vertex_x[1] - el.vertex_x[0], el.vertex_y[1] - el.vertex_y[0], el.vertex_z[1] - el.vertex_z[0]};
	//double E1[3] = {el.vertex_x[2] - el.vertex_x[0], el.vertex_y[2] - el.vertex_y[0], el.vertex_z[2] - el.vertex_z[0]};

	Eigen::Vector3d E0 = el.vertices[1] - el.vertices[0];
	Eigen::Vector3d E1 = el.vertices[2] - el.vertices[0];

	//double D[3]  = {el.vertex_x[0] - point[0], el.vertex_y[0] - point[1], el.vertex_z[0] - point[2]};

	Eigen::Vector3d D = el.vertices[0] - point;

//	double a00_ = E0[0]*E0[0] + E0[1]*E0[1] + E0[2]*E0[2];
	type a00 = E0.dot(E0);
	
//	std::cout << (a00_ - a00) << std::endl;

//	double a01_ = E0[0]*E1[0] + E0[1]*E1[1] + E0[2]*E1[2];
	type a01 = E0.dot(E1);

//	std::cout << (a01_ - a01) << std::endl;

//	double a11_ = E1[0]*E1[0] + E1[1]*E1[1] + E1[2]*E1[2];
	type a11 = E1.dot(E1);
	
//	std::cout << (a11_ - a11) << std::endl;

//	double b0_  = E0[0]*D[0]  + E0[1]*D[1]  + E0[2]*D[2];
	type b0 = E0.dot(D);

//	std::cout << (b0_ - b0) << std::endl;

//	double b1_  = E1[0]*D[0]  + E1[1]*D[1]  + E1[2]*D[2];
	type b1 = E1.dot(D);

//	std::cout << (b0_ - b0) << std::endl;

//	double c_   = D[0]*D[0]   + D[1]*D[1]   + D[2]*D[2];
	type c = D.dot(D);

//	std::cout << (c_ - c) << std::endl;

	type det = std::abs(a00*a11 - a01*a01);

	type s   = a01*b1 - a11*b0;
	type t   = a01*b0 - a00*b1;

	type sqrDistance = 0;


	if ( (s+t) <= det ) {
		if ( s < 0 ) {
			if ( t < 0 ) {
				// region 4
				if ( b0 < 0 ) {
					t = 0;
					if ( -b0 >= a00 ) {
						s = 1;
						sqrDistance = a00 + 2*b0 + c;
					} else {
						s = -b0/a00;
						sqrDistance = b0*s + c;
					}
				} else {
					s = 0;
					if ( b1 >= 0 ) {
						t = 0;
						sqrDistance = c;
					} else if ( -b1 >= a11 ) {
						t = 1;
						sqrDistance = a11 + 2*b1 + c;
					} else {
						t = -b1/a11;
						sqrDistance = b1*t + c;
					}
				}
				// end of region 4
			} else {
				// region 3
				s = 0;
				if ( b1 >= 0 ) {
					t = 0;
					sqrDistance = c;
				} else if ( -b1 >= a11 ) {
					t = 1;
					sqrDistance = a11 + 2*b1 + c;
				} else {
					t = -b1/a11;
					sqrDistance = b1*t + c;
				}
				// end of region 3
			}
		} else if ( t < 0) {
			// region 5
			t = 0;
			if ( b0 >= 0 ) {
				s = 0;
				sqrDistance = c;
			} else if ( -b0 >= a00 ) {
				s = 1;
				sqrDistance = a00 + 2*b0 + c;
			} else {
				s = -b0/a00;
				sqrDistance = b0*s + c;
			}
			// end of region 5
		} else {
			// region 0
			double invDet = 1 / det;
			s *= invDet;
			t *= invDet;
			sqrDistance = s*(a00*s + a01*t + 2*b0) + t*(a01*s + a11*t + 2*b1) + c;
			// end of region 0
		}
	} else {
		if ( s < 0 ) {
			// region 2
			double tmp0 = a01 + b0;
			double tmp1 = a11 + b1;
			if ( tmp1 > tmp0 ) { // minimum on edge s+t=1
				double numer = tmp1 - tmp0;
				double denom = a00 - 2*a01 + a11;
				if ( numer >= denom ) {
					s = 1;
					t = 0;
					sqrDistance = a00 + 2*b0 + c;
				} else {
					s = numer/denom;
					t = 1 - s;
					sqrDistance = s*(a00*s + a01*t + 2*b0) + t*(a01*s + a11*t + 2*b1) + c;
				}
			} else { // minimum on edge s=0
				s = 0;
				if ( tmp1 <= 0 ) {
					t = 1;
					sqrDistance = a11 + 2*b1 + c;
				} else if ( b1 >= 0 ) {
					t = 0;
					sqrDistance = c;
				} else {
					t = -b1/a11;
					sqrDistance = b1*t + c;
				}
			}
			// end of region 2
		} else if ( t < 0 ) {
			// region 6
			double tmp0 = a01 + b1;
			double tmp1 = a00 + b0;
			if ( tmp1 > tmp0 ) {
				double numer = tmp1 - tmp0;
				double denom = a00 - 2*a01 + a11;
				if ( numer >= denom ) {
					t = 1;
					s = 0;
					sqrDistance = a11 + 2*b1 + c;
				} else {
					t = numer/denom;
					s = 1 - t;
					sqrDistance = s*(a00*s + a01*t + 2*b0) + t*(a01*s + a11*t + 2*b1) + c;
				}
			} else {
				t = 0;
				if ( tmp1 <= 0 ) {
					s = 1;
					sqrDistance = a00 + 2*b0 + c;
				} else if ( b0 >= 0 ) {
					s = 0;
					sqrDistance = c;
				} else {
					s = -b0/a00;
					sqrDistance = b0*s + c;
				}
			}
			// end of region 6
		} else {
			// region 1
			double numer = a11 + b1 - a01 - b0;
			if ( numer <= 0 ) {
				s = 0;
				t = 1;
				sqrDistance = a11 + 2*b1 + c;
			} else {
				double denom = a00 - 2*a01 + a11;
				if ( numer >= denom ) {
					s = 1;
					t = 0;
					sqrDistance = a00 + 2*b0 + c;
				} else {
					s = numer / denom;
					t = 1 - s;
					sqrDistance = s*(a00*s + a01*t + 2*b0) + t*(a01*s + a11*t + 2*b1) + c;
				}
			}
			// end of region 1
		}
	}

	double eeps = my_eps;

	if (sqrDistance <= eeps)
		sqrDistance = 0.0;

	SignedDistance<type> dist;

	dist.dist = std::sqrt(sqrDistance);

	dist.minPoint = el.vertices[0] + s*E0 + t*E1;

//	dist.minPoint[0] = (el.vertex_x[0] + s*E0[0] + t*E1[0]);
//	dist.minPoint[1] = (el.vertex_y[0] + s*E0[1] + t*E1[1]);
//	dist.minPoint[2] = (el.vertex_z[0] + s*E0[2] + t*E1[2]);


//	double PP0[3] = { point[0] - dist.minPoint[0],
//					  point[1] - dist.minPoint[1],
//					  point[2] - dist.minPoint[2] };

//	double nrm = sqrt(PP0[0]*PP0[0] + PP0[1]*PP0[1] + PP0[2]*PP0[2]);

	// compute on which side of the triangle we are
	Eigen::Vector3d PP0 = point - dist.minPoint;
	
	if ( PP0.isZero(eeps) ) {
		// PP0 = Eigen::Vector3d::Zero();
	} else {
		PP0.normalize();
	}


	//double nrm_normal = sqrt(el.normal_x*el.normal_x + el.normal_y*el.normal_y + el.normal_z*el.normal_z);

//	double dprod = (el.normal_x/nrm_normal) * (PP0[0]/nrm) +
//				   (el.normal_y/nrm_normal) * (PP0[1]/nrm) +
//				   (el.normal_z/nrm_normal) * (PP0[2]/nrm);



	type dprod = el.normal.dot(PP0);

//	std::cout << el.normal.transpose() << " DOT " << PP0.transpose() << " = " << dprod << std::endl;

	int sgn = 0;
	if (dprod > eeps)
		sgn = 1;
	else if (dprod < -eeps)
		sgn = -1;

	dist.angle = dprod;
	dist.sign  = sgn;

	return dist;

}

template<typename type>
inline Eigen::Vector3d TriangleElement<type>::centroid() {

	Eigen::Vector3d c;

	for (auto vertex: vertices) {
		c += vertex;
	}

	return c*(1.0/3.0);

}

template<typename type> 
inline type TriangleElement<type>::computePerpendicularDistance(const TriangleElement& el, const Eigen::Vector3d& point) {
	// compute projection onto surface defined by the triangle element
	Eigen::Vector3d vec;
	vec << point - el.vertices[0];
	
	type ortDist = std::abs( vec.dot(el.normal) );
	
	return ortDist;
}


// ostream operators

template <typename type>
std::ostream& operator<<(std::ostream& os, const TriangleElement<type>& el) {
	os << "vert_0: " << el.vertices[0].transpose() << std::endl;
	os << "vert_1: " << el.vertices[1].transpose() << std::endl;
	os << "vert_2: " << el.vertices[2].transpose() << std::endl;
	os << "normal: " << el.normal.transpose();

    return os;
}





struct TriangleCompare {

	bool operator() (spatial::dimension_type dim, const TriangleElement<double>& a, const TriangleElement<double>& b) const {

		switch (dim) {
			case 0:
//				a.vertices[0]
				return (a.vertices[0][0] < b.vertices[0][0]) &&
					   (a.vertices[1][0] < b.vertices[1][0]) &&
					   (a.vertices[2][0] < b.vertices[2][0]);
				break;
			case 1:
				return (a.vertices[0][1] < b.vertices[0][1]) &&
					   (a.vertices[1][1] < b.vertices[1][1]) &&
					   (a.vertices[2][1] < b.vertices[2][1]);
//				return (a._tr[0] < b._tr[0]);
				break;
			case 2:
				return (a.vertices[0][2] < b.vertices[0][2]) &&
					   (a.vertices[1][2] < b.vertices[1][2]) &&
					   (a.vertices[2][2] < b.vertices[2][2]);
//				return (a._bl[1] < b._bl[1]);
				break;

		}

		return false;
	}

};

struct TrianglePredicate {
	const Box<double, 3>& box;

	TrianglePredicate(const Box<double, 3>& b) : box(b) {}

	spatial::relative_order
	operator() (spatial::dimension_type dim, spatial::dimension_type, const TriangleElement<double>& t ) const {

		Box<double, 3> b = t.getBoundingBox();

		switch (dim) {
			case 0:
				return (b._bl[0] > box._tr[0] && b._tr[0] > box._tr[0]) ? spatial::above :
					   (b._tr[0] < box._bl[0] && b._bl[0] < box._bl[0]) ? spatial::below :
						spatial::matching;
				break;
			case 1:
				return (b._bl[1] > box._tr[1] && b._tr[1] > box._tr[1]) ? spatial::above :
					   (b._tr[1] < box._bl[1] && b._bl[1] < box._bl[1]) ? spatial::below :
						spatial::matching;
				break;
			case 2:
				return (b._bl[2] > box._tr[2] && b._tr[2] > box._tr[2]) ? spatial::above :
					   (b._tr[2] < box._bl[2] && b._bl[2] < box._bl[2]) ? spatial::below :
						spatial::matching;
				break;
		}

	}
};

struct TriangleMetric {
    // Check that DistanceType is a fundamental floating point type
    typedef double distance_type;


    /**
     *  Compute the distance between the point of \c origin and the \c key.
     *  \return The resulting square distance.
     */
    distance_type
    distance_to_key(spatial::dimension_type rank,
                    const TriangleElement<double>& origin,
                    const TriangleElement<double>& key) const {

    	distance_type result = distance_type(0.0);

//    	std::cout << "===============================================" << std::endl;
//    	std::cout << origin << std::endl;
//    	std::cout << "--------------" << std::endl;
//    	std::cout << key << std::endl;
//    	std::cout << "--------------" << std::endl;

    	Eigen::Array3d co = (origin.vertices[0] + origin.vertices[1] + origin.vertices[2]).array() / 3.0;
    	Eigen::Array3d ck = (key.vertices[0] + key.vertices[1] + key.vertices[2]).array() / 3.0;
    	
    	// distance_type tmp = 0.0
    	
		for (int ii = 0; ii < 3; ++ii) {
			double d = std::numeric_limits<double>::max();
			for (spatial::dimension_type i = 0; i < rank; ++i) {
				double tmp = std::abs( key.vertices[i][ii] - origin.vertices[i][ii] );
				d = (tmp < d) ? tmp : d;
			}

			result += d * d;
		}
    	
//    	std::cout << "dist = " << result << std::endl;
//    	std::cout << "===============================================" << std::endl;
    	return result;

    }

    /**
     *  The distance between the point of \c origin and the closest point to
     *  the plane orthogonal to the axis of dimension \c dim and crossing \c
     *  key.
     *  \return The resulting square distance.
     */
    distance_type
    distance_to_plane(spatial::dimension_type,
    			      spatial::dimension_type dim,
                      const TriangleElement<double>& origin,
                      const TriangleElement<double>& key) const {

    	distance_type d = std::numeric_limits<double>::max();
		for (spatial::dimension_type i = 0; i < 3; ++i) {
			double tmp = std::abs( key.vertices[i][dim] - origin.vertices[i][dim] );
			d = (tmp < d) ? tmp : d;
		}
		    	
		return d * d;

    }


};











#endif /* TRIANGLEELEMENT_HPP_ */
