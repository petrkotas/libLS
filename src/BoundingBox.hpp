/*
 * Box.hpp
 *
 *  Created on: Feb 20, 2013
 *      Author: petr
 */
#pragma once
#ifndef BOUNDINGBOX_HPP_
#define BOUNDINGBOX_HPP_

#undef max
#undef min

#include <iostream>
#include <limits>
#include <Eigen/Dense>
#include "TriangleElement.hpp"



/// \brief Class representing simple bounding box.
///
/// Bounding box has one static method used for detecting intersection between two of them.
///
/// \tparam type Specifies the type of bounding box (double, int, ...)
/// \tparam dim Specifies the dimension of the bounding box (2 or 3)
template <typename type, int dim>
struct Box {

	// bottom left and top right corners
	type _bl[dim], _tr[dim];
//	double center[dim], extent[dim];

	/// \brief Default constructor. Initializes bounding box from [0,0,0,..] to [1,1,1,..]
	Box() {
		for (int i = 0; i < dim; i++) {

//			center[i] = 0.0;
//			extent[i] = -std::numeric_limits<type>::max();
			_bl[i] =  std::numeric_limits<type>::max();
			_tr[i] = -std::numeric_limits<type>::max();

		}

	}


	Box(const Box<type, dim>& b) {

		for (int i = 0; i < dim; ++i) {

//			this->center[i] = b.center[i];
//			this->extent[i] = b.extent[i];
			_bl[i] = b._bl[i];
			_tr[i] = b._tr[i];

		}

	}

	Box(const type minx[dim], const type maxx[dim]) {

		for (int i = 0; i < dim; ++i) {

			_bl[i] = minx[i];
			_tr[i] = maxx[i];

		}

	}

	Box(const type point[dim]) {

		for (int i = 0; i < dim; ++i) {

			_bl[i] = _tr[i] = point[i];

		}

	}

	Box(const Eigen::Vector3d& point) {

		for (int i = 0; i < dim; ++i) {

			_bl[i] = _tr[i] = point[i];

		}

	}


	static Box fromCenterExtent(type center[dim], type extent[dim]) {
		type _bl[dim], _tr[dim];

		for	(int i = 0; i < dim; ++i) {
			_bl[i] = center[i] - extent[i];
			_tr[i] = center[i] + extent[i];
		}

		return Box<type, dim>(_bl, _tr);
	}


	inline const type minX(int d) const {
		return _bl[d];
	}


	inline const type maxX(int d) const {
		return _tr[d];
	}

	inline double center(int d) const {
		return (_bl[d] + _tr[d]) / 2.0;
	}

	inline double extent(int d) const {
		return (_tr[d] - _bl[d]) / 2.0;
	}

	inline Eigen::Array<type, dim, 1> bl() const {
		return Eigen::Array<type, dim, 1>(_bl);
	}


	inline Eigen::Array<type, dim, 1> tr() const {
		return Eigen::Array<type, dim, 1>(_tr);
	}

	inline int dims() const {
		return dim;
	}


	inline void expand(const Box<type, dim>& bb) {
	
		for (int i = 0; i < dim; ++i) {

			if ( bb._bl[i] < _bl[i] )
				_bl[i] = bb._bl[i];

			if ( bb._tr[i] > _tr[i] )
				_tr[i] = bb._tr[i];

		}

	}


	// ratio is some positive double number. region will grow as extent * ratio
	inline void grow(const double ratio) {
		double center[dim], extent[dim];

		for (int i = 0; i < dim; ++i) {

			center[i] = double(_bl[i] + _tr[i]) / 2.0;
			extent[i] = double(_tr[i] - _bl[i]) / 2.0;

			extent[i] *= ratio;

			_bl[i] = type(center[i] - extent[i]);
			_tr[i] = type(center[i] + extent[i]);
		}
	}


	bool intersect(const Box<type, dim>& bb) const {
		double center[dim], extent[dim], bb_center[dim], bb_extent[dim];
		bool collide = true;

		for (int i = 0; i < dim; ++i) {

			center[i] = double(_bl[i] + _tr[i]) / 2.0;
			extent[i] = double(_tr[i] - _bl[i]) / 2.0;
			bb_center[i] = double(bb._bl[i] + bb._tr[i]) / 2.0;
			bb_extent[i] = double(bb._tr[i] - bb._bl[i]) / 2.0;

			collide &= std::abs( center[i] - bb_center[i] ) <= ( extent[i] + bb_extent[i] );

		}

		return collide;

	}
	

	/// \brief Determines if two bounding boxes collide
	/// \param b1 First bounding box
	/// \param b2 Second bounding box
	static bool intersects(const Box<type, dim>& b1, const Box<type, dim>& b2) {

		return b1.intersect(b2);

	}


	/// Simple test to check whether the point is inside the bounding rectangle.
	/// Mainly used as utility function in other tests.
	bool pointInside(const double point[3]) const {
		double center[dim], extent[dim];
		bool collide = true;

		for (int i = 0; i < dim; ++i) {

			center[i] = double(_bl[i] + _tr[i]) / 2.0;
			extent[i] = double(_tr[i] - _bl[i]) / 2.0;

			collide &= std::abs( point[i] - center[i] ) <= extent[i];

		}

		return collide;

	}



	void shiftCenter(const double shift[3]) {

		for (int i = 0; i < 3; ++i) {
			_bl[i] += shift[i];
			_tr[i] += shift[i];
		}

	}

	// return closest point projection onto the AABB
	Eigen::Array<type, dim, 1> closestPoint(const Eigen::Array<type, dim, 1>& pt) {

		double tmp[3];

		for (int i = 0; i < 3; ++i) {
			double v = pt[i];

			if ( v < _bl[i] ) v = _bl[i];
			if ( v > _tr[i] ) v = _tr[i];

			tmp[i] = v;
		}

		return Eigen::Array<type, dim, 1>(tmp);

	}


	double sqDistToPoint(const Eigen::Array<type, dim, 1>& pt) {
		double sqDist = 0.0;

		for (int i = 0; i < 3; ++i) {
			double v = pt[i];

			if ( v < _bl[i] ) sqDist += (_bl[i] - v) * (_bl[i] - v);
			if ( v > _tr[i] ) sqDist += (v - _tr[i]) * (v - _tr[i]);

		}

		return sqDist;
	}

	type operator() (int d) {
		switch (d) {
		case 0:
			return _bl[0];
			break;
		case 1:
			return _tr[0];
			break;
		case 2:
			return _bl[1];
			break;
		case 3:
			return _tr[1];
		}
	}

};

template <typename type, int dim, typename ownerType = void>
struct BoundingBox : public Box<type, dim> {
	const ownerType* _owner;

	BoundingBox(const type minx[dim], const type maxx[dim], const ownerType* owner) :
		Box<type, dim>(minx, maxx), _owner(owner) {}

	BoundingBox(const BoundingBox& bb) : Box<type, dim>(bb._bl, bb._tr), _owner(bb._owner) {}

	BoundingBox(const Eigen::Vector3d& point) : Box<type, dim>(point), _owner(nullptr) {}

	const ownerType* getOwner() const { return _owner; }

};


//==========================================
// ostream operators used for debug printing
//==========================================

template <typename type, int dim>
std::ostream& operator<<(std::ostream& os, const Box<type, dim>& bb) {

	os << "Bounding box : " << std::endl;

	for (int i = 0; i < dim; i++) {

		os << "| " << bb._bl[i] << " - " << bb._tr[i] << " |" << std::endl;

	}

    return os;
}

///////



struct BoundingBoxCompare {

	bool operator() (spatial::dimension_type dim, const Box<double, 3>& a, const Box<double, 3>& b) const {

		switch (dim) {
			case 0:
				return (a._bl[0] < b._bl[0]);
				break;
			case 1:
				return (a._tr[0] < b._tr[0]);
				break;
			case 2:
				return (a._bl[1] < b._bl[1]);
				break;
			case 3:
				return (a._tr[1] < b._tr[1]);
				break;
			case 4:
				return (a._bl[2] < b._bl[2]);
				break;
			case 5:
				return (a._tr[2] < b._tr[2]);
				break;
		}

		return false;
	}

};


struct BoundingBoxPredicate {
	const Box<double, 3>& box;

	BoundingBoxPredicate(const Box<double, 3>& b) : box(b) {}

	spatial::relative_order
	operator() (spatial::dimension_type dim, spatial::dimension_type, const Box<double, 3>& b ) const {

		switch (dim) {
			case 0:
			case 1:
				return (b._bl[0] > box._tr[0] && b._tr[0] > box._tr[0]) ? spatial::above :
					   (b._tr[0] < box._bl[0] && b._bl[0] < box._bl[0]) ? spatial::below :
						spatial::matching;
				break;
			case 2:
			case 3:
				return (b._bl[1] > box._tr[1] && b._tr[1] > box._tr[1]) ? spatial::above :
					   (b._tr[1] < box._bl[1] && b._bl[1] < box._bl[1]) ? spatial::below :
						spatial::matching;
				break;
			case 4:
			case 5:
				return (b._bl[2] > box._tr[2] && b._tr[2] > box._tr[2]) ? spatial::above :
					   (b._tr[2] < box._bl[2] && b._bl[2] < box._bl[2]) ? spatial::below :
						spatial::matching;
				break;
		}

	}
};



  struct BoundingBoxMetric {
    // Check that DistanceType is a fundamental floating point type
    typedef double distance_type;


    /**
     *  Compute the distance between the point of \c origin and the \c key.
     *  \return The resulting square distance.
     */
    distance_type
    distance_to_key(spatial::dimension_type rank,
                    const Box<double, 3>& origin,
                    const Box<double, 3>& key) const {

    	distance_type result = distance_type(0.0);

//    	std::cout << "===============================================" << std::endl;
//    	std::cout << origin << std::endl;
//    	std::cout << "--------------" << std::endl;
//    	std::cout << key << std::endl;
//    	std::cout << "--------------" << std::endl;

    	for (spatial::dimension_type i = 0; i < rank/2; ++i) {

    		distance_type tmp = std::abs( key.center(i) - origin.center(i) ) + key.extent(i);

    		result += tmp * tmp;

//    		if (tmp_bl < tmp_tr) result += tmp_bl*tmp_bl;
//    		else                 result += tmp_tr*tmp_tr;

//    		if ( origin._tr[i] < key._tr[i] ) result += (key._tr[i] - origin._tr[i]) * (key._tr[i] - origin._tr[i]);
//    		if ( key._bl[i] < origin._bl[i] ) result += (origin._bl[i] - key._bl[i]) * (origin._bl[i] - key._bl[i]);

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
                      const Box<double, 3>& origin,
                      const Box<double, 3>& key) const {

     	spatial::dimension_type i = dim % 3; // probably lets try it and see

    	distance_type tmp = std::abs( key.center(i) - origin.center(i) ) + key.extent(i);

		return tmp * tmp;

//    	if ( origin._tr[dim] < key._tr[dim] ) return (key._tr[dim] - origin._tr[dim]) * (key._tr[dim] - origin._tr[dim]);
//    	if ( key._bl[dim] < origin._bl[dim] ) return (origin._bl[dim] - key._bl[dim]) * (origin._bl[dim] - key._bl[dim]);

    }


  };




#endif /* BOUNDINGBOX_HPP_ */
