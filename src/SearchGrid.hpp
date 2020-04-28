/*
 * SearchGrid.hpp
 *
 *  Created on: 13. 9. 2013
 *      Author: Petr
 */
#pragma once
#ifndef SEARCHGRID_HPP_
#define SEARCHGRID_HPP_

#include <vector>
#include <list>
#include <functional>
#include <tuple>
#include <memory>
#include <cmath>
#include <Eigen/Dense>
#include "TriangleElement.hpp"
#include "BoundingBox.hpp"




template<int dim>
class SearchGrid {

public:

	std::list<TriangleElement<double> > dataFromPosition(const Eigen::Array3d& coord);


	void insertElement(const TriangleElement<double>& element);


	SearchGrid(const Eigen::Array3d& x_lo, const Eigen::Array3d& x_hi, const Eigen::Array3i& m);


	std::list<TriangleElement<double> > operator()(int i, int j, int k) const;

private:


	std::vector< std::list< TriangleElement<double> > > data; // stores data inside the grid

	Eigen::Array3d x_lo;
	Eigen::Array3d x_hi;
	Eigen::Array3d dx;
	Eigen::Array3i m;

	long int toLinear(const Eigen::Array3i& inds) const;

	Eigen::Array3i toIndices(const Eigen::Array3d& coord) const;

	Box<double, dim> cellBox(const Eigen::Array3i& inds);

	std::tuple<Eigen::Array3i, Eigen::Array3i> listOfSpannedCells(const Box<double, dim>& boundingBox);

};




template<int dim>
SearchGrid<dim>::SearchGrid(const Eigen::Array3d& x_lo, const Eigen::Array3d& x_hi, const Eigen::Array3i& m) {

	this->x_lo = x_lo;
	this->x_hi = x_hi;
	this->m    = m;

	this->dx   = (x_hi - x_lo) / (m - 1).cast<double>();

	this->data.resize( m.prod() );

}





template<int dim>
inline std::list<TriangleElement<double> > SearchGrid<dim>::operator()(int i, int j, int k) const {

	int ii, jj, kk;

	ii = (i >= m[0]) ? m[0]-1 : (i < 0) ? 0 : i;
	jj = (j >= m[1]) ? m[1]-1 : (j < 0) ? 0 : j;
	kk = (k >= m[2]) ? m[2]-1 : (k < 0) ? 0 : k;

	return data[toLinear(Eigen::Array3i(ii,jj,kk))];
}


template <int dim>
long int SearchGrid<dim>::toLinear(const Eigen::Array3i& inds) const {
	return inds[0] + inds[1]*m[0] + inds[2]*m[0]*m[1];
}


template <int dim>
Eigen::Array3i SearchGrid<dim>::toIndices(const Eigen::Array3d& coord) const{
	return ((coord - x_lo) / dx).cast<int>();
}


template <int dim>
std::list<TriangleElement<double> > SearchGrid<dim>::dataFromPosition(const Eigen::Array3d& coord) {
	using namespace Eigen;

	Eigen::Array3i inds = toIndices(coord);
	std::list<TriangleElement<double> > ret;
	int pad = 1;
	while (ret.size() == 0) {
		for (int i = inds[0]-pad; i <= inds[0]+pad; ++i) {
			for (int j = inds[1]-pad; j <= inds[1]+pad; ++j) {
				for (int k = inds[2]-pad; k <= inds[2]+pad; ++k) {
					std::list<TriangleElement<double> > tmp = this->operator ()(i,j,k);
					ret.insert(ret.begin(), tmp.begin(), tmp.end());
				}
			}
		}
		pad++;
	}

	return ret;
}


template <int dim>
Box<double, dim> SearchGrid<dim>::cellBox(const Eigen::Array3i& inds) {

	double cx[dim] = { (x_lo[0]+ 0.5*dx[0]) + inds[0]*dx[1],
					   (x_lo[1]+ 0.5*dx[1]) + inds[1]*dx[1],
					   (x_lo[2]+ 0.5*dx[2]) + inds[2]*dx[2]};
	double ex[dim] = { 0.5*dx[0], 0.5*dx[1], 0.5*dx[2] };

	return Box<double, dim>(cx, ex);

}


template <int dim>
std::tuple<Eigen::Array3i, Eigen::Array3i> SearchGrid<dim>::listOfSpannedCells(const Box<double, dim>& boundingBox) {

	using namespace Eigen;

	Eigen::Array3i bli = toIndices(boundingBox.bl());
	Eigen::Array3i tri = toIndices(boundingBox.tr());

	return std::forward_as_tuple(bli, tri);

}


template <int dim>
void SearchGrid<dim>::insertElement(const TriangleElement<double>& element) {

	auto cellSpan = listOfSpannedCells( element.getBoundingBox() );

	auto bli = std::get<0>(cellSpan);
	auto uri = std::get<1>(cellSpan);

	for (int i = bli[0]; i <= uri[0]; ++i) {
		for (int j = bli[1]; j <= uri[1]; ++j) {
			for (int k = bli[2]; k <= uri[2]; ++k) {
				if ( element.isPartlyInRegion(cellBox(Eigen::Array3i(i,j,k))) ) {
					data[toLinear(Eigen::Array3i(i,j,k))].push_back(element);
				}
			}
		}
	}

}

//==========================================
// ostream operators used for debug printing
//==========================================

template <int d>
std::ostream& operator<<(std::ostream& os, const SearchGrid<d>& sg) {

	std::cout << sg.m << std::endl;

    return os;
}



#endif















