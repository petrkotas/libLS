/*
 * kdtree.hpp
 *
 *  Created on: Jan 21, 2014
 *      Author: petr
 */
#pragma once
#ifndef KDTREE_HPP_
#define KDTREE_HPP_

#include "TriangleElement.hpp" // this should be abstract element, but for now, fuck it

class KDTreeNode {

	// yes I am going to store it here, so what
	// the triangle data is going to be stored in the node itself
	std::vector<TriangleElement<double>> element;

	// not needed, but just to store the dividing axis
	char axis;

	// childs
	KDTreeNode *left, *right;

	KDTreeNode(const TriangleElement<double>& el) {
		element = el; // should invoke copy
	}

};

class KDTree {
	int maxDepth;

};

#endif /* KDTREE_HPP_ */
