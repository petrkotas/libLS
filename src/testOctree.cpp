/*
 * testOctree.cpp
 *
 *  Created on: Aug 15, 2013
 *      Author: petr
 */

#include <iostream>
#include <Eigen/Dense>
#include <list>
#include "BoundingBox.hpp"
#include "Octree.hpp"




int main(int argc,char **argv) {

	std::cout << "Octree TEST" << std::endl;

	double center[3] = {0.5,0.5,0.5};
	double extent[3] = {0.5,0.5,0.5};

	BoundingBox<double,3> aabb(center, extent);

	TreeNode<int> tree(aabb, 4);

	// octree created now add some items




	double cx[3] = {0.5, 0.5, 0.5};
	double ex[3] = {0.2, 0.2, 0.2};

	BoundingBox<double, 3> el(cx, ex);
	tree.insertObject(el, 1);


	cx[0] = 0.26; cx[1] = 0.26; cx[2] = 0.876;
	ex[0] = 0.11; ex[1] = 0.11; ex[2] = 0.11;

	BoundingBox<double, 3> el2(cx, ex);
	tree.insertObject(el2, 2);


	cx[0] = 0.625; cx[1] = 0.625; cx[2] = 0.625;
	ex[0] = 0.0625; ex[1] = 0.0625; ex[2] = 0.0625;

	BoundingBox<double, 3> el3(cx, ex);
	tree.insertObject(el3, 3);


	cx[0] = 0.375; cx[1] = 0.375; cx[2] = 0.375;
	ex[0] = 0.06; ex[1] = 0.06; ex[2] = 0.06;

	BoundingBox<double, 3> el7(cx, ex);
	tree.insertObject(el7, 4);


	cx[0] = 0.4375; cx[1] = 0.4375; cx[2] = 0.4375;
	ex[0] = 0.05; ex[1] = 0.05; ex[2] = 0.05;

	BoundingBox<double, 3> el4(cx, ex);
	tree.insertObject(el4, 5);


	cx[0] = 0.0625; cx[1] = 0.0625; cx[2] = 0.0625;
	ex[0] = 0.02; ex[1] = 0.02; ex[2] = 0.02;

	BoundingBox<double, 3> el5(cx, ex);
	tree.insertObject(el5, 6);


	cx[0] = 0.03125; cx[1] = 0.03125; cx[2] = 0.03125;
	ex[0] = 0.01; ex[1] = 0.01; ex[2] = 0.01;

	BoundingBox<double, 3> el6(cx, ex);
	tree.insertObject(el6, 7);


	cx[0] = 0.5820; cx[1] = 0.5820; cx[2] = 0.4357;
	ex[0] = 0.02; ex[1] = 0.02; ex[2] = 0.02;

	BoundingBox<double, 3> el8(cx, ex);
	tree.insertObject(el8, 8);


	Eigen::Vector3d p;
	p << 0.5, 0.5, 0.5;

	std::list<int> data = tree.findClosest(p);


	std::cout << tree << std::endl;


//	std::list<int> data = TreeNode<int>::depthFirst(tree);

	for (auto item: data) {
		std::cout << item << ", ";
	}
	std::cout << std::endl;



	return 0;
}



