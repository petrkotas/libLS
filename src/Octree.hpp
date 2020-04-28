#pragma once
#ifndef OCTREE_HPP_
#define OCTREE_HPP_

#include <memory>
#include <vector>
#include <queue>
#include <list>
#include <iostream>
#include <exception>
#include "BoundingBox.hpp"
#include "TriangleElement.hpp"


/// Stores info about element position in space and node
template <typename storedData>
struct TreeNode {
	/// Id of parallel node
	int nodeID;
	int level;
	
	Box<double, 3>* nodeSpan;

	/// Stores indices of element possition if vector storage
	std::vector<storedData> data;
	
	TreeNode* childs[8];
	
	/// Default constructor for empty node;
	TreeNode() {
		nodeID   = -1;
		level    = -1;
		nodeSpan = nullptr;
		for (int i = 0; i < 8; ++i)
			childs[i] = nullptr;
	}


	/// Recursive constructor used for prealocation of the whole tree
	///
	/// \param center center of the node
	/// \param halfWidth extent of the node
	/// \param level set maximum level of the tree
	TreeNode(double center[3], double halfWidth, int level) {

		nodeSpan    = new Box<double, 3>(center, halfWidth);
		nodeID      = -1;
		this->level = level;

		if (level > 0 ) {

			double offset[3];
			double step = halfWidth * 0.5;

			for (int i = 0; i < 8; ++i ) {

				offset[0] = (i & 1) ? step : -step;
				offset[1] = (i & 2) ? step : -step;
				offset[2] = (i & 4) ? step : -step;

				try {
					childs[i] = new TreeNode( offset, step, (level - 1) );
				} catch(const std::exception& e) {
					std::cerr << e.what() << std::endl;
				}

			}

		} else {
			for (int i = 0; i < 8; ++i) {
				childs[i] = nullptr;
			}
		}

	}

	/// Create an empty octree, no prealocation is done!!!
	TreeNode(Box<double, 3>& span, int level) {
	
		nodeSpan    = new Box<double, 3>(span);
		nodeID      = -1;
		this->level = level;

		for (int i = 0; i < 8; ++i) {
			childs[i] = nullptr;
		}
	
	}

	/// Deep copy
	TreeNode(const TreeNode& node) {

		nodeID   = node.nodeID;
		level    = node.level;
		nodeSpan = new Box<double, 3>(*node.nodeSpan);

		data = node.data;

		if (level > 0) {
			for (int i = 0; i < 8; ++i) {
				childs[i] = new TreeNode(*node.childs[i]);
			}
		} else {
			for (int i = 0; i < 8; ++i) {
				childs[i] = nullptr;
			}
		}
	
	}

	/// Assign the deep copy to noew structure
	TreeNode operator=(const TreeNode& node) {
		return TreeNode(node);
	}


	~TreeNode() {
		delete nodeSpan;

		for (int i = 0; i < 8; ++i) {
			if (childs[i] != nullptr)
				delete childs[i];
		}

	}


	/// Insert element but create only neccesary nodes.
	void insertObject(const Box<double, 3>& objectBoundingBox, storedData objectData) {

		int index = 0, straddle = 0;
		double offset[3], half[3];

		for (int i = 0; i < 3; ++i) {
		
			double d = objectBoundingBox.center(i) - this->nodeSpan->center(i);
			if ( std::abs(d) < objectBoundingBox.extent(i) ) {
				straddle = 1;
				break;
			}

			half[i] = nodeSpan->extent(i) / 2;

			if (d > 0) {
				index |= (1 << i);
				offset[i] =  half[i];
			} else {
				offset[i] = -half[i];
			}
		
		}

		if ( !straddle && (level != 0) ) {

			if (childs[index] == nullptr) {
				Box<double, 3> childSpan(*nodeSpan);
				childSpan.shiftCenter(offset);
				childSpan.grow(0.5);
				childs[index] = new TreeNode(childSpan, level-1);
			}

			childs[index]->insertObject(objectBoundingBox, objectData);
		} else {
			data.push_back(objectData);
		}

	}


	std::vector<storedData> findClosest(Eigen::Vector3d& point) {

		// found closest matches in tree
		std::vector<storedData> closest;

		int index = 0;

		TreeNode<storedData>* currentNode = this;

		while (currentNode != nullptr) {

			for (int i = 0; i < 3; ++i) {

				double d = point[i] - currentNode->nodeSpan->center(i);

				if (d > 0) {
					index |= (1 << i);
				}

			}

			currentNode = childs[index];

		}

		return closest;

	}


	/// Loop throught the tree and return depth first order of the elements
	static std::vector<storedData> depthFirst(TreeNode& node) {

		std::vector<storedData> retList;

		retList.insert(retList.end(), node.data.begin(), node.data.end());

		for (int i = 0; i < 8; ++i) {
			if (node.childs[i] != nullptr) {
				std::vector<int> tmpRet = depthFirst(*node.childs[i]);
				retList.insert(retList.end(), tmpRet.begin(), tmpRet.end());
			}
		}

		return retList;

	}
	
	/// Loop through the tree and return the breadth first list of the elements
	static std::vector<storedData> breadthFirst(TreeNode& node) {

		std::vector<storedData> retList;

		std::queue<TreeNode* >  nodeBuffer;

		// put inputting node (take it as root node)
		nodeBuffer.push(&node);


		while ( !nodeBuffer.empty() ) {

			TreeNode* workingNode = nodeBuffer.front();
			nodeBuffer.pop();

			// process the node
			retList.insert(retList.end(), workingNode->data.begin(), workingNode->data.end());

			// check for children and put them in queue
			for (auto child: workingNode->childs) {
				if (child != nullptr) {
					nodeBuffer.push(child);
				}
			}

		}

		return retList;


	}


};






// ostream operators used for debug printing
template<typename storedData>
std::ostream& operator<<(std::ostream& os, const TreeNode<storedData>& node) {
	
	os << "Node on level (" << node.level << "), with span : ";
	os << *node.nodeSpan << "Elements: ";

	for (auto d: node.data) {
		os << d << ", ";
	}

	os << std::endl;

	if (node.level > 0) {
		for (int i = 0; i < 8; ++i)
			if (node.childs[i] != nullptr)
				os << *node.childs[i];
	}

    return os;
}


#endif // OCTREE_HPP_
