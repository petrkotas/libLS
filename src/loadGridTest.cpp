
#include <vector>
#include <iostream>
#include <algorithm>
#include <iterator>
#include <tuple>
#include <Eigen/Dense>
#include "Geometry.hpp"


#undef __FUNCT__
#define __FUNCT__ "main"
int main(int argc,char **argv) {

	bool write_out = true;

	Eigen::Array3i M;            // number of nodes in each dimension

	char           iname[120];
	char           oname[120];

	Geometry geom = Geometry::fromFile("/mnt/hgfs/Dropbox/testGeometries/dragon.stl");


	std::cout << "D O N E" << std::endl;

	return 0;


}
