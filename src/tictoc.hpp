#pragma once 

#include <petsctime.h>
#include <string>
#include <stack>

///
/// Profiling utility,  prints time into the standard output.
///
struct tictoc {
	tictoc() : prefix("#tictoc") {};
	tictoc(const std::string& pref) : prefix(pref) {};

	void tic() {
		PetscTime( &time_now );
		time.push( time_now );
	}


	void toc(int rank, int catID, const std::string& catName) {
		PetscTime(&time_now);
		std::cout << prefix << " " << rank << ", " << catID << ", " << catName << ", " << (time_now - time.top()) << std::endl;
		time.pop();
	}

private:
	std::string prefix;
	PetscLogDouble time_now;
	std::stack<PetscLogDouble> time;
};
