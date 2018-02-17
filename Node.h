#ifndef NODE_H
#define NODE_H

#include <vector>

class Node{
	public:
		vector<double> coordinates   = {0.0,0.0,0.0};
		vector<double> electricField = {0.0,0.0,0.0};
		vector<double> magneticField = {0.0,0.0,0.0};
		double         temperature   = 0.0;
		int            material      = -1;
};

#endif