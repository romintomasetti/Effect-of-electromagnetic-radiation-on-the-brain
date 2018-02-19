#ifndef NODE_H
#define NODE_H

#include <vector>

class Node{
	public:
		double coordinates[]   = {0.0,0.0,0.0};
		double electricField[] = {0.0,0.0,0.0};
		double magneticField   = {0.0,0.0,0.0};
		double temperature     = 0.0;
		unsigned char material = 0;
};

#endif