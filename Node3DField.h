#ifndef NODE3DFIELD_H
#define NODE3DFIELD_H

class Node3DField{
	public:	
		double        field[3]    = {0.0,0.0,0.0};
		unsigned char material    = 0;
		double 		  Temperature = 0;
		double        mu          = 0.0;
		double        epsilon     = 0.0;
		double        conductivity= 0.0;
};

#endif