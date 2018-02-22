/* This class calculate the electromagnetic behaviour   */
#ifndef ALGO_ELECTRO_H //
#define ALGO_ELECTRO_H


#include <cstdio>
#include <iostream>
#include <vector>
#include <map>
#include <string>

#include "Array_3D_Template.h"
#include "Node3DField.h"
#include "Materials.h"
#include "MPI_Initializer.h"

class AlgoElectro{
	private:
		int numberOfMaterials;
	public:
		// Constructor:
		AlgoElectro(int numberOfMaterials){
			this->numberOfMaterials = numberOfMaterials;
		}
		// Destructor:
		~AlgoElectro(void);

		void run(GridCreator&,MPI_Initializer&);
		void update(GridCreator&, double, double);
		// For the function communication, the class need an MPI communicator
		void communicate(GridCreator&,MPI_Initializer&);

		double Compute_dt(GridCreator &);

};

#endif
