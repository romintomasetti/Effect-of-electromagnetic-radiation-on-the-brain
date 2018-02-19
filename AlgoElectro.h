/* This class calculate the electromagnetic behaviour   */
#ifndef ALGO_ELECTRO_H //
#define ALGO_ELECTRO_H


#include <cstdio>
#include <iostream>
#include <vector>
#include <map>
#include <string>

#include "Array_3D_Template.h"
#include "Node.h"
#include "Materials.h"

class AlgoElectro{
	private:
		int numberOfMaterials;
		MPI_Initializer MPI;
	public:
		// Constructor:
		AlgoElectro(int numberOfMaterials){
			this->numberOfMaterials = numberOfMaterials;
		}
		// Destructor:
		~AlgoElectro(void);

		void run(GridCreator);
		void update(GridCreator);
		// For the function communication, the class need an MPI communicator
		void communicate(GridCreator);

}
