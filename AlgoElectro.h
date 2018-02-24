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

#include "InterfaceToParaviewer.h"

class AlgoElectro{
	private:
		unsigned long currentStep = 0;
	public:
		// Constructor:
		AlgoElectro(){}
		// Destructor:
		~AlgoElectro(void);

		void update(GridCreator&, InterfaceToParaviewer&);
		// For the function communication, the class need an MPI communicator
		void communicate(GridCreator&,MPI_Initializer&);

		double Compute_dt(GridCreator &);

};

#endif
