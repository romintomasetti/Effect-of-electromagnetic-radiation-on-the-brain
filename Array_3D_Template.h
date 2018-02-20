#ifndef ARRAY_3D_TEMPLATE_H
#define ARRAY_3D_TEMPLATE_H

#include <cstdio>
#include <iostream>
#include <vector>

using namespace std;

template <typename T>
class Array_3D_Template{
	private:
		// Data:
		T *data = NULL;
		// Sizes along x, y and z:
		size_t Nx, Ny, Nz;
		// Boolean to check if data has already been allocated:
		bool dataAlreadySet = false;
	public:
		// Default constructor:
		Array_3D_Template(){};
		// Constructor:
		Array_3D_Template(const size_t &Nx, const size_t &Ny, const size_t &Nz);
		// Copy constructor:
		Array_3D_Template(const Array_3D_Template &obj);
		// Destructor:
		~Array_3D_Template(void);
		// Accessor:
		T& operator()(const size_t i, const size_t j, const size_t k);
		// Filler (fill the data field with a given value):
		void fillIn(T val);
		// Get the size of the array:
		vector<size_t> get_size_data(void);
		// Set the size of the array:
		void set_size_data(const size_t &Nx, const size_t &Ny, const size_t &Nz);
};

#include "Array_3D_Template.tpp"

#endif