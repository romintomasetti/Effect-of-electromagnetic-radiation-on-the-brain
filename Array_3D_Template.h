#include <cstdio>
#include <iostream>

using namespace std;

template <typename T>
class Array_3D_Template{
	private:
		// Data:
		T *data = NULL;
		// Sizes along x, y and z:
		size_t Nx, Ny, Nz;
	public:
		// Constructor:
		Array_3D_Template(const size_t &Nx, const size_t &Ny, const size_t &Nz);
		// Destructor:
		~Array_3D_Template(void);
		// Accessor:
		T& operator()(const size_t &i, const size_t &j, const size_t &k);
		// Filler (fill the data field with a given value):
		void fillIn(T val);
};

#include "Array_3D_Template.tpp"
