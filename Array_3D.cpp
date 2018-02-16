#include "Array_3D.h"



string Array_3D::set_sizes(unsigned long Nx, unsigned long Ny, unsigned long Nz){
	if(this->SIZES[0] == 0 && this->SIZES[1] == 0 && this->SIZES[2] == 0){
		// The sizes were not yet set, proceed.
		this->SIZES[0] = Nx;
		this->SIZES[1] = Ny;
		this->SIZES[2] = Nz;
		return string();
	}else{
		return "Sizes already set";
	}
}

void fillVector(double *vec,unsigned long size,double value){
	for(unsigned long i = 0 ; i < size ; i ++){
		vec[i] = value;
		printf("vec[%ld] = %f.\n",i,vec[i]);
	}
}

string Array_3D::init(double val){
	if(this->values == NULL){
		// The array has never been initialized.
		this->values = new (nothrow) double[this->SIZES[0]*this->SIZES[1]*this->SIZES[2]];
		fillVector(this->values,this->SIZES[0]*this->SIZES[1]*this->SIZES[2],val);
		return string();
	}else{
		// The array has already been initialized.
		return "Already been initialized";
	}
}

string Array_3D::free(void){
	if(this->values == NULL){
		// There is nothing to free.
		return "Nothing to free";
	}else{
		delete[] this->values;
		return string();
	}
}

double Array_3D::get_value(unsigned long m, unsigned long n, unsigned long p){
	if(this->values != NULL && m < this->SIZES[0] && n < this->SIZES[1] && p < this->SIZES[2]){
		return this->values[m*this->SIZES[0]+n*this->SIZES[1]+p];
	}else{
		abort();
	}
}

bool Array_3D::set_value(unsigned long m, unsigned long n, unsigned long p, double val){
	if(m < this->SIZES[0] && n < this->SIZES[1] && p < this->SIZES[2]){
		this->values[m*this->SIZES[0]+n*this->SIZES[1]+p] = val;
		return true;
	}else{
		return false;
	}
}

string Array_3D::get_sizes(vector<unsigned long> &vect){
	vect = this->SIZES;
	return "Success";
}
