/* Constructor: */
template<typename T>
Array_3D_Template<T>::Array_3D_Template(const size_t &Nx, const size_t &Ny, const size_t &Nz){
	this->Nx = Nx;
	this->Ny = Ny;
	this->Nz = Nz;
	this->data = new T [Nx*Ny*Nz];
	#if DEBUG > 2
	cout << "ARRAY_3D_TEMPLATE::constructor::Array successfully created.\n";
	#endif
}

/* Operator () */
template<typename T>
T& Array_3D_Template<T>::operator()(const size_t &i, const size_t &j, const size_t &k)
{
	return *(this->data + i * this->Nx + j*this->Ny + k);
}

/* Destructor */
template<typename T>
Array_3D_Template<T>::~Array_3D_Template(void){
	if(this->data != NULL){
		delete[] this->data;
		this->data = NULL;
	}
	#if DEBUG > 2
	cout << "ARRAY_3D_TEMPLATE::destructor::Array successfully deleted.\n";
	#endif
}

/* Fill in the field data with a given value */
template<typename T>
void Array_3D_Template<T>::fillIn(T val){
	for(unsigned long i = 0 ; i < this->Nx*this->Ny*this->Nz; i++)
			this->data[i] = val;
}

// Get the size of the array:
template<typename T>
vector<size_t> Array_3D_Template<T>::get_size_data(void){
	vector<size_t> sizes = {this->Nx,this->Ny,this->Nz};
	return sizes;
}

// Copy constructor:
template<typename T>
Array_3D_Template<T>::Array_3D_Template(const Array_3D_Template &obj){
	#if DEBUG > 0
	cout << "Array_3D_Template<T>:: you're calling the copy constructor.\n";
	#endif
	// Copying the field 'data':
	data  = new T [obj.Nx*obj.Ny*obj.Nz];
	*data = *obj.data;
	// Copying the sizes:
	Nx = obj.Nx;
	Ny = obj.Ny;
	Nz = obj.Nz;
}