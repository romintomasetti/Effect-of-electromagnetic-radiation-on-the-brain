/* Constructor: */
template<typename T>
Array_3D_Template<T>::Array_3D_Template(const size_t &Nx, const size_t &Ny, const size_t &Nz){
	this->Nx = Nx;
	this->Ny = Ny;
	this->Nz = Nz;
	try{
		this->data = new T [Nx*Ny*Nz];
	}catch(std::bad_alloc& exc){
		printf("File %s at %d:\n\t",__FILE__,__LINE__);
		cout << "Array_3D_Template<T>::set_size_data::Failed to allocate. Aborting.\n";
		std::cerr << "bad_alloc caught: " << exc.what() << '\n';
		abort();
	}
	this->dataAlreadySet = true;
	#if DEBUG > 2
	cout << "ARRAY_3D_TEMPLATE::constructor::Array successfully created.\n";
	#endif
}

/* set_size_array */
template<typename T>
void Array_3D_Template<T>::set_size_data(const size_t &Nx, const size_t &Ny, const size_t &Nz){
	if(this->dataAlreadySet == false){
		this->Nx = Nx;
		this->Ny = Ny;
		this->Nz = Nz;
		try{
			this->data = new T [Nx*Ny*Nz];
		}catch(std::bad_alloc& exc){
			printf("File %s at %d:\n\t",__FILE__,__LINE__);
			cout << "Array_3D_Template<T>::set_size_data::Failed to allocate. Aborting.\n";
			std::cerr << "bad_alloc caught: " << exc.what() << '\n';
			abort();
		}
		this->dataAlreadySet = true;
		#if DEBUG > 2
		cout << "ARRAY_3D_TEMPLATE::set_size_data::Array successfully created.\n";
		#endif
	}
}

/* Operator () */
template<typename T>
T& Array_3D_Template<T>::operator()(const size_t i, const size_t j, const size_t k)
{
	if( (i + this->Nx * ( j + k * this->Ny )) > this->Nx*this->Ny*this->Nz){
		printf("T& Array_3D_Template<T>::operator()::ERROR\n");
		printf("Size of the array is Nx*Ny*Nz=%ld and you provide %ld (provided : [%ld,%ld,%ld]).\n",this->Nx*this->Ny*this->Nz,
							i + this->Nx * ( j + k * this->Ny ),i,j,k);
		std:abort();
	}
	return this->data[ i + this->Nx * ( j + k * this->Ny )];
}

// Accessor [] overload:
template<typename T>
T& Array_3D_Template<T>::operator[](const size_t index){
	if(index < 0 || index > this->Nx*this->Ny*this->Nz){
		printf("Array_3D_Template::ERROR::index out of bound\n");
		abort();
	}
	return this->data[index];
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
	#if DEBUG > 3
	cout << "FILLIN:: total is " << this->Nx*this->Ny*this->Nz << endl;
	#endif
	for(unsigned long i = 0 ; i < this->Nx*this->Ny*this->Nz; i++){
			this->data[i] = val;
			#if DEBUG > 3
			cout << this->data[i] << "(" << i << "),";
			#endif
	}
	#if DEBUG > 3
	cout << endl;
	#endif
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