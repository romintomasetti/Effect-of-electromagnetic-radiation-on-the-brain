/* This is a template for providing variables that can be set only once */
/* That is, when it has been set once, it can never be changed.         */
#ifndef SETONCEVARIABLE_TEMPLATE_H
#define SETONCEVARIABLE_TEMPLATE_H

#include <iostream>
#include <typeinfo>
#include "stdio.h"
#include "stdlib.h>"

template<typename T>
class SetOnceVariable_Template{
	private:
		T value;
		bool alreadySet = false;
	public:
		// Default constructor:
		SetOnceVariable_Template(){};
		// Constructor:
		SetOnceVariable_Template(T init){
			this->value = init;
			this->alreadySet = true;
		};
		// Destructor:
		~SetOnceVariable_Template(void){};
		// Set the value / overloading of the '=' operator:
		SetOnceVariable_Template<T>& operator=(const T& value){
			if(this->alreadySet == false){
				this->value = value;
				this->alreadySet = true;
			}else{
				fprintf(stderr,"In %s :: %s :: variable has already been set. Aborting.\n",
					typeid(this).name().c_str(),
					__FUNCTION__);
				fprintf(stderr,"File %s:%d\n",__FILE__,__LINE__);
				abort();
			}
			return *this;
		}
		// Get the value:
		const T& get(void){return this->value;}

		//
		bool get_alreadySet(void){
			return this->alreadySet;
		}
};

#endif