/* This is a template for providing variables that can be set only once */
/* That is, when it has been set once, it can never be changed.         */
#ifndef SETONCEVARIABLE_TEMPLATE_H
#define SETONCEVARIABLE_TEMPLATE_H

#include <iostream>

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
			std::cout << "\n\n\tCoucou ça plante !!!!!\n\n\n";
			this->alreadySet = true;
			std::cout << "ça n'a pas planté " << std::endl;
		};
		// Destructor:
		~SetOnceVariable_Template(void){};
		// Set the value / overloading of the '=' operator:
		SetOnceVariable_Template<T>& operator=(const T& value){
			if(this->alreadySet == false){
				this->value = value;
				this->alreadySet = true;
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