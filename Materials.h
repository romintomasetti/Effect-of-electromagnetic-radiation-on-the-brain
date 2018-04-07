#ifndef MATERIALS_H
#define MATERIALS_H

#include <vector>
#include <unordered_set>
#include <string>
#include <iostream>
#include <map>
#include <limits.h>
#include <iomanip>
#include <cmath>

#include "header_with_all_defines.hpp"

#include "Array_3D_Template.h"

#include "UTILS/directory_searching.hpp"

#include <boost/filesystem.hpp>

using namespace std;

typedef struct material_struct{
	std::string name;
	unsigned int ID;
	map<std::string,double> properties;
	double initial_temperature;

	inline void printf_mat(void){
		printf("\t> %s [ID %u]:\n",name.c_str(),ID);
		for(auto it = properties.cbegin(); it != properties.cend() ; it ++)
			printf("\t\t> %30s:%10lf\n",it->first.c_str(),it->second);
		printf("\n");
	}
}material_struct;

class Materials{
	private:
		bool unification_done = false;
		// Contains all the properties of all materials:
		Array_3D_Template<double> properties;
		// Number of properties:
		unsigned int numberOfProperties = 0;
		
		// Maximum number of temperature specifications:
		unsigned int maxNumberOfTemp    = 0;
		// Dictionary with the materials and the chosen unsigned char assigned to it:
		//map<string,unsigned char> materialID_FromMaterialName;
		
		
		vector<unsigned int> numberOFTempForTheMaterial;
	
		// Free the properties array (called in the destructor):
		//void   freeProperties(void);
	public:

		void printf_list_of_mat_from_dir(void);

		void get_properties_from_directory(
			boost::filesystem::path pathToData
		);

		// After the reading of the directory, unify all properties inside one field:
		std::vector<material_struct> unified_material_list;
		std::map<std::string,unsigned int> materialID_FromMaterialName_unified;
		std::map<unsigned int,std::string> materialName_FromMaterialID_unified;
		void unification_of_data_files(void);

		map<std::string,std::vector<material_struct> > list_of_mat_from_dir;
		std::vector<std::string> list_of_mat_files;
		// First, give file name, then material name:
		map<std::string,map<std::string,unsigned char> > materialID_FromMaterialName_from_dir;
		// First, give file name, then material ID:
		map<std::string,map<unsigned char,std::string> > materialName_FromMaterialID_from_dir;

		std::vector<material_struct> list_of_materials_ELECTRO;
		std::vector<std::string> list_of_properties_of_list_of_materials_ELECTRO;

		map<string,unsigned char> materialID_FromMaterialName;
		map<unsigned char,string> materialName_FromMaterialID;

		map<string,unsigned char> materialID_FromMaterialName_ELECTRO;
		map<unsigned char,string> materialName_FromMaterialID_ELECTRO;

		// Get all the properties specified in a file, and put them in a 3D array:
		void   getPropertiesFromFile(string,string);
		// Get a property for a given material at a given temperature:
		double getProperty(double, unsigned char, unsigned char,bool interpolation = false);
		// Print all the properties:
		void   printAllProperties(void);
		// Print the number of temperature lines per material:
		void   printNumberOfTempLinePerMat(void);
		// Number of materials (pour AlgoElectro.cpp)
		unsigned int numberOfMaterials  = 0;

		inline double get_init_temp(std::string const &mat_name);

		// This is given by the input parser:
		map<std::string,double> GetInitTemp_FromMaterialName;
		// Constructor:
		Materials(
			std::string dir,
			map<std::string,double> GetInitTemp_FromMaterialName
		){
			this->GetInitTemp_FromMaterialName = GetInitTemp_FromMaterialName;
			this->get_properties_from_directory(dir);
			if(this->unification_done != true)
				DISPLAY_ERROR_ABORT(
					"Unification of material data files failed."
				);
		};
		// Destructor:
		~Materials();
		// Get Dictionnary with the materials and the chosen unsigned char assigned to it:
		map<string,unsigned char> get_dictionnary_MaterialToID(void);
		map<unsigned char,string> get_dictionnary_IDToMaterial(void);

		void get_properties_from_file_ELECTRO(std::string const &filename_ELECTRO_PROPS);

		void get_properties_from_file_ALL(std::string const &filename);

		/**
		 * @brief Print all informations knows in list_of_mat_from_dir 
		 * 	about a specific material.
		 */
		void printf_all_on_one_mat_from_dir(std::string const &mat);

		/**
		 * Format all possible strings from material data files to get uniformity.
		 */
		inline void format_property_string_for_uniformity(std::string &prop);
		
		/**
		 * @brief Remove duplicated materials.
		 * Example: Bone(Cortical) and BoneCortical are the same !
		 */
		inline void remove_duplicated_material_names(std::string &mat);

		/**
		 * @brief Returns a vector without any duplicated value.
		 */
		template<typename T>
		inline void unique_vector(std::vector<T> &vec);

		/**
		 * @brief Looks inside the directory of properties for a 
		 * 	specific material property. If a property is reapeated twice (or more)
		 * 	in the files, and if they are different, take the mean and print a warning.
		 * 	If not found, returns a nan.
		 */
		double get_mean_prop_from_dir(std::string const &mat, std::string const &prop);
};

#endif
