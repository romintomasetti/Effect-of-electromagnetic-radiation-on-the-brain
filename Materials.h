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

#include "UTILS/directory_searching.hpp"

#include <boost/filesystem.hpp>

#include <climits>

using namespace std;

/**
 * @brief Structure that contains all relevant informations about a given material.
 */
typedef struct material_struct{
	std::string name;
	unsigned int ID;
	map<std::string,double> properties;
	double initial_temperature;

	inline void printf_mat(void){
		printf("\t> %s [ID %u]:\n",name.c_str(),ID);
		for(auto it = properties.cbegin(); it != properties.cend() ; it ++)
			printf("\t\t> %30s:%10lf\n",it->first.c_str(),it->second);
		printf("\t\t> Initial temperature is %.9g.\n",initial_temperature);
		printf("\n");
	}
}material_struct;


/**
 * @brief Class materials. Used to parse the material data files and to create a structure
 * 		for later use of these data.
 */
class Materials{

	/** Private objects and fields **/
	private:
		/** Rank of the MPI process. **/
		int MPI_RANK = INT_MIN;
		/** Rank of the root process **/
		int MPI_root = INT_MIN;
		/** Verbosity. **/
		unsigned int VERBOSITY = 0;
		/** Check that unification of all material files was performed sccessfully. **/
		bool unification_done = false;
		
	/** Public functions and objects and fields. **/
	public:

		/** Print the material structure generated from all data files. **/
		void printf_list_of_mat_from_dir(void);

		/** Reads data files from a directory and create a clean structure. **/
		void get_properties_from_directory(
			boost::filesystem::path pathToData
		);

		// After the reading of the directory, unify all properties inside one field:
		std::vector<material_struct> unified_material_list;
		std::map<std::string,unsigned int> materialID_FromMaterialName_unified;
		std::map<unsigned int,std::string> materialName_FromMaterialID_unified;

		void unification_of_data_files(void);

		map<std::string,std::vector<material_struct> > list_of_mat_from_dir;

		/** List of all material data files. **/
		std::vector<std::string> list_of_mat_files;

		// First, give file name, then material name:
		map<std::string,map<std::string,unsigned char> > materialID_FromMaterialName_from_dir;
		// First, give file name, then material ID:
		map<std::string,map<unsigned char,std::string> > materialName_FromMaterialID_from_dir;


		inline double get_init_temp(std::string const &mat_name);

		// This is given by the input parser:
		map<std::string,double> GetInitTemp_FromMaterialName;
		// Constructor:
		Materials(
			unsigned int VERBOSITY,
			int MPI_RANK,
			int MPI_root,
			std::string dir,
			map<std::string,double> GetInitTemp_FromMaterialName
		){
			// Set verbosity:
			this->VERBOSITY = VERBOSITY;
			// Set the MPI process rank and the root's rank:
			this->MPI_RANK = MPI_RANK;
			this->MPI_root = MPI_root;
			// Gt the initial temperatures:
			this->GetInitTemp_FromMaterialName = GetInitTemp_FromMaterialName;
			// Read the data files nd create the structure of materials:
			this->get_properties_from_directory(dir);
			// Check that all was done properly:
			if(this->unification_done != true)
				DISPLAY_ERROR_ABORT(
					"Unification of material data files failed."
				);
		};
		/** Destructor. **/
		~Materials(void);
		// Get Dictionnary with the materials and the chosen unsigned char assigned to it:
		map<string,unsigned char> get_dictionnary_MaterialToID(void);
		map<unsigned char,string> get_dictionnary_IDToMaterial(void);

		void get_properties_from_file_ELECTRO(std::string const &filename_ELECTRO_PROPS);

		void get_properties_from_file_ALL(
			std::string const &filename,
			double frequency
		);

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

		/**
		 * @brief The function reads the filename. If the frequency is specified, it
		 * 		returns the frequency;
		 */
		double find_frequency_from_mat_file_name(std::string filename);

		void printf_unified_material_list(void);
};

#endif
