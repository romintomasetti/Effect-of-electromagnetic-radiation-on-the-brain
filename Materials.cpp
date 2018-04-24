#include "Materials.h"
#include <fstream>
#include <cstring>
#include <algorithm>

#include <boost/algorithm/string.hpp>
//#include <boost/range/algorithm.hpp>

#include "CSV_parser.hpp"
#include "header_with_all_defines.hpp"

template<typename T>
inline bool is_inside(std::vector<T> &vector, T item){
	if(std::find(vector.begin(), vector.end(), item)!=vector.end()){
		return true;
	}else{
		return false;
	}
}

/**
 * @brief Comparison between two strings, not case sensitive.
 */
bool string_compare_case_insensitive (std::string str1, std::string str2)
{
  size_t i = 0;
  while ((i < str1.length()) && (i < str2.length()))
  {
    if (tolower (str1[i]) < tolower (str2[i])) return true;
    else if (tolower (str1[i]) > tolower (str2[i])) return false;
    i++;
  }

  if (str1.length() < str2.length()) return true;
  else return false;
}

inline void Materials::remove_duplicated_material_names(std::string &mat){
	if(  boost::iequals(mat,"BoneCortical")
	  || boost::iequals(mat,"Bone(Cortical)")){
		  mat = "BoneCortical";
	  }
	if(  boost::iequals(mat,"BoneCancellous")
	  || boost::iequals(mat,"Bone(Cancellous)")){
		  mat = "BoneCancellous";
	  }
	if(  boost::iequals(mat,"Lung(Deflated)")
	  || boost::iequals(mat,"LungDeflated")){
		  mat = "LungDeflated";
	  }
	if(  boost::iequals(mat,"Lung(Inflated)")
	  || boost::iequals(mat,"LungInflated")){
		  mat = "LungInflated";
	  }
	if(  boost::iequals(mat,"GallBladder")){
		mat = "GallBlader";
	}
	if(boost::iequals(mat,"CerebroSpinalFluid")){
		mat = "CerebroSpinalFluid";
	}
	if(  boost::iequals(mat,"Brain(GreyMatter)")
	  || boost::iequals(mat,"BrainGreyMatter")){
		  mat = "BrainGreyMatter";
	  }
	if(  boost::iequals(mat,"Brain(WhiteMatter)")
	  || boost::iequals(mat,"BrainWhiteMatter")){
		  mat = "BrainWhiteMatter";
	  }
}

/**
 * This function unifies all the data files.
 */
void Materials::unification_of_data_files(void){
	/// Determine the number of differet materials in all files:
	std::vector<std::string> materials_list_unique;
	std::vector<std::string> materials_prop_unique;
	/// Iterate over the files:
	for(auto it_files = this->list_of_mat_from_dir.cbegin();
			it_files != this->list_of_mat_from_dir.cend();
			it_files++)
	{
		/// Go over the materials:
		std::string curr_mat;
		std::string curr_prop;
		size_t nbr_mat = this->list_of_mat_from_dir[it_files->first].size();
		for(size_t it_mat = 0 ; it_mat < nbr_mat ; it_mat ++){
			curr_mat = it_files->second[it_mat].name;
			if(!is_inside(materials_list_unique,curr_mat)){
				// The material is not in the list yet. Add it.
				materials_list_unique.push_back(curr_mat);
			}
			for(auto it = it_files->second[it_mat].properties.cbegin();
				it != it_files->second[it_mat].properties.cend() ; 
				it ++)
			{
				curr_prop = it->first;
				if(!is_inside(materials_prop_unique,curr_prop)){
					materials_prop_unique.push_back(curr_prop);
				}
			}
		}
	}
	/// Remove the duplicates that are very similar:
	for(size_t I = 0 ; I < materials_list_unique.size() ; I++ ){
		this->remove_duplicated_material_names(materials_list_unique[I]);
	}
	this->unique_vector(materials_list_unique);
	/// Sort by alphabetical order:
    sort(materials_list_unique.begin(),materials_list_unique.end(),string_compare_case_insensitive);
	
	/*printf("There are %zu different materials.\n",materials_list_unique.size());
	for(size_t I = 0 ; I < materials_list_unique.size() ; I++)
		printf("\t>Material %zu is %s.\n",I,materials_list_unique[I].c_str());
	for(size_t I = 0 ; I < materials_prop_unique.size() ; I++)
		printf("\t>Property %zu is %s.\n",I,materials_prop_unique[I].c_str());
	abort();*/
	/**
	 * Create the final list of materials and properties.
	 */
	/// Loop over the materials:
	for(size_t mat = 0 ; mat < materials_list_unique.size() ; mat++){
		material_struct temp_mat;
		temp_mat.name               = materials_list_unique[mat];
		temp_mat.ID                 = mat;
		temp_mat.initial_temperature= this->get_init_temp(temp_mat.name);
		/// Loop over the properties:
		std::string propName;
		double      propValue;
		for(size_t prop = 0 ; prop < materials_prop_unique.size() ; prop ++){
			propName = materials_prop_unique[prop];
			/// Look for the property inside all files:
			propValue = this->get_mean_prop_from_dir(temp_mat.name,propName);
			std::pair<std::string,double> temp_pair(propName,propValue);
			temp_mat.properties.insert(temp_pair);
		}
		//temp_mat.printf_mat();
		this->unified_material_list.push_back(temp_mat);
		this->materialID_FromMaterialName_unified
			.insert(
				std::pair<std::string,unsigned int>(temp_mat.name,temp_mat.ID)
				);
		this->materialName_FromMaterialID_unified
			.insert(
				std::pair<unsigned int,std::string>(temp_mat.ID,temp_mat.name)
				);
	}
}

std::string format_vector_double(std::vector<double> const &vec){
	std::stringstream ret;
	for(size_t I = 0 ; I < vec.size() ; I ++){
		ret << std::setprecision(5) << vec[I] << ";";
	}
	return ret.str();
}

/**
 * @brief Looks inside the directory of properties for a 
 * 	specific material property. If a property is reapeated twice (or more)
 * 	in the files, and if they are different, take the mean and print a warning.
 * 	If not found, returns a nan.
 */
double Materials::get_mean_prop_from_dir(std::string const &mat, std::string const &prop){
	std::vector<double> found_values;
	/// Loop over the files:
	for(auto it_files = this->list_of_mat_from_dir.cbegin();
			it_files != this->list_of_mat_from_dir.cend();
			it_files++)
	{
		/// Go over the materials:
		std::string curr_mat;
		std::string curr_prop;
		size_t nbr_mat = this->list_of_mat_from_dir[it_files->first].size();
		for(size_t it_mat = 0 ; it_mat < nbr_mat ; it_mat ++){
			curr_mat = it_files->second[it_mat].name;
			if(curr_mat == mat){
				/// Go over the properties:
				map<std::string,double> props = it_files->second[it_mat].properties;
				for(auto it_props = props.cbegin();it_props != props.cend(); it_props++){
					if(it_props->first == prop){
						// The property is ound:
						found_values.push_back(it_props->second);
					}
				}
			}
		}
	}
	if(found_values.size() > 1){
		if ( std::adjacent_find( 
				found_values.begin(), 
				found_values.end(), 
				std::not_equal_to<double>() ) == found_values.end() )
		{
			//std::cout << ">>> All elements are equal each other" << std::endl;
			return found_values[0];
		}
		if(this->MPI_RANK == this->MPI_root && this->VERBOSITY > 1)
			DISPLAY_WARNING(
				"For material %s, property %s : I found more than one value so I take the mean."
				" Found values are [%s%s%s]",
				mat.c_str(),prop.c_str(),
				ANSI_COLOR_RED,
				format_vector_double(found_values).c_str(),
				ANSI_COLOR_RESET
			);
		return accumulate( found_values.begin(), found_values.end(), 0.0)/found_values.size();
	}else if(found_values.size() == 1){
		return found_values[0];
	}
	return nan("");

}

inline double Materials::get_init_temp(std::string const &mat_name){
	// Looks inside GetInitTemp_FromMaterialName for the initial temperature of material 'mat_name':
	std::map<std::string,double>::iterator it;
	it = this->GetInitTemp_FromMaterialName.find(mat_name);
	if(it != this->GetInitTemp_FromMaterialName.end()){
		// Material is found into initial temperature list:
		/*printf("> Material %s with init temp %lf.\n",
			mat_name.c_str(),
			this->GetInitTemp_FromMaterialName[mat_name]);*/
		return this->GetInitTemp_FromMaterialName[mat_name];
	}
	/*printf("> material %s : default temp is %lf.\n",
		mat_name.c_str(),
		this->GetInitTemp_FromMaterialName["default"]);*/
	it = this->GetInitTemp_FromMaterialName.find("default");
	if( it == this->GetInitTemp_FromMaterialName.end()){
		DISPLAY_ERROR_ABORT(
			"You must provide a default temperature in the .initTemp file."
		);
	}
	return this->GetInitTemp_FromMaterialName["default"];
}

template<typename T>
inline void Materials::unique_vector(std::vector<T> &vec){
	/**
	 * @brief This function removes duplicates.
	 */
	std::unordered_set<T> myset;
	myset.insert(vec.begin(), vec.end());
	while (!vec.empty())
    	vec.pop_back();
	for (const T& x: myset)
		vec.push_back(x);	
}

/**
 * Format all possible strings from material data files to get uniformity.
 */
inline void Materials::format_property_string_for_uniformity(std::string &prop)
{
	/// Thermal conductivity:
	if(boost::iequals(prop,"ThermalConductivity(W/m/°C)")){
		prop = "ThermalConductivity(W/m/°C)";
		boost::to_upper(prop);
		return;
	}
	/// Electrical conditivity:
	if(   boost::iequals(prop,"Elec.Cond.(S/m)")
	   || boost::iequals(prop,"Conductivity[S/m]")){
		prop = "ElectricalConductivity(S/m)";
		boost::to_upper(prop);
		return;
	}
	/// Permittivity:
	if(    boost::iequals(prop,"Permittivity")
		|| boost::iequals(prop,"Relativepermittivity")){
		prop = "RelativePermittivity";
		boost::to_upper(prop);
		return;
	}
	/// Heat capacity:
	if(boost::iequals(prop,"HeatCapacity(J/kg/°C)")){
		prop = "HeatCapacity(J/kg/°C)";
		boost::to_upper(prop);
		return;
	}
	/// Loss tangent:
	if(boost::iequals(prop,"LossTangent")){
		prop = "LossTangent";
		boost::to_upper(prop);
		return;
	}
	/// Density:
	if(boost::iequals(prop,"Density(kg/m³)")){
		prop = "Density(kg/m³)";
		boost::to_upper(prop);
		return;
	}
	/// Frequency:
	if(boost::iequals(prop,"Frequency[Hz]")){
		prop = "Frequency[Hz]";
		boost::to_upper(prop);
		return;
	}
	/// Wavelength:
	if(boost::iequals(prop,"Wavelength[m]")){
		prop = "Wavelength[m]";
		boost::to_upper(prop);
		return;
	}
	/// Properties we are not interested in:
	if(	   boost::iequals(prop,"StandardDeviation")
		|| boost::iequals(prop,"Maximum")
		|| boost::iequals(prop,"Minimum")
		|| boost::iequals(prop,"NumberOfStudies")
		|| boost::iequals(prop,"PenetrationDepth[m]")){
		prop = std::string();
		return;
	}
	DISPLAY_ERROR_ABORT(
		"No conversion for property named %s.",
		prop.c_str()
	);
}


void Materials::get_properties_from_directory(
	boost::filesystem::path pathToData
)
{
	if(!is_directory(pathToData)){
		DISPLAY_ERROR_ABORT(
			"%s is not a valid directory.",
			pathToData.string().c_str()
		);
	}
	/// Get all files from 'pathToData':
	std::vector<std::string> list_material_files;
	read_directory_for_files(
		pathToData.string(),
		list_material_files
	);

	for(size_t I = 0 ; I < list_material_files.size() ; I++){
		/// Consider only .csv files !
		if(boost::filesystem::extension(list_material_files[I]) != ".csv"){
			// Do nothing
		}else{
			this->list_of_mat_files.push_back(list_material_files[I]);
			//printf("Material file %zu is %s.\n",I,list_material_files[I].c_str());
			/// Read the file:
			this->get_properties_from_file_ALL(
				list_material_files[I]
			);
		}
	}
	
	//this->printf_list_of_mat_from_dir();
	//this->printf_all_on_one_mat_from_dir("Air");
	this->unification_of_data_files();

	this->unification_done = true;

}

/**
 * @brief Print all informations knows in list_of_mat_from_dir 
 * 	about a specific material.
 */
void Materials::printf_all_on_one_mat_from_dir(std::string const &mat)
{
	/// Iterate over the files:
	for(auto it_files = this->list_of_mat_from_dir.cbegin();
			it_files != this->list_of_mat_from_dir.cend();
			it_files++)
	{
		printf("%s> Material file %s:%s\n",
			ANSI_COLOR_GREEN,
			it_files->first.c_str(),
			ANSI_COLOR_RESET);
		size_t nbr_mat = this->list_of_mat_from_dir[it_files->first].size();
		printf("\t> Number of materials : %zu.\n",nbr_mat);
		std::string curr_mat;
		for(size_t it_mat = 0 ; it_mat < nbr_mat ; it_mat ++){
			curr_mat = it_files->second[it_mat].name;
			if(boost::iequals(curr_mat,mat)){
				this->list_of_mat_from_dir[it_files->first][it_mat].printf_mat();
			}else{
				//printf("\t> %s != %s.\n",mat.c_str(),curr_mat.c_str());
			}
		}
	}
}

/**
 * @brief Print the complicated field list_of_mat_from_dir.
 */
void Materials::printf_list_of_mat_from_dir(void){
	/// Iterate over the files:
	for(auto it_files = this->list_of_mat_from_dir.cbegin();
			it_files != this->list_of_mat_from_dir.cend();
			it_files++)
	{
		printf("%s> Material file %s:%s\n",
			ANSI_COLOR_GREEN,
			it_files->first.c_str(),
			ANSI_COLOR_RESET);
		/// Iterate over the materials:
		size_t nbr_mat = this->list_of_mat_from_dir[it_files->first].size();
		printf("\t> Number of materials : %zu.\n",nbr_mat);
		for(size_t it_mat = 0 ; it_mat < nbr_mat ; it_mat ++){
			this->list_of_mat_from_dir[it_files->first][it_mat].printf_mat();
		}
	}
}

void Materials::get_properties_from_file_ALL(
	std::string const &filename)
{
	/// Check that the file is a CSV file:
	if(filename.substr(filename.find(".")+1) != "csv"){
		DISPLAY_ERROR_ABORT("The property file %s has not a '.csv' extension.",filename.c_str());
	}
	/// Read the whole CSV file:
	std::vector<std::vector<std::string> > data;
	try{
		data = parse2DCsvFile(filename);
	}catch(...){
		DISPLAY_ERROR_ABORT(
			"There was an error in parse2DCsvFile with %s.",
			filename.c_str()
		);
	}
	/// Count the number of different materials:
	size_t nbr_mat = 0;
	std::string current_mat = data[1][0];
	std::vector<std::string> list_mat;
    for (size_t I = 1 ; I < data.size() ;  I++) {
		std::vector<std::string> l = data[I];
        if(l[0] != current_mat){
			/// Check that the material was not already encountered:
			std::vector<std::string>::iterator it;
			it = std::find (list_mat.begin(), list_mat.end(), l[0]);
			if (it == list_mat.end()){
				//std::cout << "New mat : " << l[0] << std::endl;
				std::pair<std::string,unsigned int> temp1(l[0],nbr_mat);
				std::pair<unsigned int,std::string> temp2(nbr_mat,l[0]);
				this->materialID_FromMaterialName_from_dir[filename].insert(temp1);
				this->materialName_FromMaterialID_from_dir[filename].insert(temp2);
				if(nbr_mat == 0 && l[0] != "Air"){
					DISPLAY_ERROR_ABORT(
						"The material with ID 0 MUST be air but as %s instead.",
						l[0].c_str()
					);
				}
				nbr_mat++;
				list_mat.push_back(l[0]);
				current_mat = l[0];
			}
		}
    }

	/// Check that the first column's name is something like 'material' or 'tissue':
	std::string f_col_name = data[1][0];
	boost::erase_all(f_col_name, "\"");
	if(    !boost::iequals(f_col_name,"Tissue")
		&& !boost::iequals(f_col_name,"Material")
		&& !boost::iequals(f_col_name,"Tissuename"))
	{
		DISPLAY_ERROR_ABORT(
			"Material file %s is not compliant with the required format."
			" The first column of the .csv file should have a name like"
			" material(s) or tissue(s). Has %s instead.",
			filename.c_str(),
			f_col_name.c_str()
		);
	}

	/// Count the number of properties:
	size_t nbr_prop = 0;
	/// Get the properties names (start for loop at 1 because first column is material name):
	std::vector<std::string> propNames;
	std::vector<bool>        propIgnore;
	for(size_t I = 1 ; I < data[1].size() ; I ++){
		this->format_property_string_for_uniformity(data[1][I]);
		if(data[1][I] != std::string()){
			propNames.push_back(data[1][I]);
			bool faux = false;
			propIgnore.push_back(faux);
			nbr_prop++;
		}else{
			bool vrai = true;
			propIgnore.push_back(vrai);
		}
	} 

	/// Fill in this->list_of_materials:
	std::vector<material_struct> temp(nbr_mat);
	for(size_t I = 0 ; I < nbr_mat ; I ++){
		temp[I].name = list_mat[I];
		temp[I].ID   = I;
		for(size_t J = 0 ; J < nbr_prop ; J ++){
			temp[I].properties.insert(
				std::pair<std::string,double>(propNames[J],0.0)
			);
		}
	}
	std::pair<std::string,std::vector<material_struct> > tempp(
		filename,temp
	);
	this->list_of_mat_from_dir.insert(tempp);

	for(size_t I = 1 ; I < data.size() ; I ++){
		std::vector<std::string> l = data[I];

		std::vector<std::string>::iterator it;
		it = std::find (list_mat.begin(), list_mat.end(), l[0]);

		if (it != list_mat.end()){
			std::string mat = l[0];
			size_t ID = this->materialID_FromMaterialName_from_dir[filename][mat];

			size_t counter = 0;
			//printf("Current mat %s :: Id %zu :: %zu :: %zu\n",mat.c_str(),ID,l.size(),propIgnore.size());
			for(size_t J = 1 ; J < l.size() ; J ++){
				if(propIgnore[J-1])
					continue;
				/// Identify the property:
				std::string propName = propNames[counter];
				//std::cout << propName + "=" + l[J] << std::endl;

				double propVal = 0.0;

				if(l[J] == "N/A" || l[J] == "NaN"){
					propVal = nan("");
					//printf(">>>> Has %lf...\n",propVal);
				}else{
					//std::cout << l[J] << std::endl;
					propVal = std::stod(l[J]);
				}
				
				this->list_of_mat_from_dir[filename][ID].properties[propName] = propVal;
				counter++;
			}
		}
	}
}


/* Destructor */
Materials::~Materials(void){
	if( this->VERBOSITY >= 3 )
		printf("[MPI %d] - Material destructor.\n",this->MPI_RANK);
}






