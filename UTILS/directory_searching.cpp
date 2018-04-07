#include "directory_searching.hpp"

#include <boost/range/iterator_range.hpp>

#include "header_with_all_defines.hpp"


using namespace boost::filesystem;

/**
 * REcursively searches for a directory, backward propagation.
 */
UTILS_API int search_for_directory(
  std::string const &dir,
  boost::filesystem::path &foundPath)
{
  /// Check that argument dir is a directory:
  std::string isa = what_is(dir);

  /// If it is not a directory, search for it recursively.
  if( isa != "directory"){
    /// Get current path:
    boost::filesystem::path prev(boost::filesystem::current_path());
    /// Iterate over the directories (backward propagation):
    while(true){
      prev = prev.parent_path();
      if(is_directory(prev)) {
        /// Iterate over the members of prev:
        for(auto& entry : boost::make_iterator_range(directory_iterator(prev), {})){
          if(is_directory(entry)){
            std::size_t found = entry.path().string().find_last_of("/");
            std::string f = entry.path().string().substr(found+1);
            if(f == dir){
              foundPath = entry;
              return EXIT_SUCCESS;
            }              
          }
        }
      }else{
        return EXIT_FAILURE;
      }
    }
  }else{
    return EXIT_SUCCESS;
  }
  return EXIT_FAILURE;
}

/**
 * Returns the type of name: file, directory, do not exist.
 */
UTILS_API std::string what_is(std::string const &name){
  boost::filesystem::path p (name);
  try{
    if (boost::filesystem::exists(p)){
      if (boost::filesystem::is_regular_file(p)){ 
        return "file";
      }else if (boost::filesystem::is_directory(p)){
        return "directory";
      }else{
        return "unknown";
      }
    }else{
      return "dont_exist";
    }
  }catch (const boost::filesystem::filesystem_error& ex)
  {
    std::cout << ex.what() << '\n';
  }
  return "error";
}

struct path_leaf_string
{
  std::string operator() (const boost::filesystem::directory_entry &entry) const
  {
    return entry.path().string();
  }
};

UTILS_API void read_directory_for_files(
  const std::string &name,
  std::vector<std::string> &vec)
{
  boost::filesystem::path p(name);
  if(!is_directory(p)){
    DISPLAY_ERROR_ABORT("%s is not a valid directory.",name.c_str());
  }
  boost::filesystem::directory_iterator start(p);
  boost::filesystem::directory_iterator end;
  std::transform(start,end,std::back_inserter(vec),path_leaf_string());
}
