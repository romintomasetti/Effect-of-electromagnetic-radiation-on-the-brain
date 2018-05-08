#include <fstream>
#include <iterator>
#include <string>
#include <algorithm>
#include <iostream>

bool compareFiles(const std::string& p1, const std::string& p2) {
  std::ifstream f1(p1, std::ifstream::binary|std::ifstream::ate);
  std::ifstream f2(p2, std::ifstream::binary|std::ifstream::ate);

  if (f1.fail() || f2.fail()) {
    fprintf(stderr,"Problem opening one of the files.\n");
    if( f1.fail() ){
	fprintf(stderr,"\n\t>>> File %s cannot be opened.\n",p1.c_str());
    }else{
	fprintf(stderr,"\n\t>>> File %s cannot be opened.\n",p2.c_str());
    }
    return false; //file problem
  }

  if (f1.tellg() != f2.tellg()) {
    fprintf(stderr,"Size mismatch.\n");
    return false; //size mismatch
  }

  //seek back to beginning and use std::equal to compare contents
  f1.seekg(0, std::ifstream::beg);
  f2.seekg(0, std::ifstream::beg);
  return std::equal(std::istreambuf_iterator<char>(f1.rdbuf()),
                    std::istreambuf_iterator<char>(),
                    std::istreambuf_iterator<char>(f2.rdbuf()));
}

int main(int argc, char *argv[]){

	for(int i = 0 ; i < argc ; i ++){
		printf("argv[%d] = %s\n",i,argv[i]);
	}

	if( argc < 3){
		fprintf(stderr,"Less than 3 args.\n");
		return EXIT_FAILURE;
	}else{
		fprintf(stderr,"There are %d args.\n",argc);
	}

	std::string file1(argv[1]);
	std::string file2(argv[2]);

	std::cout << file1 << std::endl;
	std::cout << file2 << std::endl;

	if(compareFiles(file1,file2) == true){
		fprintf(stderr,"The files are identical (returning %d).\n",EXIT_SUCCESS);
		return EXIT_SUCCESS;
	}else{
		fprintf(stderr,"Files %s and %s are different (returning %d).\n",
			file1.c_str(),file2.c_str(),EXIT_FAILURE);
		return EXIT_FAILURE;
	}

}
