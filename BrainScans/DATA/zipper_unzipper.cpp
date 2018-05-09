#include <cstdlib>
#include <stdio.h>
#include <string>
#include <fstream>
#include <iostream>
#include <cstring>

inline bool fileExists(const std::string& name) {
    std::ifstream f(name.c_str());
    return f.good();
}

inline std::string getFilenameWithoutExtension(std::string const& name){
	size_t pos = name.find(".");
	std::string res = name.substr(0,pos);
	return res;
}

void compressFile(const std::string& filename){
	std::string command  = "tar -zcvf ";
	if( fileExists(filename) ){
		printf("Compresing %s...\n",filename.c_str());
		std::string filenameNoExt = getFilenameWithoutExtension(filename);
        filenameNoExt.append(".tar");
		command.append(filenameNoExt);
        command.append(" ");
        command.append(filename);
        printf(">>> Command is %s\n",command.c_str());
		system(command.c_str());
		printf("File %s compressed successfully !\n",filename.c_str());
	}else{
		printf(">>> File %s was not found.\n",filename.c_str());
		exit(-1);
	}
}

void uncompressFile(const std::string& filename){
	std::string command  = "tar -zxvf ";
	if( fileExists(filename) ){
		printf("Uncompressing %s...\n",filename.c_str());
		command.append(filename);
		system(command.c_str());
		printf("File %s uncompressed successfully !\n",filename.c_str());
	}else{
		printf(">>> File %s was not found.\n",filename.c_str());
		exit(-1);
	}
}

int main(int argc, char *argv[]){

	/* Number of arguments: executable name, file name, compress or uncompress. */
	/* To uncompress, give -u. To compress, give -c. */

	if( argc != 3 ){
		printf(">>> Wrong command. Example: zipper_unzipper -u file.tar.gz\n");
		exit(-1);
	}
	
	std::string filename(argv[2]);
	if( ! fileExists(filename) ){
		printf(">>> File %s was not found.\n",filename.c_str());
		exit(-1);
	}

	if( strcmp(argv[1], "-c") == 0 ){
		compressFile(filename);
		return EXIT_SUCCESS;
	}else if( strcmp(argv[1], "-u") == 0){
		uncompressFile(filename);
		return EXIT_SUCCESS;
	}else{
		printf(">>> Unknown command %s.\n",argv[1]);
		return EXIT_FAILURE;
	}

	return EXIT_SUCCESS;
}
