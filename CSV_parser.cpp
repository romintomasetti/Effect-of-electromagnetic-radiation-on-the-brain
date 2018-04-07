#include "CSV_parser.hpp"

#include <string>
#include <iostream>
#include <fstream>
#include <map>
#include <vector>
#include <sstream>
#include <stdio.h>
#include <stdlib.h>
#include <algorithm>
#include <cctype>

#include <boost/algorithm/string.hpp>

#include <cstring>

#include <sys/stat.h>

#include <cmath>

#include <dirent.h>
#include <fstream>

#include <time.h>

inline void RemoveAnyBlankSpaceInStr(std::string &str){
	str.erase(std::remove_if(str.begin(),
			 str.end(), [](unsigned char x){return std::isspace(x);}),
			 str.end());
}

std::vector<std::vector<std::string> > parse2DCsvFile(std::string inputFileName) {
 
    std::vector<std::vector<std::string> > data;
    std::ifstream inputFile(inputFileName);
    int l = 0;
 
    while (inputFile) {
        l++;
        std::string s;
        if (!getline(inputFile, s)) break;
        if (s[0] != '#') {
            std::istringstream ss(s);
            std::vector<std::string> record;
 
            while (ss) {
                std::string line;
                if (!getline(ss, line, ','))
                    break;
                try {
                    RemoveAnyBlankSpaceInStr(line);
                    boost::erase_all(line, "\"");
                    record.push_back(line);
                }
                catch (const std::invalid_argument &e) {
                    std::cout << "NaN found in file " << inputFileName << " line " << l
                         << std::endl;
                    e.what();
                }
            }
 
            data.push_back( record );
        }
    }
 
    if (!inputFile.eof()) {
        std::cerr << "Could not read file " << inputFileName << "\n";
        throw std::invalid_argument("File not found.");
    }
 
    return data;
}