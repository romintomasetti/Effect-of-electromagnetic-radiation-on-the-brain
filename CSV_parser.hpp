#ifndef CSV_PARSER_HPP
#define CSV_PARSER_HPP

#include <string>
#include <vector>
#include <sstream>
#include <iostream>
#include <fstream>
/**
 * Reads csv file into table, exported as a vector of vector of doubles.
 * @param inputFileName input file name (full path).
 * @return data as vector of vector of doubles.
 */ 
std::vector<std::vector<std::string> > parse2DCsvFile(std::string inputFileName);

void RemoveAnyBlankSpaceInStr(std::string &str);

#endif
