#ifndef PROFILINGCLASS_H
#define PROFILINGCLASS_H

#include <string>
#include <cstring>
#include <stdio.h>
#include <map>
#include <stdlib.h>
#include <iostream>
#include <fstream>

class ProfilingClass{
    private:
        // Used memory:
        double usedMemoryInMegaBytes = 0.0;
        
        // Dictionnary for arbitrary time input:
        std::map<std::string,double> time_taken_for;

        // File for output:
        std::ofstream outputFile;
        std::string   outputFileName = std::string();

    public:
        // Constructor:
        ProfilingClass(void){}
        // Destructor:
        ~ProfilingClass(void);

        // Adding memory usage:
        void addMemoryUsage(std::string,double);

        // Add a timing input:
        void addTimingInputToDictionnary(std::string);

        // Increment a timing input:
        void incrementTimingInput(std::string, double);

        // Set the output file's name:
        void setOutputFileName(std::string);

};

#endif
