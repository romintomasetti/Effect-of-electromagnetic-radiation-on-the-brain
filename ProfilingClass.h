#ifndef PROFILINGCLASS_H
#define PROFILINGCLASS_H

#include <string>
#include <cstring>
#include <stdio.h>
#include <map>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <ctime>
#include <iostream>
#include "SetOnceVariable_Template.h"

class ProfilingClass{
    private:
        // Used memory:
        double usedMemoryInMegaBytes = 0.0;
        double peakUsedMemoryInMegaBytes = 0.0;
        
        // Dictionnary for arbitrary time input:
        std::map<std::string,double> time_taken_for;

        // File for output:
        std::ofstream outputFile;
        std::string   outputFileName = std::string();

        // Program starting time:
        SetOnceVariable_Template<std::time_t> program_starting_time;

    public:
        // Constructor:
        ProfilingClass(void){this->set_program_starting_time();}
        // Destructor:
        ~ProfilingClass(void);

        // Set program starting time:
        void set_program_starting_time();

        // Adding memory usage:
        void addMemoryUsage(std::string,double);

        // Add a timing input:
        void addTimingInputToDictionnary(std::string);

        // Increment a timing input:
        void incrementTimingInput(std::string, double);

        // Set the output file's name:
        void setOutputFileName(std::string);

        // Remove some memory usage:
        void removeMemoryUsage(std::string type, double mem,
                        std::string senderMessage = std::string());

};

#endif
