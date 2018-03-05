#ifndef PROFILINGCLASS_H
#define PROFILINGCLASS_H

#include <string>
#include <cstring>
#include <stdio.h>
#include <map>
#include <stdlib.h>

class ProfilingClass{
    private:
        // Used memory:
        double usedMemoryInMegaBytes = 0.0;
        
        // Dictionnary for arbitrary time input:
        std::map<std::string,double> time_taken_for;

    public:
        // Constructor:
        ProfilingClass(void){}
        // Destructor:
        ~ProfilingClass(void){}

        // Adding memory usage:
        void addMemoryUsage(std::string,double);

        // Add a timing input:
        void addTimingInputToDictionnary(std::string);

        // Increment a timing input:
        void incrementTimingInput(std::string, double);

};

#endif
