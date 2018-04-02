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

#include "MemoryUsage.h"

#include "header_with_all_defines.hpp"

class ProfilingClass{
    private:

        size_t mem_usage_peak_rss_mega_bytes = 0;
        size_t total_mem_usage_peak_rss_mega_bytes_all_mpi = 0;
        
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

        // Add a timing input:
        void addTimingInputToDictionnary(
            std::string,
            bool try_and_do_nothing_if_exist = true
        );

        // Increment a timing input:
        void incrementTimingInput(std::string, double);

        // Get the time of an input:
        double getTimingInput(std::string);

        // Set the output file's name:
        void setOutputFileName(std::string);

        void gatherMemoryUsageMPI(void);

        void storePeakMemoryUsage(size_t MEM_IN_MEGA_BYTES);

        void probeMaxRSS(void){
            this->mem_usage_peak_rss_mega_bytes = getPeakRSS()/1E6;
        }

};

#endif
