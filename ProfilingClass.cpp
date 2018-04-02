#include "ProfilingClass.h"

#include <iostream>
#include <cstring>

#include <iomanip>
#include <sstream>

#include "mpi.h"



// Set program starting time:
void ProfilingClass::set_program_starting_time(){
    this->program_starting_time = std::time(nullptr);
}

/**
 * @brief Add a timing input.
 * 
 * Arguments:
 *      1) name of the new input
 *      2) if true, do nothing if input already exists. If false, abort if
 *          input already exists.
 */
void ProfilingClass::addTimingInputToDictionnary(
    std::string newTimingInput,
    bool try_and_do_nothing_if_exist // = true by default
){

    // Determine if the field already exists.
    if ( this->time_taken_for.find(newTimingInput) == this->time_taken_for.end() ) {
        // The timing input doesn't exist. Adding it.
        this->time_taken_for.insert(std::pair<std::string,double>(newTimingInput,0.0));
        #ifndef NDEBUG
            printf("ProfilingClass::addTimingInputToDictionnary:: timing input %s successfully added.\n",
                newTimingInput.c_str());
        #endif
    }else if(try_and_do_nothing_if_exist == true){
        // Do nothing
    } else {
        // The timing input already exists. Aborting.
        fprintf(stderr,"ProfilingClass::addTimingInputToDictionnary::ERROR\n");
        fprintf(stderr,"\t>>> At file %s:%d\n",__FILE__,__LINE__);
        fprintf(stderr,"\t>>> The timing input '%s' already exists.\n",newTimingInput.c_str());
        abort();
    }
}

// Get the time of an input:
double ProfilingClass::getTimingInput(std::string timingInput){
    if ( this->time_taken_for.find(timingInput) == this->time_taken_for.end() ) {
        // The timing input doesn't exist. Aborting.
        fprintf(stderr,"ProfilingClass::%s::ERROR\n",__FUNCTION__);
        fprintf(stderr,"\t>>> At file %s:%d\n",__FILE__,__LINE__);
        fprintf(stderr,"\t>>> The timing input '%s' doesn't exist.\n",timingInput.c_str());
        
        std::map<std::string,double>::iterator it;
        std::cout << "\t>>> Timing inputs are:\n";
        for (it=this->time_taken_for.begin(); it!=this->time_taken_for.end(); ++it)
            std::cout << "\t\t>>> " << it->first << " => " << it->second << "\n";
        abort();
    } else {
        // The timing input exists.
        return this->time_taken_for[timingInput];
    }
}

// Increment a timing input:
void ProfilingClass::incrementTimingInput(std::string timingInput, double incr){

    if(incr < 0){
        fprintf(stderr,"ProfilingClass::incrementTimingInput::ERROR\n");
        fprintf(stderr,">>> incr should be positive ut has %f.\n",incr);
        fprintf(stderr,">>> Aborting. File %s:%d.\n",__FILE__,__LINE__);
        abort();
    }

    if ( this->time_taken_for.find(timingInput) == this->time_taken_for.end() ) {
        // The timing input doesn't exist. Aborting.
        fprintf(stderr,"ProfilingClass::addTimingInputToDictionnary::ERROR\n");
        fprintf(stderr,"\t>>> At file %s:%d\n",__FILE__,__LINE__);
        fprintf(stderr,"\t>>> The timing input '%s' doesn't exist.\n",timingInput.c_str());
        
        std::map<std::string,double>::iterator it;
        std::cout << "\t>>> Timing inputs are:\n";
        for (it=this->time_taken_for.begin(); it!=this->time_taken_for.end(); ++it)
            std::cout << "\t\t>>> " << it->first << " => " << it->second << "\n";
        abort();
    } else {
        // The timing input exists.
        this->time_taken_for[timingInput] += incr;
        #ifndef NDEBUG
            printf("Time for %s is now %f seconds.\n",
                        timingInput.c_str(),this->time_taken_for[timingInput]);
        #endif
    }
}

/**
 * @brief Writes everything inside the outputfile.
 */
void ProfilingClass::writeToOutputFile(void){

    if(this->outputFileName == std::string()){
        // The output file's name has not been set. Aborting.
        DISPLAY_ERROR_ABORT(
            "The output file's name has not been set. Aborting.\n"
        );
    }
    /* Write the timings and memory usage inside the output file */
    this->outputFile.open(this->outputFileName.c_str(),std::fstream::app);
    if(this->outputFile.is_open()){
        // Write inside file.
        #ifndef NDEBUG
            std::cout << "ProfilingClass::~ProfilingClass::WRITING" << std::endl;
        #endif
        // Timing inputs:
        std::map<std::string,double>::iterator it;
        for (it=this->time_taken_for.begin(); it!=this->time_taken_for.end();++it ){
            this->outputFile << ">>> " << it->first;
            this->outputFile << " => " ;
            this->outputFile << std::fixed << std::setprecision(10) << it->second;
            this->outputFile << " seconds.\n";
        }
        
        if(this->mem_usage_peak_rss_mega_bytes == 0){
            DISPLAY_WARNING(
                "You never called a RSS prober. Do not trust the output is %s about RSS usage.",
                this->outputFileName.c_str()
            );
        }
        this->gatherMemoryUsageMPI();

        this->outputFile << "Peak RSS for this MPI is " << this->mem_usage_peak_rss_mega_bytes;
        this->outputFile << " MBytes" << std::endl;
        this->outputFile << "Peak RSS for all MPI  is " << this->total_mem_usage_peak_rss_mega_bytes_all_mpi;
        this->outputFile << " MBytes" << std::endl;
        this->outputFile.close();
    }else{
        // File could not be opened.
        DISPLAY_ERROR_ABORT(
            "Cannot open the file %s.",this->outputFileName.c_str()
        );
    }
}

// Set the output file's name:
void ProfilingClass::setOutputFileName(std::string str){
    if(this->outputFileName == std::string()){
        // The output file's name was not set. Set it.
        this->outputFileName = str;
        // Open it once to delete all:
        this->outputFile.open(str.c_str());
        if(this->outputFile.is_open()){
            // Program starting time:
            this->outputFile << std::asctime(std::localtime(&this->program_starting_time.get()));
        }
        this->outputFile.close();
    }else{
        // The output file's nam was not set. Print a warning.
        printf("ProfilingClass::setOutputFileName::WARNING\n");
        printf("Output file's name was already set to %s.\n",this->outputFileName.c_str());
        printf("I do nothing and keep this name.\n");
    }
}

void ProfilingClass::storePeakMemoryUsage(size_t MEM_IN_MEGA_BYTES){
    this->mem_usage_peak_rss_mega_bytes += MEM_IN_MEGA_BYTES;
}

void ProfilingClass::gatherMemoryUsageMPI(void){
    #ifdef MPI_COMM_WORLD
        MPI_Allreduce(
            &this->mem_usage_peak_rss_mega_bytes,
            &this->total_mem_usage_peak_rss_mega_bytes_all_mpi,
            1,
            my_MPI_SIZE_T,
            MPI_SUM,
            MPI_COMM_WORLD);
    #else
        DISPLAY_ERROR_ABORT(
            "MPI_COMM_WORLD is not defined."
        );
    #endif
}