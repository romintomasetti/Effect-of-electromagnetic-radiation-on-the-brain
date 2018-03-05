#include "ProfilingClass.h"

#include <iostream>

// Adding memory usage:
void ProfilingClass::addMemoryUsage(std::string type,double mem){
    /* type is the type of memory the user adds: in bytes, kiloBytes, megaBytes, etc.*/
    if(std::strcmp(type.c_str(),"BYTES") == 0){
        this->usedMemoryInMegaBytes += mem*1E-6;
    }else if(std::strcmp(type.c_str(),"KILOBYTES") == 0){
        this->usedMemoryInMegaBytes += mem*1E-3;
    }else if(std::strcmp(type.c_str(),"MEGABYTES") == 0){
        this->usedMemoryInMegaBytes += mem;
    }else if(std::strcmp(type.c_str(),"GIGABYTES") == 0){
        this->usedMemoryInMegaBytes += mem*1E3;
    }else{
        fprintf(stderr,"ProfilingClass::addMemoryUsage::ERROR\n");
        fprintf(stderr,"\n>>> std::string type has %s, which is not known.\nAborting.\n",type.c_str());
        fprintf(stderr,"\n>>> In file %s:%d\n",__FILE__,__LINE__);
        abort();
    }
    printf(">>> Memory in use is %f MBytes.\n",this->usedMemoryInMegaBytes);
}

// Add a timing input:
void ProfilingClass::addTimingInputToDictionnary(std::string newTimingInput){
    // Determine if the field already exists.
    if ( this->time_taken_for.find(newTimingInput) == this->time_taken_for.end() ) {
        // The timing input doesn't exist. Adding it.
        this->time_taken_for.insert(std::pair<std::string,double>(newTimingInput,0.0));
        printf("ProfilingClass::addTimingInputToDictionnary:: timing input %s successfully added.\n",
            newTimingInput.c_str());
    } else {
        // The timing input already exists. Aborting.
        fprintf(stderr,"ProfilingClass::addTimingInputToDictionnary::ERROR\n");
        fprintf(stderr,"\t>>> At file %s:%d\n",__FILE__,__LINE__);
        fprintf(stderr,"\t>>> The timing input '%s' already exists.\n",newTimingInput.c_str());
        abort();
    }
}

// Increment a timing input:
void ProfilingClass::incrementTimingInput(std::string timingInput, double incr){
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
        printf("Time for %s is now %f seconds.\n",timingInput.c_str(),this->time_taken_for[timingInput]);
    }
}