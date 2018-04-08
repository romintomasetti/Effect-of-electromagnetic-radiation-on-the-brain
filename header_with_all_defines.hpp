#ifndef HEADER_WITH_ALL_DEFINES_HPP
#define HEADER_WITH_ALL_DEFINES_HPP
/**
 * Define some colors:
 */
#define ANSI_COLOR_RED     "\x1b[31m"
#define ANSI_COLOR_GREEN   "\x1b[32m"
#define ANSI_COLOR_YELLOW  "\x1b[33m"
#define ANSI_COLOR_BLUE    "\x1b[34m"
#define ANSI_COLOR_MAGENTA "\x1b[35m"
#define ANSI_COLOR_CYAN    "\x1b[36m"
#define ANSI_COLOR_RESET   "\x1b[0m"


#include <cmath>
/**
 * PHYSICAL VARIABLES:
 */
#define VACUUM_PERMITTIVITY 8.8541878176E-12
#define VACUUM_PERMEABILITY 4*M_PI*1E-7

/**
 * @brief Custom assert function.
 * 
 * CMake adds -DNDEBUG to the CMAKE_C_FLAGS_{RELEASE, MINSIZEREL} 
 * by default. So, you can use #ifndef NDEBUG.
 */

#define KRED  "\x1B[31m"
#define KNRM  "\x1B[0m"

/**
 * CUSTOM ASSERT WHEN CMAKE IS IN DEBUG MODE
 */
#ifndef NDEBUG
    #define ASSERT(left,operator,right) { \
        if(!((left) operator (right))){ \
            printf("%sASSERT FAILED: " #left "(%lf) " #operator " " #right "(%lf) @ %s:%d.%s\n",\
                    KRED,(double)(left),(double)(right),__FILE__,__LINE__,KNRM); abort(); }\
        }
#else
    #define ASSERT(left,operator,right) 
#endif


/**
 * CUSTOM ABORT FUNCTION
 */
#ifdef MPI_COMM_WORLD
        #define ABORT_MPI(ARG){\
			int done_already;\
			MPI_Initialized(&done_already);\
			if (!done_already)\
				abort();\
			else\
				MPI_Abort(MPI_COMM_WORLD,ARG);\
		}
#else
        #define ABORT_MPI(ARG) abort();
#endif

/**
 * Custom FPRINTF.
 */
#ifdef MPI_COMM_WORLD
#define FPRINTF(...){\
        int ID_MPI_Process = INT_MIN;\
		int done_already;\
		MPI_Initialized(&done_already);\
		if(!done_already){\
			ID_MPI_Process = INT_MIN;\
		}else{\
			MPI_Comm_rank( MPI_COMM_WORLD, &ID_MPI_Process);\
		}\
        fprintf(stderr,"[MPI %d] - ",ID_MPI_Process);\
        fprintf(__VA_ARGS__);\
}
#else
#define FPRINTF(...)\
        fprintf(__VA_ARGS__);
#endif

/**
 * DISPLAY A MESSAGE AND ABORT
 */
#define DISPLAY_ERROR_ABORT(...){              \
        FPRINTF(stderr,"%sIn %s :: ERROR :: %s",      \
                    ANSI_COLOR_RED,                   \
                    __FUNCTION__,                     \
                    ANSI_COLOR_GREEN);                \
        fprintf(stderr,__VA_ARGS__);           \
        fprintf(stderr,"%s Aborting.%s\nIn %s:%d\n",   \
                    ANSI_COLOR_YELLOW,                \
                    ANSI_COLOR_RESET,                 \
                    __FILE__,                         \
                    __LINE__);                        \
        ABORT_MPI(-1);                                \
        }

/**
 * DISPLAY A WARNING. DO NOT ABORT.
 */
#define DISPLAY_WARNING(...){                       \
        FPRINTF(stderr,"%sIn %s :: WARNING :: %s",  \
                    ANSI_COLOR_YELLOW,              \
                    __FUNCTION__,                   \
                    ANSI_COLOR_GREEN);              \
        fprintf(stderr,__VA_ARGS__);                \
        fprintf(stderr,"\n%sIn %s:%d\n",            \
                    ANSI_COLOR_RESET,               \
                    __FILE__,                       \
                    __LINE__);                      \
        }

/**
 * CUSTOM SIZE_T DATATYPE FOR MPI COMMUNICATION.
 */
#include <stdint.h>
#include <limits.h>

#if SIZE_MAX == UCHAR_MAX
   #define my_MPI_SIZE_T MPI_UNSIGNED_CHAR
#elif SIZE_MAX == USHRT_MAX
   #define my_MPI_SIZE_T MPI_UNSIGNED_SHORT
#elif SIZE_MAX == UINT_MAX
   #define my_MPI_SIZE_T MPI_UNSIGNED
#elif SIZE_MAX == ULONG_MAX
   #define my_MPI_SIZE_T MPI_UNSIGNED_LONG
#elif SIZE_MAX == ULLONG_MAX
   #define my_MPI_SIZE_T MPI_UNSIGNED_LONG_LONG
#else
   #error "what is happening here?"
#endif

#endif
