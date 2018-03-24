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

/**
 * @brief Custom assert function.
 * 
 * CMake adds -DNDEBUG to the CMAKE_C_FLAGS_{RELEASE, MINSIZEREL} 
 * by default. So, you can use #ifndef NDEBUG.
 */

#define KRED  "\x1B[31m"
#define KNRM  "\x1B[0m"
#ifndef NDEBUG
    #define ASSERT(left,operator,right) { \
        if(!((left) operator (right))){ \
            printf("%sASSERT FAILED: " #left "(%lf) " #operator " " #right "(%lf) @ %s:%d.%s\n",\
                    KRED,(double)(left),(double)(right),__FILE__,__LINE__,KNRM); abort(); }\
        }
#else
    #define ASSERT(left,operator,right) 
#endif