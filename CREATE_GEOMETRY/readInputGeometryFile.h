#ifndef READINPUTGEOMETRYFILE_H
#define READINPUTGEOMETRYFILE_H

#if defined(WIN32)
#ifdef discr_integr_EXPORTS
#define CREATE_GEOMETRY_API __declspec(dllexport)
#else
#define CREATE_GEOMETRY_API __declspec(dllimport)
#endif
#else
#define CREATE_GEOMETRY_API
#endif

#include <stdio.h>
#include <stdlib.h>

#include <string>
#include <iostream>
#include <vector>
#include <fstream>

#include <cctype>

#include <sstream>

using namespace std;

CREATE_GEOMETRY_API unsigned int* read_input_geometry_file(std::string filename, size_t *size_read);

#endif
