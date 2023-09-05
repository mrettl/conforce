/**
 * DO NOT INCLUDE OR RUN THIS FILE. THIS IS JUST A TEMPLATE USED BY THE CODE GENERATION!
 * Compile this file with:
 *
 * "gcc -shared -O3 -fPIC -std=c99 -o {library name}.{library extension} {this file name}.c"
 *
 * The {library extension} is "dll" on Windows and "so" on Linux
 *
 */

#ifdef _WIN32
#    define CF_API __declspec(dllexport)
#else
#    define CF_API
#endif

#include <math.h>  // contains power
#include <stddef.h>  // contains size_t
