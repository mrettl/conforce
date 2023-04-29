#ifdef _WIN32
#    define CF_API __declspec(dllexport)
#else
#    define CF_API
#endif

#include <math.h>  // contains power
#include <stddef.h>  // contains size_t
