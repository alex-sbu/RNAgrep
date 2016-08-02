// Debugging 

// To add debugging statements, add a new flag to 
// _debug_flags, set the flag:
//  DEBUGGING_FLAGS |= flag;
// and call debug as needed:
// debug(flag, "Mr. Hooper is: %x", 57005);

#ifndef DEB
#define DEB

#include <stdarg.h>

enum _debug_flags {
    DPARALLEL = 1 << 0,
    DHELIX = 1 << 1,
    DHELIXa = 1 << 2
};

#define debug(flag, ...) (debug_print(flag, __VA_ARGS__))

int DEBUGGING_FLAGS = 0;

void debug_print(int flag, char *format, ...)
{
    if (DEBUGGING_FLAGS & flag) {
        va_list args;
        va_start(args, format);
        vfprintf(OPTION.outfile, format, args);
        va_end(args);
    }
}

#define ifdebug(x) if (DEBUGGING_FLAGS & x)

//#define debug(flag, format, ...) 
#endif
