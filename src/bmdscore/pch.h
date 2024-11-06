// pch.h: This is a precompiled header file.
// Files listed below are compiled only once, improving build performance for future builds.
// This also affects IntelliSense performance, including code completion and many code browsing
// features. However, files listed here are ALL re-compiled if any one of them is updated between
// builds. Do not add files here that you will be updating frequently as this negates the
// performance advantage.

#ifndef PCH_H
#define PCH_H

// add headers that you want to pre-compile here
#pragma once
// #pragma message ("Entered pch.h")

#include <stdio.h>

#include <limits>
#define NEAR_ZERO 1e-8
// #pragma message ("defining _USE_MATH_DEFINES in pch.h")
#define _USE_MATH_DEFINES

#ifndef R_COMPILATION
#  if defined WIN32 || defined _WINDOWS
#    include <math.h>
#    include <tchar.h>
#    include <windows.h>

#    include "targetver.h"
#  else  // i.e., not WIN32 or _WINDOWS
#    include <cfloat>
#  endif  // defined WIN32 || defined _WINDOWS
#endif    // ifndef R_COMPILATION

// Custom DEBUG and LOG macros
// - Do not define DEBUGLOG when building R because it will break the build
#ifdef DEBUGLOG
#  include <fstream>
#  define DEBUG_OPEN_LOG(fname, fvar)                           \
    ofstream(fvar);                                             \
    (fvar).open(fname, fstream::app);                           \
    (fvar) << __FUNCTION__ << " at line: " << __LINE__ << endl; \
    /*(fvar) << std::scientific << setprecision(17);*/
#  define DEBUG_LOG(f, s) \
    (f) << s << endl;     \
    flush((f));
#  define DEBUG_CLOSE_LOG(fvar) (fvar).close();
#else
#  define DEBUG_OPEN_LOG(fname, fvar)
#  define DEBUG_LOG(f, s)
#  define DEBUG_CLOSE_LOG(fvar)
#endif  // DEBUGLOG
#include "framework.h"

#endif  // PCH_H
