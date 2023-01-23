/*
 * This file is part of
 * Surface Maps via Adaptive Triangulations
 * (https://github.com/patr-schm/surface-maps-via-adaptive-triangulations)
 * and is released under the MIT license.
 *
 * Authors: Janis Born
 */
#pragma once

namespace SurfaceMaps
{

//! Prints the current call stack.
//! Automatically tries to demangle symbols of C++ function names.
void print_stack_trace();

//! Call this at the beginning of your program to enable automatic stack traces on a crash.
//! A stack trace will be printed on SIGSEGV, SIGILL, SIGABRT, and SIGFPE signals.
void register_segfault_handler();

}
