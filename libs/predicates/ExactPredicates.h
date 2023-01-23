#pragma once

#ifdef __cplusplus
extern "C" {
#endif

#include <cassert>

void exactinit();
double orient2d(const double *pa, const double *pb, const double *pc);
double orient3d(const double *pa, const double *pb, const double *pc, const double *pd);

#ifdef __cplusplus
} // extern "C"
#endif
