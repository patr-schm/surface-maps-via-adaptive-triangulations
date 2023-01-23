/*
 * This file is part of
 * Surface Maps via Adaptive Triangulations
 * (https://github.com/patr-schm/surface-maps-via-adaptive-triangulations)
 * and is released under the MIT license.
 *
 * Authors: Patrick Schmidt, Janis Born
 */
#pragma once

#include "AnsiColorCodes.hh"
#include <iostream>
#include <sstream>
#include <math.h>

#define ISM_INFO_wait(str) \
        {std::cout << SurfaceMaps::ANSI_FG_GREEN << str << SurfaceMaps::ANSI_RESET; \
        std::cout.flush();}

#define ISM_INFO(str) \
        {std::cout << SurfaceMaps::ANSI_FG_GREEN << str << SurfaceMaps::ANSI_RESET << std::endl; \
        std::cout.flush();}

#define ISM_HIGHLIGHT(str) \
        {std::cout << SurfaceMaps::ANSI_FG_GREEN << SurfaceMaps::ANSI_BG_MAGENTA << str << SurfaceMaps::ANSI_RESET << std::endl; \
        std::cout.flush();}

#define ISM_DEBUG_OUT(str) \
        {std::cout << SurfaceMaps::ANSI_FG_MAGENTA \
                   << "[DEBUG] " \
                   << str \
                   << SurfaceMaps::ANSI_RESET << std::endl; \
        std::cout.flush();}

#define ISM_DEBUG_VAR(var) \
        {ISM_DEBUG_OUT(#var << " = " << var)}

#define ISM_WARNING(str) \
        {std::cout << SurfaceMaps::ANSI_FG_YELLOW \
                  << "[WARNING] " \
                  << str \
                  << SurfaceMaps::ANSI_RESET \
                  << " (in function " << __FUNCTION__ << ":" << __LINE__ \
                  << " in file " << __FILE__ << ")" \
                  << std::endl; \
        std::cout.flush();}

#define ISM_ERROR(str) \
        std::cout << SurfaceMaps::ANSI_FG_RED \
                  << "[ERROR] " \
                  << str \
                  << SurfaceMaps::ANSI_RESET \
                  << " (in function " << __FUNCTION__ << ":" << __LINE__ \
                  << " in file " << __FILE__ << ")" \
                  << std::endl

#define ISM_ERROR_throw(st) \
        {ISM_ERROR(st); \
        std::stringstream str_strm; \
        str_strm << "[ERROR] " << st; \
        throw std::runtime_error(str_strm.str());}

// Assertions

#define ISM_ASSERT(exp) \
        {if(!(exp)) ISM_ERROR_throw("Assertion failed: " << (#exp));}

#define ISM_ASSERT_MSG(exp, msg) \
        {if(!(exp)) ISM_ERROR_throw("Assertion failed: " << (#exp) << ". Message: " << msg);}

#define ISM_ASSERT_EQ(a, b) \
        {if((a) != (b)) ISM_ERROR_throw("Assertion failed: " << (a) << " == " << (b));}

#define ISM_ASSERT_NEQ(a, b) \
        {if((a) == (b)) ISM_ERROR_throw("Assertion failed: " << (a) << " != " << (b));}

#define ISM_ASSERT_G(a, b) \
        {if((a) <= (b)) ISM_ERROR_throw("Assertion failed: " << (a) << " > " << (b));}

#define ISM_ASSERT_GEQ(a, b) \
        {if((a) < (b)) ISM_ERROR_throw("Assertion failed: " << (a) << " >= " << (b));}

#define ISM_ASSERT_L(a, b) \
        {if((a) >= (b)) ISM_ERROR_throw("Assertion failed: " << (a) << " < " << (b));}

#define ISM_ASSERT_LEQ(a, b) \
        {if((a) > (b)) ISM_ERROR_throw("Assertion failed: " << (a) << " <= " << (b));}

#define ISM_ASSERT_FINITE(a) \
        {ISM_ASSERT(isfinite(a));}

#define ISM_ASSERT_FINITE_MAT(A) \
{ \
    auto A_eval = A; \
    for (int i = 0; i < (A_eval).rows(); ++i) \
    { \
        for (int j = 0; j < (A_eval).cols(); ++j) \
        { \
            if (!isfinite((A_eval)(i, j))) \
                ISM_ERROR_throw("Assertion failed: Not finite " << (A_eval)); \
        } \
    } \
}

#define ISM_ASSERT_EPS(a, b, eps) \
        {if(std::abs((a) - (b)) >= eps) ISM_ERROR_throw("Assertion failed: |" << (a) << " - " << (b) << "| < " << eps);}

#define ISM_ASSERT_REL(a, b, eps) \
        {ISM_ASSERT_NEQ(b, 0.0);\
         const double rel_error = std::abs(a - b) / std::abs(b);\
         if(rel_error >= eps) ISM_ERROR_throw("Assertion failed: relative error between " << a << " and " << b << " is " << rel_error << " and exceeds " << eps);}

#define ISM_ASSERT_ELEMENT_EPS(a, b, eps) \
        {ISM_ASSERT_EQ(a.size(), b.size()); \
         for (int i = 0; i < a.size(); ++i) \
         { \
             if(std::abs((a[i]) - (b[i])) >= eps) \
                 ISM_ERROR_throw("Assertion failed for element " << i << ": |" << (a[i]) << " - " << (b[i]) << "| < " << eps); \
         }}

#define ISM_ASSERT_NORM_EPS(a, b, eps) \
        {const auto error = ((a) - (b)).norm(); if (error >= eps) ISM_ERROR_throw("Assertion failed: ||[" << (a) << "] - [" << (b) << "]|| = " << error << " !< " << eps);}

#define ISM_ASSERT_IS_NAN(exp) \
        {if((exp) == (exp)) ISM_ERROR_throw("Assertion failed, Should be NAN: " << (#exp));}

#define ISM_ASSERT_NOT_NAN(exp) \
        {if((exp) != (exp)) ISM_ERROR_throw("Assertion failed, NAN: " << (#exp));}

#define ISM_ASSERT_NOT_NAN2(exp) \
        {if((exp[0]) != (exp[0]) || (exp[1]) != (exp[1])) ISM_ERROR_throw("Assertion failed, NAN: " << (#exp));}

#define ISM_ASSERT_NOT_NAN3(exp) \
        {if((exp[0]) != (exp[0]) || (exp[1]) != (exp[1]) || (exp[2]) != (exp[2])) ISM_ERROR_throw("Assertion failed, NAN: " << (#exp));}

#define ISM_ASSERT_NOT_NAN_SPARSE(A) \
{ \
for (uint i = 0; i < (A).outerSize(); ++i) \
{ \
    for (SparseMatrix::InnerIterator it((A), i); it; ++it) \
        ISM_ASSERT_NOT_NAN(it.value()) \
} \
}

#define ISM_ASSERT_SYMMETRIC(A, eps) \
{ \
    if (((A) - (A).transpose()).array().abs().maxCoeff() > eps) ISM_ERROR_throw("Matrix not symmetric"); \
}

// Warnings

#define ISM_EXPECT(exp) \
        {if(!(exp)) ISM_WARNING("Expectation not met: " << (#exp));}

#define ISM_EXPECT_MSG(exp, msg) \
        {if(!(exp)) ISM_WARNING("Expectation not met: " << (#exp) << ". Message: " << msg);}

#define ISM_EXPECT_EQ(a, b) \
        {if((a) != (b)) ISM_WARNING("Expectation not met: " << (a) << " == " << (b));}

#define ISM_EXPECT_NEQ(a, b) \
        {if((a) == (b)) ISM_WARNING("Expectation not met: " << (a) << " != " << (b));}

#define ISM_EXPECT_G(a, b) \
        {if((a) <= (b)) ISM_WARNING("Expectation not met: " << (a) << " > " << (b));}

#define ISM_EXPECT_GEQ(a, b) \
        {if((a) < (b)) ISM_WARNING("Expectation not met: " << (a) << " >= " << (b));}

#define ISM_EXPECT_L(a, b) \
        {if((a) >= (b)) ISM_WARNING("Expectation not met: " << (a) << " < " << (b));}

#define ISM_EXPECT_LEQ(a, b) \
        {if((a) > (b)) ISM_WARNING("Expectation not met: " << (a) << " <= " << (b));}

#define ISM_EXPECT_FINITE(a) \
        {ISM_EXPECT(isfinite(a));}

#define ISM_EXPECT_EPS(a, b, eps) \
        {if(std::abs((a) - (b)) >= eps) ISM_WARNING("Expectation not met: |" << (a) << " - " << (b) << "| < " << eps);}

#define ISM_EXPECT_REL(a, b, eps) \
        {ISM_ASSERT_NEQ(b, 0.0);\
         const double rel_error = std::abs(a - b) / std::abs(b);\
         if(rel_error >= eps) ISM_WARNING("Expectation not met: relative error between " << a << " and " << b << " is " << rel_error << " and exceeds " << eps);}

#define ISM_EXPECT_ELEMENT_EPS(a, b, eps) \
        {ISM_ASSERT_EQ(a.size(), b.size()); \
         for (int i = 0; i < a.size(); ++i) \
         { \
             if(std::abs((a[i]) - (b[i])) >= eps) \
                 ISM_WARNING("Expectation not met for element " << i << ": |" << (a[i]) << " - " << (b[i]) << "| < " << eps); \
         }}

#define ISM_EXPECT_NORM_EPS(a, b, eps) \
        {const auto error = ((a) - (b)).norm(); if (error >= eps) ISM_WARNING("Expectation not met: ||[" << (a) << "] - [" << (b) << "]|| = " << error << " !< " << eps);}

#define ISM_EXPECT_NAN(exp) \
        {if((exp) == (exp)) ISM_WARNING("Expectation not met, Should be NAN: " << (#exp));}

#define ISM_EXPECT_NOT_NAN(exp) \
        {if((exp) != (exp)) ISM_WARNING("Expectation not met, NAN: " << (#exp));}

#define ISM_EXPECT_NOT_NAN2(exp) \
        {if((exp[0]) != (exp[0]) || (exp[1]) != (exp[1])) ISM_WARNING("Expectation not met, NAN: " << (#exp));}
