/*===========================================================================*\
 *                                                                           *
 *                              OpenFlipper                                  *
 *           Copyright (c) 2001-2015, RWTH-Aachen University                 *
 *           Department of Computer Graphics and Multimedia                  *
 *                          All rights reserved.                             *
 *                            www.openflipper.org                            *
 *                                                                           *
 *---------------------------------------------------------------------------*
 * This file is part of OpenFlipper.                                         *
 *---------------------------------------------------------------------------*
 *                                                                           *
 * Redistribution and use in source and binary forms, with or without        *
 * modification, are permitted provided that the following conditions        *
 * are met:                                                                  *
 *                                                                           *
 * 1. Redistributions of source code must retain the above copyright notice, *
 *    this list of conditions and the following disclaimer.                  *
 *                                                                           *
 * 2. Redistributions in binary form must reproduce the above copyright      *
 *    notice, this list of conditions and the following disclaimer in the    *
 *    documentation and/or other materials provided with the distribution.   *
 *                                                                           *
 * 3. Neither the name of the copyright holder nor the names of its          *
 *    contributors may be used to endorse or promote products derived from   *
 *    this software without specific prior written permission.               *
 *                                                                           *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS       *
 * "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED *
 * TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A           *
 * PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER *
 * OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,  *
 * EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,       *
 * PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR        *
 * PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF    *
 * LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING      *
 * NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS        *
 * SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.              *
 *                                                                           *
\*===========================================================================*/



//=============================================================================
//
//  CLASS QuadricT
//
//=============================================================================

#ifndef ACG_QUADRIC_HH
#define ACG_QUADRIC_HH


//== INCLUDES =================================================================


#include "../../Math/VectorT.hh"


//== NAMESPACE ================================================================


namespace ACG {
namespace Geometry {


//== CLASS DEFINITION =========================================================


/** /class QuadricT QuadricT.hh <ACG/Geometry/Types/QuadricT.hh>

    Stores a quadric as a 4x4 symmetrix matrix. Used by the
    error quadric based mesh decimation algorithms.
**/

template <class Scalar>
class ACGDLLEXPORT QuadricT
{
public:


  /// construct with upper triangle of symmetrix 4x4 matrix
  QuadricT(Scalar _a, Scalar _b, Scalar _c, Scalar _d,
	              Scalar _e, Scalar _f, Scalar _g,
	                         Scalar _h, Scalar _i,
	                                    Scalar _j)
    : a(_a), b(_b), c(_c), d(_d),
             e(_e), f(_f), g(_g),
                    h(_h), i(_i),
                           j(_j)
  {}


  /// constructor from given plane equation: ax+by+cz+d=0
  QuadricT( Scalar _a=0.0, Scalar _b=0.0, Scalar _c=0.0, Scalar _d=0.0 )
    :  a(_a*_a), b(_a*_b),  c(_a*_c),  d(_a*_d),
                 e(_b*_b),  f(_b*_c),  g(_b*_d),
                            h(_c*_c),  i(_c*_d),
                                       j(_d*_d)
  {}


  /// set all entries to zero
  void clear()  { a = b = c = d = e = f = g = h = i = j = 0.0; }


  /// add quadrics
  QuadricT<Scalar>& operator+=( const QuadricT<Scalar>& _q )
  {
    a += _q.a;  b += _q.b;  c += _q.c;  d += _q.d;
                e += _q.e;  f += _q.f;  g += _q.g;
                            h += _q.h;  i += _q.i;
			                j += _q.j;
    return *this;
  }


  /// multiply by scalar
  QuadricT<Scalar>& operator*=( Scalar _s)
  {
    a *= _s;  b *= _s;  c *= _s;  d *= _s;
              e *= _s;  f *= _s;  g *= _s;
                        h *= _s;  i *= _s;
                                  j *= _s;
    return *this;
  }


  /// multiply 4D vector from right: Q*v
  VectorT<Scalar,4> operator*(const VectorT<Scalar,4>& _v) const
  {
    return VectorT<Scalar,4>(_v[0]*a + _v[1]*b + _v[2]*c + _v[3]*d,
			     _v[0]*b + _v[1]*e + _v[2]*f + _v[3]*g,
			     _v[0]*c + _v[1]*f + _v[2]*h + _v[3]*i,
			     _v[0]*d + _v[1]*g + _v[2]*i + _v[3]*j);
  }


  /// evaluate quadric Q at vector v: v*Q*v
  Scalar operator()(const VectorT<Scalar,3> _v) const
  {
    Scalar x(_v[0]), y(_v[1]), z(_v[2]);
    return a*x*x + 2.0*b*x*y + 2.0*c*x*z + 2.0*d*x
                 +     e*y*y + 2.0*f*y*z + 2.0*g*y
	                     +     h*z*z + 2.0*i*z
                                         +     j;
  }


  /// evaluate quadric Q at vector v: v*Q*v
  Scalar operator()(const VectorT<Scalar,4> _v) const
  {
    Scalar x(_v[0]), y(_v[1]), z(_v[2]), w(_v[3]);
    return a*x*x + 2.0*b*x*y + 2.0*c*x*z + 2.0*d*x*w
                 +     e*y*y + 2.0*f*y*z + 2.0*g*y*w
	                     +     h*z*z + 2.0*i*z*w
                                         +     j*w*w;
  }


private:

  Scalar a, b, c, d,
            e, f, g,
               h, i,
                  j;
};


/// Quadric using floats
typedef QuadricT<float> Quadricf;

/// Quadric using double
typedef QuadricT<double> Quadricd;


//=============================================================================
} // namespace Geometry
} // namespace ACG
//=============================================================================
#endif // ACG_QUADRIC_HH defined
//=============================================================================
