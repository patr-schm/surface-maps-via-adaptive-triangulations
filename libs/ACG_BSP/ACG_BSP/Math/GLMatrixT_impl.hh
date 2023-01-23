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
//  CLASS GLMatrixT - IMPLEMENTATION
//
//=============================================================================


#define ACG_GLMATRIX_C


//== INCLUDES =================================================================


#include "GLMatrixT.hh"
#include <cmath>


//== IMPLEMENTATION ========================================================== 


namespace ACG {


//-----------------------------------------------------------------------------


template<typename Scalar> 
void
GLMatrixT<Scalar>::
scale( Scalar _x, Scalar _y, Scalar _z, 
       MultiplyFrom _mult_from )
{
  GLMatrixT<Scalar> m;
  m.identity();

  m(0,0) = _x;
  m(1,1) = _y;
  m(2,2) = _z;
  
  if (_mult_from == MULT_FROM_RIGHT)  *this *= m;
  else                                leftMult(m);
}


//-----------------------------------------------------------------------------


template<typename Scalar> 
void
GLMatrixT<Scalar>::
translate( Scalar _x, Scalar _y, Scalar _z, 
	   MultiplyFrom _mult_from )
{
  GLMatrixT<Scalar> m;
  m.identity();

  m(0,3) = _x;
  m(1,3) = _y;
  m(2,3) = _z;

  if (_mult_from == MULT_FROM_RIGHT)  *this *= m;
  else                                leftMult(m);
}


//-----------------------------------------------------------------------------


template<typename Scalar> 
void
GLMatrixT<Scalar>::
rotateXYZ( Axis _axis, Scalar _angle, 
	   MultiplyFrom _mult_from )
{
  GLMatrixT<Scalar> m;
  m.identity();

  Scalar ca = cos(_angle * (M_PI/180.0));
  Scalar sa = sin(_angle * (M_PI/180.0));
    
  switch (_axis) 
  {
    case X:
      m(1,1) = m(2,2) = ca;  m(2,1) = sa;  m(1,2) = -sa;
      break;

    case Y:
      m(0,0) = m(2,2) = ca;  m(0,2) = sa;  m(2,0) = -sa;
      break;

    case Z:
      m(0,0) = m(1,1) = ca;  m(1,0) = sa;  m(0,1) = -sa;
      break;
  }


  if (_mult_from == MULT_FROM_RIGHT)  *this *= m;
  else                                leftMult(m);
}


//-----------------------------------------------------------------------------


/* Rotation matrix (taken from Mesa3.1)
   original function contributed by Erich Boleyn (erich@uruk.org) */
template <typename Scalar> 
void
GLMatrixT<Scalar>::
rotate( Scalar angle, Scalar x, Scalar y, Scalar z,
	MultiplyFrom _mult_from )
{
  GLMatrixT<Scalar> m;
  Scalar  mag, s, c;
  Scalar  xx, yy, zz, xy, yz, zx, xs, ys, zs, one_c;

  mag = sqrt( x*x + y*y + z*z );
  if(mag == 0.) {
    return;
  }
   
  s = sin( angle * ( Scalar(M_PI) / Scalar(180.) ) );
  c = cos( angle * ( Scalar(M_PI) / Scalar(180.) ) );

  x /= mag;
  y /= mag;
  z /= mag;
  
  xx = x * x;
  yy = y * y;
  zz = z * z;
  xy = x * y;
  yz = y * z;
  zx = z * x;
  xs = x * s;
  ys = y * s;
  zs = z * s;
  one_c = static_cast<Scalar>(1.0) - c;

  m(0,0) = (one_c * xx) + c;
  m(0,1) = (one_c * xy) - zs;
  m(0,2) = (one_c * zx) + ys;
  m(0,3) = static_cast<Scalar>(0.0);

  m(1,0) = (one_c * xy) + zs;
  m(1,1) = (one_c * yy) + c;
  m(1,2) = (one_c * yz) - xs;
  m(1,3) = static_cast<Scalar>(0.0);
  
  m(2,0) = (one_c * zx) - ys;
  m(2,1) = (one_c * yz) + xs;
  m(2,2) = (one_c * zz) + c;
  m(2,3) = static_cast<Scalar>(0.0);
  
  m(3,0) = static_cast<Scalar>(0.0);
  m(3,1) = static_cast<Scalar>(0.0);
  m(3,2) = static_cast<Scalar>(0.0);
  m(3,3) = static_cast<Scalar>(1.0);

  if (_mult_from == MULT_FROM_RIGHT)  *this *= m;
  else                                leftMult(m);
}


//-----------------------------------------------------------------------------


template<typename Scalar> 
void
GLMatrixT<Scalar>::
lookAt(const VectorT<Scalar,3>& eye,
       const VectorT<Scalar,3>& center,
       const VectorT<Scalar,3>& up)
{
  VectorT<Scalar,3> z(eye);
  z -= center;
  z.normalize();

  VectorT<Scalar,3> x(up % z);
  x.normalize();

  VectorT<Scalar,3> y(z % x);
  y.normalize();

  GLMatrixT<Scalar> m;
  m(0,0)=x[0]; m(0,1)=x[1]; m(0,2)=x[2]; m(0,3)=0.0;
  m(1,0)=y[0]; m(1,1)=y[1]; m(1,2)=y[2]; m(1,3)=0.0;
  m(2,0)=z[0]; m(2,1)=z[1]; m(2,2)=z[2]; m(2,3)=0.0;
  m(3,0)=static_cast<Scalar>(0.0);  m(3,1)=static_cast<Scalar>(0.0);  m(3,2)=static_cast<Scalar>(0.0);  m(3,3)=static_cast<Scalar>(1.0);

  *this *= m;
  translate(-eye[0], -eye[1], -eye[2]);
}


//-----------------------------------------------------------------------------


template<typename Scalar> 
void
GLMatrixT<Scalar>::
inverse_lookAt(const VectorT<Scalar,3>& eye,
	       const VectorT<Scalar,3>& center,
	       const VectorT<Scalar,3>& up)
{
  VectorT<Scalar,3> z(eye);
  z -= center;
  z.normalize();

  VectorT<Scalar,3> x(up % z);
  x.normalize();

  VectorT<Scalar,3> y(z % x);
  y.normalize();

  GLMatrixT<Scalar> m;
  m(0,0)=x[0]; m(0,1)=y[0]; m(0,2)=z[0]; m(0,3)=eye[0];
  m(1,0)=x[1]; m(1,1)=y[1]; m(1,2)=z[1]; m(1,3)=eye[1];
  m(2,0)=x[2]; m(2,1)=y[2]; m(2,2)=z[2]; m(2,3)=eye[2];
  m(3,0)=static_cast<Scalar>(0.0);  m(3,1)=static_cast<Scalar>(0.0);  m(3,2)=static_cast<Scalar>(0.0);  m(3,3)=static_cast<Scalar>(1.0);

  leftMult(m);
}


//-----------------------------------------------------------------------------


template<typename Scalar> 
void
GLMatrixT<Scalar>::
perspective(Scalar fovY, Scalar aspect, 
	    Scalar near_plane, Scalar far_plane)
{
  Scalar ymax = near_plane * tan( fovY * static_cast<Scalar>(M_PI) / static_cast<Scalar>(360.0) );
  Scalar ymin = -ymax;
  Scalar xmin = ymin * aspect;
  Scalar xmax = ymax * aspect;
  frustum(xmin, xmax, ymin, ymax, near_plane, far_plane);
}


//-----------------------------------------------------------------------------


template<typename Scalar> 
void
GLMatrixT<Scalar>::
inverse_perspective(Scalar fovY, Scalar aspect, 
		    Scalar near_plane, Scalar far_plane)
{
  Scalar ymax = near_plane * tan( fovY *  static_cast<Scalar>(M_PI) / static_cast<Scalar>(360.0)  );
  Scalar ymin = -ymax;
  Scalar xmin = ymin * aspect;
  Scalar xmax = ymax * aspect;
  inverse_frustum(xmin, xmax, ymin, ymax, near_plane, far_plane);
}


//-----------------------------------------------------------------------------


template<typename Scalar> 
void
GLMatrixT<Scalar>::
frustum( Scalar left,Scalar right, 
	 Scalar bottom, Scalar top, 
	 Scalar near_plane, Scalar far_plane )
{
  assert(near_plane > static_cast<Scalar>(0.0) && far_plane > near_plane);

  GLMatrixT<Scalar> m;
  Scalar x, y, a, b, c, d;

  x = (static_cast<Scalar>(2.0)*near_plane) / (right-left);
  y = (static_cast<Scalar>(2.0)*near_plane) / (top-bottom);
  a = (right+left) / (right-left);
  b = (top+bottom) / (top-bottom);
  c = -(far_plane+near_plane) / ( far_plane-near_plane);
  d = -(static_cast<Scalar>(2.0)*far_plane*near_plane) / (far_plane-near_plane);

  m(0,0) = x;                        m(0,1) = static_cast<Scalar>(0.0); m(0,2) = a;                         m(0,3) = static_cast<Scalar>(0.0);
  m(1,0) = static_cast<Scalar>(0.0); m(1,1) = y;                        m(1,2) = b;                         m(1,3) = static_cast<Scalar>(0.0);
  m(2,0) = static_cast<Scalar>(0.0); m(2,1) = static_cast<Scalar>(0.0); m(2,2) = c;                         m(2,3) = d;
  m(3,0) = static_cast<Scalar>(0.0); m(3,1) = static_cast<Scalar>(0.0); m(3,2) = static_cast<Scalar>(-1.0); m(3,3) = static_cast<Scalar>(0.0);

  *this *= m;
}


//-----------------------------------------------------------------------------


template<typename Scalar> 
void
GLMatrixT<Scalar>::
inverse_frustum(Scalar left,Scalar right,
		Scalar bottom, Scalar top,
		Scalar near_plane, Scalar far_plane)
{
  assert(near_plane > 0.0 && far_plane > near_plane);

  GLMatrixT<Scalar> m;
  Scalar x, y, a, b, c, d;

  x = (right-left) / (static_cast<Scalar>(2.0)*near_plane);
  y = (top-bottom) / (static_cast<Scalar>(2.0)*near_plane);
  a = (right+left) / (static_cast<Scalar>(2.0)*near_plane);
  b = (top+bottom) / (static_cast<Scalar>(2.0)*near_plane);
  c = (far_plane+near_plane) / (static_cast<Scalar>(2.0)*far_plane*near_plane);
  d = (near_plane-far_plane) / (static_cast<Scalar>(2.0)*far_plane*near_plane);

  m(0,0)=x;                         m(0,1)= static_cast<Scalar>(0.0);  m(0,2)= static_cast<Scalar>(0.0);  m(0,3)= a;
  m(1,0)=static_cast<Scalar>(0.0);  m(1,1)= y;                         m(1,2)= static_cast<Scalar>(0.0);  m(1,3)= b;
  m(2,0)=static_cast<Scalar>(0.0);  m(2,1)= static_cast<Scalar>(0.0);  m(2,2)= static_cast<Scalar>(0.0);  m(2,3)= static_cast<Scalar>(-1.0);
  m(3,0)=static_cast<Scalar>(0.0);  m(3,1)= static_cast<Scalar>(0.0);  m(3,2)= d;                         m(3,3)= c;

  leftMult(m);
}


//-----------------------------------------------------------------------------


template<typename Scalar> 
void
GLMatrixT<Scalar>::
ortho(Scalar left, Scalar right,
      Scalar bottom, Scalar top,
      Scalar near_plane, Scalar far_plane)
{
  GLMatrixT<Scalar> m;

  Scalar r_l = right-left;
  Scalar t_b = top - bottom;
  Scalar f_n = far_plane - near_plane;

//   assert(r_l > 0.0 && t_b > 0.0  && f_n > 0.0);

  m(1,0) = m(0,1) = m(0,2) =
  m(2,0) = m(2,1) = m(1,2) =
  m(3,0) = m(3,1) = m(3,2) = static_cast<Scalar>(0.0);
  
  m(0,0) = static_cast<Scalar>( 2.0) / r_l;
  m(1,1) = static_cast<Scalar>( 2.0) / t_b;
  m(2,2) = static_cast<Scalar>(-2.0) / f_n;
  
  m(0,3) = -(right + left) / r_l;
  m(1,3) = -(top + bottom) / t_b;
  m(2,3) = -(far_plane + near_plane) / f_n;

  m(3,3) = static_cast<Scalar>(1.0);

  *this *= m;
}


//-----------------------------------------------------------------------------

template<typename Scalar> 
void
GLMatrixT<Scalar>::
inverse_ortho(Scalar left, Scalar right,
	      Scalar bottom, Scalar top,
	      Scalar near_plane, Scalar far_plane)
{
  GLMatrixT<Scalar> m;

  m(1,0) = m(0,1) = m(0,2) =
  m(2,0) = m(2,1) = m(1,2) =
  m(3,0) = m(3,1) = m(3,2) = static_cast<Scalar>(0.0);
  

  m(0,0) = static_cast<Scalar>( 0.5) * (right - left);
  m(1,1) = static_cast<Scalar>( 0.5) * (top - bottom);
  m(2,2) = static_cast<Scalar>(-0.5) * (far_plane - near_plane);
  
  m(0,3) = static_cast<Scalar>( 0.5) * (right + left);
  m(1,3) = static_cast<Scalar>( 0.5) * (top + bottom);
  m(2,3) = static_cast<Scalar>(-0.5) * (far_plane + near_plane);

  m(3,3) = static_cast<Scalar>(1.0);

  leftMult(m);
}


//-----------------------------------------------------------------------------

template<typename Scalar> 
VectorT<Scalar, 2>
  GLMatrixT<Scalar>::
  extract_planes_perspective() const
{
  Scalar a = (*this)(2,2);
  Scalar b = (*this)(2,3);

  return VectorT<Scalar, 2>(b / (static_cast<Scalar>(-1.0) + a), b / (static_cast<Scalar>(1.0) + a));
}

//-----------------------------------------------------------------------------

template<typename Scalar> 
VectorT<Scalar, 2>
  GLMatrixT<Scalar>::
  extract_planes_ortho() const
{
  Scalar a = (*this)(2,2);
  Scalar b = (*this)(2,3);

  return VectorT<Scalar, 2>((1.0 + b) / a, (-1.0 + b) / a);
}

//-----------------------------------------------------------------------------

template<typename Scalar> 
bool
  GLMatrixT<Scalar>::
  isPerspective() const
{
  // check if matrix matches the pattern of frustum()

  // expect at least 5 nonzeros
  const int nnz = 5;

  // nonzero entries
  const int nze[] = 
  {
    0,0,
    1,1,
    2,2,
    2,3,
    3,2
  };

  // expect at least 9 zeros
  const int nz = 9;

  // zero entries
  const int ze[] =
  {
    0,1,
    0,3,
    1,0,
    1,3,
    2,0,
    2,1,
    3,0,
    3,1,
    3,3,
  };
  
  for (int i = 0; i < nnz; ++i)
  {
    if ( checkEpsilon( (*this)(nze[i*2],nze[i*2+1]) ) )
      return false;
  }

  // expect -1 at (3,2)
  if ( !checkEpsilon( (*this)(3,2) + 1.0 ) )
    return false;

  // check rest for zero
  for (int i = 0; i < nz; ++i)
  {
    if ( !checkEpsilon( (*this)(ze[i*2],ze[i*2+1]) ) )
      return false;
  }

  return true; 
}

//-----------------------------------------------------------------------------

template<typename Scalar> 
bool
  GLMatrixT<Scalar>::
  isOrtho() const
{
  // check if matrix matches the pattern of ortho()

  // expect at least 5 nonzeros (diagonal) and m(2,3)

  for (int i = 0; i < 3; ++i)
  {
    if ( checkEpsilon( (*this)(i,i) ) )
      return false;
  }

  // expect 1 at (3,3)
  if ( !checkEpsilon( (*this)(3,3) - 1.0 ) )
    return false;

  // m(2,3) must be nonzero
  if ( checkEpsilon( (*this)(2,3) ) )
    return false;


  // expect at least 9 zeros
  const int nz = 9;

  // zero entries
  const int ze[] =
  {
    0,1,
    0,2,
    1,0,
    1,2,
    2,0,
    2,1,
    3,0,
    3,1,
    3,2,
  };

  // check rest for zero (except last column)
  for (int i = 0; i < nz; ++i)
  {
    if ( !checkEpsilon( (*this)(ze[i*2],ze[i*2+1]) ) )
      return false;
  }

  return true; 
}

//-----------------------------------------------------------------------------

template<typename Scalar> 
VectorT<Scalar, 2>
  GLMatrixT<Scalar>::
  extract_planes() const
{
  if (isPerspective())
    return extract_planes_perspective();

  if (isOrtho())
    return extract_planes_ortho();

  // return invalid planes
  return VectorT<Scalar,2>(1.0, -1.0);
}

//-----------------------------------------------------------------------------


#undef MAT
#undef M


//=============================================================================
} // namespace ACG
//=============================================================================
