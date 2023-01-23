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
//  CLASS GLMatrixT
//
//=============================================================================


#ifndef ACG_GLMATRIX_HH
#define ACG_GLMATRIX_HH


//== INCLUDES =================================================================


#include "Matrix4x4T.hh"
#include "../Config/ACGDefines.hh"
#include <cmath>


namespace ACG {

	      
//== CLASS DEFINITION =========================================================


/** enum to chose whether to multiply new matrices from right
    (default) or left to the current matrix */
enum MultiplyFrom { MULT_FROM_RIGHT, MULT_FROM_LEFT };



/// 4x4 matrix implementing OpenGL commands.
template <class Scalar>
class GLMatrixT : public Matrix4x4T<Scalar>
{
public:
   
  typedef VectorT<Scalar, 3> Vec3;


  /// constructor: uninitialized values
  GLMatrixT() {}
 
  /// construct from other matrix type
  template <class OtherScalar>
  inline GLMatrixT(const Matrix4x4T<OtherScalar>& _rhs)
    : Matrix4x4T<Scalar>(_rhs)
  {}

  /// copy constructor
  inline GLMatrixT(const GLMatrixT<Scalar>& _rhs) = default;

  /** setup matrix using an array of N*N scalar values.
      elements are ordered 'column first' (like OpenGL) */
  inline GLMatrixT(const Scalar _array[16]) : Matrix4x4T<Scalar>(_array) {}

  /**
   * Initialize an affine matrix from column vectors.
   */
  inline GLMatrixT(const Vec3 &col1, const Vec3 &col2, const Vec3 &col3) {
      /*
       * Don't try to optimize anything by hand, here. gcc -O2 does the right thing:
       *
       * mov    %rax,-0x70(%rsp)
       * mov    %rdx,-0x68(%rsp)
       * mov    %rsi,-0x50(%rsp)
       * lea    -0x68(%rsp),%rdx
       * mov    %rax,0x8(%rsp)
       * mov    %rcx,-0x60(%rsp)
       * lea    0x10(%rsp),%rsi
       * movq   $0x0,-0x58(%rsp)
       * mov    %rdi,-0x48(%rsp)
       * xor    %eax,%eax
       * mov    %r8,-0x40(%rsp)
       * movq   $0x0,-0x38(%rsp)
       * mov    %r9,-0x30(%rsp)
       * mov    %r10,-0x28(%rsp)
       * mov    %r11,-0x20(%rsp)
       * movq   $0x0,-0x18(%rsp)
       * movq   $0x0,-0x10(%rsp)
       * movq   $0x0,-0x8(%rsp)
       * movq   $0x0,(%rsp)
       *
       */
      memcpy(this->mat_ + 0, col1.data(), sizeof(Scalar) * 3);
      this->mat_[3] = 0;
      memcpy(this->mat_ + 4, col2.data(), sizeof(Scalar) * 3);
      this->mat_[7] = 0;
      memcpy(this->mat_ + 8, col3.data(), sizeof(Scalar) * 3);
      for (int i = 11; i < 15; ++i) this->mat_[i] = 0;
      this->mat_[15] = 1;
  }


  /// destructor
  ~GLMatrixT() {}


  /// assignement from other matrix type
  template<typename otherScalar>
  inline GLMatrixT<Scalar>& operator=(const GLMatrixT<otherScalar>& _rhs)
  { Matrix4x4T<Scalar>::operator=(_rhs); return *this; }

  /// assignement from other matrix type
  template<typename otherScalar>
  inline GLMatrixT<Scalar>& operator=(const Matrix4x4T<otherScalar>& _rhs)
  { Matrix4x4T<Scalar>::operator=(_rhs); return *this; }


      
  /// multiply self with scaling matrix (x,y,z)
  inline void scale( Scalar _x, Scalar _y, Scalar _z, 
		     MultiplyFrom _mult_from = MULT_FROM_RIGHT );
  /// multiply self with scaling matrix (x,y,z)
  inline void scale( const Vec3& _v,
		     MultiplyFrom _mult_from = MULT_FROM_RIGHT ) {
    scale(_v[0], _v[1], _v[2], _mult_from);
  }


  /// multiply self with translation matrix (x,y,z)
  inline void translate( Scalar _x, Scalar _y, Scalar _z,
			 MultiplyFrom _mult_from = MULT_FROM_RIGHT );
  /// multiply self with translation matrix (x,y,z)
  inline void translate( const Vec3& _v,
			 MultiplyFrom _mult_from = MULT_FROM_RIGHT ) {
    translate(_v[0], _v[1], _v[2], _mult_from);
  }


  /** multiply self with a rotation matrix
      (angle in degree, arbitrary axis given by xyz) */
  void rotate( Scalar angle, Scalar x, Scalar y, Scalar z,
	       MultiplyFrom _mult_from = MULT_FROM_RIGHT );
  /** multiply self with a rotation matrix
      (angle in degree, arbitrary axis given by xyz) */
  void rotate( Scalar _angle, const Vec3& _axis,
	       MultiplyFrom _mult_from = MULT_FROM_RIGHT ) {
    rotate(_angle, _axis[0], _axis[1], _axis[2], _mult_from);
  }



  /// multiply self with a rotation matrix (angle in degree, x-axis) 
  inline void rotateX( Scalar _angle,
		       MultiplyFrom _mult_from = MULT_FROM_RIGHT ) {
    rotateXYZ( X, _angle, _mult_from );
  }

  /// multiply self with a rotation matrix (angle in degree, y-axis) 
  inline void rotateY( Scalar _angle,
		       MultiplyFrom _mult_from = MULT_FROM_RIGHT ) {
    rotateXYZ( Y, _angle, _mult_from );
  }

  /// multiply self with a rotation matrix (angle in degree, z-axis) 
  inline void rotateZ( Scalar _angle,
		       MultiplyFrom _mult_from = MULT_FROM_RIGHT ) {
    rotateXYZ( Z, _angle, _mult_from );
  }




  /** multiply self with a viewing transformation given by 
      eye position, reference point (center) and an up vector.
      (similar to gluLookAt) */
  void lookAt(const Vec3& eye,
	      const Vec3& center,
	      const Vec3& up);

  /// multiply self from left with inverse lookAt matrix
  void inverse_lookAt(const Vec3& eye,
		      const Vec3& center,
		      const Vec3& up);


  /** \brief multiply self with a perspective projection matrix
   *
   * The specified perspective projection will be multiplied with the
   * matrix already stored in this matix.
   *
   * @param fovY        Half of the Field of View in y direction angle (Degree)
   * @param aspect      aspect ratio: x = y * aspect
   * @param near_plane  Distance to near plane ( > 0 )
   * @param far_plane   Distance to far plane ( > near_plane > 0 )
   */
  void perspective(Scalar fovY, Scalar aspect, 
		   Scalar near_plane, Scalar far_plane);

  /// multiply self from left with inverse of perspective projection matrix 
  void inverse_perspective(Scalar fovY, Scalar aspect,
			   Scalar near_plane,Scalar far_plane);

  /// multiply self with a perspective projection matrix 
  void frustum(Scalar left, Scalar right,
	       Scalar bottom, Scalar top,
	       Scalar near_plane, Scalar far_plane);

  /// multiply self from left with inverse of perspective projection matrix 
  void inverse_frustum(Scalar left,Scalar right,
		       Scalar bottom, Scalar top,
		       Scalar near_plane, Scalar far_plane);

  /// multiply self with orthographic projection matrix
  void ortho(Scalar left, Scalar right,
	     Scalar bottom, Scalar top,
	     Scalar near_plane, Scalar far_plane);

  /// multiply self from left with inverse orthographic projection matrix
  void inverse_ortho(Scalar left, Scalar right,
		     Scalar bottom, Scalar top,
		     Scalar near_plane, Scalar far_plane);



  /// extract near and far clipping planes from a perspective projection matrix
  VectorT<Scalar, 2> extract_planes_perspective() const;

  /// extract near and far clipping planes from an orthographic projection matrix
  VectorT<Scalar, 2> extract_planes_ortho() const;

  /// check if the matrix is a perspective projection matrix
  bool isPerspective() const;

  /// check if the matrix is an orthographic projection matrix
  bool isOrtho() const;

  /// detect type of projection matrix and extract near and far clipping planes
  VectorT<Scalar, 2> extract_planes() const;

  //----------------------------------------------------- overloaded operators 

  GLMatrixT& operator+= ( const Matrix4x4T<Scalar>& _rhs) {
    Matrix4x4T<Scalar>::operator+=(_rhs); return *this;
  }
  GLMatrixT& operator-= ( const Matrix4x4T<Scalar>& _rhs) {
    Matrix4x4T<Scalar>::operator-=(_rhs); return *this;
  }
  GLMatrixT& operator*= ( const Matrix4x4T<Scalar>& _rhs) {
    Matrix4x4T<Scalar>::operator*=(_rhs); return *this;
  }
  GLMatrixT& leftMult(const Matrix4x4T<Scalar>& _rhs) {
    Matrix4x4T<Scalar>::leftMult(_rhs); return *this;
  }

  GLMatrixT operator+ (const Matrix4x4T<Scalar>& _rhs) const {
    return GLMatrixT<Scalar>(*this) += _rhs;
  }    
  GLMatrixT operator- (const Matrix4x4T<Scalar>& _rhs) const {
    return GLMatrixT<Scalar>(*this) -= _rhs;
  }
  GLMatrixT operator*(const Matrix4x4T<Scalar>& _rhs) const {
    return GLMatrixT<Scalar>(*this) *= _rhs;
  }

  template <typename T>
  inline VectorT<T,4> operator*(const VectorT<T,4>& _v) const {
    return Matrix4x4T<Scalar>::operator*(_v);
  }



private:

  enum Axis { X, Y, Z };
  void rotateXYZ( Axis _axis, Scalar _angle, MultiplyFrom _mult_from );
};




//=============================================================================


/// typedef
typedef GLMatrixT<float>  GLMatrixf;
/// typedef
typedef GLMatrixT<double> GLMatrixd;


//=============================================================================
} // namespace ACG
//=============================================================================
#if defined(INCLUDE_TEMPLATES) && !defined(ACG_GLMATRIX_C)
#define ACG_GLMATRIX_TEMPLATES
#include "GLMatrixT_impl.hh"
#endif
//=============================================================================
#endif // ACG_GLMATRIX_HH defined
//=============================================================================

