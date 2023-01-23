#pragma once
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

#include <ACG_BSP/Math/Matrix4x4T.hh>
#include <ACG_BSP/Math/VectorT.hh>
#include <ACG_BSP/Config/ACGDefines.hh>

namespace ACG {
namespace Geometry {

class ACGDLLEXPORT Plane {

public:
    Plane() {}
    Plane(const ACG::Vec3d& _p, const ACG::Vec3d& _n, const ACG::Vec3d& _x, const ACG::Vec3d& _y) :
        position(_p), normal(_n), xDirection(_x), yDirection(_y) {
    }

  /** \brief Set plane
   *
   * @param _position   One point on the plane. Will be used as corner point point for rendering in the PlaneNode
   * @param _xDirection Vector pointing in planes x direction
   * @param _yDirection Vector pointing in planes y direction
   */
  void setPlane(const ACG::Vec3d& _position, const ACG::Vec3d& _xDirection, const ACG::Vec3d& );

   /** \brief Set plane with given normal and one point
    *
    * @param _position One point on the plane. Will be used as corner point for rendering in the PlaneNode
    * @param _normal   Plane normal
    */
  void setPlane(const ACG::Vec3d & _position, const ACG::Vec3d & _normal);

  /** \brief Set plane size
   *
   * Scales the plane such that the x and y direction vectors have the given lengths
   *
   * @param _xDirection Size in x direction
   * @param _yDirection Size in y direction
   */
  void setSize(double _xDirection, double _yDirection);

  /** \brief Transform the plane with given matrix
   *
   *
   * @param _mat Transformation matrix.
   */
  void transform(const ACG::Matrix4x4d & _mat);

public:

  ACG::Vec3d position;
  ACG::Vec3d normal;
  ACG::Vec3d xDirection;
  ACG::Vec3d yDirection;

};

}} // namespaces
