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

#include "PlaneType.hh"

namespace ACG {
namespace Geometry {


void Plane::setPlane(const ACG::Vec3d& _position, const ACG::Vec3d& _xDirection, const ACG::Vec3d& _yDirection)
{
  position   = _position;
  xDirection = _xDirection;
  yDirection = _yDirection;
  normal     = (_xDirection % _yDirection).normalize();
}

//----------------------------------------------------------------

void Plane::setPlane(const ACG::Vec3d& _position, const ACG::Vec3d& _normal)
{

  //find a non zero component
  int comp = -1;
  for (int i=0; i < 3; i++)
    if ( _normal[i] != 0.0 ){
      comp = i;
      break;
    }

  if (comp == -1){
    std::cerr << "PlaneNode: normal is invalid!" << std::endl;
    return;
  }

  //compute orthogonal vectors in the plane
  xDirection[comp] = (-_normal[ (comp + 1) % 3 ] - _normal[(comp + 2) % 3]) / _normal[comp];
  xDirection[ (comp + 1) % 3 ] = 1;
  xDirection[ (comp + 2) % 3 ] = 1;
  xDirection = xDirection.normalize();

  yDirection = _normal % xDirection;
  yDirection = yDirection.normalize();

  position = _position;
  normal   = _normal;
}

//----------------------------------------------------------------

void Plane::transform(const ACG::Matrix4x4d& _mat)
{
  position    = _mat.transform_point(position);
  xDirection  = _mat.transform_vector(xDirection);
  yDirection  = _mat.transform_vector(yDirection);

  normal      = (xDirection % yDirection).normalize();
}

//----------------------------------------------------------------

void Plane::setSize(double _xDirection, double _yDirection)
{
  xDirection = xDirection.normalize() * _xDirection;
  yDirection = yDirection.normalize() * _yDirection;
}

}} // namespaces

//=============================================================================
