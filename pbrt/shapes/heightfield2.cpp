
/*
  pbrt source code Copyright(c) 1998-2012 Matt Pharr and Greg Humphreys.

  This file is part of pbrt.

  Redistribution and use in source and binary forms, with or without
  modification, are permitted provided that the following conditions are
  met:

  - Redistributions of source code must retain the above copyright
    notice, this list of conditions and the following disclaimer.

  - Redistributions in binary form must reproduce the above copyright
    notice, this list of conditions and the following disclaimer in the
    documentation and/or other materials provided with the distribution.

  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS
  IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED
  TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
  PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
  HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
  SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
  LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
  DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
  THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
  (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
  OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

 */

// shapes/heightfield2.cpp*
#include "stdafx.h"
#include "shapes/heightfield2.h"
#include "shapes/trianglemesh.h"
#include "paramset.h"

#include <stdexcept>

// Heightfield2 Method Definitions
Heightfield2::Heightfield2(const Transform *o2w, const Transform *w2o,
    bool ro, int x, int y, const float *zs)
  : Shape(o2w, w2o, ro) {
  nx = x;
  ny = y;
  z = new float[nx*ny];
  memcpy(z, zs, nx*ny*sizeof(float));
}


Heightfield2::~Heightfield2() {
  delete[] z;
}


BBox Heightfield2::ObjectBound() const {
  float minz = z[0], maxz = z[0];
  for (int i = 1; i < nx*ny; ++i) {
    if (z[i] < minz) minz = z[i];
    if (z[i] > maxz) maxz = z[i];
  }
  return BBox(Point(0,0,minz), Point(1,1,maxz));
}


void Heightfield2::Refine(vector<Reference<Shape> > &) const { throw std::runtime_error("Heightfield2::Refine"); }

bool Heightfield2::Intersect(
  const Ray &ray,
  float *tHit,
  float *rayEpsilon,
  DifferentialGeometry *dg) const
{
  return false;
}
bool Heightfield2::IntersectP(const Ray &ray) const {
  return false;
}


Heightfield2 *CreateHeightfieldShape2(
  const Transform *o2w,
  const Transform *w2o,
  bool reverseOrientation,
  const ParamSet &params)
{
  int nu = params.FindOneInt("nu", -1);
  int nv = params.FindOneInt("nv", -1);
  int nitems;
  const float *Pz = params.FindFloat("Pz", &nitems);
  Assert(nitems == nu*nv);
  Assert(nu != -1 && nv != -1 && Pz != NULL);
  return new Heightfield2(o2w, w2o, reverseOrientation, nu, nv, Pz);
}


