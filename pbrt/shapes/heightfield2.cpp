
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
#include "core/geometry.h"
#include "core/transform.h"
#include "paramset.h"

#include <stdexcept>
#include <limits>
#include <vector>

using std::vector;

struct Heightfield2_impl {
  float *z;
  int nx, ny;
  int www;
};

// Heightfield2 Method Definitions
Heightfield2::Heightfield2(const Transform *o2w, const Transform *w2o,
    bool ro, int x, int y, const float *zs)
  : Shape(o2w, w2o, ro) {
  impl = new Heightfield2_impl;
  impl->nx = x;
  impl->ny = y;
  impl->z = new float[impl->nx*impl->ny];
  memcpy(impl->z, zs, impl->nx*impl->ny*sizeof(float));
}


Heightfield2::~Heightfield2() {
  delete[] impl->z;
  delete impl;
}


BBox Heightfield2::ObjectBound() const {
  float minz = impl->z[0], maxz = impl->z[0];
  for (int i = 1; i < impl->nx*impl->ny; ++i) {
    if (impl->z[i] < minz) minz = impl->z[i];
    if (impl->z[i] > maxz) maxz = impl->z[i];
  }
  return BBox(Point(0,0,minz), Point(1,1,maxz));
}


void Heightfield2::Refine(vector<Reference<Shape>> &) const { throw std::runtime_error("Heightfield2::Refine"); }

static const float EPS = 1e-9;

inline static bool jiao_tri(
  const Ray &ray,
  const Point& pt0,
  const Point& pt1,
  const Point& pt2,
  float* tHit)
{
  const Vector u {pt1 - pt0}, w {pt2 - pt0};
  const Vector n {Cross(u, w)};
  const float deno = Dot(n, ray.d);

  if (-EPS<deno && deno<EPS)
    return false;

  const float nume = Dot(n, pt0-ray.o);
  const float t = nume/deno;
  if (t < ray.mint || t > ray.maxt)
    return false;

  const Vector v {ray(t) - pt0};
  float alpha, beta, A[2][2] {{u.x, w.x},{u.y, w.y}}, B[2] {v.x, v.y};
  if (!SolveLinearSystem2x2(A, B, &alpha, &beta))
    return false;

  *tHit = t;

  return 0<=alpha && alpha<=1 && 0<=beta && beta<=1;
}

bool Heightfield2::Intersect(
  const Ray &r,
  float *tHit,
  float *rayEpsilon,
  DifferentialGeometry *dg) const
{
  vector<vector<Point>> pt;
  for (int y = 0; y < impl->ny-1; ++y) {
    pt.push_back({});
    for (int x = 0; x < impl->nx-1; ++x) {
      pt.back().emplace_back( static_cast<float>(x)/(impl->nx-1)
                            , static_cast<float>(y)/(impl->ny-1)
                            , impl->z[y*impl->ny + x] );
    }
  }

  Ray ray;
  (*WorldToObject)(r, &ray);

  bool hit = false;
  float tmin = std::numeric_limits<float>::max()  ;
  int argtri[3][2];
  for (int y = 0; y < impl->ny-1; ++y) {
    for (int x = 0; x < impl->nx-1; ++x) {
      float t;
      if (jiao_tri(ray, pt[y][x], pt[y][x+1], pt[y+1][x+1], &t)) {
        hit = true;
        if (t < tmin) {
          tmin = t;
          argtri[0][0] = x;   argtri[0][1] = y;
          argtri[1][0] = x+1; argtri[1][1] = y;
          argtri[2][0] = x+1; argtri[2][1] = y+1;
        }
      }
      if (jiao_tri(ray, pt[y][x], pt[y+1][x+1], pt[y+1][x], &t)) {
        hit = true;
        if (t < tmin) {
          tmin = t;
          argtri[0][0] = x;   argtri[0][1] = y;
          argtri[1][0] = x+1; argtri[1][1] = y+1;
          argtri[2][0] = x; argtri[2][1] = y+1;
        }
      }
    }
  }
  if (!hit) return false;

  const Point p {ray(tmin)};
  const Vector dpdu{1,0,0}, dpdv{0,1,0};
  const Vector n {Normalize(Cross(dpdu,dpdv))};
  Normal dndu, dndv;

  const Transform &o2w = *ObjectToWorld;
  *dg = { o2w(p), o2w(dpdu), o2w(dpdv), o2w(dndu), o2w(dndv), p.x, p.y, this };
  dg->dudx = 1;
  dg->dvdy = 1;

  *tHit = tmin;
  *rayEpsilon = 1e-3f * tmin;
  return true;
}


bool Heightfield2::IntersectP(const Ray &ray) const {
  float tHit, rayEpsilon;
  DifferentialGeometry dg;
  return this->Intersect(ray, &tHit, &rayEpsilon, &dg);
}

void Heightfield2::GetShadingGeometry(
  const Transform &obj2world,
  const DifferentialGeometry &dg,
  DifferentialGeometry *dgShading) const
{
  *dgShading = dg;
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


