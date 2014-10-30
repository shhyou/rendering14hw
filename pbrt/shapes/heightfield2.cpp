
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
using std::array;
using std::pair;

struct Heightfield2_impl {
  float *z;
  int nx, ny;
  vector<Point> pt;
  vector<vector<vector<pair<int,int>>>> tri;
};

#define COORD(y,x) ((y)*impl->nx + (x))

// Heightfield2 Method Definitions
Heightfield2::Heightfield2(const Transform *o2w, const Transform *w2o,
    bool ro, int x, int y, const float *zs)
  : Shape(o2w, w2o, ro) {
  impl = new Heightfield2_impl;
  impl->nx = x;
  impl->ny = y;
  impl->z = new float[impl->nx*impl->ny];
  memcpy(impl->z, zs, impl->nx*impl->ny*sizeof(float));
  for (int y = 0; y < impl->ny; ++y) {
    for (int x = 0; x < impl->nx; ++x) {
      impl->pt.emplace_back( static_cast<float>(x)/(impl->nx-1)
                           , static_cast<float>(y)/(impl->ny-1)
                           , impl->z[COORD(y,x)] );
    }
  }
  for (int y = 0; y < impl->ny; ++y) {
    impl->tri.push_back({});
    for (int x = 0; x < impl->nx; ++x) {
      impl->tri.back().push_back({{y,x}, {y,x+1},   {y+1,x+1}});
      impl->tri.back().push_back({{y,x}, {y+1,x+1}, {y+1,x}});
    }
  }
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

inline static bool jiao_plane(
  const Ray &ray,
  const Point& pt0,
  const Vector& n,
  float *tHit)
{
  const float deno = Dot(n, ray.d);

  if (-EPS<deno && deno<EPS)
    return false;

  const float nume = Dot(n, pt0-ray.o);
  const float t = nume/deno;
  if (t < ray.mint || t > ray.maxt)
    return false;

  *tHit = t;
  return true;
}

inline static bool jiao_tri(
  const Heightfield2_impl *impl,
  const Ray &ray,
  const vector<pair<int,int>>& _pt,
  float* tHit)
{
  const Point &pt0 = impl->pt[COORD(_pt[0].first, _pt[0].second)]
            , &pt1 = impl->pt[COORD(_pt[1].first, _pt[1].second)]
            , &pt2 = impl->pt[COORD(_pt[2].first, _pt[2].second)];
  const Vector u {pt1 - pt0}, w {pt2 - pt0};
  const Vector n {Cross(u, w)};
  float t;
  if (!jiao_plane(ray, pt0, n, &t))
    return false;

  const Vector v {ray(t) - pt0};
  float alpha, beta;
  const float A[2][2] {{u.x, w.x},{u.y, w.y}}, B[2] {v.x, v.y};
  if (!SolveLinearSystem2x2(A, B, &alpha, &beta))
    return false;

  *tHit = t;

  return 0<=alpha && alpha<=1 && 0<=beta && beta<=1 && alpha+beta<=1;
}

bool Heightfield2::Intersect(
  const Ray &r,
  float *tHit,
  float *rayEpsilon,
  DifferentialGeometry *dg) const
{
//  puts("A");
  vector<vector<vector<pair<int,int>>>> &tri = impl->tri;
  Ray ray;
  (*WorldToObject)(r, &ray);

//  puts("B");
  bool hit = false;
  float tmin = std::numeric_limits<float>::max();
  int argtri[2] {0,0};

  const float dy = 1.0/(impl->ny - 1);
  if (-EPS<ray.d.y && ray.d.y<EPS) { // horizontal
    float y = 0;
    int yidx;
    for (yidx = 0; yidx < impl->ny-1; ++yidx, y += dy)
      if (y<=ray.o.y && ray.o.y<=y+dy) break;
    for (int xidx = 0; xidx < 2*(impl->nx-1); ++xidx) {
      float t;
      if (jiao_tri(this->impl, ray, tri[yidx][xidx], &t)) {
        hit = true;
        if (t < tmin) {
          tmin = t;
          argtri[0] = yidx;
          argtri[1] = xidx;
        }
      }
    }
  } else {
    float t0, t1;
    if (!jiao_plane(ray, Point{0,0,0}, Vector{0,1,0}, &t0)
       || !jiao_plane(ray, Point{0,dy,0}, Vector{0,1,0}, &t1))
      return false;
    const float dt = t1 - t0;

    float dx, x, tray = t0;
    int ix, dix;
    if (ray.d.x >= 0) {
      dx = 1.0/(impl->nx-1);
      x = 0;
      ix = 0;
      dix = 1;
    } else {
      dx = -1.0/(impl->nx-1);
      x = 1 + dx;
      ix = impl->nx-2;
      dix = -1;
    }
    const float absdx = fabs(dx);
    for (int yidx = 0; yidx < impl->ny-1; ++yidx, tray += dt) {
      const float xtop = ray(tray+dt).x;
      while (0 <= ix && ix < impl->nx-1) {
        for (int k = 0; k < 2; ++k) {
          float t;
          if (jiao_tri(this->impl, ray, tri[yidx][2*ix+k], &t)) {
            hit = true;
            if (t < tmin) {
              tmin = t;
              argtri[0] = yidx;
              argtri[1] = 2*ix+k;
            }
          }
        }
        if (x<=xtop && xtop<=x+absdx) break;
        x += dx;
        ix += dix;
      }
      if (ix<0 || ix>=impl->nx-1) break;
    }
  }

//  printf("C: hit=%d",hit);
  if (!hit) return false;

  const vector<pair<int,int>> &_pt = impl->tri[argtri[0]][argtri[1]];
//  puts("D");
  const Point p {ray(tmin)}
            , &pt0 = impl->pt[COORD(_pt[0].first, _pt[0].second)]
            , &pt1 = impl->pt[COORD(_pt[1].first, _pt[1].second)]
            , &pt2 = impl->pt[COORD(_pt[2].first, _pt[2].second)];
  const Vector ntri {Cross(pt1-pt0, pt2-pt0)};
  const Vector dpdu {1,0,-ntri.x/ntri.z}, dpdv {0,1,-ntri.y/ntri.z};
  const Vector n {Normalize(Cross(dpdu,dpdv))};
  const Normal dndu {0,0,0}, dndv {0,0,0};

//  puts("E");
  const Transform &o2w = *ObjectToWorld;
  *dg = { o2w(p), o2w(dpdu), o2w(dpdv), o2w(dndu), o2w(dndv), p.x, p.y, this };
  dg->dudx = 1;
  dg->dvdy = 1;

//  puts("F");
  *tHit = tmin;
  *rayEpsilon = 1e-3f * tmin;
//  puts("G");
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

