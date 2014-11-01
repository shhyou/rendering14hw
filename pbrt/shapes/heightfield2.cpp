
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

#include <cassert>
#include <stdexcept>
#include <limits>
#include <vector>

using std::vector;
using std::array;
using std::pair;

struct Heightfield2_impl {
  const int nx, ny;
  const float dx, dy;
  float * const z, minz, maxz;
  BBox bbox;
  vector<Point> pt;
  vector<vector<vector<pair<int,int>>>> tri;
};

#define COORD(y,x) ((y)*impl->nx + (x))

// Heightfield2 Method Definitions
Heightfield2::Heightfield2(const Transform *o2w, const Transform *w2o,
    bool ro, int x, int y, const float *zs)
  : Shape(o2w, w2o, ro),
    impl(new Heightfield2_impl {x, y, 1.0f/(x - 1), 1.0f/(y - 1),
                                new float[x*y]})
{
  memcpy(impl->z, zs, impl->nx*impl->ny*sizeof(float));
  ///////////////////////////////////////////////////////////////////
  for (int y = 0; y < impl->ny; ++y)
    for (int x = 0; x < impl->nx; ++x)
      impl->pt.emplace_back( x*impl->dx, y*impl->dy, impl->z[COORD(y,x)] );
  for (int y = 0; y < impl->ny-1; ++y) {
    impl->tri.push_back({});
    for (int x = 0; x < impl->nx-1; ++x) {
      impl->tri.back().push_back({{y,x}, {y,x+1},   {y+1,x+1}});
      impl->tri.back().push_back({{y,x}, {y+1,x+1}, {y+1,x}});
    }
  }
  ///////////////////////////////////////////////////////////////////
  float minz = impl->z[0], maxz = impl->z[0];
  for (int i = 1; i < impl->nx*impl->ny; ++i) {
    if (impl->z[i] < minz) minz = impl->z[i];
    if (impl->z[i] > maxz) maxz = impl->z[i];
  }
  impl->minz = minz;
  impl->maxz = maxz;
  impl->bbox = BBox(Point(0,0,impl->minz), Point(1,1,impl->maxz));
}


Heightfield2::~Heightfield2() {
  delete[] impl->z;
  delete impl;
}


BBox Heightfield2::ObjectBound() const {
  return impl->bbox;
}


void Heightfield2::Refine(vector<Reference<Shape>> &) const { throw std::runtime_error("Heightfield2::Refine"); }

static const float EPS = 1e-10;

inline static bool jiao_tri(
  const Heightfield2_impl *impl,
  const Ray &ray,
  const vector<pair<int,int>>& _pt,
  float* tHit)
{
  const Point &pt0 = impl->pt[COORD(_pt[0].first, _pt[0].second)]
            , &pt1 = impl->pt[COORD(_pt[1].first, _pt[1].second)]
            , &pt2 = impl->pt[COORD(_pt[2].first, _pt[2].second)];
  const Vector e1 = pt1 - pt0;
  const Vector e2 = pt2 - pt0;
  const Vector s1 = Cross(ray.d, e2);
  const float divisor = Dot(s1, e1);

  if (divisor == 0.)
      return false;
  const float invDivisor = 1.f / divisor;

  // Compute first barycentric coordinate
  const Vector s = ray.o - pt0;
  float b1 = Dot(s, s1) * invDivisor;
  if (b1 < 0. || b1 > 1.)
      return false;

  // Compute second barycentric coordinate
  Vector s2 = Cross(s, e1);
  float b2 = Dot(ray.d, s2) * invDivisor;
  if (b2 < 0. || b1 + b2 > 1.)
      return false;

  // Compute _t_ to intersection point
  float t = Dot(e2, s2) * invDivisor;
  if (t < ray.mint || t > ray.maxt)
    return false;

  *tHit = t;
  return true;
}

bool Heightfield2::Intersect(
  const Ray &r,
  float *tHit,
  float *rayEpsilon,
  DifferentialGeometry *dg) const
{
  Ray ray;
  (*WorldToObject)(r, &ray);

  Point pst;
  const float d_inv[3] { 1.0f/ray.d.x, 1.0f/ray.d.y, 1.0f/ray.d.z };
  if (impl->bbox.Inside(ray(ray.mint))) {
//    printf("\nintersect: inside\n");
    pst = ray(ray.mint);
  } else {
    bool hit = false;
    float t0 = std::numeric_limits<float>::max();
    for (int coord = 0; coord != 2; ++coord) {
      for (int k = 0; k < 2; ++k) {
        const float border = (k==0? impl->dx:impl->dy)/2;
        const float t = ((coord? (1.0 - border) : border) - ray.o[k])*d_inv[k];
        if (impl->bbox.Inside(ray(t)) && t<t0) {
          t0 = t;
          hit = true;
        }
      }
    }
    if (!hit)
      return false;
    pst = ray(t0);
//    printf("\nintersect: at %f (%.8f,%.8f,%.8f)\n", t0, pst.x, pst.y, pst.z);
  }

//  printf("direction: (%f,%f,%f)\n", ray.d.x, ray.d.y, ray.d.z);

  int coord[2] { static_cast<int>(pst.y * (impl->ny-1) + EPS)
               , static_cast<int>(pst.x * (impl->nx-1) + EPS) };

  int argtri[2] {0,0};

  const float t_step[2] { static_cast<float>(fabs(impl->dy * d_inv[1]))
                        , static_cast<float>(fabs(impl->dx * d_inv[0])) };
  float t_nxt[2];
  int coord_step[2];

  for (int k = 0; k != 2; ++k) {
    const float d = k==0? impl->dy : impl->dx;
    if (ray.d[1-k] >= 0) {
      t_nxt[k] = ((coord[k]+1)*d - ray.o[1-k]) * d_inv[1-k];
      coord_step[k] = 1;
    } else {
      t_nxt[k] = (coord[k]*d - ray.o[1-k]) * d_inv[1-k];
      coord_step[k] = -1;
    }
  }

  bool hit = false;
  float tmin = std::numeric_limits<float>::max();

  do {
//    printf("at (%d,%d) in (%d,%d)\n", coord[0], coord[1], impl->ny-1, impl->nx-1);
    for (int k = 0; k < 2; ++k) {
      float t;
      if (jiao_tri(impl, ray, impl->tri[coord[0]][coord[1]*2 + k], &t)) {
        hit = true;
        if (t < tmin) {
          tmin = t;
          argtri[0] = coord[0];
          argtri[1] = coord[1]*2+k;
        }
      }
    }
    if (hit || (t_nxt[0]>ray.maxt && t_nxt[1]>ray.maxt)) break;
    if (t_nxt[0] <= t_nxt[1]) {
      coord[0] += coord_step[0];
      t_nxt[0] += t_step[0];
    } else {
      coord[1] += coord_step[1];
      t_nxt[1] += t_step[1];
    }
  } while (0<=coord[0] && coord[0]<impl->ny-1
        && 0<=coord[1] && coord[1]<impl->nx-1);

  if (!hit) return false;

  const vector<pair<int,int>> &_pt = impl->tri[argtri[0]][argtri[1]];
  const Point p {ray(tmin)}
            , &pt0 = impl->pt[COORD(_pt[0].first, _pt[0].second)]
            , &pt1 = impl->pt[COORD(_pt[1].first, _pt[1].second)]
            , &pt2 = impl->pt[COORD(_pt[2].first, _pt[2].second)];
  const Vector ntri {Cross(pt1-pt0, pt2-pt0)};
  const float inv_z = 1.0/ntri.z;
  const Vector dpdu {1,0,-ntri.x/inv_z}, dpdv {0,1,-ntri.y/inv_z};
  const Vector n {Cross(dpdu,dpdv)};
  const static Normal dndu {0,0,0}, dndv {0,0,0};

  const Transform &o2w = *ObjectToWorld;
  *dg = { o2w(p), o2w(dpdu), o2w(dpdv), o2w(dndu), o2w(dndv), p.x, p.y, this };
  dg->dudx = 1;
  dg->dvdy = 1;

  *tHit = tmin;
  *rayEpsilon = 1e-3f * tmin;
  return true;
}


#if 1

bool Heightfield2::IntersectP(const Ray &ray) const {
  float tHit, rayEpsilon;
  DifferentialGeometry dg;
  return this->Intersect(ray, &tHit, &rayEpsilon, &dg);
}

#else

bool Heightfield2::IntersectP(const Ray &r) const {
  vector<vector<vector<pair<int,int>>>> &tri = impl->tri;
  Ray ray;
  (*WorldToObject)(r, &ray);

  bool hit = false;
  float tmin = std::numeric_limits<float>::max();

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
            }
          }
        }
        if (dt>0 && hit) break;
        if (x<=xtop && xtop<=x+absdx) break;
        x += dx;
        ix += dix;
      }
      if (dt>0 && hit) break;
      if (ix<0 || ix>=impl->nx-1) break;
    }
  }
  return hit;
}

#endif

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

