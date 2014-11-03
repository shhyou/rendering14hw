
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
#include <cstdio>
#include <stdexcept>
#include <limits>
#include <vector>

#define PRRT 0

using std::vector;
using std::array;
using std::pair;

struct Heightfield2_impl {
  const int nx, ny;
  const float dx, dy;
  float * const z, minz, maxz;
  BBox bbox;
  vector<Point> pt;
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
  ///////////////////////////////////////////////////////////////////
  float minz = impl->z[0], maxz = impl->z[0];
  for (int i = 1; i < impl->nx*impl->ny; ++i) {
    if (impl->z[i] < minz) minz = impl->z[i];
    if (impl->z[i] > maxz) maxz = impl->z[i];
  }
  impl->minz = minz;
  impl->maxz = maxz;
  impl->bbox = BBox(Point(0,0,impl->minz), Point(1,1,impl->maxz));
  ///////////////////////////////////////////////////////////////////
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
  const int y1, const int x1,
  const int y2, const int x2,
  const int y3, const int x3,
  float* tHit)
{
  const Point &pt0 = impl->pt[COORD(y1, x1)]
            , &pt1 = impl->pt[COORD(y2, x2)]
            , &pt2 = impl->pt[COORD(y3, x3)];
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

  float t_init;
  const float d_inv[3] { 1.0f/ray.d.x, 1.0f/ray.d.y, 1.0f/ray.d.z };
  if (impl->bbox.Inside(ray(ray.mint))) {
    t_init = ray.mint;
  } else if (!impl->bbox.IntersectP(ray, &t_init)) {
    return false;
  }

#if PRRT
  printf("\nintersect: t=%f, (%f,%f,%f);\n           ray=(%f,%f,%f)+t(%f,%f,%f)\n",
    t_init, ray(t_init)[0],ray(t_init)[1],ray(t_init)[2],
    ray.o[0],ray.o[1],ray.o[2],ray.d[0],ray.d[1],ray.d[2]);
  printf("           nx=%d,ny=%d,dx=%f,dy=%f,minz=%f,maxz=%f\n",impl->nx,impl->ny,impl->dx,impl->dy,impl->minz,impl->maxz);
#endif

  int coord[3] = { static_cast<int>((ray.o.x + t_init*ray.d.x)*(impl->nx-1) + EPS)
                 , static_cast<int>((ray.o.y + t_init*ray.d.y)*(impl->ny-1) + EPS)
                 , 0 };

  if (coord[0] == impl->nx-1) --coord[0];
  if (coord[1] == impl->ny-1) --coord[1];

  float t_nxt[3], t_step[3];
  int coord_step[3], boundary[3];

  // handle z
  if (ray.d[2] >= 0) {
    t_nxt[2] = (impl->maxz - ray.o.z)*d_inv[2];
    t_step[2] = (impl->maxz-impl->minz)*d_inv[2]; // unimportant
    coord_step[2] = 1;
    boundary[2] = 1;
  } else {
    t_nxt[2] = (impl->minz - ray.o.z)*d_inv[2];
    t_step[2] = -(impl->maxz-impl->minz)*d_inv[2]; // unimportant
    coord_step[2] = -1;
    boundary[2] = -1;
  }

  // handle y
  if (ray.d[1] >= 0) {
    t_nxt[1] = ((coord[1]+1)*impl->dy - ray.o.y)*d_inv[1];
    t_step[1] = impl->dy * d_inv[1];
    coord_step[1] = 1;
    boundary[1] = impl->ny - 1;
  } else {
    t_nxt[1] = (coord[1]*impl->dy - ray.o.y)*d_inv[1];
    t_step[1] = - impl->dy * d_inv[1];
    coord_step[1] = -1;
    boundary[1] = -1;
  }

#if PRRT
  printf("           %f v.s. %f\n",
    ((coord[0]+1)*impl->dx - ray.o.x)*d_inv[0],
    (coord[0]*impl->dx - ray.o.x)*d_inv[0]);
#endif

  // handle x
  if (ray.d[0] >= 0) {
    t_nxt[0] = ((coord[0]+1)*impl->dx - ray.o.x)*d_inv[0];
    t_step[0] = impl->dx * d_inv[0];
    coord_step[0] = 1;
    boundary[0] = impl->nx - 1;
  } else {
    t_nxt[0] = (coord[0]*impl->dx - ray.o.x)*d_inv[0];
    t_step[0] = - impl->dx * d_inv[0];
    coord_step[0] = -1;
    boundary[0] = -1;
  }

  bool hit = false;
  int argtri[3] {0,0,0};
  float tmin = std::numeric_limits<float>::max();

  for (;;) {
#if PRRT
    printf("@(%d,%d,%d) t_nxt={%f,%f,%f} pt_nxt=(%f,%f,%f)\n",
      coord[0],coord[1],coord[2],t_nxt[0],t_nxt[1],t_nxt[2],
      ray(t_nxt[t_nxt[1]>=t_nxt[0]])[0],
      ray(t_nxt[t_nxt[1]>=t_nxt[0]])[1],
      ray(t_nxt[t_nxt[1]>=t_nxt[0]])[2]);
#endif

    for (int k = 0; k < 2; ++k) {
      float t;
      const int&y = coord[1], &x = coord[0];
      if (jiao_tri(impl, ray, y, x, y, x+1, y+1, x+1, &t)) {
        hit = true;
        if (t < tmin) {
          tmin = t;
          argtri[0] = COORD(y,x);
          argtri[1] = COORD(y,x+1);
          argtri[2] = COORD(y+1,x+1);
        }
      }
      if (jiao_tri(impl, ray, y, x, y+1, x+1, y+1, x, &t)) {
        hit = true;
        if (t < tmin) {
          tmin = t;
          argtri[0] = COORD(y,x);
          argtri[1] = COORD(y+1,x+1);
          argtri[2] = COORD(y+1,x);
        }
      }
    }
    if (hit) break;
    int cmp_res = ((t_nxt[0]>=t_nxt[1])<<2)
                | ((t_nxt[0]>=t_nxt[2])<<1)
                | (t_nxt[1]>=t_nxt[2]);
    static const int lookup_cmp[8] = { 0           // 0<1   0<2   1<2
                                     , 0           // 0<1   0<2   2<=1
                                     , 2147483647  // 0<1   2<=0  1<2
                                     , 2           // 0<1   2<=0  2<=1
                                     , 1           // 1<=0  0<2   1<2
                                     , 2147483647  // 1<=0  0<2   2<=1
                                     , 1           // 1<=0  2<=0  1<2
                                     , 2 };        // 1<=0  2<=0  2<=1
    int k = lookup_cmp[cmp_res];
    if (t_nxt[k] > ray.maxt) break;
    coord[k] += coord_step[k];
    if (coord[k] == boundary[k]) break;
    t_nxt[k] += t_step[k];
  }

#if PRRT
  printf(hit?"hit=#t,tmin=%f,pt=(%f,%f,%f),argtri=(%d,%d)\n":"hit=#f\n",tmin,
    ray(tmin)[0],ray(tmin)[1],ray(tmin)[2],
    argtri[0],argtri[1]);
  for (int y = 0; y < impl->ny-1; ++y) {
    for (int x = 0; x < impl->nx-1; ++x) {
      for (int k = 0; k < 2; ++k) {
        float t;
        if (jiao_tri(impl, ray, impl->tri[y][x*2+k], &t)) {
          if (t<tmin) {
            printf("\nQQ:jiao at (%d,%d), t=%.10f; tmin=%.10f; arg=(%d,%d)\n",
              y,x*2+k,t,tmin,argtri[0],argtri[1]);
            std::getchar();
          } else if (t<=tmin+1e-4 && (y!=argtri[0]||x*2+k!=argtri[1])) {
            printf("\nXD:jiao at (%d,%d), t=%.10f; tmin=%.10f; arg=(%d,%d)\n",
              y,x*2+k,t,tmin,argtri[0],argtri[1]);
          }
        }
      }
    }
  }
#endif

  if (!hit) return false;

  const Transform &real_o2w = *ObjectToWorld;
  const static Normal dndu {0,0,0}, dndv {0,0,0};
  const Point &pt0 = impl->pt[argtri[0]]
            , &pt1 = impl->pt[argtri[1]]
            , &pt2 = impl->pt[argtri[2]];
  const Vector ntri {Cross(pt1-pt0, pt2-pt0)};
  const float inv_z = 1.0/ntri.z;
  const Vector dpdu {1,0,-ntri.x*inv_z}, dpdv {0,1,-ntri.y*inv_z};
  const Point pt = ray(tmin);
  *dg = {real_o2w(pt), real_o2w(dpdu), real_o2w(dpdv),
                                       real_o2w(dndu), real_o2w(dndv),
                                       pt.x, pt.y, this};
  dg->dudx = 1;
  dg->dvdy = 1;
  *tHit = tmin;
  *rayEpsilon = 1e-3f * tmin;
  return true;
}


bool Heightfield2::IntersectP(const Ray &r) const {
  Ray ray;
  (*WorldToObject)(r, &ray);

  float t_init;
  const float d_inv[3] { 1.0f/ray.d.x, 1.0f/ray.d.y, 1.0f/ray.d.z };
  if (impl->bbox.Inside(ray(ray.mint))) {
    t_init = ray.mint;
  } else if (!impl->bbox.IntersectP(ray, &t_init)) {
    return false;
  }

  int coord[3] = { static_cast<int>((ray.o.x + t_init*ray.d.x)*(impl->nx-1) + EPS)
                 , static_cast<int>((ray.o.y + t_init*ray.d.y)*(impl->ny-1) + EPS)
                 , 0 };

  if (coord[0] == impl->nx-1) --coord[0];
  if (coord[1] == impl->ny-1) --coord[1];

  float t_nxt[3], t_step[3];
  int coord_step[3], boundary[3];

  // handle z
  if (ray.d[2] >= 0) {
    t_nxt[2] = (impl->maxz - ray.o.z)*d_inv[2];
    t_step[2] = (impl->maxz-impl->minz)*d_inv[2]; // unimportant
    coord_step[2] = 1;
    boundary[2] = 1;
  } else {
    t_nxt[2] = (impl->minz - ray.o.z)*d_inv[2];
    t_step[2] = -(impl->maxz-impl->minz)*d_inv[2]; // unimportant
    coord_step[2] = -1;
    boundary[2] = -1;
  }

  // handle y
  if (ray.d[1] >= 0) {
    t_nxt[1] = ((coord[1]+1)*impl->dy - ray.o.y)*d_inv[1];
    t_step[1] = impl->dy * d_inv[1];
    coord_step[1] = 1;
    boundary[1] = impl->ny - 1;
  } else {
    t_nxt[1] = (coord[1]*impl->dy - ray.o.y)*d_inv[1];
    t_step[1] = - impl->dy * d_inv[1];
    coord_step[1] = -1;
    boundary[1] = -1;
  }

  // handle x
  if (ray.d[0] >= 0) {
    t_nxt[0] = ((coord[0]+1)*impl->dx - ray.o.x)*d_inv[0];
    t_step[0] = impl->dx * d_inv[0];
    coord_step[0] = 1;
    boundary[0] = impl->nx - 1;
  } else {
    t_nxt[0] = (coord[0]*impl->dx - ray.o.x)*d_inv[0];
    t_step[0] = - impl->dx * d_inv[0];
    coord_step[0] = -1;
    boundary[0] = -1;
  }

  for (;;) {
    for (int k = 0; k < 2; ++k) {
      float t;
      const int&y = coord[1], &x = coord[0];
      if (jiao_tri(impl, ray, y, x, y, x+1, y+1, x+1, &t)
       || jiao_tri(impl, ray, y, x, y+1, x+1, y+1, x, &t))
        return true;
    }
    int cmp_res = ((t_nxt[0]>=t_nxt[1])<<2)
                | ((t_nxt[0]>=t_nxt[2])<<1)
                | (t_nxt[1]>=t_nxt[2]);
    static const int lookup_cmp[8] = { 0           // 0<1   0<2   1<2
                                     , 0           // 0<1   0<2   2<=1
                                     , 2147483647  // 0<1   2<=0  1<2
                                     , 2           // 0<1   2<=0  2<=1
                                     , 1           // 1<=0  0<2   1<2
                                     , 2147483647  // 1<=0  0<2   2<=1
                                     , 1           // 1<=0  2<=0  1<2
                                     , 2 };        // 1<=0  2<=0  2<=1
    int k = lookup_cmp[cmp_res];
    if (t_nxt[k] > ray.maxt) break;
    coord[k] += coord_step[k];
    if (coord[k] == boundary[k]) break;
    t_nxt[k] += t_step[k];
  }
  return false;
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

