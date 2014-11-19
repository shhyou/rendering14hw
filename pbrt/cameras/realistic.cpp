
#include "stdafx.h"
#include "film.h"
#include "montecarlo.h"
#include "sampler.h"
#include "cameras/realistic.h"

#include <cassert>
#include <cmath>
#include <cstdio>
#include <fstream>
#include <string>
#include <vector>
#include <limits>
#include <stdexcept>
#include <algorithm>

#define VERBOSE_ 0

static constexpr float EPS = 1e-6, PI = acos(-1.f);
static constexpr bool ifz(float f) { return -EPS<f && f<EPS; }

struct lens_t {
  float radius;
  float axpos;
  float n;
  float aperture;
};

struct RealisticCamer_Impl {
  const Film * const film;
  const float film_diag, aper_diam;
  float R, film_z;
  vector<lens_t> lenses;
};

RealisticCamera::RealisticCamera(
  const AnimatedTransform &cam2world,
  float hither, float yon, 
  float sopen, float sclose, 
  float filmdistance, float aperture_diameter, string specfile, 
  float filmdiag, Film *f)
  : Camera(cam2world, sopen, sclose, f),
    impl(new RealisticCamer_Impl {f, filmdiag, aperture_diameter})
{
  {
    std::ifstream fspec(specfile);
    if (!fspec.good())
      throw std::runtime_error("Cannot open input file '" + specfile + "'");
    string s;
    while (getline(fspec, s)) {
      if (s[0] == '#') continue;
      float radius, axpos, n, aperture;
      sscanf(s.c_str(), "%f%f%f%f", &radius, &axpos, &n, &aperture);
      impl->lenses.push_back({radius, axpos, n, aperture});
      if (ifz(radius) && ifz(n) && aperture_diameter<aperture)
        aperture = aperture_diameter;
#if 0
      if (ifz(radius) && ifz(n))
        impl->R = aperture/2.0f;
#else
      impl->R = aperture/2.0f; /* sample on the last lens */
#endif
    }
    fspec.close();
  }

  assert(ifz(impl->lenses.back().axpos));
  impl->lenses.back().axpos = filmdistance;
  impl->film_z = 0;
  for (const lens_t& lens : impl->lenses)
    impl->film_z -= lens.axpos;

  fprintf(stderr, "[+] Read %u lens(es).\n", impl->lenses.size());
  fprintf(stderr, "[ ] dist = %f, z = %f, diag = %f, R = %f\n", filmdistance, impl->film_z,impl->film_diag, impl->R);
  fprintf(stderr, "[ ] resolution: %d * %d\n", impl->film->xResolution, impl->film->yResolution);
}

static inline bool sphere_intersect(const float& rad, const float& zcenter, const Ray &r, float* t) {
  const float A = r.d.x*r.d.x + r.d.y*r.d.y + r.d.z*r.d.z;
  const float B = 2*(r.d.x*r.o.x + r.d.y*r.o.y + r.d.z*(r.o.z - zcenter));
  const float C = r.o.x*r.o.x + r.o.y*r.o.y + (r.o.z-zcenter)*(r.o.z-zcenter)
                  - rad*rad;
  float ts[2];
  if (!Quadratic(A, B, C, &ts[0], &ts[1]))
    return false;

  if (rad >= 0)
    std::swap(ts[0], ts[1]);

  *t = ts[0];
  return r.mint<=ts[0] && ts[0]<=r.maxt;
}

static inline bool refraction(
  const Ray& r,
  const Point& pt,
  const Vector& _normal, /* normal vector of the plane */
  const float n_prv,
  const float n,
  Ray *r_new)
{
  assert(!ifz(n) && !ifz(n_prv));
  const Vector d {Normalize(r.d)};
  const Vector normal {Dot(d, _normal)/Dot(_normal, _normal)*_normal};
  const Vector dir {d - normal};
#if VERBOSE_
  fprintf(stderr, "        d=(%f,%f,%f) normal=(%f,%f,%f)\n",
    d.x,d.y,d.z,
    normal.x,normal.y,normal.z);
  fprintf(stderr, "        dir=(%f,%f,%f)\n", dir.x,dir.y,dir.z);
#endif
  const float cos_thetai = Dot(d, normal)/sqrt(Dot(normal, normal));
  const float sin_thetai = sqrt(1 - cos_thetai*cos_thetai);
  const float sin_thetao = sin_thetai * n_prv / n;
  if (sin_thetao > 1.f)
    return false;
  const float cos_thetao = sqrt(1 - sin_thetao*sin_thetao);
#if VERBOSE_
  fprintf(stderr, "        cos_i=%f,sin_i=%f,cos_o=%f,sin_o=%f\n",cos_thetai,sin_thetai,cos_thetao,sin_thetao);
#endif
  const Vector dir_out {n_prv*cos_thetai/(n*cos_thetao)*dir};
  *r_new = {pt, normal + dir_out, 0.0f, std::numeric_limits<float>::max()};
  return true;
}

float RealisticCamera::GenerateRay(const CameraSample &sample, Ray *ray) const {
  const float reso_diag {static_cast<float>(sqrt(
      static_cast<float>(impl->film->xResolution)*impl->film->xResolution
    + static_cast<float>(impl->film->yResolution)*impl->film->yResolution))};

  const float reso_scale {impl->film_diag / reso_diag};

  const Point pimg { -(sample.imageX - impl->film->xResolution/2.0f)*reso_scale
                   , (sample.imageY - impl->film->yResolution/2.0f)*reso_scale
                   , impl->film_z };

  float lensU, lensV;
#if 1
  UniformSampleDisk(sample.lensU, sample.lensV, &lensU, &lensV);
#else
  ConcentricSampleDisk(sample.lensU, sample.lensV, &lensU, &lensV);
#endif
  lensU *= impl->R;
  lensV *= impl->R;

  const Point psample {lensU, lensV, impl->film_z + impl->lenses.back().axpos};
  float z {impl->film_z}, n_cur, n_nxt;
  const Vector d0 {Normalize(psample - pimg)};

  Ray r {pimg, d0, 0., std::numeric_limits<float>::max()};

#if VERBOSE_
  fprintf(stderr, "resolution scale: %f; p0: (%f, %f, %f)\n", reso_scale, pimg.x, pimg.y, pimg.z);
  fprintf(stderr, "(u,v) = (%f,%f); p1 = (%f,%f,%f)\n", sample.lensU, sample.lensV, psample.x, psample.y, psample.z);
  fprintf(stderr, "ray => (%f, %f, %f)\n", r.d.x, r.d.y, r.d.z);
#endif

  for (auto it = impl->lenses.rbegin(); it != impl->lenses.rend(); ++it) {
    z += it->axpos;

    // TODO: fixme: handle aperture
    if (ifz(it->radius)) {
      /* (r.o + t*r.d)*n = z */
      const float t = (z - r.o.z)/r.d.z;
      const Point pt {r(t)};
      const float dist = pt.x*pt.x + pt.y*pt.y, apt = it->aperture*it->aperture/4.0f;
      if (dist > apt)
        return 0.0;
      continue;
    }

    const float zcenter = z - it->radius;

#if VERBOSE_
    fprintf(stderr, "    z = %f, zcenter = %f\n", z, zcenter);
#endif

    float t;
    if (!sphere_intersect(it->radius, zcenter, r, &t)) {
#if VERBOSE_
      fprintf(stderr, "    No intersection\n");
      getchar();
#endif
      return 0.0;
    }

    const Point pt {r(t)};

#if VERBOSE_
    fprintf(stderr, "    t = %f, (%f,%f,%f)\n", t, pt.x, pt.y, pt.z);
#endif

#if 1
    const float dist = pt.x*pt.x + pt.y*pt.y, apt = it->aperture*it->aperture/4.0f;
    if (dist > apt) {
#if VERBOSE_
      fprintf(stderr, "    Out of aperture\n");
      getchar();
#endif
      return 0.0;
    }
#endif

    Ray r_new;
    n_cur = it->n;
    n_nxt = it+1==impl->lenses.rend()? 1.0 : (it+1)->n;
    if (ifz(n_cur)) n_cur = 1;
    if (ifz(n_nxt)) n_nxt = 1;
    if (!refraction(r, pt, pt - Point {0,0,zcenter}, n_cur, n_nxt, &r_new)) {
#if VERBOSE_
      fprintf(stderr, "    Refraction\n");
      getchar();
#endif
      return 0.0;
    }

#if VERBOSE_
    fprintf(stderr, "    ray => (%f,%f,%f)\n", r_new.d.x, r_new.d.y, r_new.d.z);
#endif
    r = r_new;
    assert(Dot(r.d, Vector {0,0,1}) > 0);
  }

  *ray = CameraToWorld(r);

  const float cos_theta = Dot(d0, Vector {0,0,1});
#if VERBOSE_
  getchar();
#endif
  return PI * impl->R * impl->R * pow(cos_theta,4.f)
       / (impl->lenses.back().axpos * impl->lenses.back().axpos);
}


RealisticCamera *CreateRealisticCamera(const ParamSet &params,
        const AnimatedTransform &cam2world, Film *film) {
  // Extract common camera parameters from \use{ParamSet}
  float hither = params.FindOneFloat("hither", -1);
  float yon = params.FindOneFloat("yon", -1);
  float shutteropen = params.FindOneFloat("shutteropen", -1);
  float shutterclose = params.FindOneFloat("shutterclose", -1);

  // Realistic camera-specific parameters
  string specfile = params.FindOneString("specfile", "");
  float filmdistance = params.FindOneFloat("filmdistance", 70.0); // about 70 mm default to film
   float fstop = params.FindOneFloat("aperture_diameter", 1.0);  
  float filmdiag = params.FindOneFloat("filmdiag", 35.0);

  Assert(hither != -1 && yon != -1 && shutteropen != -1 &&
    shutterclose != -1 && filmdistance!= -1);
  if (specfile == "") {
      Severe( "No lens spec file supplied!\n" );
  }
  return new RealisticCamera(cam2world, hither, yon,
           shutteropen, shutterclose, filmdistance, fstop, 
           specfile, filmdiag, film);
}
