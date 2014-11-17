
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
#if 0
      if (ifz(radius) && ifz(n))
        impl->R = aperture;
#else
      impl->R = aperture/2.0f; /* sample on the last lens */
#endif
    }
    fspec.close();
  }
  impl->lenses.back().axpos = filmdistance;
  impl->film_z = 0;
  for (const lens_t& lens : impl->lenses)
    impl->film_z -= lens.axpos;

  fprintf(stderr, "[+] Read %u lens(es).\n", impl->lenses.size());
  fprintf(stderr, "[ ] dist = %f, z = %f, diag = %f, R = %f\n", filmdistance, impl->film_z,impl->film_diag, impl->R);
}

static inline bool sphere_intersect(const float& rad, const float& zcenter, const Ray &r, float* t) {
  const float A = r.d.x*r.d.x + r.d.y*r.d.y + r.d.z*r.d.z;
  const float B = 2*(r.d.x*r.o.x + r.d.y*r.o.y + r.d.z*(r.o.z - zcenter));
  const float C = r.o.x*r.o.x + r.o.y*r.o.y + (r.o.z-zcenter)*(r.o.z-zcenter)
                  - rad*rad;
  float t0, t1;
  if (!Quadratic(A, B, C, &t0, &t1))
    return false;

  *t = rad >= 0? t1 : t0;
//  return *t >= 0;
  return true;
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
  const float cos_thetai = Dot(d, normal)/Dot(normal, normal);
  const float sin_thetai = sqrt(1 - cos_thetai*cos_thetai);
  const float sin_thetao = sin_thetai * n_prv / n;
  const float cos_thetao = sqrt(1 - sin_thetao*sin_thetao);
  const Vector dir_out {sin_thetao*cos_thetai/(sin_thetai*cos_thetao)*dir};
  *r_new = {pt, normal + dir_out, 0.0f, std::numeric_limits<float>::max()};
  return true;
}

float RealisticCamera::GenerateRay(const CameraSample &sample, Ray *ray) const {
  const float reso_diag {static_cast<float>(sqrt(
      static_cast<float>(impl->film->xResolution)*impl->film->xResolution
    + static_cast<float>(impl->film->yResolution)*impl->film->yResolution))};

  const float reso_scale {impl->film_diag / reso_diag};

  const Point pimg { (sample.imageX - impl->film->xResolution/2.0f)*reso_scale
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
  float z {impl->film_z}, n_prv {1.0f};
  const Vector d0 {Normalize(psample - pimg)};

//  fprintf(stderr, "generate ray: ");
  Ray r {pimg, d0, 0., std::numeric_limits<float>::max()};
  for (auto it = impl->lenses.rbegin(); it != impl->lenses.rend(); ++it) {
    z += it->axpos;
    float t;

    // TODO: fixme: handle aperture
    if (ifz(it->radius)) {
      n_prv = 1.0f;
      continue;
    }

    const float zcenter = z - it->radius;
    if (!sphere_intersect(it->radius, zcenter, r, &t)) {
//      fprintf(stderr, "No intersection\n");
      return 0.0;
    }

    const Point pt {r(t)};
    if (pt.z >= it->aperture/2.0f || pt.z <= -it->aperture/2.0f) {
//      fprintf(stderr, "Out of aperture\n");
      return 0.0;
    }

    Ray r_new;
    if (!refraction(r, pt, Normalize(pt - Point {0,0,zcenter}), n_prv, it->n, &r_new)) {
//      fprintf(stderr, "Impossible: refraction\n");
      return 0.0;
    }

    r = r_new;
    n_prv = it->n;
  }
//  fprintf(stderr, "\n");

  *ray = CameraToWorld(r);

  const float cos_theta = Dot(d0, Vector  {0,0,1});
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
