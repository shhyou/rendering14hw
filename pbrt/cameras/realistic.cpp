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

static constexpr float EPS = 1e-6;
static constexpr bool ifz(float f) { return -EPS<f && f<EPS; }

struct lens_t {
  float radius;
  float axpos;
  float n_ratio;
  float R2; /* (aperture/2)^2, square of the length to the boundary */
};

struct RealisticCamer_Impl {
  const Film * const film;
  const float raster_scale;
  float R, film_z, sample_area_over_Z2;
  vector<lens_t> lenses;
};

RealisticCamera::RealisticCamera(
  const AnimatedTransform &cam2world,
  float hither, float yon, 
  float sopen, float sclose, 
  float filmdistance, float aperture_diameter, string specfile, 
  float filmdiag, Film *f)
  : Camera(cam2world, sopen, sclose, f),
    impl(new RealisticCamer_Impl
          { f
          , static_cast<float>(filmdiag
                               /
                               sqrt(static_cast<float>(f->xResolution)*f->xResolution
                                  + static_cast<float>(f->yResolution)*f->yResolution))})
{
  {
    std::ifstream fspec(specfile);
    if (!fspec.good())
      throw std::runtime_error("Cannot open input file '" + specfile + "'");
    string s;
    vector<float> ns;
    ns.push_back(1.f);
    while (getline(fspec, s)) {
      if (s[0] == '#') continue;
      float radius, axpos, n, aperture;
      sscanf(s.c_str(), "%f%f%f%f", &radius, &axpos, &n, &aperture);
      if (ifz(radius) && ifz(n) && aperture_diameter<aperture) /* is aperture */
        aperture = aperture_diameter;
      impl->lenses.push_back({radius, axpos, 0., (aperture/2.f)*(aperture/2.f)});
      ns.push_back(ifz(n)? 1.f : n);
#if 0
      if (ifz(radius) && ifz(n))  /* sample on aperture */
        impl->R = aperture/2.0f;
#else
      impl->R = aperture/2.0f;    /* sample on the last lens */
#endif
    }
    fspec.close();

    for (size_t i = 1; i != ns.size(); ++i)
      impl->lenses[i-1].n_ratio = ns[i]/ns[i-1];
  }

  impl->sample_area_over_Z2 = acos(-1.f)*impl->R*impl->R / (filmdistance*filmdistance);

  assert(ifz(impl->lenses.back().axpos));
  impl->lenses.back().axpos = filmdistance;
  impl->film_z = 0;
  for (const lens_t& lens : impl->lenses)
    impl->film_z -= lens.axpos;

  fprintf(stderr, "[+] Read %u lens(es).\n", impl->lenses.size());
  fprintf(stderr, "[ ] dist = %f, z = %f, raster_scale = %f, R = %f\n", filmdistance, impl->film_z, impl->raster_scale, impl->R);
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
  const float n_ratio,
  Ray *r_new)
{
  const Vector d {Normalize(r.d)};
  const Vector n {Dot(d, _normal)/Dot(_normal, _normal)*_normal};
  const Vector dir {d - n};
  const float dot = Dot(d, n);
  const float cos_thetai2 = dot*dot/Dot(n, n);
  const float sin_thetai2 = 1 - cos_thetai2;
  const float sin_thetao2 = sin_thetai2 * n_ratio * n_ratio;
  if (sin_thetao2 > 1.f)
    return false;
  const float cos_thetao2 = 1 - sin_thetao2;
  const Vector dir_out {n_ratio*sqrt(cos_thetai2/cos_thetao2)*dir};
  *r_new = {pt, n + dir_out, 0.0f, std::numeric_limits<float>::max()};
  return true;
}

float RealisticCamera::GenerateRay(const CameraSample &sample, Ray *ray) const {
  const Point pimg { -(sample.imageX - impl->film->xResolution/2.0f)*impl->raster_scale
                   ,  (sample.imageY - impl->film->yResolution/2.0f)*impl->raster_scale
                   , impl->film_z };

  float lensU, lensV;
#if 0
  UniformSampleDisk(sample.lensU, sample.lensV, &lensU, &lensV);
#else
  ConcentricSampleDisk(sample.lensU, sample.lensV, &lensU, &lensV);
#endif
  lensU *= impl->R;
  lensV *= impl->R;

  const Point psample {lensU, lensV, impl->film_z + impl->lenses.back().axpos};
  float z {impl->film_z};
  const Vector d0 {Normalize(psample - pimg)};

  Ray r {pimg, d0, 0., std::numeric_limits<float>::max()};
  Vector normal;
  Point pt;

  for (auto it = impl->lenses.rbegin(); it != impl->lenses.rend(); ++it) {
    z += it->axpos;

    if (ifz(it->radius)) {
      /* (r.o + t*r.d)*n = z */
      const float t = (z - r.o.z)/r.d.z;
      pt = {r(t)};
      const float dist = pt.x*pt.x + pt.y*pt.y;
      if (dist > it->R2)
        return 0.0;
      normal = {0,0,1};
    } else {
      const float zcenter = z - it->radius;

      float t;
      if (!sphere_intersect(it->radius, zcenter, r, &t))
        return 0.0;

      pt = {r(t)};
      const float dist = pt.x*pt.x + pt.y*pt.y;
      if (dist > it->R2)
        return 0.0;
      normal = pt - Point {0,0,zcenter};
    }
    Ray r_new;
    if (!refraction(r, pt, normal, it->n_ratio, &r_new))
      return 0.0;
    r = r_new;
  }

  *ray = CameraToWorld(r);

  const float cos_theta = Dot(d0, Vector {0,0,1});
  return impl->sample_area_over_Z2 * pow(cos_theta,4.f);
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
