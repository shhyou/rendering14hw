
#include "stdafx.h"
#include "film.h"
#include "montecarlo.h"
#include "sampler.h"
#include "cameras/realistic.h"

#include <cmath>
#include <cstdio>
#include <fstream>
#include <string>
#include <vector>
#include <limits>
#include <stdexcept>

static constexpr float EPS = 1e-6;
static constexpr bool ifz(float f) { return -EPS<f && f<EPS; }

struct lens_t {
  float radius;
  float axpos;
  float n;
  float aperture;
};

struct RealisticCamer_Impl {
  const Film * const film;
  float film_dist, film_diag, aper_diam;
  float R;
  vector<lens_t> lenses;
};

RealisticCamera::RealisticCamera(
  const AnimatedTransform &cam2world,
  float hither, float yon, 
  float sopen, float sclose, 
  float filmdistance, float aperture_diameter, string specfile, 
  float filmdiag, Film *f)
  : Camera(cam2world, sopen, sclose, f),
    impl(new RealisticCamer_Impl {f, filmdistance, filmdiag, aperture_diameter})
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
      impl->lenses.push_back({radius, axpos,n, aperture});
#if 0
      if (ifz(radius) && ifz(n))
        impl->R = aperture;
#else
      impl->R = aperture/2.0f; /* sample on the last lens */
#endif
    }
    fspec.close();
    fprintf(stderr, "[+] Read %u lens(es).\n", impl->lenses.size());
  }
}

bool sphere_intersect(const float& rad, const float& zcenter, const Ray &r, float* t) {
  const float A = r.d.x*r.d.x + r.d.y*r.d.y + r.d.z*r.d.z;
  const float B = 2*(r.d.x*r.o.x + r.d.y*r.o.y + r.d.z*(r.o.z - zcenter));
  const float C = r.o.x*r.o.x + r.o.y*r.o.y + (r.o.z-zcenter)*(r.o.z-zcenter)
                  - rad*rad;
  float t0, t1;
  if (!Quadratic(A, B, C, &t0, &t1))
    return false;

  *t = rad >= 0? t1 : t0;
  return true;
}

float RealisticCamera::GenerateRay(const CameraSample &sample, Ray *ray) const {
  const float reso_diag {static_cast<float>(sqrt(
      static_cast<float>(impl->film->xResolution)*impl->film->xResolution
    + static_cast<float>(impl->film->yResolution)*impl->film->yResolution))};

  const float reso_scale {impl->film_diag / reso_diag};

  const Point pimg { (sample.imageX - impl->film->xResolution/2.0f)*reso_scale
                   , (sample.imageY - impl->film->yResolution/2.0f)*reso_scale
                   , -impl->film_dist };

  float lensU, lensV;
#if 1
  UniformSampleDisk(sample.lensU, sample.lensV, &lensU, &lensV);
#else
  ConcentricSampleDisk(sample.lensU, sample.lensV, &lensU, &lensV);
#endif
  lensU *= impl->R;
  lensV *= impl->R;

  const Point psample {lensU, lensV, 0.0};
  float z {0}; 

  Ray r {pimg, Normalize(psample - pimg), 0., std::numeric_limits<float>::max()};
  for (auto it = impl->lenses.rbegin(); it != impl->lenses.rend(); ++it) {
    z += it->axpos;
    float t;

    if (!sphere_intersect(it->radius, z - it->radius, r, &t))
      return 0.0;

    const Point pt = r(t);
  }

  *ray = CameraToWorld(r);
  return 1.0;
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
