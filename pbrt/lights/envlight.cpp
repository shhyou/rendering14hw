
/*
  pbrt source code Copyright(c) 1998-2012 Matt Pharr and Greg Humphreys.

  This file is part of pbrt.

  Reimpl->distribution and use in source and binary forms, with or without
  modification, are permitted provided that the following conditions are
  met:

  - Reimpl->distributions of source code must retain the above copyright
    notice, this list of conditions and the following disclaimer.

  - Reimpl->distributions in binary form must reproduce the above copyright
    notice, this list of conditions and the following disclaimer in the
    documentation and/or other materials provided with the impl->distribution.

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


// lights/environment.cpp*
#include "stdafx.h"
#include "lights/envlight.h"
#include "lights/distant.h"
#include "sh.h"
#include "montecarlo.h"
#include "paramset.h"
#include "imageio.h"
#include "mipmap.h"

// MedianCutEnvironmentLight private implementation
struct MedCutEnvImpl {
  // MedianCutEnvironmentLight Private Data
  MIPMap<RGBSpectrum> *radianceMap;
  Distribution2D *distribution;
};

// MedianCutEnvironmentLight Utility Classes
struct MedianCutEnvironmentCube {
  // MedianCutEnvironmentCube Public Methods
  MedianCutEnvironmentCube(const MedianCutEnvironmentLight *l, const Scene *s,
             float t, bool cv, float pe)
    : light(l), scene(s), time(t), pEpsilon(pe), computeVis(cv) { }
  Spectrum operator()(int, int, const Point &p, const Vector &w) {
    Ray ray(p, w, pEpsilon, INFINITY, time);
    if (!computeVis || !scene->IntersectP(ray))
      return light->Le(RayDifferential(ray));
    return 0.f;
  }
  const MedianCutEnvironmentLight *light;
  const Scene *scene;
  float time, pEpsilon;
  bool computeVis;
};

// MedianCutEnvironmentLight Method Definitions
MedianCutEnvironmentLight::~MedianCutEnvironmentLight() {
  delete impl->distribution;
  delete impl->radianceMap;
}


MedianCutEnvironmentLight::MedianCutEnvironmentLight(const Transform &light2world,
    const Spectrum &L, int ns, const string &texmap)
  : Light(light2world, ns) {
  int width = 0, height = 0;
  RGBSpectrum *texels = NULL;
  // Read texel data from _texmap_ into _texels_
  if (texmap != "") {
    texels = ReadImage(texmap, &width, &height);
    if (texels)
      for (int i = 0; i < width * height; ++i)
        texels[i] *= L.ToRGBSpectrum();
  }
  if (!texels) {
    width = height = 1;
    texels = new RGBSpectrum[1];
    texels[0] = L.ToRGBSpectrum();
  }
  impl->radianceMap = new MIPMap<RGBSpectrum>(width, height, texels);
  delete[] texels;
  // Initialize sampling PDFs for environment area light

  // Compute scalar-valued image _img_ from environment map
  float filter = 1.f / max(width, height);
  float *img = new float[width*height];
  for (int v = 0; v < height; ++v) {
    float vp = (float)v / (float)height;
    float sinTheta = sinf(M_PI * float(v+.5f)/float(height));
    for (int u = 0; u < width; ++u) {
      float up = (float)u / (float)width;
      img[u+v*width] = impl->radianceMap->Lookup(up, vp, filter).y();
      img[u+v*width] *= sinTheta;
    }
  }

  // Compute sampling impl->distributions for rows and columns of image
  impl->distribution = new Distribution2D(img, width, height);
  delete[] img;
}


Spectrum MedianCutEnvironmentLight::Power(const Scene *scene) const {
  Point worldCenter;
  float worldRadius;
  scene->WorldBound().BoundingSphere(&worldCenter, &worldRadius);
  return M_PI * worldRadius * worldRadius *
    Spectrum(impl->radianceMap->Lookup(.5f, .5f, .5f), SPECTRUM_ILLUMINANT);
}


Spectrum MedianCutEnvironmentLight::Le(const RayDifferential &r) const {
  Vector wh = Normalize(WorldToLight(r.d));
  float s = SphericalPhi(wh) * INV_TWOPI;
  float t = SphericalTheta(wh) * INV_PI;
  return Spectrum(impl->radianceMap->Lookup(s, t), SPECTRUM_ILLUMINANT);
}


void MedianCutEnvironmentLight::SHProject(const Point &p, float pEpsilon,
    int lmax, const Scene *scene, bool computeLightVis,
    float time, RNG &rng, Spectrum *coeffs) const {
  // Project _MedianCutEnvironmentLight_ to SH using Monte Carlo if visibility needed
  if (computeLightVis) {
    Light::SHProject(p, pEpsilon, lmax, scene, computeLightVis,
             time, rng, coeffs);
    return;
  }
  for (int i = 0; i < SHTerms(lmax); ++i)
    coeffs[i] = 0.f;
  int ntheta = impl->radianceMap->Height(), nphi = impl->radianceMap->Width();
  if (min(ntheta, nphi) > 50) {
    // Project _MedianCutEnvironmentLight_ to SH from lat-long representation

    // Precompute $\theta$ and $\phi$ values for lat-long map projection
    float *buf = new float[2*ntheta + 2*nphi];
    float *bufp = buf;
    float *sintheta = bufp;  bufp += ntheta;
    float *costheta = bufp;  bufp += ntheta;
    float *sinphi = bufp;  bufp += nphi;
    float *cosphi = bufp;
    for (int theta = 0; theta < ntheta; ++theta) {
      sintheta[theta] = sinf((theta + .5f)/ntheta * M_PI);
      costheta[theta] = cosf((theta + .5f)/ntheta * M_PI);
    }
    for (int phi = 0; phi < nphi; ++phi) {
      sinphi[phi] = sinf((phi + .5f)/nphi * 2.f * M_PI);
      cosphi[phi] = cosf((phi + .5f)/nphi * 2.f * M_PI);
    }
    float *Ylm = ALLOCA(float, SHTerms(lmax));
    for (int theta = 0; theta < ntheta; ++theta) {
      for (int phi = 0; phi < nphi; ++phi) {
        // Add _MedianCutEnvironmentLight_ texel's contribution to SH coefficients
        Vector w = Vector(sintheta[theta] * cosphi[phi],
                  sintheta[theta] * sinphi[phi],
                  costheta[theta]);
        w = Normalize(LightToWorld(w));
        Spectrum Le = Spectrum(impl->radianceMap->Texel(0, phi, theta),
                     SPECTRUM_ILLUMINANT);
        SHEvaluate(w, lmax, Ylm);
        for (int i = 0; i < SHTerms(lmax); ++i)
          coeffs[i] += Le * Ylm[i] * sintheta[theta] *
            (M_PI / ntheta) * (2.f * M_PI / nphi);
      }
    }

    // Free memory used for lat-long theta and phi values
    delete[] buf;
  }
  else {
    // Project _MedianCutEnvironmentLight_ to SH from cube map sampling
    SHProjectCube(MedianCutEnvironmentCube(this, scene, time, computeLightVis,
                     pEpsilon),
            p, 200, lmax, coeffs);
  }
}


MedianCutEnvironmentLight *CreateMedianCutEnvironmentLight(const Transform &light2world,
    const ParamSet &paramSet) {
  Spectrum L = paramSet.FindOneSpectrum("L", Spectrum(1.0));
  Spectrum sc = paramSet.FindOneSpectrum("scale", Spectrum(1.0));
  string texmap = paramSet.FindOneFilename("mapname", "");
  int nSamples = paramSet.FindOneInt("nsamples", 1);
  if (PbrtOptions.quickRender) nSamples = max(1, nSamples / 4);
  return new MedianCutEnvironmentLight(light2world, L * sc, nSamples, texmap);
}


Spectrum MedianCutEnvironmentLight::Sample_L(const Point &p, float pEpsilon,
    const LightSample &ls, float time, Vector *wi, float *pdf,
    VisibilityTester *visibility) const {
  // Find $(u,v)$ sample coordinates in environment light texture
  float uv[2], mapPdf;
  impl->distribution->SampleContinuous(ls.uPos[0], ls.uPos[1], uv, &mapPdf);
  if (mapPdf == 0.f) return 0.f;

  // Convert environment light sample point to direction
  float theta = uv[1] * M_PI, phi = uv[0] * 2.f * M_PI;
  float costheta = cosf(theta), sintheta = sinf(theta);
  float sinphi = sinf(phi), cosphi = cosf(phi);
  *wi = LightToWorld(Vector(sintheta * cosphi, sintheta * sinphi,
                costheta));

  // Compute PDF for sampled environment light direction
  *pdf = mapPdf / (2.f * M_PI * M_PI * sintheta);
  if (sintheta == 0.f) *pdf = 0.f;

  // Return radiance value for environment light direction
  visibility->SetRay(p, pEpsilon, *wi, time);
  Spectrum Ls = Spectrum(impl->radianceMap->Lookup(uv[0], uv[1]),
               SPECTRUM_ILLUMINANT);
  return Ls;
}


float MedianCutEnvironmentLight::Pdf(const Point &, const Vector &w) const {
  Vector wi = WorldToLight(w);
  float theta = SphericalTheta(wi), phi = SphericalPhi(wi);
  float sintheta = sinf(theta);
  if (sintheta == 0.f) return 0.f;
  float p = impl->distribution->Pdf(phi * INV_TWOPI, theta * INV_PI) /
       (2.f * M_PI * M_PI * sintheta);
  return p;
}


Spectrum MedianCutEnvironmentLight::Sample_L(const Scene *scene,
    const LightSample &ls, float u1, float u2, float time,
    Ray *ray, Normal *Ns, float *pdf) const {
  // Compute direction for environment light sample ray

  // Find $(u,v)$ sample coordinates in environment light texture
  float uv[2], mapPdf;
  impl->distribution->SampleContinuous(ls.uPos[0], ls.uPos[1], uv, &mapPdf);
  if (mapPdf == 0.f) return Spectrum(0.f);

  float theta = uv[1] * M_PI, phi = uv[0] * 2.f * M_PI;
  float costheta = cosf(theta), sintheta = sinf(theta);
  float sinphi = sinf(phi), cosphi = cosf(phi);
  Vector d = -LightToWorld(Vector(sintheta * cosphi, sintheta * sinphi,
                  costheta));
  *Ns = (Normal)d;

  // Compute origin for environment light sample ray
  Point worldCenter;
  float worldRadius;
  scene->WorldBound().BoundingSphere(&worldCenter, &worldRadius);
  Vector v1, v2;
  CoordinateSystem(-d, &v1, &v2);
  float d1, d2;
  ConcentricSampleDisk(u1, u2, &d1, &d2);
  Point Pdisk = worldCenter + worldRadius * (d1 * v1 + d2 * v2);
  *ray = Ray(Pdisk + worldRadius * -d, d, 0., INFINITY, time);

  // Compute _MedianCutEnvironmentLight_ ray PDF
  float directionPdf = mapPdf / (2.f * M_PI * M_PI * sintheta);
  float areaPdf = 1.f / (M_PI * worldRadius * worldRadius);
  *pdf = directionPdf * areaPdf;
  if (sintheta == 0.f) *pdf = 0.f;
  Spectrum Ls = (impl->radianceMap->Lookup(uv[0], uv[1]), SPECTRUM_ILLUMINANT);
  return Ls;
}


