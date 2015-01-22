
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

#define SHOW_SAMPLE 0
#define PAUSE_LOOKUP 0 // 1 for display, 2 for pause
#define RECON_FILTER 1

#include <vector>
#include <iterator>
#include <algorithm>

// renderers/reconrenderer.cpp*
#include "stdafx.h"
#include "renderers/reconrenderer.h"
#include "samplers/lowdiscrepancy.h"
#include "filters/gaussian.h"
#include "scene.h"
#include "film.h"
#include "volume.h"
#include "sampler.h"
#include "integrator.h"
#include "progressreporter.h"
#include "camera.h"
#include "intersection.h"

using std::vector;
using std::back_insert_iterator;
using std::back_inserter;

struct ReconSample_t {
  float distz, depth;
  float imageX, imageY, lensU, lensV;
  Spectrum L;
  float reproj_x(float u) const {
    return imageX + (u - lensU)*depth;
  }
  float reproj_y(float v) const {
    return imageY + (v - lensV)*depth;
  }
};

const int search_t_d {4};

struct search_t {
  int max_x, max_y;
  const float dist2;
  vector<vector<vector<ReconSample_t*>>> sampls;
  search_t(vector<ReconSample_t>& sampls_, int nsamp)
    : dist2(4.f/nsamp)
  {
    max_x = 0, max_y = 0;
    for (ReconSample_t& sampl : sampls_) {
      max_x = std::max(max_x, static_cast<int>(sampl.imageX));
      max_y = std::max(max_y, static_cast<int>(sampl.imageY));
    }
    ++max_y;
    ++max_x;
    sampls.resize(max_y, vector<vector<ReconSample_t*>>(max_x));
    for (ReconSample_t& sampl : sampls_) {
      int y {static_cast<int>(sampl.imageY)}
        , x {static_cast<int>(sampl.imageX)};
      sampls[y][x].push_back(&sampl);
    }
  }
  void lookup(float x, float y, float u, float v,
              vector<ReconSample_t>& res) const
  {
    int y0 {static_cast<int>(y)}, x0 {static_cast<int>(x)};
#if PAUSE_LOOKUP
printf("search (%f,%f) (%f,%f)\n", x, y, u, v);
#endif
    for (int dy = -search_t_d; dy <= search_t_d; ++dy) {
      for (int dx = -search_t_d; dx <= search_t_d; ++dx) {
        if (0<=y0+dy && y0+dy<max_y && 0<=x0+dx && x0+dx<max_x) {
          for (auto psampl : sampls[y0+dy][x0+dx]) {
            const float disty {y - psampl->reproj_y(v)}
                      , distx {x - psampl->reproj_x(u)};
            if (distx*distx + disty*disty <= dist2) {
#if PAUSE_LOOKUP
printf("    (%f,%f) (%f,%f)\n",
  psampl->imageX, psampl->imageY,
  psampl->lensU, psampl->lensV);
printf("    reproj (%f,%f) z=%f\n",
  psampl->reproj_x(u), psampl->reproj_y(v), psampl->distz);
#endif
              res.push_back(*psampl);
            }
          }
        }
      }
    }
#if PAUSE_LOOKUP > 1
getchar();
#elif PAUSE_LOOKUP
putchar('\n');
#endif
    std::sort(res.begin(), res.end(), [](const ReconSample_t& r1, const ReconSample_t& r2) {
      return r1.distz < r2.distz;
    });
  }
};

struct ReconRenderer_Impl {
  Sampler *sampler;
  Transform *WorldToCamera;
  vector<ReconSample_t> sampls;
  search_t *searcher;
};

// ReconRendererInit Declarations
class ReconRendererInit : public Task {
public:
  // ReconRendererInit Public Methods
  ReconRendererInit(const Scene *sc, Renderer *ren, Camera *c,
                    ProgressReporter &pr, Sampler *ms, Sample *sam,
                    int nsamp_, int tn, int tc,
                    Transform* WorldToCamera_,
                    back_insert_iterator<vector<ReconSample_t>>& samplit_)
    : reporter(pr), WorldToCamera(WorldToCamera_), samplit(samplit_)
  {
    scene = sc; renderer = ren; camera = c; mainSampler = ms;
    origSample = sam; nsamp = nsamp_; taskNum = tn; taskCount = tc;
  }
  ~ReconRendererInit() = default;
  void Run();
private:
  // ReconRendererInit Private Data
  ProgressReporter &reporter;
  Transform *WorldToCamera;
  back_insert_iterator<vector<ReconSample_t>>& samplit;
  const Scene *scene;
  const Renderer *renderer;
  Camera *camera;
  Sampler *mainSampler;
  Sample *origSample;
  int nsamp;
  int taskNum, taskCount;
};

// ReconRendererInit Definitions
void ReconRendererInit::Run() {
  PBRT_STARTED_RENDERTASK(taskNum);
  // Get sub-_Sampler_ for _ReconRendererTask_
  Sampler *sampler = mainSampler->GetSubSampler(taskNum, taskCount);
  if (!sampler) {
    reporter.Update();
    PBRT_FINISHED_RENDERTASK(taskNum);
    return;
  }

  // Declare local variables used for rendering loop
  MemoryArena arena;
  RNG rng(taskNum);

  // Allocate space for samples and intersections
  const int& maxSamples = sampler->MaximumSampleCount();
  Sample *samples = origSample->Duplicate(maxSamples);
  RayDifferential *rays = new RayDifferential[maxSamples];
  Spectrum *Ls_a = new Spectrum[maxSamples];
  Spectrum *Ts_a = new Spectrum[maxSamples];
  Intersection *isects_a = new Intersection[maxSamples];

  // Get samples from _Sampler_ and update image
  int sampleCount;
  while ((sampleCount = sampler->GetMoreSamples(samples, rng)) > 0) {
    // Generate camera rays and compute radiance along rays
    for (int i = 0; i < sampleCount; ++i) {
      // Find camera ray for _sample[i]_
      PBRT_STARTED_GENERATING_CAMERA_RAY(&samples[i]);
      float rayWeight = camera->GenerateRayDifferential(samples[i], &rays[i]);
      rays[i].ScaleDifferentials(1.f / sqrtf(sampler->samplesPerPixel));
      PBRT_FINISHED_GENERATING_CAMERA_RAY(&samples[i], &rays[i], rayWeight);

      // Evaluate radiance along camera ray
      PBRT_STARTED_CAMERA_RAY_INTEGRATION(&rays[i], &samples[i]);
      if (rayWeight > 0.f) {
        Ls_a[i] = rayWeight * renderer->Li(scene, rays[i], &samples[i], rng,
                                           arena, &isects_a[i], &Ts_a[i]);
      } else {
        Ls_a[i] = 0.f;
        Ts_a[i] = 1.f;
      }

      // Issue warning if unexpected radiance value returned
      if (Ls_a[i].HasNaNs() || Ls_a[i].y() < -1e-5 || std::isinf(Ls_a[i].y())) {
        Error("NaN/negative/infinite radiance value returned "
              "for image sample.  Setting to black.");
        Ls_a[i] = Spectrum(0.f);
      }
      PBRT_FINISHED_CAMERA_RAY_INTEGRATION(&rays[i], &samples[i], &Ls_a[i]);

      const Point p {(*WorldToCamera)(isects_a[i].dg.p)};
      *samplit++ = {p.z, 1.f/p.z, samples[i].imageX, samples[i].imageY,
                    samples[i].lensU, samples[i].lensV, Ls_a[i]};
    }

#if SHOW_SAMPLE
    // Report sample results to _Sampler_, add contributions to image
    if (sampler->ReportResults(samples, rays, Ls_a, isects_a, sampleCount)) {
      for (int i = 0; i < sampleCount; ++i) {
        PBRT_STARTED_ADDING_IMAGE_SAMPLE(&samples[i], &rays[i], &Ls_a[i], &Ts_a[i]);
        camera->film->AddSample(samples[i], Ls_a[i]);
        PBRT_FINISHED_ADDING_IMAGE_SAMPLE();
      }
    }
#endif

    // Free _MemoryArena_ memory from computing image sample values
    arena.FreeAll();
  }

  // Clean up after _ReconRendererTask_ is done with its image region
  camera->film->UpdateDisplay(sampler->xPixelStart,
                              sampler->yPixelStart, sampler->xPixelEnd+1, sampler->yPixelEnd+1);
  delete sampler;
  delete[] samples;
  delete[] rays;
  delete[] Ls_a;
  delete[] Ts_a;
  delete[] isects_a;
  reporter.Update();
  PBRT_FINISHED_RENDERTASK(taskNum);
}




// ReconRendererTask Declarations
class ReconRendererTask : public Task {
public:
  // ReconRendererTask Public Methods
  ReconRendererTask(const Scene *sc, Renderer *ren, Camera *c,
                    ProgressReporter &pr, Sampler *ms, Sample *sam,
                    int nsamp_, float delta_, int tn, int tc,
                    const search_t& searcher_)
    : reporter(pr), searcher(searcher_), delta(delta_)
  {
    scene = sc; renderer = ren; camera = c; mainSampler = ms;
    origSample = sam; nsamp = nsamp_; taskNum = tn; taskCount = tc;
  }
  bool SameSurface(const ReconSample_t& s1, const ReconSample_t& s2,
                   float u, float v);
  ~ReconRendererTask() = default;
  void Run();
private:
  // ReconRendererTask Private Data
  ProgressReporter &reporter;
  const search_t& searcher;
  float delta;
  const Scene *scene;
  const Renderer *renderer;
  Camera *camera;
  Sampler *mainSampler;
  Sample *origSample;
  int nsamp;
  int taskNum, taskCount;

  // construct surfaces at a given location
  void ConstructSurfaces(float x, float y, float u, float v,
                         vector<vector<ReconSample_t>>& surfaces);
};

bool ReconRendererTask::SameSurface(
  const ReconSample_t& s1,
  const ReconSample_t& s2,
  float u,
  float v)
{
  const float x1[2] = {
    s1.reproj_x(u + delta),
    s1.reproj_x(u - delta),
  };
  const float x2[2] = {
    s2.reproj_x(u + delta),
    s2.reproj_x(u - delta),
  };
  const float y1[2] = {
    s1.reproj_y(v + delta),
    s1.reproj_y(v - delta),
  };
  const float y2[2] = {
    s2.reproj_y(v + delta),
    s2.reproj_y(v - delta),
  };
  const unsigned int xres =
       (x1[0]-x2[0] > 0)
    | ((x1[1]-x2[1] > 0) << 1);
  const unsigned int yres =
       (y1[0]-y2[0] > 0)
    | ((y1[1]-y2[1] > 0) << 1);
  static const bool consistent[4] = {
    true,   // <=   <=
    false,  // <=    >
    false,  //  >   <=
    true    //  >    >
  };
  return consistent[xres] && consistent[yres];
}

static inline float Cross(float x0, float y0, float x1, float y1) {
  return x0*y1 - x1*y0;
}

void ReconRendererTask::ConstructSurfaces(float x, float y, float u, float v,
                                          vector<vector<ReconSample_t>>& surfaces)
{
  vector<ReconSample_t> N;
  searcher.lookup(x, y, u, v, N); /* output is sorted */
  vector<vector<ReconSample_t>>().swap(surfaces);

  if (N.empty())
    return;

  auto it = begin(N);
  do {
    surfaces.push_back({});
    vector<ReconSample_t>& S = surfaces.back();
    while (it!=end(N) && (S.size()<3 || this->SameSurface(S.back(), *it, u, v))) {
      S.push_back(*it);
      ++it;
    }
  } while (it != end(N));
}

// ReconRendererTask Definitions
void ReconRendererTask::Run() {
  PBRT_STARTED_RENDERTASK(taskNum);
  // Get sub-_Sampler_ for _ReconRendererTask_
  Sampler *sampler = mainSampler->GetSubSampler(taskNum, taskCount);
  if (!sampler) {
    reporter.Update();
    PBRT_FINISHED_RENDERTASK(taskNum);
    return;
  }

  // Declare local variables used for rendering loop
  RNG rng(taskNum);

  // Allocate space for samples and intersections
  const int& maxSamples = sampler->MaximumSampleCount();
  Sample *samples = origSample->Duplicate(maxSamples);
  RayDifferential *rays = new RayDifferential[maxSamples];
  Spectrum *Ls_a = new Spectrum[maxSamples];
  Spectrum *Ts_a = new Spectrum[maxSamples];
  Intersection *isects_a = new Intersection[maxSamples];

  const GaussianFilter filt(2.f/sqrt(nsamp), 2.f/sqrt(nsamp), 100.f);

  // Get samples from _Sampler_ and update image
  int sampleCount;
  while ((sampleCount = sampler->GetMoreSamples(samples, rng)) > 0) {
    // Generate camera rays and compute radiance along rays
    for (int i = 0; i < sampleCount; ++i) {
      // Find camera ray for _sample[i]_
      PBRT_STARTED_GENERATING_CAMERA_RAY(&samples[i]);
      float rayWeight = camera->GenerateRayDifferential(samples[i], &rays[i]);
      rays[i].ScaleDifferentials(1.f / sqrtf(sampler->samplesPerPixel));
      PBRT_FINISHED_GENERATING_CAMERA_RAY(&samples[i], &rays[i], rayWeight);

      // Evaluate radiance along camera ray
      PBRT_STARTED_CAMERA_RAY_INTEGRATION(&rays[i], &samples[i]);
      if (rayWeight > 0.f) {
        vector<vector<ReconSample_t>> surfaces;
        this->ConstructSurfaces(samples[i].imageX, samples[i].imageY,
                                samples[i].lensU,  samples[i].lensV,
                                surfaces);
        auto it = begin(surfaces);
#if 0
printf("\n(%f,%f) (%f,%f): %u surfaces\n",
  samples[i].imageX, samples[i].imageY,
  samples[i].lensU, samples[i].lensV,
  surfaces.size());
#endif
        if (surfaces.empty())
          continue;
        bool found = false;
        for (; it != end(surfaces); ++it) {
          const vector<ReconSample_t>& S = *it;
#if 0
for (const ReconSample_t& s : S) {
  printf("    sample (%f,%f) (%f,%f)\n", s.imageX, s.imageY, s.lensU, s.lensV);
}
printf("\n");
#endif
          for (size_t a = 2; !found && a < S.size(); ++a) {
            const float x0 = S[a].reproj_x(samples[i].lensU)
                      , y0 = S[a].reproj_y(samples[i].lensV);
            for (size_t b = 1; !found && b < a; ++b) {
              const float x1 = S[b].reproj_x(samples[i].lensU)
                        , y1 = S[b].reproj_y(samples[i].lensV);
              for (size_t c = 0; !found && c < b; ++c) {
                const float x2 = S[c].reproj_x(samples[i].lensU)
                          , y2 = S[c].reproj_y(samples[i].lensV);
                const int sgn0 = Cross(x1-x0, y1-y0, samples[i].imageX-x0, samples[i].imageY-y0) > 0
                        , sgn1 = Cross(x2-x1, y2-y1, samples[i].imageX-x1, samples[i].imageY-y1) > 0
                        , sgn2 = Cross(x0-x2, y0-y2, samples[i].imageX-x2, samples[i].imageY-y2) > 0;
                if (sgn0==sgn1 && sgn1==sgn2) {
                  found = true;
                  break;
                }
              }
              if (found) break;
            }
            if (found) break;
          }
          if (found) break;
        }
//getchar();
        const vector<ReconSample_t>& S = it!=end(surfaces)? *it : surfaces.front();
        Ls_a[i] = 0.f; // XXX TODO FIXME: add filter
        Ts_a[i] = 0.f;
        for (const ReconSample_t& s : S) {
#if RECON_FILTER
          Ls_a[i] += s.L * filt.Evaluate(s.reproj_x(samples[i].lensU) - samples[i].imageX,
                                         s.reproj_y(samples[i].lensV) - samples[i].imageY);
#else
          Ls_a[i] += s.L;
#endif
        }
        Ls_a[i] *= rayWeight;
      } else {
        Ls_a[i] = 0.f;
        Ts_a[i] = 1.f;
      }

      // Issue warning if unexpected radiance value returned
      if (Ls_a[i].HasNaNs() || Ls_a[i].y() < -1e-5 || std::isinf(Ls_a[i].y())) {
        Error("NaN/negative/infinite recons radiance value' returned "
              "for image sample.  Setting to black.");
        Ls_a[i] = Spectrum(0.f);
      }
      PBRT_FINISHED_CAMERA_RAY_INTEGRATION(&rays[i], &samples[i], &Ls_a[i]);
    }

    // Report sample results to _Sampler_, add contributions to image
    if (sampler->ReportResults(samples, rays, Ls_a, isects_a, sampleCount)) {
      for (int i = 0; i < sampleCount; ++i) {
        PBRT_STARTED_ADDING_IMAGE_SAMPLE(&samples[i], &rays[i], &Ls_a[i], &Ts_a[i]);
        camera->film->AddSample(samples[i], Ls_a[i]);
        PBRT_FINISHED_ADDING_IMAGE_SAMPLE();
      }
    }
  }

  // Clean up after _ReconRendererTask_ is done with its image region
  camera->film->UpdateDisplay(sampler->xPixelStart,
                              sampler->yPixelStart, sampler->xPixelEnd+1, sampler->yPixelEnd+1);
  delete sampler;
  delete[] samples;
  delete[] rays;
  delete[] Ls_a;
  delete[] Ts_a;
  delete[] isects_a;
  reporter.Update();
  PBRT_FINISHED_RENDERTASK(taskNum);
}




// ReconRenderer Method Definitions
ReconRenderer::ReconRenderer(Sampler *s, Camera *c,
                             SurfaceIntegrator *si, VolumeIntegrator *vi,
                             int nsamp_, Transform* WorldToCamera_)
  : impl(new ReconRenderer_Impl)
{
  sampler = s;
  camera = c;
  impl->WorldToCamera = WorldToCamera_;
  surfaceIntegrator = si;
  volumeIntegrator = vi;
  nsamp = nsamp_;
  {
    ParamSet param;
    param.AddInt("pixelsamples", &nsamp_, 1);
    impl->sampler = CreateLowDiscrepancySampler(param, c->film, c);
  }
  impl->searcher = nullptr;
}


ReconRenderer::~ReconRenderer() {
  delete impl->searcher;
  delete impl->sampler;
  delete impl;
  delete sampler;
  delete camera;
  delete surfaceIntegrator;
  delete volumeIntegrator;
}


void ReconRenderer::Render(const Scene *scene) {
  PBRT_FINISHED_PARSING();
  // Allow integrators to do preprocessing for the scene
  PBRT_STARTED_PREPROCESSING();
  surfaceIntegrator->Preprocess(scene, camera, this);
  volumeIntegrator->Preprocess(scene, camera, this);
  PBRT_FINISHED_PREPROCESSING();
  PBRT_STARTED_RENDERING();

  // Create and launch _ReconRendererTask_s for rendering image

  ProgressReporter reporter(1+RoundUpPow2(camera->film->xResolution * camera->film->yResolution / (16*16)), "Rendering"); // XXX TODO FIXME
  back_insert_iterator<vector<ReconSample_t>> samplit = back_inserter(impl->sampls);
  {
    // No parallelization exists; we do the whole-image processing
    const int nTasks = RoundUpPow2(1);

    // Allocate and initialize _sample_
    Sample *sample = new Sample(impl->sampler, surfaceIntegrator,
                                volumeIntegrator, scene);
    vector<Task *> renderTasks;
    for (int i = 0; i < nTasks; ++i) {
      renderTasks.push_back(new ReconRendererInit(scene, this, camera,
                                                  reporter, impl->sampler, sample,
                                                  nsamp,
                                                  nTasks-1-i, nTasks,
                                                  impl->WorldToCamera,
                                                  samplit));
    }
    EnqueueTasks(renderTasks);
    WaitForAllTasks();
    for (uint32_t i = 0; i < renderTasks.size(); ++i)
      delete renderTasks[i];
    delete sample;
  }
  impl->searcher = new search_t(impl->sampls, nsamp);
  {
    // Compute number of _SamplerRendererTask_s to create for rendering
    const int nPixels = camera->film->xResolution * camera->film->yResolution;
    const int nTasks = RoundUpPow2(max(32 * NumSystemCores(), nPixels / (16*16)));
    float delta = 1.0f/sqrt(nsamp);

    // Allocate and initialize _sample_
    Sample *sample = new Sample(sampler, surfaceIntegrator,
                                volumeIntegrator, scene);
    vector<Task *> renderTasks;
    for (int i = 0; i < nTasks; ++i) {
      renderTasks.push_back(new ReconRendererTask(scene, this, camera,
                                                  reporter, sampler, sample, 
                                                  nsamp, delta,
                                                  nTasks-1-i, nTasks,
                                                  *impl->searcher));
    }
#if !SHOW_SAMPLE
    EnqueueTasks(renderTasks);
#endif
    WaitForAllTasks();
    for (uint32_t i = 0; i < renderTasks.size(); ++i)
      delete renderTasks[i];
    delete sample;
  }
  reporter.Done();
  PBRT_FINISHED_RENDERING();
  // Clean up after rendering and store final image
  camera->film->WriteImage();
}


Spectrum ReconRenderer::Li(const Scene *scene,
                           const RayDifferential &ray, const Sample *sample, RNG &rng,
                           MemoryArena &arena, Intersection *isect, Spectrum *T) const
{
  Assert(ray.time == sample->time);
  Assert(!ray.HasNaNs());
  // Allocate local variables for _isect_ and _T_ if needed
  Spectrum localT;
  if (!T) T = &localT;
  Intersection localIsect;
  if (!isect) isect = &localIsect;
  Spectrum Li = 0.f;
  if (scene->Intersect(ray, isect)) {
    Li = surfaceIntegrator->Li(scene, this, ray, *isect, sample,
                   rng, arena);
  } else {
    // Handle ray that doesn't intersect any geometry
    for (uint32_t i = 0; i < scene->lights.size(); ++i)
       Li += scene->lights[i]->Le(ray);
  }
  Spectrum Lvi = volumeIntegrator->Li(scene, this, ray, sample, rng,
                    T, arena);
  return *T * Li + Lvi;
}


Spectrum ReconRenderer::Transmittance(const Scene *scene,
                                      const RayDifferential &ray, const Sample *sample, RNG &rng,
                                      MemoryArena &arena) const
{
  return volumeIntegrator->Transmittance(scene, this, ray, sample,
                                         rng, arena);
}
