
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

#include <vector>
#include <iterator>

// renderers/reconrenderer.cpp*
#include "stdafx.h"
#include "renderers/reconrenderer.h"
#include "samplers/lowdiscrepancy.h"
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
  Sample sampl;
  Spectrum L, T;
  Intersection isect;
};

struct ReconRenderer_Impl {
  Sampler *sampler;
  vector<ReconSample_t> sampls;
};

// ReconRendererInit Declarations
class ReconRendererInit : public Task {
public:
  // ReconRendererInit Public Methods
  ReconRendererInit(const Scene *sc, Renderer *ren, Camera *c,
                    ProgressReporter &pr, Sampler *ms, Sample *sam,
                    int nsamp_, int tn, int tc,
                    back_insert_iterator<vector<ReconSample_t>>& samplit_)
    : reporter(pr), samplit(samplit_)
  {
    scene = sc; renderer = ren; camera = c; mainSampler = ms;
    origSample = sam; nsamp = nsamp_; taskNum = tn; taskCount = tc;
  }
  ~ReconRendererInit() = default;
  void Run();
private:
  // ReconRendererInit Private Data
  ProgressReporter &reporter;
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

      *samplit++ = {samples[i], Ls_a[i], Ts_a[i], isects_a[i]};
    }

#if 0
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
                    int nsamp_, int tn, int tc,
                    const vector<ReconSample_t>& sampls_)
    : reporter(pr), sampls(sampls_)
  {
    scene = sc; renderer = ren; camera = c; mainSampler = ms;
    origSample = sam; nsamp = nsamp_; taskNum = tn; taskCount = tc;
  }
  ~ReconRendererTask() = default;
  void Run();
private:
  // ReconRendererTask Private Data
  ProgressReporter &reporter;
  const vector<ReconSample_t>& sampls;
  const Scene *scene;
  const Renderer *renderer;
  Camera *camera;
  Sampler *mainSampler;
  Sample *origSample;
  int nsamp;
  int taskNum, taskCount;
};

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
//        Ls_a[i] = rayWeight * renderer->Li(scene, rays[i], &samples[i], rng,
//                                         arena, &isects_a[i], &Ts_a[i]);
        Ls_a[i] = -1.f; // XXX TODO FIXME
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
                             int nsamp_)
  : impl(new ReconRenderer_Impl)
{
  sampler = s;
  camera = c;
  surfaceIntegrator = si;
  volumeIntegrator = vi;
  nsamp = nsamp_;
  {
    ParamSet param;
    param.AddInt("pixelsamples", &nsamp_, 1);
    impl->sampler = CreateLowDiscrepancySampler(param, c->film, c);
  }
}


ReconRenderer::~ReconRenderer() {
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

  ProgressReporter reporter(1*2, "Rendering"); // XXX TODO FIXME
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
                                                  samplit));
    }
    EnqueueTasks(renderTasks);
    WaitForAllTasks();
    for (uint32_t i = 0; i < renderTasks.size(); ++i)
      delete renderTasks[i];
    delete sample;
  }
  {
    // Compute number of _SamplerRendererTask_s to create for rendering
    const int nPixels = camera->film->xResolution * camera->film->yResolution;
    const int nTasks = RoundUpPow2(max(32 * NumSystemCores(), nPixels / (16*16)));

    // Allocate and initialize _sample_
    Sample *sample = new Sample(sampler, surfaceIntegrator,
                                volumeIntegrator, scene);
    vector<Task *> renderTasks;
    for (int i = 0; i < nTasks; ++i) {
      renderTasks.push_back(new ReconRendererTask(scene, this, camera,
                                                  reporter, sampler, sample, 
                                                  nsamp,
                                                  nTasks-1-i, nTasks,
                                                  impl->sampls));
    }
    EnqueueTasks(renderTasks);
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
