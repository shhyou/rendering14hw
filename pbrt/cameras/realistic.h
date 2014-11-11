#ifndef PBRT_CAMERAS_REALISTIC_H
#define PBRT_CAMERAS_REALISTIC_H

#include "camera.h"
#include "paramset.h"
#include "film.h"

struct RealisticCamer_Impl;

// RealisticCamera Declarations
class RealisticCamera : public Camera {
public:
  // RealisticCamera Public Methods
  RealisticCamera(
    const AnimatedTransform &cam2world,
    float hither, float yon, float sopen,
    float sclose, float filmdistance, float aperture_diameter, string specfile,
    float filmdiag, Film *film);
  float GenerateRay(const CameraSample &sample, Ray *) const;

private:
  RealisticCamer_Impl * const impl;
};


RealisticCamera *CreateRealisticCamera(
  const ParamSet &params,
  const AnimatedTransform &cam2world, Film *film);


#endif  // PBRT_CAMERAS_REALISTIC_H
