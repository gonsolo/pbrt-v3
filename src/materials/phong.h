#if defined(_MSC_VER)
#define NOMINMAX
#pragma once
#endif

#ifndef PBRT_MATERIALS_PHONG_H
#define PBRT_MATERIALS_PHONG_H

#include "material.h"
#include "pbrt.h"

namespace pbrt {

class PhongMaterial : public Material {
  public:
    PhongMaterial(const std::shared_ptr<Texture<Spectrum>> &color,
                   const std::shared_ptr<Texture<Float>> &phongExp)
        : color(color),
          phongExp(phongExp) {}

    void ComputeScatteringFunctions(SurfaceInteraction *si, MemoryArena &arena,
                                    TransportMode mode,
                                    bool allowMultipleLobes) const;

  private:
    std::shared_ptr<Texture<Spectrum>> color;
    std::shared_ptr<Texture<Float>> phongExp;
};

PhongMaterial *CreatePhongMaterial(const TextureParams &mp);

}  // namespace pbrt

#endif  // PBRT_MATERIALS_PHONG_H
