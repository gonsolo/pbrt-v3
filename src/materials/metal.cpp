
/*
    pbrt source code is Copyright(c) 1998-2016
                        Matt Pharr, Greg Humphreys, and Wenzel Jakob.

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


// materials/metal.cpp*
#include "materials/metal.h"
#include "reflection.h"
#include "paramset.h"
#include "texture.h"
#include "interaction.h"

namespace pbrt {

// MetalMaterial Method Definitions
MetalMaterial::MetalMaterial(const std::shared_ptr<Texture<Spectrum>> &eta,
                             const std::shared_ptr<Texture<Spectrum>> &k,
                             const std::shared_ptr<Texture<Float>> &roughness,
                             const std::shared_ptr<Texture<Float>> &uRoughness,
                             const std::shared_ptr<Texture<Float>> &vRoughness,
                             const std::shared_ptr<Texture<Float>> &bumpMap,
                             bool remapRoughness)
    : eta(eta),
      k(k),
      roughness(roughness),
      uRoughness(uRoughness),
      vRoughness(vRoughness),
      bumpMap(bumpMap),
      remapRoughness(remapRoughness) {}

void MetalMaterial::ComputeScatteringFunctions(SurfaceInteraction *si,
                                               MemoryArena &arena,
                                               TransportMode mode,
                                               bool allowMultipleLobes) const {
    // Perform bump mapping with _bumpMap_, if present
    if (bumpMap) Bump(bumpMap, si);
    si->bsdf = ARENA_ALLOC(arena, BSDF)(*si);

    Float uRough =
        uRoughness ? uRoughness->Evaluate(*si) : roughness->Evaluate(*si);
    Float vRough =
        vRoughness ? vRoughness->Evaluate(*si) : roughness->Evaluate(*si);
    if (remapRoughness) {
        uRough = TrowbridgeReitzDistribution::RoughnessToAlpha(uRough);
        vRough = TrowbridgeReitzDistribution::RoughnessToAlpha(vRough);
    }
    Fresnel *frMf = ARENA_ALLOC(arena, FresnelConductor)(1., eta->Evaluate(*si),
                                                         k->Evaluate(*si));
    MicrofacetDistribution *distrib =
        ARENA_ALLOC(arena, TrowbridgeReitzDistribution)(uRough, vRough);
    si->bsdf->Add(ARENA_ALLOC(arena, MicrofacetReflection)(1., distrib, frMf));
}

const int CopperSamples = 56;
const Float CopperWavelengths[CopperSamples] = {
    298.7570554f, 302.4004341f, 306.1337728f, 309.960445f,  313.8839949f,
    317.9081487f, 322.036826f,  326.2741526f, 330.6244747f, 335.092373f,
    339.6826795f, 344.4004944f, 349.2512056f, 354.2405086f, 359.374429f,
    364.6593471f, 370.1020239f, 375.7096303f, 381.4897785f, 387.4505563f,
    393.6005651f, 399.9489613f, 406.5055016f, 413.2805933f, 420.2853492f,
    427.5316483f, 435.0322035f, 442.8006357f, 450.8515564f, 459.2006593f,
    467.8648226f, 476.8622231f, 486.2124627f, 495.936712f,  506.0578694f,
    516.6007417f, 527.5922468f, 539.0616435f, 551.0407911f, 563.5644455f,
    576.6705953f, 590.4008476f, 604.8008683f, 619.92089f,   635.8162974f,
    652.5483053f, 670.1847459f, 688.8009889f, 708.4810171f, 729.3186941f,
    751.4192606f, 774.9011125f, 799.8979226f, 826.5611867f, 855.0632966f,
    885.6012714f};

const Float CopperN[CopperSamples] = {
    1.400313f, 1.38f,  1.358438f, 1.34f,  1.329063f, 1.325f, 1.3325f,   1.34f,
    1.334375f, 1.325f, 1.317812f, 1.31f,  1.300313f, 1.29f,  1.281563f, 1.27f,
    1.249062f, 1.225f, 1.2f,      1.18f,  1.174375f, 1.175f, 1.1775f,   1.18f,
    1.178125f, 1.175f, 1.172812f, 1.17f,  1.165312f, 1.16f,  1.155312f, 1.15f,
    1.142812f, 1.135f, 1.131562f, 1.12f,  1.092437f, 1.04f,  0.950375f, 0.826f,
    0.645875f, 0.468f, 0.35125f,  0.272f, 0.230813f, 0.214f, 0.20925f,  0.213f,
    0.21625f,  0.223f, 0.2365f,   0.25f,  0.254188f, 0.26f,  0.28f,     0.3f};

const Float CopperK[CopperSamples] = {
    1.662125f, 1.687f, 1.703313f, 1.72f,  1.744563f, 1.77f,  1.791625f, 1.81f,
    1.822125f, 1.834f, 1.85175f,  1.872f, 1.89425f,  1.916f, 1.931688f, 1.95f,
    1.972438f, 2.015f, 2.121562f, 2.21f,  2.177188f, 2.13f,  2.160063f, 2.21f,
    2.249938f, 2.289f, 2.326f,    2.362f, 2.397625f, 2.433f, 2.469187f, 2.504f,
    2.535875f, 2.564f, 2.589625f, 2.605f, 2.595562f, 2.583f, 2.5765f,   2.599f,
    2.678062f, 2.809f, 3.01075f,  3.24f,  3.458187f, 3.67f,  3.863125f, 4.05f,
    4.239563f, 4.43f,  4.619563f, 4.817f, 5.034125f, 5.26f,  5.485625f, 5.717f};

MetalMaterial *CreateMetalMaterial(const TextureParams &mp) {
    static Spectrum copperN =
        Spectrum::FromSampled(CopperWavelengths, CopperN, CopperSamples);
    std::shared_ptr<Texture<Spectrum>> eta =
        mp.GetSpectrumTexture("eta", copperN);
    static Spectrum copperK =
        Spectrum::FromSampled(CopperWavelengths, CopperK, CopperSamples);
    std::shared_ptr<Texture<Spectrum>> k = mp.GetSpectrumTexture("k", copperK);
    std::shared_ptr<Texture<Float>> roughness =
        mp.GetFloatTexture("roughness", .01f);
    std::shared_ptr<Texture<Float>> uRoughness =
        mp.GetFloatTextureOrNull("uroughness");
    std::shared_ptr<Texture<Float>> vRoughness =
        mp.GetFloatTextureOrNull("vroughness");
    std::shared_ptr<Texture<Float>> bumpMap =
        mp.GetFloatTextureOrNull("bumpmap");
    bool remapRoughness = mp.FindBool("remaproughness", true);
    return new MetalMaterial(eta, k, roughness, uRoughness, vRoughness, bumpMap,
                             remapRoughness);
}

}  // namespace pbrt
