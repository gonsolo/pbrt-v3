
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

#if defined(_MSC_VER)
#define NOMINMAX
#pragma once
#endif

#ifndef PBRT_CORE_SAMPLING_H
#define PBRT_CORE_SAMPLING_H

// core/sampling.h*
#include "pbrt.h"
#include "geometry.h"
#include "rng.h"
#include <algorithm>

namespace pbrt {

// Sampling Declarations
void StratifiedSample1D(Float *samples, int nsamples, RNG &rng,
                        bool jitter = true);
void StratifiedSample2D(Point2f *samples, int nx, int ny, RNG &rng,
                        bool jitter = true);
void LatinHypercube(Float *samples, int nSamples, int nDim, RNG &rng);
struct Distribution1D {
    // Distribution1D Public Methods
    Distribution1D(const Float *f, size_t n) : func(f, f + n), cdf(n + 1) {
        // Compute integral of step function at $x_i$
        cdf[0] = 0;
        for (int i = 1; i < n + 1; ++i) cdf[i] = cdf[i - 1] + func[i - 1] / n;

        // Transform step function integral into CDF
        funcInt = cdf[n];
        if (funcInt == 0) {
            for (int i = 1; i < n + 1; ++i) cdf[i] = Float(i) / Float(n);
        } else {
            for (int i = 1; i < n + 1; ++i) cdf[i] /= funcInt;
        }
    }
    size_t Count() const { return func.size(); }
    Float SampleContinuous(Float u, Float *pdf, size_t *off = nullptr) const {
        // Find surrounding CDF segments and _offset_
        size_t offset = FindInterval(cdf.size(),
                                  [&](size_t index) { return cdf[index] <= u; });
        if (off) *off = offset;
        // Compute offset along CDF segment
        Float du = u - cdf[offset];
        if ((cdf[offset + 1] - cdf[offset]) > 0) {
            CHECK_GT(cdf[offset + 1], cdf[offset]);
            du /= (cdf[offset + 1] - cdf[offset]);
        }
        DCHECK(!std::isnan(du));

        // Compute PDF for sampled offset
        if (pdf) *pdf = (funcInt > 0) ? func[offset] / funcInt : 0;

        // Return $x\in{}[0,1)$ corresponding to sample
        return (offset + du) / Count();
    }
    size_t SampleDiscrete(Float u, Float *pdf = nullptr,
                       Float *uRemapped = nullptr) const {
        // Find surrounding CDF segments and _offset_
        size_t offset = FindInterval(cdf.size(),
                                  [&](size_t index) { return cdf[index] <= u; });
        if (pdf) *pdf = (funcInt > 0) ? func[offset] / (funcInt * Count()) : 0;
        if (uRemapped)
            *uRemapped = (u - cdf[offset]) / (cdf[offset + 1] - cdf[offset]);
        if (uRemapped) CHECK(*uRemapped >= 0.f && *uRemapped <= 1.f);
        return offset;
    }
    Float DiscretePDF(size_t index) const {
        CHECK(index >= 0 && index < Count());
        return func[index] / (funcInt * Count());
    }

    // Distribution1D Public Data
    std::vector<Float> func, cdf;
    Float funcInt;
};

Point2f RejectionSampleDisk(RNG &rng);
Vector3f UniformSampleHemisphere(const Point2f &u);
Float UniformHemispherePdf();
Vector3f UniformSampleSphere(const Point2f &u);
Float UniformSpherePdf();
Vector3f UniformSampleCone(const Point2f &u, Float thetamax);
Vector3f UniformSampleCone(const Point2f &u, Float thetamax, const Vector3f &x,
                           const Vector3f &y, const Vector3f &z);
Float UniformConePdf(Float thetamax);
Point2f UniformSampleDisk(const Point2f &u);
Point2f ConcentricSampleDisk(const Point2f &u);
Point2f UniformSampleTriangle(const Point2f &u);
class Distribution2D {
  public:
    // Distribution2D Public Methods
    Distribution2D(const Float *data, int nu, int nv);
    Point2f SampleContinuous(const Point2f &u, Float *pdf) const {
        Float pdfs[2];
        size_t v;
        Float d1 = pMarginal->SampleContinuous(u[1], &pdfs[1], &v);
        Float d0 = pConditionalV[v]->SampleContinuous(u[0], &pdfs[0]);
        *pdf = pdfs[0] * pdfs[1];
        return Point2f(d0, d1);
    }
    Float Pdf(const Point2f &p) const {
        size_t iu = Clamp(size_t(p[0] * pConditionalV[0]->Count()), 0,
                       pConditionalV[0]->Count() - 1);
        size_t iv =
            Clamp(size_t(p[1] * pMarginal->Count()), 0, pMarginal->Count() - 1);
        return pConditionalV[iv]->func[iu] / pMarginal->funcInt;
    }

  private:
    // Distribution2D Private Data
    std::vector<std::unique_ptr<Distribution1D>> pConditionalV;
    std::unique_ptr<Distribution1D> pMarginal;
};

// Sampling Inline Functions
template <typename T>
void Shuffle(T *samp, size_t count, size_t nDimensions, RNG &rng) {
    for (int i = 0; i < count; ++i) {
        int other = i + rng.UniformUInt32(static_cast<uint32_t>(count - i));
        for (int j = 0; j < nDimensions; ++j)
            std::swap(samp[nDimensions * i + j], samp[nDimensions * other + j]);
    }
}

inline Vector3f CosineSampleHemisphere(const Point2f &u) {
    Point2f d = ConcentricSampleDisk(u);
    Float z = std::sqrt(std::max((Float)0, 1 - d.x * d.x - d.y * d.y));
    return Vector3f(d.x, d.y, z);
}

inline Float CosineHemispherePdf(Float cosTheta) { return cosTheta * InvPi; }

inline Float BalanceHeuristic(int nf, Float fPdf, int ng, Float gPdf) {
    return (nf * fPdf) / (nf * fPdf + ng * gPdf);
}

inline Float PowerHeuristic(int nf, Float fPdf, int ng, Float gPdf) {
    Float f = nf * fPdf, g = ng * gPdf;
    return (f * f) / (f * f + g * g);
}

}  // namespace pbrt

#endif  // PBRT_CORE_SAMPLING_H
