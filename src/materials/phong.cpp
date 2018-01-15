#include "materials/phong.h"
#include "bssrdf.h"
#include "interaction.h"
#include "paramset.h"
#include "reflection.h"
#include "stats.h"
#include "stringprint.h"
#include "texture.h"

namespace pbrt {

class PhongBRDF : public BxDF {
  public:
    PhongBRDF(const Spectrum &R, Float n)
        : BxDF(BxDFType(BSDF_REFLECTION | BSDF_GLOSSY)), R(R), n(n) {}
    Spectrum f(const Vector3f &wo, const Vector3f &wi) const;
	Spectrum f_analytical(const Vector3f &wo, const Vector3f &wi, Float cosThetaLight, Float cosNormalLight) const;

    Spectrum Sample_f(const Vector3f &wo, Vector3f *wi, const Point2f &u,
                      Float *pdf, BxDFType *sampledType) const;
    Float Pdf(const Vector3f &wo, const Vector3f &wi) const;
	
	//Spectrum rho(const Vector3f &, int, const Point2f *) const { return R; }
    //Spectrum rho(int, const Point2f *, const Point2f *) const { return R; }
    std::string ToString() const;

  private:
    Spectrum R;
	Float n;
};

Spectrum PhongBRDF::f(const Vector3f &wo, const Vector3f &wi) const {
	//std::cout << "Phong: " << std::endl;
	Vector3f r = Reflect(wi, Vector3f(0,0,1));
	//std::cout << r << std::endl;
	Float d = AbsDot(wo, r);
	//std::cout << d << std::endl;
	Float dn = powf(d, n);
	//std::cout << dn << std::endl;
	Spectrum result = R * dn;
	//std::cout << result << std::endl;
	return result;
}

Spectrum PhongBRDF::f_analytical(const Vector3f &wo, const Vector3f &wi, Float cosThetaLight, Float cosNormalLight) const {

	Spectrum result;
	// TODO
	return result;
}

std::string PhongBRDF::ToString() const {
    return StringPrintf("[ PhongBRDF R: %s ]", R.ToString().c_str());
}

Spectrum PhongBRDF::Sample_f(const Vector3f &wo, Vector3f *wi,
                                   const Point2f &u, Float *pdf,
                                   BxDFType *sampledType) const {
	Spectrum s;
	// TODO
    return s;
}

Float PhongBRDF::Pdf(const Vector3f &wo, const Vector3f &wi) const {
	// TODO
	return 1.f;
}

void PhongMaterial::ComputeScatteringFunctions(SurfaceInteraction *si,
                                                MemoryArena &arena,
                                                TransportMode mode,
                                                bool allowMultipleLobes) const {

    // Evaluate textures for _DisneyMaterial_ material and allocate BRDF
    si->bsdf = ARENA_ALLOC(arena, BSDF)(*si);

    // Diffuse
    Spectrum c = color->Evaluate(*si).Clamp();
    Float n = phongExp->Evaluate(*si);
	si->bsdf->Add(ARENA_ALLOC(arena, PhongBRDF)(c, n));
}

PhongMaterial *CreatePhongMaterial(const TextureParams &mp) {
    std::shared_ptr<Texture<Spectrum>> color = mp.GetSpectrumTexture("color", Spectrum(0.5f));
    std::shared_ptr<Texture<Float>> phongExp = mp.GetFloatTexture("phongexp", 1.f);
    return new PhongMaterial(color, phongExp);
}

}  // namespace pbrt
