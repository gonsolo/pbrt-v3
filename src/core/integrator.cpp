
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

#include <algorithm>
 // core/integrator.cpp*
#include "integrator.h"
#include "scene.h"
#include "geometry.h"
#include "interaction.h"
#include "sampling.h"
#include "parallel.h"
#include "film.h"
#include "sampler.h"
#include "integrator.h"
#include "progressreporter.h"
#include "camera.h"
#include "stats.h"

#include "lights/diffuse.h"
#include "shapes/sphere.h"

namespace pbrt {

	STAT_COUNTER("Integrator/Camera rays traced", nCameraRays);

	// Integrator Method Definitions
	Integrator::~Integrator() {}

	Float GonzoCombinedCosine;

	Sphere* getSphere(const Light& light) {
		const Light* lightPointer = &light;
		const DiffuseAreaLight* dal = dynamic_cast<const DiffuseAreaLight*>(lightPointer);
		auto shape = dal->shape;
		const std::type_info& info = typeid(*shape);
		Shape* tmp = shape.get();
		Sphere* sphere = dynamic_cast<Sphere*>(tmp);
		return sphere;
	}

	Vector3f getIncoming(const Light& light, const Interaction& it) {
		Sphere* sphere = getSphere(light);
		Point3f sphereCenter = (*sphere->ObjectToWorld)(Point3f(0, 0, 0));
		Vector3f wiu = sphereCenter - it.p;
		//Vector3f wi = Normalize(wiu);
		return wiu;
	}

	Float getRadius(const Light& light) {
		Sphere* sphere = getSphere(light);
		return sphere->radius;
	}

	// Integrator Utility Functions
	Spectrum UniformSampleAllLights(const Interaction &it, const Scene &scene,
		MemoryArena &arena, Sampler &sampler,
		const std::vector<int> &nLightSamples,
		bool handleMedia) {
		ProfilePhase p(Prof::DirectLighting);
		Spectrum L(0.f);
		for (size_t j = 0; j < scene.lights.size(); ++j) {
			// Accumulate contribution of _j_th light to _L_
			const std::shared_ptr<Light> &light = scene.lights[j];
			int nSamples = nLightSamples[j];
			const Point2f *uLightArray = sampler.Get2DArray(nSamples);
			const Point2f *uScatteringArray = sampler.Get2DArray(nSamples);
			if (!uLightArray || !uScatteringArray) {
				// Use a single sample for illumination from _light_
				Point2f uLight = sampler.Get2D();
				Point2f uScattering = sampler.Get2D();
				L += EstimateDirect(it, uScattering, *light, uLight, scene, sampler,
					arena, handleMedia);
			}
			else {
				// Estimate direct lighting using sample arrays
				Spectrum Ld(0.f);
				GonzoCombinedCosine = 0;
				for (int k = 0; k < nSamples; ++k)
					Ld += EstimateDirect(it, uScatteringArray[k], *light,
						uLightArray[k], scene, sampler, arena,
						handleMedia);

				Vector3f wiu = getIncoming(*light, it);
				Vector3f wi = Normalize(wiu);
				Float wiLength = wiu.Length();
				const SurfaceInteraction &isect = (const SurfaceInteraction &)it;
				Float radius = getRadius(*light);
				//std::cout << "wi, wiLength, n, radius: " << wi << " " << wiLength << " " << isect.shading.n << " " << radius << std::endl;
				Float Cosine = AbsDot(wi, isect.shading.n);
				Float w = std::acos(Cosine);
				Float a = std::asin(radius / wiLength);
				//std::cout << "w, a: " << w << " " << a << std::endl;
				//std::cout << "Cosine, Combined Cosine: " << Cosine << " " << GonzoCombinedCosine / nSamples << std::endl;


				L += Ld / nSamples;
			}
		}
		return L;
	}

	Spectrum UniformSampleOneLight(const Interaction &it, const Scene &scene,
		MemoryArena &arena, Sampler &sampler,
		bool handleMedia, const Distribution1D *lightDistrib) {
		ProfilePhase p(Prof::DirectLighting);
		// Randomly choose a single light to sample, _light_
		int nLights = int(scene.lights.size());
		if (nLights == 0) return Spectrum(0.f);
		size_t lightNum;
		Float lightPdf;
		if (lightDistrib) {
			lightNum = lightDistrib->SampleDiscrete(sampler.Get1D(), &lightPdf);
			if (lightPdf == 0) return Spectrum(0.f);
		}
		else {
			lightNum = std::min((int)(sampler.Get1D() * nLights), nLights - 1);
			lightPdf = Float(1) / nLights;
		}
		const std::shared_ptr<Light> &light = scene.lights[lightNum];
		Point2f uLight = sampler.Get2D();
		Point2f uScattering = sampler.Get2D();
		return EstimateDirect(it, uScattering, *light, uLight,
			scene, sampler, arena, handleMedia) / lightPdf;
	}

	// ThetaLight, ThetaNormal, result

	const int CircleCosineLength = 256;
	double CircleCosine[CircleCosineLength][3] = {
		{ 0.01, 0., -1.0 },
		{ 0.01, 0.10, 0.9949547518 },
		{ 0.01, 0.20, 0.9800189546 },
		{ 0.01, 0.30, 0.9552919707 },
		{ 0.01, 0.40, 0.9210211231 },
		{ 0.01, 0.50, 0.8775493389 },
		{ 0.01, 0.60, 0.8253119966 },
		{ 0.01, 0.70, 0.7648334270 },
		{ 0.01, 0.80, 0.6967250360 },
		{ 0.01, 0.90, 0.6217027544 },
		{ 0.01, 1.00, 0.5402315347 },
		{ 0.01, 1.10, 0.4534176941 },
		{ 0.01, 1.20, 0.3622548919 },
		{ 0.01, 1.30, 0.2674249375 },
		{ 0.01, 1.40, 0.1699107401 },
		{ 0.01, 1.50, 0.07069376491 },
		{ 0.11, 0., -1.0 },
		{ 0.11, 0.10, 0.9890315520 },
		{ 0.11, 0.20, 0.9743119047 },
		{ 0.11, 0.30, 0.9499601222 },
		{ 0.11, 0.40, 0.9162521794 },
		{ 0.11, 0.50, 0.8735885242 },
		{ 0.11, 0.60, 0.8225275782 },
		{ 0.11, 0.70, 0.7638985663 },
		{ 0.11, 0.80, 0.6993212855 },
		{ 0.11, 0.90, 0.6382264489 },
		{ 0.11, 1.00, 0.5317751532 },
		{ 0.11, 1.10, 0.4276738447 },
		{ 0.11, 1.20, 0.3496837318 },
		{ 0.11, 1.30, 0.2585130451 },
		{ 0.11, 1.40, 0.1631317988 },
		{ 0.11, 1.50, 0.06547993100 },
		{ 0.21, 0., -1.0 },
		{ 0.21, 0.10, 0.9732988662 },
		{ 0.21, 0.20, 0.9591688393 },
		{ 0.21, 0.30, 0.9358450300 },
		{ 0.21, 0.40, 0.9036920685 },
		{ 0.21, 0.50, 0.8632954214 },
		{ 0.21, 0.60, 0.8156385315 },
		{ 0.21, 0.70, 0.7627776948 },
		{ 0.21, 0.80, 0.7141070868 },
		{ 0.21, 0.90, 0.6650239738 },
		{ 0.21, 1.00, 0.5095711657 },
		{ 0.21, 1.10, 0.3766895794 },
		{ 0.21, 1.20, 0.3110120846 },
		{ 0.21, 1.30, 0.2341696925 },
		{ 0.21, 1.40, 0.1449247920 },
		{ 0.21, 1.50, 0.05155460747 },
		{ 0.31, 0., -1.0 },
		{ 0.31, 0.10, 0.9479280482 },
		{ 0.31, 0.20, 0.9348023897 },
		{ 0.31, 0.30, 0.9132467726 },
		{ 0.31, 0.40, 0.8838193790 },
		{ 0.31, 0.50, 0.8475465766 },
		{ 0.31, 0.60, 0.8066252658 },
		{ 0.31, 0.70, 0.7700803058 },
		{ 0.31, 0.80, 0.7388300073 },
		{ 0.31, 0.90, 0.6507309135 },
		{ 0.31, 1.00, 0.4745366098 },
		{ 0.31, 1.10, 0.3335644127 },
		{ 0.31, 1.20, 0.2526718372 },
		{ 0.31, 1.30, 0.1896937596 },
		{ 0.31, 1.40, 0.1146687291 },
		{ 0.31, 1.50, 0.02879184928 },
		{ 0.41, 0., -1.0 },
		{ 0.41, 0.10, 0.9132012537 },
		{ 0.41, 0.20, 0.9015821467 },
		{ 0.41, 0.30, 0.8827288599 },
		{ 0.41, 0.40, 0.8576296277 },
		{ 0.41, 0.50, 0.8284913337 },
		{ 0.41, 0.60, 0.8040406272 },
		{ 0.41, 0.70, 0.7860477978 },
		{ 0.41, 0.80, 0.7323245121 },
		{ 0.41, 0.90, 0.6028873158 },
		{ 0.41, 1.00, 0.4280793723 },
		{ 0.41, 1.10, 0.2900799029 },
		{ 0.41, 1.20, 0.1986043018 },
		{ 0.41, 1.30, 0.1296618940 },
		{ 0.41, 1.40, 0.06832386241 },
		{ 0.41, 1.50, -0.003350604343 },
		{ 0.51, 0., -1.0 },
		{ 0.51, 0.10, 0.8695213448 },
		{ 0.51, 0.20, 0.8600932218 },
		{ 0.51, 0.30, 0.8453173012 },
		{ 0.51, 0.40, 0.8273324552 },
		{ 0.51, 0.50, 0.8147074580 },
		{ 0.51, 0.60, 0.8091092904 },
		{ 0.51, 0.70, 0.7754078958 },
		{ 0.51, 0.80, 0.6841238566 },
		{ 0.51, 0.90, 0.5332428740 },
		{ 0.51, 1.00, 0.3719928144 },
		{ 0.51, 1.10, 0.2425416479 },
		{ 0.51, 1.20, 0.1471146677 },
		{ 0.51, 1.30, 0.07291460683 },
		{ 0.51, 1.40, 0.009609746331 },
		{ 0.51, 1.50, -0.04845457201 },
		{ 0.61, 0., -1.0 },
		{ 0.61, 0.10, 0.8174416407 },
		{ 0.61, 0.20, 0.8113093032 },
		{ 0.61, 0.30, 0.8032024209 },
		{ 0.61, 0.40, 0.8016784453 },
		{ 0.61, 0.50, 0.8078927458 },
		{ 0.61, 0.60, 0.7883045041 },
		{ 0.61, 0.70, 0.7207791769 },
		{ 0.61, 0.80, 0.5996571630 },
		{ 0.61, 0.90, 0.4505179826 },
		{ 0.61, 1.00, 0.3083305766 },
		{ 0.61, 1.10, 0.1902629735 },
		{ 0.61, 1.20, 0.09607441538 },
		{ 0.61, 1.30, 0.01943568549 },
		{ 0.61, 1.40, -0.04566080339 },
		{ 0.61, 1.50, -0.1033488520 },
		{ 0.71, 0., -1.0 },
		{ 0.71, 0.10, 0.7577687682 },
		{ 0.71, 0.20, 0.7572605146 },
		{ 0.71, 0.30, 0.7654450383 },
		{ 0.71, 0.40, 0.7824147201 },
		{ 0.71, 0.50, 0.7725163911 },
		{ 0.71, 0.60, 0.7188967782 },
		{ 0.71, 0.70, 0.6201058076 },
		{ 0.71, 0.80, 0.4931801267 },
		{ 0.71, 0.90, 0.3610136041 },
		{ 0.71, 1.00, 0.2392762634 },
		{ 0.71, 1.10, 0.1338209737 },
		{ 0.71, 1.20, 0.04438558057 },
		{ 0.71, 1.30, -0.03161046296 },
		{ 0.71, 1.40, -0.09701804835 },
		{ 0.71, 1.50, -0.1542070699 },
		{ 0.81, 0., -1.0 },
		{ 0.81, 0.10, 0.6920509508 },
		{ 0.81, 0.20, 0.7073734447 },
		{ 0.81, 0.30, 0.7327525111 },
		{ 0.81, 0.40, 0.7261479363 },
		{ 0.81, 0.50, 0.6781596682 },
		{ 0.81, 0.60, 0.5950914704 },
		{ 0.81, 0.70, 0.4908879154 },
		{ 0.81, 0.80, 0.3791315379 },
		{ 0.81, 0.90, 0.2693892591 },
		{ 0.81, 1.00, 0.1670211982 },
		{ 0.81, 1.10, 0.07434058266 },
		{ 0.81, 1.20, -0.008153908290 },
		{ 0.81, 1.30, -0.08088834190 },
		{ 0.81, 1.40, -0.1446736025 },
		{ 0.81, 1.50, -0.2004130566 },
		{ 0.91, 0., -1.0 },
		{ 0.91, 0.10, 0.6296462073 },
		{ 0.91, 0.20, 0.6556467307 },
		{ 0.91, 0.30, 0.6396356318 },
		{ 0.91, 0.40, 0.5919606400 },
		{ 0.91, 0.50, 0.5241420114 },
		{ 0.91, 0.60, 0.4439872118 },
		{ 0.91, 0.70, 0.3570923257 },
		{ 0.91, 0.80, 0.2677058477 },
		{ 0.91, 0.90, 0.1790828028 },
		{ 0.91, 1.00, 0.09366012246 },
		{ 0.91, 1.10, 0.01318557425 },
		{ 0.91, 1.20, -0.06116260629 },
		{ 0.91, 1.30, -0.1286614516 },
		{ 0.91, 1.40, -0.1889387646 },
		{ 0.91, 1.50, -0.2418825480 },
		{ 1.01, 0., -1.0 },
		{ 1.01, 0.10, 0.5140095783 },
		{ 1.01, 0.20, 0.4935586053 },
		{ 1.01, 0.30, 0.4607801394 },
		{ 1.01, 0.40, 0.4167614120 },
		{ 1.01, 0.50, 0.3631349782 },
		{ 1.01, 0.60, 0.3018200441 },
		{ 1.01, 0.70, 0.2348898361 },
		{ 1.01, 0.80, 0.1644547669 },
		{ 1.01, 0.90, 0.09256123971 },
		{ 1.01, 1.00, 0.02111062588 },
		{ 1.01, 1.10, -0.04819927190 },
		{ 1.01, 1.20, -0.1139106026 },
		{ 1.01, 1.30, -0.1748220433 },
		{ 1.01, 1.40, -0.2299890991 },
		{ 1.01, 1.50, -0.2787120709 },
		{ 1.11, 0., -1.0 },
		{ 1.11, 0.10, 0.4266365118 },
		{ 1.11, 0.20, 0.3737334932 },
		{ 1.11, 0.30, 0.3279944848 },
		{ 1.11, 0.40, 0.2837738214 },
		{ 1.11, 0.50, 0.2366003340 },
		{ 1.11, 0.60, 0.1852245346 },
		{ 1.11, 0.70, 0.1299064109 },
		{ 1.11, 0.80, 0.07158546895 },
		{ 1.11, 0.90, 0.01151129321 },
		{ 1.11, 1.00, -0.04894186727 },
		{ 1.11, 1.10, -0.1083829856 },
		{ 1.11, 1.20, -0.1654745460 },
		{ 1.11, 1.30, -0.2189809332 },
		{ 1.11, 1.40, -0.2678016220 },
		{ 1.11, 1.50, -0.3109924344 },
		{ 1.21, 0., -1.0 },
		{ 1.21, 0.10, 0.3429958701 },
		{ 1.21, 0.20, 0.3102842292 },
		{ 1.21, 0.30, 0.2517553714 },
		{ 1.21, 0.40, 0.1963726543 },
		{ 1.21, 0.50, 0.1440924759 },
		{ 1.21, 0.60, 0.09285574459 },
		{ 1.21, 0.70, 0.04140715926 },
		{ 1.21, 0.80, -0.01062507720 },
		{ 1.21, 0.90, -0.06299550930 },
		{ 1.21, 1.00, -0.1150748959 },
		{ 1.21, 1.10, -0.1660116276 },
		{ 1.21, 1.20, -0.2148400509 },
		{ 1.21, 1.30, -0.2605548432 },
		{ 1.21, 1.40, -0.3021636223 },
		{ 1.21, 1.50, -0.3387246654 },
		{ 1.31, 0., -1.0 },
		{ 1.31, 0.10, 0.2506426848 },
		{ 1.31, 0.20, 0.2286164050 },
		{ 1.31, 0.30, 0.1891086570 },
		{ 1.31, 0.40, 0.1292396488 },
		{ 1.31, 0.50, 0.07176934369 },
		{ 1.31, 0.60, 0.01773024369 },
		{ 1.31, 0.70, -0.03357475112 },
		{ 1.31, 0.80, -0.08282637562 },
		{ 1.31, 0.90, -0.1303576135 },
		{ 1.31, 1.00, -0.1761389772 },
		{ 1.31, 1.10, -0.2198519954 },
		{ 1.31, 1.20, -0.2609717758 },
		{ 1.31, 1.30, -0.2988370554 },
		{ 1.31, 1.40, -0.3327047476 },
		{ 1.31, 1.50, -0.3617907109 },
		{ 1.41, 0., -1.0 },
		{ 1.41, 0.10, 0.1546023605 },
		{ 1.41, 0.20, 0.1380044297 },
		{ 1.41, 0.30, 0.1098491621 },
		{ 1.41, 0.40, 0.06782404636 },
		{ 1.41, 0.50, 0.009369397113 },
		{ 1.41, 0.60, -0.04635433327 },
		{ 1.41, 0.70, -0.09811444887 },
		{ 1.41, 0.80, -0.1460069604 },
		{ 1.41, 0.90, -0.1903196880 },
		{ 1.41, 1.00, -0.2312513025 },
		{ 1.41, 1.10, -0.2688210246 },
		{ 1.41, 1.20, -0.3028629786 },
		{ 1.41, 1.30, -0.3330512988 },
		{ 1.41, 1.40, -0.3589316634 },
		{ 1.41, 1.50, -0.3799490090 },
		{ 1.51, 0., -1.0 },
		{ 1.51, 0.10, 0.05652971095 },
		{ 1.51, 0.20, 0.04382887442 },
		{ 1.51, 0.30, 0.02256692978 },
		{ 1.51, 0.40, -0.007656015433 },
		{ 1.51, 0.50, -0.04889745963 },
		{ 1.51, 0.60, -0.1035055507 },
		{ 1.51, 0.70, -0.1546792402 },
		{ 1.51, 0.80, -0.2011928308 },
		{ 1.51, 0.90, -0.2428715110 },
		{ 1.51, 1.00, -0.2797777350 },
		{ 1.51, 1.10, -0.3120036002 },
		{ 1.51, 1.20, -0.3395709254 },
		{ 1.51, 1.30, -0.3623916228 },
		{ 1.51, 1.40, -0.3802570566 },
		{ 1.51, 1.50, -0.3928375137 },
	};



	static bool verbose = true;

	Spectrum EstimateIncomingAnalytical(const DiffuseAreaLight* dal,
			const Interaction& it,
			Float* pdf,
			VisibilityTester* vis,
			Float* cosThetaLight,
			Vector3f* wi) {

		Spectrum Ld;
		auto shape = dal->shape;
		const std::type_info& info = typeid(*shape);
		Shape* tmp = shape.get();
		Sphere* sphere = dynamic_cast<Sphere*>(tmp);
		Point3f sphereCenter = (*sphere->ObjectToWorld)(Point3f(0, 0, 0));

		Float radius = sphere->radius;
		auto Lemit = dal->Lemit;

		// Uniform cone sampling from PBRT always gives a noise free image since
		// the computation is independent from distance and angle of the sample
		// point. This is because we don't have to convert from area measure to
		// solid angle measure.
		Float radiusSquared = Square(radius);
		Float distanceSquared = DistanceSquared(sphereCenter, it.p);
		Float sinThetaLight2 = radiusSquared / distanceSquared;
		*cosThetaLight = std::max(0.f, std::sqrt(1.f - sinThetaLight2));
		Spectrum Li = Lemit * 2 * Pi * (1 - *cosThetaLight);
		*pdf = 1; // TODO: Higher?
		Interaction lightit;
		lightit.gonzoSphericalAreaLight = true;
		lightit.gonzoRadius = radius;
		lightit.p = sphereCenter;
	    *vis = VisibilityTester(it, lightit);
		*wi = Normalize(sphereCenter - it.p);
		return Li; // Only incoming light
	}

	extern Float GonzoCombinedCosine;


	Spectrum EstimateDirect(const Interaction &it, const Point2f &uScattering,
		const Light &light, const Point2f &uLight,
		const Scene &scene, Sampler &sampler,
		MemoryArena &arena, bool handleMedia, bool specular) {
		BxDFType bsdfFlags =
			specular ? BSDF_ALL : BxDFType(BSDF_ALL & ~BSDF_SPECULAR);
		Spectrum Ld(0.f);
		// Sample light source with multiple importance sampling
		Vector3f wi;
		Float lightPdf = 0, scatteringPdf = 0;
		VisibilityTester visibility;

		const Light* lightPointer = &light;
		const DiffuseAreaLight* dal = dynamic_cast<const DiffuseAreaLight*>(lightPointer);

		Spectrum Li;
		Float cosThetaLight;
		if (dal && dal->analytical) {
			Li = EstimateIncomingAnalytical(dal, it, &lightPdf, &visibility, &cosThetaLight, &wi);
			//std::cout << "Analytical Li: " << Li << std::endl;
		}
		else {
			Li = light.Sample_Li(it, uLight, &wi, &lightPdf, &visibility);
		}

		VLOG(2) << "EstimateDirect uLight:" << uLight << " -> Li: " << Li << ", wi: "
			<< wi << ", pdf: " << lightPdf;
		if (lightPdf > 0 && !Li.IsBlack()) {
			// Compute BSDF or phase function's value for light sample
			Spectrum f;
			if (it.IsSurfaceInteraction()) {
				// Evaluate BSDF for light sampling strategy
				const SurfaceInteraction &isect = (const SurfaceInteraction &)it;
				if (dal && dal->analytical) {
					Float cosNormalLight = AbsDot(wi, isect.shading.n);
					f = isect.bsdf->f_analytical(isect.wo, wi, cosThetaLight, cosNormalLight, bsdfFlags);
					//std::cout << "Analytical BSDF: " << f << std::endl;
					scatteringPdf = 1.f; // GONZO: TODO
				}
				else {
					f = isect.bsdf->f(isect.wo, wi, bsdfFlags) * AbsDot(wi, isect.shading.n);
					scatteringPdf = isect.bsdf->Pdf(isect.wo, wi, bsdfFlags);
				}
				VLOG(2) << "  surf f*dot :" << f << ", scatteringPdf: " << scatteringPdf;
			}
			else {
				// Evaluate phase function for light sampling strategy
				const MediumInteraction &mi = (const MediumInteraction &)it;
				Float p = mi.phase->p(mi.wo, wi);
				f = Spectrum(p);
				scatteringPdf = p;
				VLOG(2) << "  medium p: " << p;
			}
			if (!f.IsBlack()) {
				// Compute effect of visibility for light source sample
				if (handleMedia) {
					Li *= visibility.Tr(scene, sampler);
					VLOG(2) << "  after Tr, Li: " << Li;
				}
				else {
					if (!visibility.Unoccluded(scene)) {
						VLOG(2) << "  shadow ray blocked";
						Li = Spectrum(0.f);
					}
					else
						VLOG(2) << "  shadow ray unoccluded";
				}

				// Add light's contribution to reflected radiance
				if (!Li.IsBlack()) {
					if (IsDeltaLight(light.flags))
						Ld += f * Li / lightPdf;
					else {
						Float weight =
							PowerHeuristic(1, lightPdf, 1, scatteringPdf);
						Ld += f * Li * weight / lightPdf;
					}
				}
			}
		}

		// GONZO: No BSDF sampling for now
		return Ld;

		// Sample BSDF with multiple importance sampling
		if (!IsDeltaLight(light.flags)) {
			Spectrum f;
			bool sampledSpecular = false;
			if (it.IsSurfaceInteraction()) {
				// Sample scattered direction for surface interactions
				BxDFType sampledType;
				const SurfaceInteraction &isect = (const SurfaceInteraction &)it;
				f = isect.bsdf->Sample_f(isect.wo, &wi, uScattering, &scatteringPdf,
					bsdfFlags, &sampledType);
				f *= AbsDot(wi, isect.shading.n);
				sampledSpecular = (sampledType & BSDF_SPECULAR) != 0;
			}
			else {
				// Sample scattered direction for medium interactions
				const MediumInteraction &mi = (const MediumInteraction &)it;
				Float p = mi.phase->Sample_p(mi.wo, &wi, uScattering);
				f = Spectrum(p);
				scatteringPdf = p;
			}
			VLOG(2) << "  BSDF / phase sampling f: " << f << ", scatteringPdf: " <<
				scatteringPdf;
			if (!f.IsBlack() && scatteringPdf > 0) {
				// Account for light contributions along sampled direction _wi_
				Float weight = 1;
				if (!sampledSpecular) {
					lightPdf = light.Pdf_Li(it, wi);
					if (lightPdf == 0) return Ld;
					weight = PowerHeuristic(1, scatteringPdf, 1, lightPdf);
				}

				// Find intersection and compute transmittance
				SurfaceInteraction lightIsect;
				Ray ray = it.SpawnRay(wi);
				Spectrum Tr(1.f);
				bool foundSurfaceInteraction =
					handleMedia ? scene.IntersectTr(ray, sampler, &lightIsect, &Tr)
					: scene.Intersect(ray, &lightIsect);

				// Add light contribution from material sampling
				Spectrum Li(0.f);
				if (foundSurfaceInteraction) {
					if (lightIsect.primitive->GetAreaLight() == &light)
						Li = lightIsect.Le(-wi);
				}
				else
					Li = light.Le(ray);
				if (!Li.IsBlack()) Ld += f * Li * Tr * weight / scatteringPdf;
			}
		}
		return Ld;
	}

	std::unique_ptr<Distribution1D> ComputeLightPowerDistribution(
		const Scene &scene) {
		if (scene.lights.empty()) return nullptr;
		std::vector<Float> lightPower;
		for (const auto &light : scene.lights)
			lightPower.push_back(light->Power().y());
		return std::unique_ptr<Distribution1D>(
			new Distribution1D(&lightPower[0], lightPower.size()));
	}

	// SamplerIntegrator Method Definitions
	void SamplerIntegrator::Render(const Scene &scene) {
		Preprocess(scene, *sampler);
		// Render image tiles in parallel

		// Compute number of tiles, _nTiles_, to use for parallel rendering
		Bounds2i sampleBounds = camera->film->GetSampleBounds();
		Vector2i sampleExtent = sampleBounds.Diagonal();
		const Int tileSize = 16;
		Point2i nTiles((sampleExtent.x + tileSize - 1) / tileSize,
			(sampleExtent.y + tileSize - 1) / tileSize);
		ProgressReporter reporter(nTiles.x * nTiles.y, "Rendering");
		{
			ParallelFor2D([&](Point2i tile) {
				// Render section of image corresponding to _tile_

				// Allocate _MemoryArena_ for tile
				MemoryArena arena;

				// Get sampler instance for tile
				int seed = tile.y * nTiles.x + tile.x;
				std::unique_ptr<Sampler> tileSampler = sampler->Clone(seed);

				// Compute sample bounds for tile
				Int x0 = sampleBounds.pMin.x + tile.x * tileSize;
				Int x1 = std::min(x0 + tileSize, sampleBounds.pMax.x);
				Int y0 = sampleBounds.pMin.y + tile.y * tileSize;
				Int y1 = std::min(y0 + tileSize, sampleBounds.pMax.y);
				Bounds2i tileBounds(Point2i(x0, y0), Point2i(x1, y1));
				LOG(INFO) << "Starting image tile " << tileBounds;

				// Get _FilmTile_ for tile
				std::unique_ptr<FilmTile> filmTile =
					camera->film->GetFilmTile(tileBounds);

				// Loop over pixels in tile to render them
				for (Point2i pixel : tileBounds) {
					{
						ProfilePhase pp(Prof::StartPixel);
						tileSampler->StartPixel(pixel);
					}

					// Do this check after the StartPixel() call; this keeps
					// the usage of RNG values from (most) Samplers that use
					// RNGs consistent, which improves reproducability /
					// debugging.
					if (!InsideExclusive(pixel, pixelBounds))
						continue;

					do {
						// Initialize _CameraSample_ for current sample
						CameraSample cameraSample =
							tileSampler->GetCameraSample(pixel);

						// Generate camera ray for current sample
						RayDifferential ray;
						Float rayWeight =
							camera->GenerateRayDifferential(cameraSample, &ray);
						ray.ScaleDifferentials(
							1 / std::sqrt((Float)tileSampler->samplesPerPixel));
						++nCameraRays;

						// Evaluate radiance along camera ray
						Spectrum L(0.f);
						if (rayWeight > 0) L = Li(ray, scene, *tileSampler, arena);

						// Issue warning if unexpected radiance value returned
						if (L.HasNaNs()) {
							LOG(ERROR) << StringPrintf(
								"Not-a-number radiance value returned "
								"for pixel (%d, %d), sample %d. Setting to black.",
								pixel.x, pixel.y,
								(int)tileSampler->CurrentSampleNumber());
							L = Spectrum(0.f);
						}
						else if (L.y() < -1e-5) {
							LOG(ERROR) << StringPrintf(
								"Negative luminance value, %f, returned "
								"for pixel (%d, %d), sample %d. Setting to black.",
								L.y(), pixel.x, pixel.y,
								(int)tileSampler->CurrentSampleNumber());
							L = Spectrum(0.f);
						}
						else if (std::isinf(L.y())) {
							LOG(ERROR) << StringPrintf(
								"Infinite luminance value returned "
								"for pixel (%d, %d), sample %d. Setting to black.",
								pixel.x, pixel.y,
								(int)tileSampler->CurrentSampleNumber());
							L = Spectrum(0.f);
						}
						VLOG(1) << "Camera sample: " << cameraSample << " -> ray: " <<
							ray << " -> L = " << L;

						// Add camera ray's contribution to image
						filmTile->AddSample(cameraSample.pFilm, L, rayWeight);

						// Free _MemoryArena_ memory from computing image sample
						// value
						arena.Reset();
					} while (tileSampler->StartNextSample());
				}
				LOG(INFO) << "Finished image tile " << tileBounds;

				// Merge image tile into _Film_
				camera->film->MergeFilmTile(std::move(filmTile));
				reporter.Update();
			}, nTiles);
			reporter.Done();
		}
		LOG(INFO) << "Rendering finished";

		// Save final image after rendering
		camera->film->WriteImage();
	}

	Spectrum SamplerIntegrator::SpecularReflect(
		const RayDifferential &ray, const SurfaceInteraction &isect,
		const Scene &scene, Sampler &sampler, MemoryArena &arena, int depth) const {
		// Compute specular reflection direction _wi_ and BSDF value
		Vector3f wo = isect.wo, wi;
		Float pdf;
		BxDFType type = BxDFType(BSDF_REFLECTION | BSDF_SPECULAR);
		Spectrum f = isect.bsdf->Sample_f(wo, &wi, sampler.Get2D(), &pdf, type);

		// Return contribution of specular reflection
		const Normal3f &ns = isect.shading.n;
		if (pdf > 0.f && !f.IsBlack() && AbsDot(wi, ns) != 0.f) {
			// Compute ray differential _rd_ for specular reflection
			RayDifferential rd = isect.SpawnRay(wi);
			if (ray.hasDifferentials) {
				rd.hasDifferentials = true;
				rd.rxOrigin = isect.p + isect.dpdx;
				rd.ryOrigin = isect.p + isect.dpdy;
				// Compute differential reflected directions
				Normal3f dndx = isect.shading.dndu * isect.dudx +
					isect.shading.dndv * isect.dvdx;
				Normal3f dndy = isect.shading.dndu * isect.dudy +
					isect.shading.dndv * isect.dvdy;
				Vector3f dwodx = -ray.rxDirection - wo,
					dwody = -ray.ryDirection - wo;
				Float dDNdx = Dot(dwodx, ns) + Dot(wo, dndx);
				Float dDNdy = Dot(dwody, ns) + Dot(wo, dndy);
				rd.rxDirection =
					wi - dwodx + 2.f * Vector3f(Dot(wo, ns) * dndx + dDNdx * ns);
				rd.ryDirection =
					wi - dwody + 2.f * Vector3f(Dot(wo, ns) * dndy + dDNdy * ns);
			}
			return f * Li(rd, scene, sampler, arena, depth + 1) * AbsDot(wi, ns) /
				pdf;
		}
		else
			return Spectrum(0.f);
	}

	Spectrum SamplerIntegrator::SpecularTransmit(
		const RayDifferential &ray, const SurfaceInteraction &isect,
		const Scene &scene, Sampler &sampler, MemoryArena &arena, int depth) const {
		Vector3f wo = isect.wo, wi;
		Float pdf;
		const Point3f &p = isect.p;
		const Normal3f &ns = isect.shading.n;
		const BSDF &bsdf = *isect.bsdf;
		Spectrum f = bsdf.Sample_f(wo, &wi, sampler.Get2D(), &pdf,
			BxDFType(BSDF_TRANSMISSION | BSDF_SPECULAR));
		Spectrum L = Spectrum(0.f);
		if (pdf > 0.f && !f.IsBlack() && AbsDot(wi, ns) != 0.f) {
			// Compute ray differential _rd_ for specular transmission
			RayDifferential rd = isect.SpawnRay(wi);
			if (ray.hasDifferentials) {
				rd.hasDifferentials = true;
				rd.rxOrigin = p + isect.dpdx;
				rd.ryOrigin = p + isect.dpdy;

				Float eta = bsdf.eta;
				Vector3f w = -wo;
				if (Dot(wo, ns) < 0) eta = 1.f / eta;

				Normal3f dndx = isect.shading.dndu * isect.dudx +
					isect.shading.dndv * isect.dvdx;
				Normal3f dndy = isect.shading.dndu * isect.dudy +
					isect.shading.dndv * isect.dvdy;

				Vector3f dwodx = -ray.rxDirection - wo,
					dwody = -ray.ryDirection - wo;
				Float dDNdx = Dot(dwodx, ns) + Dot(wo, dndx);
				Float dDNdy = Dot(dwody, ns) + Dot(wo, dndy);

				Float mu = eta * Dot(w, ns) - Dot(wi, ns);
				Float dmudx =
					(eta - (eta * eta * Dot(w, ns)) / Dot(wi, ns)) * dDNdx;
				Float dmudy =
					(eta - (eta * eta * Dot(w, ns)) / Dot(wi, ns)) * dDNdy;

				rd.rxDirection =
					wi + eta * dwodx - Vector3f(mu * dndx + dmudx * ns);
				rd.ryDirection =
					wi + eta * dwody - Vector3f(mu * dndy + dmudy * ns);
			}
			L = f * Li(rd, scene, sampler, arena, depth + 1) * AbsDot(wi, ns) / pdf;
		}
		return L;
	}

}  // namespace pbrt
