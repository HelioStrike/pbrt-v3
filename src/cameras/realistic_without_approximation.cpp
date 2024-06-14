#include "realistic_without_approximation.h"

#include <camera.h>
#include <floatfile.h>
#include <paramset.h>
#include <pbrt.h>
#include <reflection.h>

#include <vector>

namespace pbrt {

RealisticCameraWithoutApproximation::RealisticCameraWithoutApproximation(
    const AnimatedTransform &cam2world, Film *film, const Medium *medium,
    Float shutter_open, Float shutter_close, Float film_distance,
    Float aperture_diameter, Float film_diag, Float hither, Float yon,
    std::vector<LensElement> lens_elements)
    : Camera(cam2world, shutterOpen, shutterClose, film, medium),
      film_distance_(film_distance),
      lens_elements_(lens_elements) {
    int nSamples = 64;
    exit_pupil_bounds_.resize(nSamples);
    ParallelFor(
        [&](int i) {
            Float r0 = (Float)i / nSamples * film->diagonal / 2;
            Float r1 = (Float)(i + 1) / nSamples * film->diagonal / 2;
            exit_pupil_bounds_[i] = BoundExitPupil(r0, r1);
        },
        nSamples);
}

Float RealisticCameraWithoutApproximation::GenerateRay(
    const CameraSample &sample, Ray *ray) const {
    Point2f s(sample.pFilm.x / film->fullResolution.x,
              sample.pFilm.y / film->fullResolution.y);
    Point2f pFilm2 = film->GetPhysicalExtent().Lerp(s);
    Point3f pFilm(-pFilm2.x, pFilm2.y, LensFrontZ());

    Float exitPupilBoundsArea;
    Point3f pRear = SampleExitPupil(Point2f(pFilm.x, pFilm.y), sample.pLens,
                                    &exitPupilBoundsArea);
    Ray rFilm(pFilm, pRear - pFilm, Infinity,
              Lerp(sample.time, shutterOpen, shutterClose));
    if (!TraceLensesFromFilm(rFilm, ray)) return 0;

    *ray = CameraToWorld(*ray);
    ray->d = Normalize(ray->d);
    ray->medium = medium;

    Float cosTheta = Normalize(rFilm.d).z;
    Float cos4Theta = (cosTheta * cosTheta) * (cosTheta * cosTheta);
    return (shutterClose - shutterOpen) * (cos4Theta * exitPupilBoundsArea) /
           (LensRearZ() * LensRearZ());
}

bool RealisticCameraWithoutApproximation::TraceLensesFromFilm(
    const Ray &rCamera, Ray *rOut) const {
    Float elementZ = 0;

    static const Transform CameraToLens = Scale(1, 1, -1);
    Ray rLens = CameraToLens(rCamera);

    for (int i = lens_elements_.size() - 1; i >= 0; i--) {
        const LensElement &element = lens_elements_[i];
        elementZ -= element.thickness;

        Float t;
        Normal3f n;
        bool isStop = (element.radius == 0);
        if (isStop) {
            t = (elementZ - rLens.o.z) / rLens.d.z;
        } else {
            Float radius = element.radius;
            Float zCenter = elementZ + element.radius;
            if (!IntersectSphericalElement(radius, zCenter, rLens, &t, &n))
                return false;
        }
        Point3f pHit = rLens(t);
        float r2 = pHit.x * pHit.x + pHit.y * pHit.y;
        if (r2 > element.aperture * element.aperture) return false;
        rLens.o = pHit;

        if (!isStop) {
            Vector3f w;
            Float ni = element.N;
            Float nt = (i > 0 && lens_elements_[i - 1].N != 0)
                           ? lens_elements_[i - 1].N
                           : 1;
            if (!Refract(Normalize(-rLens.d), n, ni / nt, &w)) return false;
            rLens.d = w;
        }
    }
    if (rOut != nullptr) {
        static const Transform LensToCamera = Scale(1, 1, -1);
        *rOut = LensToCamera(rLens);
    }
    return true;
}

RealisticCameraWithoutApproximation *CreateRealisticCameraWithoutApproximation(
    const ParamSet &params, const AnimatedTransform &cam2world, Film *film,
    const Medium *medium) {
    std::string spec_file = params.FindOneFilename("specfile", "");
    if (spec_file.empty()) {
        Error("No lens description file supplied!");
        return nullptr;
    }
    Float film_distance = params.FindOneFloat("filmdistance", 10.f);
    Float aperture_diameter = params.FindOneFloat("aperture_diameter", 0.f);
    Float film_diag = params.FindOneFloat("filmdiag", 50.f);
    Float hither = params.FindOneFloat("hither", 0.1f);
    Float yon = params.FindOneFloat("yon", 10000000.0f);
    Float shutter_open = params.FindOneFloat("shutteropen", 0.f);
    Float shutter_close = params.FindOneFloat("shutterclose", 1.f);
    if (shutter_close < shutter_open) {
        Warning("Shutter close time [%f] < shutter open [%f].  Swapping them.",
                shutter_close, shutter_open);
        std::swap(shutter_close, shutter_open);
    }

    // Load element data from lens description file
    // `element_data` has contains data in the sequence (`radius`, `thickness`,
    // `N`, `aperture`).
    std::vector<Float> element_data;
    if (!ReadFloatFile(spec_file.c_str(), &element_data)) {
        Error("Error reading lens specification file \"%s\".",
              spec_file.c_str());
        return nullptr;
    }
    if (element_data.size() % 4 != 0) {
        Error(
            "Excess values in lens specification file \"%s\"; "
            "must be multiple-of-four values, read %d.",
            spec_file.c_str(), (int)element_data.size());
        return nullptr;
    }
    int number_of_elements = element_data.size() / 4;
    std::vector<LensElement> lens_elements;
    for (int i = 0; i < number_of_elements; i++) {
        lens_elements.push_back({element_data[4 * i], element_data[4 * i + 1],
                                 element_data[4 * i + 2],
                                 element_data[4 * i + 3]});
    }

    return new RealisticCameraWithoutApproximation(
        cam2world, film, medium, shutter_open, shutter_close, film_distance,
        aperture_diameter, film_diag, hither, yon, lens_elements);
}

}  // namespace pbrt
