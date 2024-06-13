#ifndef PBRT_CAMERAS_REALISTIC_WITHOUT_APPROXIMATION_H
#define PBRT_CAMERAS_REALISTIC_WITHOUT_APPROXIMATION_H

// cameras/realistic_without_approximation.h*
#include "camera.h"
#include "cameras/realistic_without_approximation.h"
#include "film.h"
#include "lowdiscrepancy.h"
#include "pbrt.h"

namespace pbrt {

struct LensElement {
    float radius, thickness, N, aperture;
};

// RealisticCamera Declarations
class RealisticCameraWithoutApproximation : public Camera {
  public:
    // RealisticCamera Public Methods
    RealisticCameraWithoutApproximation(
        const AnimatedTransform &cam2world, Film *film, const Medium *medium,
        Float shutter_open, Float shutter_close, Float film_distance,
        Float aperture_diameter, Float film_diag, Float hither, Float yon,
        std::vector<LensElement> lensElements);
    Float GenerateRay(const CameraSample &sample, Ray *ray) const;
    bool TraceLensesFromFilm(const Ray &rCamera, Ray *rOut) const;

  private:
    Float LensFrontZ() const {
        Float zSum = 0;
        for (auto &lens : lens_elements_) {
            zSum += lens.thickness;
        }
        return zSum;
    }

    Float LensRearZ() const { return lens_elements_.back().thickness; }

    Float RearElementRadius() const { return lens_elements_.back().radius; }

    bool IntersectSphericalElement(Float radius, Float zCenter, const Ray &ray,
                                   Float *t, Normal3f *n) const {
        // Compute _t0_ and _t1_ for ray--element intersection
        Point3f o = ray.o - Vector3f(0, 0, zCenter);
        Float A = ray.d.x * ray.d.x + ray.d.y * ray.d.y + ray.d.z * ray.d.z;
        Float B = 2 * (ray.d.x * o.x + ray.d.y * o.y + ray.d.z * o.z);
        Float C = o.x * o.x + o.y * o.y + o.z * o.z - radius * radius;
        Float t0, t1;
        if (!Quadratic(A, B, C, &t0, &t1)) return false;

        // Select intersection $t$ based on ray direction and element curvature
        bool useCloserT = (ray.d.z > 0) ^ (radius < 0);
        *t = useCloserT ? std::min(t0, t1) : std::max(t0, t1);
        if (*t < 0) return false;

        // Compute surface normal of element at ray intersection point
        *n = Normal3f(Vector3f(o + *t * ray.d));
        *n = Faceforward(Normalize(*n), -ray.d);
        return true;
    }

    Bounds2f BoundExitPupil(Float pFilmX0, Float pFilmX1) {
        Bounds2f pupilBounds;
        const int nSamples = 1024 * 1024;
        int nExitingRays = 0;

        Float rearRadius = RearElementRadius();
        Bounds2f projRearBounds(Point2f(-1.5f * rearRadius, -1.5f * rearRadius),
                                Point2f(1.5f * rearRadius, 1.5f * rearRadius));

        for (int i = 0; i < nSamples; i++) {
            Point3f pFilm(Lerp((i + 0.5f) / nSamples, pFilmX0, pFilmX1), 0, 0);
            Float u[2] = {RadicalInverse(0, i), RadicalInverse(1, i)};
            Point3f pRear(
                Lerp(u[0], projRearBounds.pMin.x, projRearBounds.pMax.x),
                Lerp(u[1], projRearBounds.pMin.y, projRearBounds.pMax.y),
                LensRearZ());

            if (Inside(Point2f(pRear.x, pRear.y), pupilBounds) ||
                TraceLensesFromFilm(Ray(pFilm, pRear - pFilm), nullptr)) {
                pupilBounds = Union(pupilBounds, Point2f(pRear.x, pRear.y));
                ++nExitingRays;
            }
        }

        if (nExitingRays == 0) return projRearBounds;

        pupilBounds =
            Expand(pupilBounds, 2 * projRearBounds.Diagonal().Length() /
                                    std::sqrt(nSamples));
        return pupilBounds;
    }

    Point3f SampleExitPupil(const Point2f &pFilm, const Point2f &lensSample,
                            Float *sampleBoundsArea) const {
        Float rFilm = std::sqrt(pFilm.x * pFilm.x + pFilm.y * pFilm.y);
        int rIndex = rFilm / (film->diagonal / 2) * exit_pupil_bounds_.size();
        rIndex = std::min((int)exit_pupil_bounds_.size() - 1, rIndex);
        Bounds2f pupilBounds = exit_pupil_bounds_[rIndex];
        if (sampleBoundsArea) *sampleBoundsArea = pupilBounds.Area();

        Point2f pLens = pupilBounds.Lerp(lensSample);

        Float sinTheta = (rFilm != 0) ? pFilm.y / rFilm : 0;
        Float cosTheta = (rFilm != 0) ? pFilm.x / rFilm : 1;
        return Point3f(cosTheta * pLens.x - sinTheta * pLens.y,
                       sinTheta * pLens.x + cosTheta * pLens.y, LensRearZ());
    }

    Float film_distance_;
    std::vector<LensElement> lens_elements_;
    std::vector<Bounds2f> exit_pupil_bounds_;
};

RealisticCameraWithoutApproximation *CreateRealisticCameraWithoutApproximation(
    const ParamSet &params, const AnimatedTransform &cam2world, Film *film,
    const Medium *medium);

}  // namespace pbrt

#endif  // PBRT_CAMERAS_REALISTIC_WITHOUT_APPROXIMATION_H
