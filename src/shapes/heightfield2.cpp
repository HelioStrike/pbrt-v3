// shapes/heightfield2.cpp*
#include "shapes/heightfield2.h"

#include <memory>

#include "paramset.h"
#include "shapes/triangle.h"

namespace pbrt {

HeightField2::HeightField2(const Transform *ObjectToWorld,
                           const Transform *WorldToObject,
                           bool reverseOrientation, int nx, int ny,
                           const Float *z)
    : Shape(ObjectToWorld, WorldToObject, reverseOrientation),
      nx_(nx),
      ny_(ny),
      z_(z),
      ntris_(2 * (nx - 1) * (ny - 1)),
      nverts_(nx * ny),
      indices_(std::unique_ptr<int[]>(new int[3 * ntris_])),
      P_(std::unique_ptr<Point3f[]>(new Point3f[nx * ny])),
      uvs_(std::unique_ptr<Point2f[]>(new Point2f[nx * ny])) {
    // Compute heightfield2 vertex positions
    int pos = 0;
    for (int y = 0; y < ny; ++y) {
        for (int x = 0; x < nx; ++x) {
            P_[pos].x = uvs_[pos].x = (float)x / (float)(nx - 1);
            P_[pos].y = uvs_[pos].y = (float)y / (float)(ny - 1);
            P_[pos].z = z[pos];
            bounds_ = Union(bounds_, P_[pos]);
            ++pos;
        }
    }
    // bounds_ = Union(bounds_, Point3f(-1000, -1000, -1000));
    // bounds_ = Union(bounds_, Point3f(1000, 1000, 1000));
    // std::cout << "BRUHBRUH: " << bounds_ << '\n';

    // Fill in heightfield2 vertex offset array
    int *vp = indices_.get();
    for (int y = 0; y < ny - 1; ++y) {
        for (int x = 0; x < nx - 1; ++x) {
#define VERT(x, y) ((x) + (y) * nx)
            *vp++ = VERT(x, y);
            *vp++ = VERT(x + 1, y);
            *vp++ = VERT(x + 1, y + 1);

            *vp++ = VERT(x, y);
            *vp++ = VERT(x + 1, y + 1);
            *vp++ = VERT(x, y + 1);
        }
#undef VERT
    }
}

Bounds3f HeightField2::ObjectBound() const { return bounds_; }

// Checks whether `p` is above/on the plane formed by po, p1, and p2.
bool AbovePlane(Point3f po, Point3f p1, Point3f p2, Point3f p) {
    auto v1 = p1 - po;
    auto v2 = p2 - po;
    auto n = Normalize(Cross(v1, v2));
    auto d = -Dot(po, n);
    return Dot(p, n) + d >= 0;
}

bool HeightField2::Intersect(const Ray &ray, Float *tHit,
                             SurfaceInteraction *isect,
                             bool testAlphaTexture) const {
    auto ray_to_object = (*WorldToObject)(ray);
    for (int i = 0; i < ntris_; i++) {
        auto i0 = indices_[3 * i];
        auto i1 = indices_[3 * i + 1];
        auto i2 = indices_[3 * i + 2];
        auto p0 = P_[i0];
        auto p1 = P_[i1];
        auto p2 = P_[i2];

        auto normal = Normalize(Cross(p1 - p0, p2 - p0));
        auto d = -Dot(p0, normal);
        auto time_hit = Float(-(Dot(ray_to_object.o, normal) + d) /
                              Dot(ray_to_object.d, normal));
        auto intersection_point = ray_to_object.o + ray_to_object.d * time_hit;

        auto ap0 = AbovePlane(ray_to_object.o, p0, p1, intersection_point);
        auto ap1 = AbovePlane(ray_to_object.o, p1, p2, intersection_point);
        auto ap2 = AbovePlane(ray_to_object.o, p2, p0, intersection_point);
        if (ap0 == ap1 && ap1 == ap2) {
            Vector2f duv02 = uvs_[i0] - uvs_[i2], duv12 = uvs_[i1] - uvs_[i2];
            Vector3f dp02 = p0 - p2, dp12 = p1 - p2;
            Float determinant = duv02[0] * duv12[1] - duv02[1] * duv12[0];
            bool degenerateUV = std::abs(determinant) < 1e-8;

            Vector3f dpdu, dpdv;

            if (degenerateUV || Cross(dpdu, dpdv).LengthSquared() == 0) {
                if (normal.LengthSquared() == 0) return false;
                CoordinateSystem(normal, &dpdu, &dpdv);
            } else {
                Float invdet = 1 / determinant;
                dpdu = (duv12[1] * dp02 - duv02[1] * dp12) * invdet;
                dpdv = (-duv12[0] * dp02 + duv02[0] * dp12) * invdet;
            }

            if (tHit != nullptr) {
                *tHit = time_hit;
            }
            if (isect != nullptr) {
                *isect = SurfaceInteraction(
                    /*p=*/(*ObjectToWorld)(intersection_point),
                    /*pError=*/Vector3f(0, 0, 0),
                    /*uv=*/Point2f(intersection_point.x, intersection_point.y),
                    /*wo=*/-ray.d,
                    /*dpdu=*/(*ObjectToWorld)(dpdu),
                    /*dpdv=*/(*ObjectToWorld)(dpdv),
                    /*dndu=*/Normal3f(0, 0, 0),
                    /*dndv=*/Normal3f(0, 0, 0),
                    /*time=*/time_hit,
                    /*shape=*/this, /*faceIndex=*/0);
            }
            return true;
        }
    }
    return false;
}

Float HeightField2::Area() const {
    Float a = 0.0f;
    for (int i = 0; i < nverts_; i++) {
        const Point3f &p0 = P_[indices_[3 * i]];
        const Point3f &p1 = P_[indices_[3 * i + 1]];
        const Point3f &p2 = P_[indices_[3 * i + 2]];
        a += 0.5 * Cross(p1 - p0, p2 - p0).Length();
    }
    return a;
}

Interaction HeightField2::Sample(const Point2f &u, Float *pdf) const {
    return Interaction();
}

// Heightfield2 Definitions
std::shared_ptr<Shape> CreateHeightfield2(const Transform *ObjectToWorld,
                                          const Transform *WorldToObject,
                                          bool reverseOrientation,
                                          const ParamSet &params) {
    int nx = params.FindOneInt("nu", -1);
    int ny = params.FindOneInt("nv", -1);
    int nitems;
    const Float *z = params.FindFloat("Pz", &nitems);
    CHECK_EQ(nitems, nx * ny);
    CHECK(nx != -1 && ny != -1 && z != nullptr);

    return std::make_shared<HeightField2>(ObjectToWorld, WorldToObject,
                                          reverseOrientation, nx, ny, z);
}

}  // namespace pbrt
