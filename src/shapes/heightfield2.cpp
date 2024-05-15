// shapes/heightfield2.cpp*
#include "shapes/heightfield2.h"

#include <math.h>

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
      FN_(std::unique_ptr<Normal3f[]>(new Normal3f[ntris_])),
      PN_(std::unique_ptr<Normal3f[]>(new Normal3f[nx * ny])),
      PN_div_(std::unique_ptr<int[]>(new int[nx * ny])),
      uvs_(std::unique_ptr<Point2f[]>(new Point2f[nx * ny])) {
    // Compute heightfield2 vertex positions
    int pos = 0;
    for (int y = 0; y < ny; ++y) {
        for (int x = 0; x < nx; ++x) {
            P_[pos].x = uvs_[pos].x = (float)x / (float)(nx - 1);
            P_[pos].y = uvs_[pos].y = (float)y / (float)(ny - 1);
            P_[pos].z = z[pos];
            zMin = std::min(zMin, z[pos]);
            zMax = std::max(zMax, z[pos]);
            ++pos;
        }
    }

    // Fill in heightfield2 vertex offset array
    pos = 0;
    for (int y = 0; y < ny - 1; ++y) {
        for (int x = 0; x < nx - 1; ++x) {
#define VERT(x, y) ((x) + (y) * nx)
            indices_[3 * pos] = VERT(x, y);
            indices_[3 * pos + 1] = VERT(x + 1, y);
            indices_[3 * pos + 2] = VERT(x + 1, y + 1);
            FN_[pos] = Normal3f(Normalize(
                Cross(P_[indices_[3 * pos + 2]] - P_[indices_[3 * pos]],
                      P_[indices_[3 * pos + 1]] - P_[indices_[3 * pos]])));
            PN_[VERT(x, y)] += FN_[pos];
            PN_[VERT(x + 1, y)] += FN_[pos];
            PN_[VERT(x + 1, y + 1)] += FN_[pos];
            PN_div_[VERT(x, y)]++, PN_div_[VERT(x + 1, y)]++,
                PN_div_[VERT(x + 1, y + 1)]++;
            pos++;

            indices_[3 * pos] = VERT(x, y);
            indices_[3 * pos + 1] = VERT(x + 1, y + 1);
            indices_[3 * pos + 2] = VERT(x, y + 1);
            FN_[pos] = Normal3f(Normalize(
                Cross(P_[indices_[3 * pos + 2]] - P_[indices_[3 * pos]],
                      P_[indices_[3 * pos + 1]] - P_[indices_[3 * pos]])));
            PN_[VERT(x, y)] += FN_[pos];
            PN_[VERT(x + 1, y + 1)] += FN_[pos];
            PN_[VERT(x, y + 1)] += FN_[pos];
            PN_div_[VERT(x, y)]++, PN_div_[VERT(x, y + 1)]++,
                PN_div_[VERT(x + 1, y + 1)]++;
            pos++;
        }
#undef VERT
    }

    pos = 0;
    for (int y = 0; y < ny; ++y) {
        for (int x = 0; x < nx; ++x) {
            PN_[pos] /= PN_div_[pos];
            ++pos;
        }
    }

    nVoxels[0] = nx_;
    nVoxels[1] = ny_;
    nVoxels[2] = 1;
    width[2] = zMax - zMin;
    invWidth[2] = 1 / (zMax - zMin);
}

Bounds3f HeightField2::ObjectBound() const {
    return Bounds3f(Point3f(0, 0, zMin), Point3f(1, 1, zMax));
}

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
    auto rayToObject = (*WorldToObject)(ray);

    float rayT;
    if (Inside(rayToObject(0), bounds_))
        rayT = 0;
    else if (!bounds_.IntersectP(rayToObject, &rayT)) {
        return false;
    }
    Point3f gridIntersect = rayToObject(rayT);

    // Set up 3D DDA for ray
    float NextCrossingT[3], DeltaT[3];
    int Step[3], Out[3], Pos[3];
    int nVoxels[3] = {nx_, ny_, 1};
    for (int axis = 0; axis < 3; ++axis) {
        if (rayToObject.d[axis] == -0.f) rayToObject.d[axis] = 0.f;
        // Compute current voxel for axis
        Pos[axis] = posToVoxel(gridIntersect, axis);
        if (rayToObject.d[axis] >= 0) {
            // Handle ray with positive direction for voxel stepping
            NextCrossingT[axis] =
                rayT + (voxelToPos(Pos[axis] + 1, axis) - gridIntersect[axis]) /
                           rayToObject.d[axis];
            DeltaT[axis] = 1 / rayToObject.d[axis];
            Step[axis] = 1;
            Out[axis] = nVoxels[axis];
        } else {
            // Handle ray with negative direction for voxel stepping
            NextCrossingT[axis] =
                rayT + (voxelToPos(Pos[axis], axis) - gridIntersect[axis]) /
                           rayToObject.d[axis];
            DeltaT[axis] = -width[axis] / rayToObject.d[axis];
            Step[axis] = -1;
            Out[axis] = -1;
        }
    }

    // Walk ray through voxel grid
    bool hitSomething = false;
    for (;;) {
    }
    return hitSomething;

    for (int i = 0; i < ntris_; i++) {
        auto i0 = indices_[3 * i];
        auto i1 = indices_[3 * i + 1];
        auto i2 = indices_[3 * i + 2];

        // Get triangle vertices in _p0_, _p1_, and _p2_
        const Point3f &p0 = P_[i0];
        const Point3f &p1 = P_[i1];
        const Point3f &p2 = P_[i2];

        const Normal3f &n0 = PN_[i0];
        const Normal3f &n1 = PN_[i1];
        const Normal3f &n2 = PN_[i2];

        // Perform ray--triangle intersection test

        // Transform triangle vertices to ray coordinate space

        // Translate vertices based on ray origin
        Point3f p0t = p0 - Vector3f(rayToObject.o);
        Point3f p1t = p1 - Vector3f(rayToObject.o);
        Point3f p2t = p2 - Vector3f(rayToObject.o);

        // Permute components of triangle vertices and ray direction
        int kz = MaxDimension(Abs(rayToObject.d));
        int kx = kz + 1;
        if (kx == 3) kx = 0;
        int ky = kx + 1;
        if (ky == 3) ky = 0;
        Vector3f d = Permute(rayToObject.d, kx, ky, kz);
        p0t = Permute(p0t, kx, ky, kz);
        p1t = Permute(p1t, kx, ky, kz);
        p2t = Permute(p2t, kx, ky, kz);

        // Apply shear transformation to translated vertex positions
        Float Sx = -d.x / d.z;
        Float Sy = -d.y / d.z;
        Float Sz = 1.f / d.z;
        p0t.x += Sx * p0t.z;
        p0t.y += Sy * p0t.z;
        p1t.x += Sx * p1t.z;
        p1t.y += Sy * p1t.z;
        p2t.x += Sx * p2t.z;
        p2t.y += Sy * p2t.z;

        // Compute edge function coefficients _e0_, _e1_, and _e2_
        Float e0 = p1t.x * p2t.y - p1t.y * p2t.x;
        Float e1 = p2t.x * p0t.y - p2t.y * p0t.x;
        Float e2 = p0t.x * p1t.y - p0t.y * p1t.x;

        // Fall back to double precision test at triangle edges
        if (sizeof(Float) == sizeof(float) &&
            (e0 == 0.0f || e1 == 0.0f || e2 == 0.0f)) {
            double p2txp1ty = (double)p2t.x * (double)p1t.y;
            double p2typ1tx = (double)p2t.y * (double)p1t.x;
            e0 = (float)(p2typ1tx - p2txp1ty);
            double p0txp2ty = (double)p0t.x * (double)p2t.y;
            double p0typ2tx = (double)p0t.y * (double)p2t.x;
            e1 = (float)(p0typ2tx - p0txp2ty);
            double p1txp0ty = (double)p1t.x * (double)p0t.y;
            double p1typ0tx = (double)p1t.y * (double)p0t.x;
            e2 = (float)(p1typ0tx - p1txp0ty);
        }

        if ((e0 < 0 || e1 < 0 || e2 < 0) && (e0 > 0 || e1 > 0 || e2 > 0))
            continue;
        Float det = e0 + e1 + e2;
        if (det == 0) continue;

        // Compute scaled hit distance to triangle and test against ray $t$
        // range
        p0t.z *= Sz;
        p1t.z *= Sz;
        p2t.z *= Sz;
        Float tScaled = e0 * p0t.z + e1 * p1t.z + e2 * p2t.z;
        if (det < 0 && (tScaled >= 0 || tScaled < rayToObject.tMax * det))
            continue;
        else if (det > 0 && (tScaled <= 0 || tScaled > rayToObject.tMax * det))
            continue;

        // Compute barycentric coordinates and $t$ value for triangle
        // intersection
        Float invDet = 1 / det;
        Float b0 = e0 * invDet;
        Float b1 = e1 * invDet;
        Float b2 = e2 * invDet;
        Float t = tScaled * invDet;

        // Ensure that computed triangle $t$ is conservatively greater than
        // zero

        // Compute $\delta_z$ term for triangle $t$ error bounds
        Float maxZt = MaxComponent(Abs(Vector3f(p0t.z, p1t.z, p2t.z)));
        Float deltaZ = gamma(3) * maxZt;

        // Compute $\delta_x$ and $\delta_y$ terms for triangle $t$ error
        // bounds
        Float maxXt = MaxComponent(Abs(Vector3f(p0t.x, p1t.x, p2t.x)));
        Float maxYt = MaxComponent(Abs(Vector3f(p0t.y, p1t.y, p2t.y)));
        Float deltaX = gamma(5) * (maxXt + maxZt);
        Float deltaY = gamma(5) * (maxYt + maxZt);

        // Compute $\delta_e$ term for triangle $t$ error bounds
        Float deltaE =
            2 * (gamma(2) * maxXt * maxYt + deltaY * maxXt + deltaX * maxYt);

        // Compute $\delta_t$ term for triangle $t$ error bounds and check
        // _t_
        Float maxE = MaxComponent(Abs(Vector3f(e0, e1, e2)));
        Float deltaT =
            3 * (gamma(3) * maxE * maxZt + deltaE * maxZt + deltaZ * maxE) *
            std::abs(invDet);
        if (t <= deltaT) continue;

        // Compute triangle partial derivatives
        Vector3f dpdu, dpdv;

        // Compute deltas for triangle partial derivatives
        Vector2f duv02 = uvs_[0] - uvs_[2], duv12 = uvs_[1] - uvs_[2];
        Vector3f dp02 = p0 - p2, dp12 = p1 - p2;
        Float determinant = duv02[0] * duv12[1] - duv02[1] * duv12[0];
        bool degenerateUV = std::abs(determinant) < 1e-8;
        if (!degenerateUV) {
            Float invdet = 1 / determinant;
            dpdu = (duv12[1] * dp02 - duv02[1] * dp12) * invdet;
            dpdv = (-duv12[0] * dp02 + duv02[0] * dp12) * invdet;
        }
        if (degenerateUV || Cross(dpdu, dpdv).LengthSquared() == 0) {
            // Handle zero determinant for triangle partial derivative
            // matrix
            Vector3f ng = Cross(p2 - p0, p1 - p0);
            if (ng.LengthSquared() == 0)
                // The triangle is actually degenerate; the intersection is
                // bogus.
                continue;

            CoordinateSystem(Normalize(ng), &dpdu, &dpdv);
        }

        // Compute error bounds for triangle intersection
        Float xAbsSum =
            (std::abs(b0 * p0.x) + std::abs(b1 * p1.x) + std::abs(b2 * p2.x));
        Float yAbsSum =
            (std::abs(b0 * p0.y) + std::abs(b1 * p1.y) + std::abs(b2 * p2.y));
        Float zAbsSum =
            (std::abs(b0 * p0.z) + std::abs(b1 * p1.z) + std::abs(b2 * p2.z));
        Vector3f pError = gamma(7) * Vector3f(xAbsSum, yAbsSum, zAbsSum);

        // Interpolate $(u,v)$ parametric coordinates and hit point
        Point3f pHit = b0 * p0 + b1 * p1 + b2 * p2;
        Point2f uvHit = b0 * uvs_[0] + b1 * uvs_[1] + b2 * uvs_[2];
        Normal3f normalInterpolated = b0 * n0 + b1 * n1 + b2 * n2;

        // Fill in _SurfaceInteraction_ from triangle hit
        if (isect != nullptr) {
            *isect = (*ObjectToWorld)(
                SurfaceInteraction(pHit, pError, uvHit, -rayToObject.d, dpdu,
                                   dpdv, Normal3f(0, 0, 0), Normal3f(0, 0, 0),
                                   rayToObject.time, this, i));
            isect->n = isect->shading.n = (*ObjectToWorld)(normalInterpolated);
        }
        if (tHit != nullptr) {
            *tHit = t;
        }
        return true;
    }
    return false;
}

Float HeightField2::Area() const {
    Float a = 0.0f;
    for (int i = 0; i < ntris_; i++) {
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
