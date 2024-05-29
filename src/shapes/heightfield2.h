#ifndef PBRT_SHAPES_HEIGHTFIELD2_H
#define PBRT_SHAPES_HEIGHTFIELD2_H

// shapes/heightfield2.h*
#include "pbrt.h"
#include "shape.h"

namespace pbrt {

// HeightField2 Declarations
class HeightField2 : public Shape {
  public:
    // HeightField2 Public Methods
    HeightField2(const Transform *ObjectToWorld, const Transform *WorldToObject,
                 bool reverseOrientation, int nx, int ny, const Float *z);
    Bounds3f ObjectBound() const;
    bool Intersect(const Ray &ray, Float *tHit, SurfaceInteraction *isect,
                   bool testAlphaTexture) const;
    Float Area() const;
    Interaction Sample(const Point2f &u, Float *pdf) const;

  private:
    bool IntersectWithTriangle(int index, Ray rayToObject, Float *tHit,
                               SurfaceInteraction *isect) const;
    int posToVoxel(const Point3f &P, int axis) const {
        int v = int((P[axis] - bounds_.pMin[axis]) * invWidth[axis]);
        return Clamp(v, 0, nVoxels[axis] - 1);
    }
    float voxelToPos(int p, int axis) const {
        return bounds_.pMin[axis] + p * width[axis];
    }
    inline int offset(int x, int y, int z) const {
        return z * nVoxels[0] * nVoxels[1] + y * nVoxels[0] + x;
    }

    int nx_, ny_, ntris_, nverts_;
    int nVoxels[3];
    int *verts_;
    // Width per voxel.
    Float width[3] = {1, 1, 1}, invWidth[3] = {1, 1, 1};
    Float zMin = 10000000, zMax = -10000000;
    const Float *z_;
    std::unique_ptr<int[]> indices_;
    std::unique_ptr<Point3f[]> P_;
    std::unique_ptr<Normal3f[]> FN_;
    std::unique_ptr<Normal3f[]> PN_;
    std::unique_ptr<Point2f[]> uvs_;
    Bounds3f bounds_;
};

// Heightfield2 Declarations
std::shared_ptr<Shape> CreateHeightfield2(const Transform *o2w,
                                          const Transform *w2o, bool ro,
                                          const ParamSet &params);

}  // namespace pbrt

#endif  // PBRT_SHAPES_HEIGHTFIELD2_H
