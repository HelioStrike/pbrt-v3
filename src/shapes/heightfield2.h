#ifndef PBRT_SHAPES_HEIGHTFIELD2_H
#define PBRT_SHAPES_HEIGHTFIELD2_H

// shapes/heightfield2.h*
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
    int nx_, ny_, ntris_, nverts_;
    int *verts_;
    const Float *z_;
    std::unique_ptr<int[]> indices_;
    std::unique_ptr<Point3f[]> P_;
    std::unique_ptr<Point2f[]> uvs_;
    Bounds3f bounds_;
};

// Heightfield2 Declarations
std::shared_ptr<Shape> CreateHeightfield2(const Transform *o2w,
                                          const Transform *w2o, bool ro,
                                          const ParamSet &params);

}  // namespace pbrt

#endif  // PBRT_SHAPES_HEIGHTFIELD2_H
