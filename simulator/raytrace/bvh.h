#ifndef __BVH_H__
#define __BVH_H__

#include <vector>
#include "../primitive.h"
#include "../boundingbox.h"

/// Bounding Volume Hierarchy.
/// This is like a KD tree, but the data we store in it
/// are primitives that have a volume, not single points.
/// And theoretically we split on objects rather than space, but eh.
class BVH {
public:
    
    // CONSTRUCTOR & DESTRUCTOR
    static BVH* Build(std::vector<Primitive*>& primitives);
    friend BVH* BuildHelper(std::vector<Primitive*>& primitives, uint start, uint end, int height);

    ~BVH();

    // ACCESSORS
    const BoundingBox& getBoundingBox() const { return bbox; }


    bool isLeaf() const {
        if (child1 == NULL && child2 == NULL) return true;
        assert (child1 != NULL && child2);
        return false;
    }

    const BVH* getChild1() const { return child1; }
    const BVH* getChild2() const { return child2; }
    BVH* getChild1() { return child1; }
    BVH* getChild2() { return child2; }

    const std::vector<Primitive*>& getPrimitives() const {
        return primitives;
    }

    // count the number of primitives we store
    uint primitiveCount() const;

    // Computes the closest intersection.
    bool CastRay(const Ray& ray, Hit& hit) const;

    // for OpenGL rendering
    int triCount() const;

    void packMesh(float*& current) const;


    // Do some invariant checking.
    // Return the number of violations.
    int checkRepr() const;


private:

    BVH(const BoundingBox &_bbox, int _depth=0) : 
        bbox(_bbox), depth(_depth),
        child1(nullptr), child2(nullptr),
        primitives() {} 

    BoundingBox bbox;
    int depth;

    // We own our children
    BVH* child1;
    BVH* child2;

    // We don't own the primitives, and in general they might not be valid
    // After the scene geometry changes
    std::vector<Primitive*> primitives;
};

#endif
