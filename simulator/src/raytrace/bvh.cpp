#include "bvh.h"
#include "hit.h"

BVH::~BVH() {
    delete child1;
    delete child2;

    // we do not own the primitives !
}

BVH* BuildHelper(std::vector<Primitive*>& primitives, uint start, uint end, int depth) {
    BoundingBox bbox;
    for (uint i = start; i < end; i++) {
        bbox.Extend(primitives[i]->getBoundingBox());
    }
    BVH* node = new BVH(bbox, depth);

    // make a leaf node
    if (end - start <= 2) {
        for (uint i = start; i < end; i++) {
            node->primitives.push_back(primitives[i]);
        }
        return node;
    }

    int split_axis;
    Vec3f size = bbox.getMax() - bbox.getMin();
    if (size.x() >= size.y() && size.x() >= size.z()) {
        split_axis = 0;
    } else if (size.y() >= size.x() && size.y() >= size.z()) {
        split_axis = 1;
    } else {
        split_axis = 2;
    }

    // sort by desired axis
    std::sort(primitives.begin() + start, primitives.begin() + end, [split_axis](Primitive* a, Primitive* b) {
        float a_coord = a->getBoundingBox().getCenter()[split_axis];
        float b_coord = b->getBoundingBox().getCenter()[split_axis];
        return a_coord < b_coord;
    });

    // half of the primitives go to one child, half to the other
    uint mid = (start + end) / 2;
    node->child1 = BuildHelper(primitives, start, mid, depth + 1);
    node->child2 = BuildHelper(primitives, mid, end, depth + 1);

    return node;
}

bool BVH::CastRay(const Ray& ray, Hit& hit) const {
    Hit hit1;

    if (!bbox.intersect(ray, hit1)) {
        // no intersection, we can prune this node
        return false; 
    }
    
    if (isLeaf()) {
        bool ans = false;
        // simply check each primitive in our possession
        // std::cout << "Leaf node: intersecting " << primitives.size() << " primitives\n";
        //std::cout<<"1"<<std::endl;
        for (Primitive* p : primitives) {

            if (p->intersect(ray, hit)) {
                ans = true;
            }
        }
        
        // std::cout << "result: " << (ans? "found hit" : "nothing") << "\n";
        return ans; 
    }
    
    assert (child1 != nullptr);
    assert (child2 != nullptr);

    bool ans = false;

    if (child1->CastRay(ray, hit)) {
        ans = true;
    }
    if (child2->CastRay(ray, hit)) {
        ans = true;
    }
    
    return ans;
}

BVH* BVH::Build(std::vector<Primitive*>& primitives) {
    uint actual_count = primitives.size();
    BVH* bvh = BuildHelper(primitives, 0, primitives.size(), 0);

    assert(primitives.size() == actual_count);
    uint count = bvh->primitiveCount();
    assert(count == actual_count);

    int errs = bvh->checkRepr();
    assert(errs == 0);

    return bvh;
}

uint BVH::primitiveCount() const {
    if (isLeaf()) {
        return primitives.size();
    }
    assert (child1);
    assert (child2);
    return child1->primitiveCount() + child2->primitiveCount();
}
        
    

int BVH::triCount() const {
    int count = BoundingBox::triCount();
    if (child1) {
        count += child1->triCount();
    }
    if (child2) {
        count += child2->triCount();
    }
    return count;
}

void BVH::packMesh(float*& current) const {
    bbox.packMesh(current);
    if (child1) {
        child1->packMesh(current);
    }
    if (child2) {
        child2->packMesh(current);
    }
}

int BVH::checkRepr() const {
    int count = 0;
    if (isLeaf()) {
        for (Primitive* p : primitives) {
            if (!p->getBoundingBox().isSubset(bbox)) {
                std::cerr << "WARNING: primitive is not contained within leaf bbox\n";
                count += 1;
            }
        }
        return count;
    }
    if (!primitives.empty()) {
        std::cerr << "WARNING: inner node contains " << primitives.size() << "primitives\n";
        count += 1;
    }
    if (!child1->bbox.isSubset(bbox)) {
        std::cerr << "WARNING: child bounding box is not a subset of parent\n";
        count += 1;
    }
    if (!child2->bbox.isSubset(bbox)) {
        std::cerr << "WARNING: cihld bounding box is not a subset of parent\n";
        count += 1;
    }
    count += child1->checkRepr();
    count += child2->checkRepr();
    return count;
}
