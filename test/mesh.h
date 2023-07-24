#ifndef MESH_H
#define MESH_H

#include <vector>
#include "hash.h"
#include "material.h"
#include "vertex.h"
#include "face.h"

class Vertex;
class Edge;
class BoundingBox;
class Face;
class Primitive;
class ArgParser;
class Ray;
class Hit;
class Camera;

enum FACE_TYPE { FACE_TYPE_ORIGINAL, FACE_TYPE_RASTERIZED, FACE_TYPE_SUBDIVIDED, FACE_TYPE_PERSON};

// ======================================================================
// ======================================================================
// A class to store all objects in the scene.  The quad faces of the
// mesh can be subdivided to improve the resolution of the radiosity
// solution.  The original mesh is maintained for efficient occlusion
// testing.

class Mesh {

public:

    // ===============================
    // CONSTRUCTOR & DESTRUCTOR & LOAD
    Mesh(ArgParser *_args) {
        args = _args; 
        bbox = NULL;
    }
    virtual ~Mesh();

    /// Parses model data from the given file, stopping after consuming the first '}' or EOF
    void Load(std::ifstream& file, Material* active_material, const Matrix& transform);

    /// Parses material data from the given file
    Material* LoadMaterial(const std::string& path);

    // ========
    // VERTICES
    int numVertices(bool pro) const {
        if (pro) {
            return person_vertices.size();
        }else {
            return vertices.size();
        }
    }
    Vertex* addVertex(const Vec3f &pos, int index);
    Vertex* Mesh::addPersonVertex(const Vec3f& position, int index);
    // look up vertex by index from original .obj file
    Vertex* getVertex(int i) const {
        assert (i >= 0 && i < numVertices(false));
        return vertices[i];
    }
    Vertex* getPersonVertex(int i) const {
        assert(i >= 0 && i < numVertices(true));
        return person_vertices[i];
    }
    // this creates a relationship between 3 vertices (2 parents, 1 child)
    void setParentsChild(Vertex *p1, Vertex *p2, Vertex *child);
    // this accessor will find a child vertex (if it exists) when given
    // two parent vertices
    Vertex* getChildVertex(Vertex *p1, Vertex *p2) const;

    // =====
    // EDGES
    int numEdges() const { return edges.size(); }
    // this efficiently looks for an edge with the given vertices, using a hash table
    Edge* getEdge(Vertex *a, Vertex *b) const;
    const edgeshashtype& getEdges() const { return edges; }

    // =================
    // ACCESS THE LIGHTS
    std::vector<Face*>& getLights() { return original_lights; }

    // ==================================
    // ACCESS THE QUADS (for ray tracing)
    int numOriginalQuads() const { return original_quads.size(); }
    Face* getOriginalQuad(int i) const {
        assert (i < numOriginalQuads());
        return original_quads[i]; }

    // =======================================
    // ACCESS THE PRIMITIVES (for ray tracing)
    int numPrimitives() const { return primitives.size(); }
    Primitive* getPrimitive(int i) const {
        assert (i >= 0 && i < numPrimitives()); 
        return primitives[i]; }
    // ACCESS THE PRIMITIVES (for radiosity)
    int numRasterizedPrimitiveFaces() const { return rasterized_primitive_faces.size(); }
    Face* getRasterizedPrimitiveFace(int i) const {
        assert (i >= 0 && i < numRasterizedPrimitiveFaces());
        return rasterized_primitive_faces[i]; }

    // ==============================================================
    // ACCESS THE SUBDIVIDED QUADS + RASTERIZED FACES (for radiosity)
    int numFaces() const { return subdivided_quads.size() + rasterized_primitive_faces.size()+person_quads.size(); }
    Face* getFace(int i) const {
        int num_faces = numFaces();
        assert (i >= 0 && i < num_faces);
        if (i < (int)subdivided_quads.size()) return subdivided_quads[i];
        else if (i < (int)(person_quads.size() + subdivided_quads.size())) return(person_quads[(i - subdivided_quads.size())]);
        else return getRasterizedPrimitiveFace(i-subdivided_quads.size()); }

    int triCount() const; 
    void packMesh(float*& current);

    void move_person(Vec3f dir) {
        double totle = abs(dir.x()) + abs(dir.z());
        Vec3f transfer;
        transfer.setx(dir.x() / totle * 0.01);
        transfer.setz(dir.z() / totle * 0.01);
        for (int i = 0; i < person_vertices.size(); i++) {
            person_vertices[i]->move_vertex(transfer);
        }
        if (GLOBAL_args->camera->first_person) {
            GLOBAL_args->camera->camera_position += transfer;
            GLOBAL_args->camera->point_of_interest += transfer;
        }
    }

    // ============================
    // CREATE OR SUBDIVIDE GEOMETRY
    Face* addRasterizedPrimitiveFace(Vertex *a, Vertex *b, Vertex *c, Vertex *d, Material *material) {
        return addFace(a,b,c,d,material,FACE_TYPE_RASTERIZED);
    }
    Face* addpersonQuad(Vertex* a, Vertex* b, Vertex* c, Vertex* d, Material* material) {
        return addFace(a, b, c, d, material, FACE_TYPE_PERSON);
    }
    Face* addOriginalQuad(Vertex *a, Vertex *b, Vertex *c, Vertex *d, Material *material) {
        return addFace(a,b,c,d,material,FACE_TYPE_ORIGINAL);
    }
    Face* addSubdividedQuad(Vertex *a, Vertex *b, Vertex *c, Vertex *d, Material *material) {
        return addFace(a,b,c,d,material,FACE_TYPE_SUBDIVIDED);
    }

    // ===============
    // OTHER ACCESSORS
    BoundingBox* getBoundingBox() const { return bbox; }

    // ===============
    // OTHER FUNCTIONS
    void Subdivision();

private:

    // ==================================================
    // HELPER FUNCTIONS FOR CREATING/SUBDIVIDING GEOMETRY
    Vertex* AddEdgeVertex(Vertex *a, Vertex *b);
    Vertex* AddMidVertex(Vertex *a, Vertex *b, Vertex *c, Vertex *d);
    Face* addFace(Vertex *a, Vertex *b, Vertex *c, Vertex *d, Material *material, enum FACE_TYPE face_type);
    void removeFaceEdges(Face *f);
    void addPrimitive(Primitive *p); 
    bool person;

    // ==============
    // REPRESENTATION
    ArgParser *args;
public:
    std::vector<Material*> materials;
private:

    // the bounding box of all rasterized faces in the scene
    BoundingBox *bbox; 
    BoundingBox* P_bbox;

    // the vertices & edges used by all quads (including rasterized primitives)
    std::vector<Vertex*> vertices;  
    std::vector<Vertex*> person_vertices;
    edgeshashtype edges;
    vphashtype vertex_parents;

    // the quads from the .obj file (need to move)
    std::vector<Face*> person_quads;
    // the quads from the .obj file (before subdivision)
    std::vector<Face*> original_quads;
    // the quads from the .obj file that have non-zero emission value
    std::vector<Face*> original_lights; 
    // all primitives (spheres, rings, fluids, etc.)
    std::vector<Primitive*> primitives;
    // the primitives converted to quads
    std::vector<Face*> rasterized_primitive_faces;
    // the quads from the .obj file after subdivision
    std::vector<Face*> subdivided_quads;
};

// ======================================================================
// ======================================================================

#endif




