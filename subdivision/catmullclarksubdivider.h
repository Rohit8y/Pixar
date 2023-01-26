#ifndef CATMULL_CLARK_SUBDIVIDER_H
#define CATMULL_CLARK_SUBDIVIDER_H

#include "mesh/mesh.h"
#include "subdivider.h"

/**
 * @brief The CatmullClarkSubdivider class is a subdivider class that performs
 * Catmull-Clark subdivision on meshes.
 */
class CatmullClarkSubdivider : public Subdivider {
    public:
        CatmullClarkSubdivider();
        Mesh subdivide(Mesh& mesh) const override;

    private:
        void reserveSizes(Mesh& mesh, Mesh& newMesh) const;
        void geometryRefinement(Mesh& mesh, Mesh& newMesh) const;
        void topologyRefinement(Mesh& mesh, Mesh& newMesh) const;

        void setHalfEdgeData(Mesh& newMesh, int h, int edgeIdx, int vertIdx,
                             int twinIdx, double sharpness) const;

        int getAdjacentSharpEdgesValence(const Vertex &vertex) const;
        QVector<QVector3D> getAdjacentSharpEdgesMidpoints(const Vertex &vertex) const;
        float getAvgSharpness(const Vertex &vertex) const;

        QVector3D facePoint(const Face& face) const;

        QVector3D edgePoint(const HalfEdge& edge) const;
        QVector3D sharpEdgePoint(const HalfEdge &edge) const;
        QVector3D edgeBlend(const HalfEdge &edge) const;
        QVector3D boundaryEdgePoint(const HalfEdge& edge) const;

        QVector3D vertexPoint(const Vertex& vertex) const;
        QVector3D vertexCorner(const Vertex& vertex) const;
        QVector3D vertexCrease(const Vertex &vertex, QVector<QVector3D> sharpEdges_midpoints) const;
        QVector3D vertexBlend(const QVector3D &corner_vertex, const QVector3D &crease_vertex ) const;
        QVector3D boundaryVertexPoint(const Vertex& vertex) const;
        void updateSharpnessOfTwinEdges(Mesh &newMesh) const;
};

#endif  // CATMULL_CLARK_SUBDIVIDER_H
