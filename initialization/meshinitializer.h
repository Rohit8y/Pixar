#ifndef MESH_INITIALIZER_H
#define MESH_INITIALIZER_H

#include "../mesh/mesh.h"
#include "objfile.h"

/**
 * @brief The MeshInitializer class initializes half-edge meshes from OBJFiles.
 */
class MeshInitializer {
    public:
        MeshInitializer();
        Mesh constructHalfEdgeMesh(const OBJFile& loadedOBJFile);

    private:
        void initGeometry(Mesh& mesh, int numVertices,
                          const QVector<QVector3D>& vertexCoords);
        void initTopology(Mesh& mesh, int numFaces,
                          const QVector<QVector<int>>& faceCoordInd,
                          const QVector<int>& sharpHeInd,
                          const QVector<double>& sharpnessValue);
        void addHalfEdge(Mesh& mesh, int h, Face* face,
                         const QVector<int>& faceIndices, int i, double sharpness);
        void setTwins(Mesh& mesh, int h, int vertIdx1, int vertIdx2);

        void updateSharpnessOfTwinEdges(Mesh &mesh) const;


        QList<QPair<int, int>> edgeList;
        QVector<int> edgeIndices;
};

#endif  // MESH_INITIALIZER_H
