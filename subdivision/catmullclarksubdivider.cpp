#include "catmullclarksubdivider.h"

#include <QDebug>

/**
 * @brief CatmullClarkSubdivider::CatmullClarkSubdivider Creates a new empty
 * Catmull Clark subdivider.
 */
CatmullClarkSubdivider::CatmullClarkSubdivider() {}

/**
 * @brief CatmullClarkSubdivider::subdivide Subdivides the provided control mesh
 * and returns the subdivided mesh. Performs just a single subdivision step. The
 * subdivision follows the indexing rules of this paper:
 * https://diglib.eg.org/bitstream/handle/10.1111/cgf14381/v40i8pp057-070.pdf?sequence=1&isAllowed=y
 * The descriptions of the face, edge and vertex points are also quoted from
 * this paper.
 * @param controlMesh The mesh to be subdivided.
 * @return The mesh resulting of applying a single subdivision step on the
 * control mesh.
 */
Mesh CatmullClarkSubdivider::subdivide(Mesh &mesh) const {
    Mesh newMesh;
    reserveSizes(mesh, newMesh);
    geometryRefinement(mesh, newMesh);
    topologyRefinement(mesh, newMesh);
    return newMesh;
}

/**
 * @brief CatmullClarkSubdivider::reserveSizes Resizes the vertex, half-edge and
 * face vectors. Aslo recalculates the edge count.
 * @param controlMesh The control mesh.
 * @param newMesh The new mesh. At this point, the mesh is fully empty.
 */
void CatmullClarkSubdivider::reserveSizes(Mesh &controlMesh,
                                            Mesh &newMesh) const {
    int newNumEdges = 2 * controlMesh.numEdges() + controlMesh.numHalfEdges();
    int newNumFaces = controlMesh.numHalfEdges();
    int newNumHalfEdges = controlMesh.numHalfEdges() * 4;
    int newNumVerts =
        controlMesh.numVerts() + controlMesh.numFaces() + controlMesh.numEdges();

    newMesh.getVertices().resize(newNumVerts);
    newMesh.getHalfEdges().resize(newNumHalfEdges);
    newMesh.getFaces().resize(newNumFaces);
    newMesh.edgeCount = newNumEdges;
}

/**
 * @brief CatmullClarkSubdivider::geometryRefinement Performs the geometry
 * refinement. In other words, it calculates the coordinates of the vertex, edge
 * and face points. Note that this also sets the valences.
 * The valence of a new face point is equal to the valence of the face.
 * The valence of a new edge point is always 4, unless the edge point lies on a
 * boundary, in which case the valence will be 3.
 * The valence of a new vertex point is equal to the valence of the original
 * vertex point.
 * @param controlMesh The control mesh.
 * @param newMesh The new mesh. At the start of this function, the only
 * guarantee you have of this newMesh is that the vertex, half-edge and face
 * vectors have the correct sizes.
 */
void CatmullClarkSubdivider::geometryRefinement(Mesh &controlMesh,
                                                Mesh &newMesh) const {
    QVector<Vertex> &newVertices = newMesh.getVertices();
    QVector<Vertex> &vertices = controlMesh.getVertices();
    QVector<Face> &faces = controlMesh.getFaces();

    // Face Points
    for (int f = 0; f < controlMesh.numFaces(); f++) {
        QVector3D coords = facePoint(faces[f]);
        int i = controlMesh.numVerts() + faces[f].index;
        // Face points always inherit the valence of the face
        Vertex facePoint(coords, nullptr, faces[f].valence, i);
        newVertices[i] = facePoint;
    }

    // Edge Points
    QVector<HalfEdge> &halfEdges = controlMesh.getHalfEdges();
    for (int h = 0; h < controlMesh.numHalfEdges(); h++) {
        HalfEdge currentEdge = halfEdges[h];
        // Only create a new vertex per set of halfEdges (i.e,. once per undirected edge)
        if (h > currentEdge.twinIdx()) {
            int v = controlMesh.numVerts() + controlMesh.numFaces() + currentEdge.edgeIdx();
            int valence;
            QVector3D coords;
            if (currentEdge.sharpness == 0) {
                // Apply edge smooth rule
                if (currentEdge.isBoundaryEdge()) {
                    coords = boundaryEdgePoint(currentEdge);
                    valence = 3;
                } else {
                    coords = edgePoint(currentEdge);
                    valence = 4;
                }
            } else if (currentEdge.sharpness >= 1) {
                // Apply edge sharp rule
                coords = sharpEdgePoint(currentEdge);
                valence = 4;
            } else if (currentEdge.sharpness < 1) {
                // Apply blend edge sharp rule
                coords = edgeBlend(currentEdge);
                valence = 4;
            }
            newVertices[v] = Vertex(coords, nullptr, valence, v);
        }
    }

    // Vertex Points
    for (int v = 0; v < controlMesh.numVerts(); v++) {
        QVector3D coords;
        // Get the number of adjacent sharp edges to the current vertex
        int sharpEdges_num = getAdjacentSharpEdgesValence(vertices[v]);
        if (sharpEdges_num < 2) {
            // Apply vertex smooth rule
            if (vertices[v].isBoundaryVertex()) {
                coords = boundaryVertexPoint(vertices[v]);
            } else {
                coords = vertexPoint(vertices[v]);
            }
        } else if (sharpEdges_num > 2) {
            // Apply corner vertex sharp rule
            coords = vertexCorner(vertices[v]);
        } else if (sharpEdges_num == 2) {
            // Get the average sharpness value of the two adjacent sharp edges to the current vertex
            float avg_sharpness = getAvgSharpness(vertices[v]);
            // Get the midpoints of the two adjacent sharp edges to the current vertex
            QVector<QVector3D> sharpEdges_midpoints = getAdjacentSharpEdgesMidpoints(vertices[v]);
            if (avg_sharpness >= 1.0) {
                // Apply crease vertex sharp rule
                coords = vertexCrease(vertices[v], sharpEdges_midpoints);
            } else {
                // Apply corner vertex sharp rule
                QVector3D corner_coords = vertexCorner(vertices[v]);
                // Apply crease vertex sharp rule
                QVector3D crease_coords = vertexCrease(vertices[v], sharpEdges_midpoints);
                // Apply blend vertex sharp rule
                coords = vertexBlend(corner_coords, crease_coords);
            }
        }
        newVertices[v] = Vertex(coords, nullptr, vertices[v].valence, v);
    }
}

/**
 * @brief CatmullClarkSubdivider::vertexPoint Calculates the new position of the
 * provided vertex. It does so according to the formula for smooth vertex
 * points:
 *
 * Q/n + 2R/n + S(n-3)/n
 *
 * where
 *
 * Q = the average of the new face points of all faces adjacent to the old
 * vertex point.
 * R = the average of the midpoints of all edges incident to the
 * old vertex point.
 * S = old vertex point.
 * n = valence of the vertex.
 *
 * @param vertex The vertex to calculate the new position of. Note that this
 * vertex is the vertex from the control mesh.
 * @return The coordinates of the new vertex point.
 */
QVector3D CatmullClarkSubdivider::vertexPoint(const Vertex &vertex) const {
    HalfEdge *edge = vertex.out;
    QVector3D R;  // average of edge mid points
    QVector3D Q;  // average of face points
    for (int i = 0; i < vertex.valence; i++) {
        R += (edge->origin->coords + edge->next->origin->coords) / 2.0;
        Q += facePoint(*edge->face);
        edge = edge->prev->twin;
    }
    float n = float(vertex.valence);
    Q /= n;
    R /= n;
    // See Equation 1 of the aforementioned paper
    return (Q + 2 * R + (vertex.coords * (n - 3.0f))) / n;
}

/**
 * @brief CatmullClarkSubdivider::vertexCorner Calculates the new position of the
 * provided vertex if it is a corner vertex. A vertex that is considered as corner
 * remains unchanged, i.e., does not move during subdivision.
 *
 * @param vertex The vertex to calculate the new position of. Note that this
 * vertex is the vertex from the control mesh.
 * @return The coordinates of the new vertex point.
 */
QVector3D CatmullClarkSubdivider::vertexCorner(const Vertex &vertex) const {
    return vertex.coords;
}

/**
 * @brief CatmullClarkSubdivider::vertexCrease Calculates the new position of a
 * vertex if it is a crease vertex. It does so according to the formula for crease vertex
 * points:
 *
 * (e_1 + 6S + e_2) / 8
 *
 * where
 *
 * e_1 = the midpoint of the first sharp edge adjacent to the provided vertex
 * e_2 = the midpoint of the second sharp edge adjacent to the provided vertex
 * S = old vertex point
 *
 * @param vertex The vertex to calculate the new position of. Note that this
 * vertex is the vertex from the control mesh.
 * @param sharpEdges_midpoints The midpoints of the two sharp edges adjacent to the
 * provided vertex.
 * @return The coordinates of the new vertex point.
 */
QVector3D CatmullClarkSubdivider::vertexCrease(const Vertex &vertex, QVector<QVector3D> sharpEdges_midpoints) const {
    QVector3D edgeMidpoint1 = sharpEdges_midpoints[0];
    QVector3D edgeMidpoint2 = sharpEdges_midpoints[1];

    return (edgeMidpoint1 + 6.0 * vertex.coords + edgeMidpoint2) / 8.0;
}

/**
 * @brief CatmullClarkSubdivider::vertexBlend Calculates the new position of a
 * vertex if it is a blend vertex. It does so according to the formula for blend vertex
 * points:
 *
 * (Vcrease + Vcorner) / 2
 *
 * where
 *
 * Vcrease = the coordinates of the vertex after applying the crease vertex rule
 * Vcorner = the coordinates of the vertex after applying the corner vertex rule
 *
 * @param vertex The vertex to calculate the new position of. Note that this
 * vertex is the vertex from the control mesh.
 * @return The coordinates of the new vertex point.
 */
QVector3D CatmullClarkSubdivider::vertexBlend(const QVector3D &cornerVertex, const QVector3D &creaseVertex) const {
    return (cornerVertex + creaseVertex) / 2.0;
}

/**
 * @brief CatmullClarkSubdivider::boundaryVertexPoint Calculates the new
 * position of the provided vertex. It does so according to the formula for
 * boundary vertex points:
 *
 * (R + S) / 2
 *
 * where
 *
 * R = the average of the midpoints of all edges incident to the
 * old vertex point.
 * S = old vertex point.
 *
 * @param vertex The vertex to calculate the new position of. Note that this
 * vertex is the vertex from the control mesh.
 * @return The coordinates of the new boundary vertex point.
 */
QVector3D CatmullClarkSubdivider::boundaryVertexPoint(const Vertex &vertex) const {
    QVector3D boundPoint = vertex.coords * 2;
    boundPoint += boundaryEdgePoint(*vertex.nextBoundaryHalfEdge());
    boundPoint += boundaryEdgePoint(*vertex.prevBoundaryHalfEdge());

    return boundPoint / 4.0;
}

/**
 * @brief CatmullClarkSubdivider::edgePoint Calculates the position of the edge
 * point according to the formula for smooth edge points:
 *
 * (M + Q) / 2
 *
 * where
 *
 * Q = the average of the new face points of the two faces adjacent to the old
 * edge.
 * M = the midpoint of the edge.
 *
 * @param edge One of the half-edges that lives on the edge to calculate
 * the edge point. Note that this half-edge is the half-edge from the control
 * mesh.
 * @return The coordinates of the new edge point.
 */
QVector3D CatmullClarkSubdivider::edgePoint(const HalfEdge &edge) const {
    QVector3D edgePt = boundaryEdgePoint(edge);
    edgePt += (facePoint(*edge.face) + facePoint(*edge.twin->face)) / 2.0;
    return edgePt /= 2.0;
}

QVector3D CatmullClarkSubdivider::sharpEdgePoint(const HalfEdge &edge) const {
    return boundaryEdgePoint(edge);
}

QVector3D CatmullClarkSubdivider::edgeBlend(const HalfEdge &edge) const {
    float sharpness = edge.sharpness;
    QVector3D smoothCoord;
    QVector3D sharpCoord;

    if (edge.isBoundaryEdge()) {
        smoothCoord = boundaryEdgePoint(edge);
    } else {
        smoothCoord = edgePoint(edge);
    }
    sharpCoord = sharpEdgePoint(edge);

    return (1 - sharpness) * smoothCoord + (sharpness) * sharpCoord;
}

/**
 * @brief CatmullClarkSubdivider::boundaryEdgePoint Calculates the position of
 * the boundary edge point by taking the midpoint of the edge.
 * @param edge One of the half-edges that lives on the edge to calculate
 * the edge point. Note that this half-edge is the half-edge from the control
 * mesh.
 * @return The coordinates of the new boundary edge point.
 */
QVector3D CatmullClarkSubdivider::boundaryEdgePoint(const HalfEdge &edge) const {

    return (edge.origin->coords + edge.next->origin->coords) / 2.0f;
}

/**
 * @brief CatmullClarkSubdivider::facePoint Calculates the position of the face
 * point by averaging the positions of all vertices adjacent to the face.
 * @param face The face to calculate the face point of. Note that this face is
 * the face from the control mesh.
 * @return The coordinates of the new face point.
 */
QVector3D CatmullClarkSubdivider::facePoint(const Face &face) const {
    QVector3D edgePt;
    HalfEdge *edge = face.side;
    for (int side = 0; side < face.valence; side++) {
        edgePt += edge->origin->coords;
        edge = edge->next;
    }
    return edgePt / face.valence;
}

/**
 * @brief CatmullClarkSubdivider::getAdjacentSharpEdgesValence Calculates the number of
 * sharp edges adjacent to the provided vertex.
 * @param Vertex The vertex to calculate the number of its adjacent sharp edges.
 * @return The number of adjacent sharp edges to the input vertex.
 */
/*int CatmullClarkSubdivider::getAdjacentSharpEdgesValence(const Vertex &vertex) const {
    int sharpCounter = 0;
    HalfEdge *currentEdge = vertex.out;
    for (int i = 0; i < vertex.valence; i++) {
        if (currentEdge->sharpness > 0) {
            sharpCounter++;
        }
        currentEdge = currentEdge->twin->next;
    }
    return sharpCounter;
}*/

int CatmullClarkSubdivider::getAdjacentSharpEdgesValence(const Vertex &vertex) const {
    int sharpCounter = 0;
    HalfEdge *currentEdge = vertex.out;
    for (int i = 0; i < vertex.valence; i++) {
        if (!currentEdge->isBoundaryEdge()){
            currentEdge = currentEdge->twin;
        }
        if (currentEdge->sharpness > 0) {
            sharpCounter++;
        }
        currentEdge = currentEdge->prev;
        if (currentEdge->sharpness > 0) {
            sharpCounter++;
        }
    }

    return sharpCounter;
}

/**
 * @brief CatmullClarkSubdivider::getAdjacentSharpEdgesMidpoints Extracts the midpoints of
 * all the adjacent sharp edges to the provided vertex.
 * @param Vertex The vertex to get its adjacent sharp edges.
 * @return The midpoints of the sharp edges adjacent to the provided vertex.
 */
QVector<QVector3D> CatmullClarkSubdivider::getAdjacentSharpEdgesMidpoints(const Vertex &vertex) const {
    QVector<QVector3D> sharpEdges_midpoints;
    HalfEdge *currentEdge = vertex.out;
    QVector3D sharpEdge_midpoint;
    for (int i = 0; i < vertex.valence; i++) {
        if (currentEdge->sharpness > 0) {
            sharpEdge_midpoint = (currentEdge->origin->coords + currentEdge->next->origin->coords) / 2.0;
            sharpEdges_midpoints.push_back(sharpEdge_midpoint);
        }
        currentEdge = currentEdge->twin->next;
    }
    return sharpEdges_midpoints;
}

/**
 * @brief CatmullClarkSubdivider::getAvgSharpness Calculates the average sharpness value
 * of the two adjacent sharp edges to the provided vertex.
 * @param Vertex The vertex to calculate the average sharpness value of its two adjacent sharp edges.
 * @return The average sharpness value of the sharp edges adjacent to the provided vertex.
 */
float CatmullClarkSubdivider::getAvgSharpness(const Vertex &vertex) const {
    HalfEdge *currentEdge = vertex.out;
    float avg_sharpness = 0.0;
    for (int i = 0; i < vertex.valence; i++) {
        if (currentEdge->sharpness > 0) {
            avg_sharpness += currentEdge->sharpness;
        }
        currentEdge = currentEdge->twin->next;
    }

    return avg_sharpness /= 2.0;
}

/**
 * @brief CatmullClarkSubdivider::topologyRefinement Performs the topology
 * refinement. Every face is split into n new faces, where n is the valence of
 * the original face. Newly generated faces are all quads.
 * @param controlMesh The control mesh.
 * @param newMesh The new mesh.
 */
void CatmullClarkSubdivider::topologyRefinement(Mesh &controlMesh,
                                                Mesh &newMesh) const {
    for (int f = 0; f < newMesh.numFaces(); ++f) {
        newMesh.faces[f].index = f;
        newMesh.faces[f].valence = 4;
    }

    // Split halfedges
    for (int h = 0; h < controlMesh.numHalfEdges(); ++h) {
        HalfEdge *edge = &controlMesh.halfEdges[h];

        int h1 = 4 * h;
        int h2 = 4 * h + 1;
        int h3 = 4 * h + 2;
        int h4 = 4 * h + 3;

        int twinIdx1 = edge->twinIdx() < 0 ? -1 : 4 * edge->twin->next->index + 3;
        int twinIdx2 = 4 * edge->next->index + 2;
        int twinIdx3 = 4 * edge->prev->index + 1;
        int twinIdx4 = 4 * edge->prev->twinIdx();

        int vertIdx1 = edge->origin->index;
        int vertIdx2 =
            controlMesh.numVerts() + controlMesh.numFaces() + edge->edgeIndex;
        int vertIdx3 = controlMesh.numVerts() + edge->faceIdx();
        int vertIdx4 =
            controlMesh.numVerts() + controlMesh.numFaces() + edge->prev->edgeIndex;

        int edgeIdx1 = 2 * edge->edgeIndex + (h > edge->twinIdx() ? 0 : 1);
        int edgeIdx2 = 2 * controlMesh.numEdges() + h;
        int edgeIdx3 = 2 * controlMesh.numEdges() + edge->prev->index;
        int edgeIdx4 = 2 * edge->prev->edgeIndex +
                       (edge->prevIdx() > edge->prev->twinIdx() ? 1 : 0);

        setHalfEdgeData(newMesh, h1, edgeIdx1, vertIdx1, twinIdx1);
        setHalfEdgeData(newMesh, h2, edgeIdx2, vertIdx2, twinIdx2);
        setHalfEdgeData(newMesh, h3, edgeIdx3, vertIdx3, twinIdx3);
        setHalfEdgeData(newMesh, h4, edgeIdx4, vertIdx4, twinIdx4);
    }
}

/**
 * @brief LoopSubdivider::setHalfEdgeData Sets the data of a single half-edge
 * (and the corresponding vertex and face)
 * @param newMesh The new mesh this half-edge will live in.
 * @param h Index of the half-edge.
 * @param edgeIdx Index of the (undirected) edge this half-edge will belong to.
 * @param vertIdx Index of the vertex that this half-edge will originate from.
 * @param twinIdx Index of the twin of this half-edge. -1 if the half-edge lies
 * on a boundary.
 */
void CatmullClarkSubdivider::setHalfEdgeData(Mesh &newMesh, int h, int edgeIdx,
                                                int vertIdx, int twinIdx) const {
    HalfEdge *halfEdge = &newMesh.halfEdges[h];

    halfEdge->edgeIndex = edgeIdx;
    halfEdge->index = h;
    halfEdge->origin = &newMesh.vertices[vertIdx];
    halfEdge->face = &newMesh.faces[halfEdge->faceIdx()];
    halfEdge->next = &newMesh.halfEdges[halfEdge->nextIdx()];
    halfEdge->prev = &newMesh.halfEdges[halfEdge->prevIdx()];
    halfEdge->twin = twinIdx < 0 ? nullptr : &newMesh.halfEdges[twinIdx];

    halfEdge->origin->out = halfEdge;
    halfEdge->origin->index = vertIdx;
    halfEdge->face->side = halfEdge;
}
