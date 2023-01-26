#ifndef SETTINGS_H
#define SETTINGS_H

#include <QMatrix4x4>

#include "shadertypes.h"
#include "mesh/mesh.h"

/**
 * Struct that contains all the settings of the program. Initialised with a
 * number of default values.
 */
typedef struct Settings {
    bool modelLoaded = false;
    bool wireframeMode = true;
    bool tesselationMode = false;
    bool showCpuMesh = true;

    int subDivValue = 0;
    bool isEdgeSelected = false;

    float FoV = 80;
    float dispRatio = 16.0f / 9.0f;
    float rotAngle = 0.0f;

    bool uniformUpdateRequired = true;

    ShaderType currentShader = ShaderType::PHONG;

    QMatrix4x4 modelViewMatrix, projectionMatrix;
    QMatrix3x3 normalMatrix;

    QVector<Mesh> meshes;
    QVector<unsigned int> edgeSlected;
    HalfEdge* selectedHE;
    double intSharpnessSelectedHE;
    double decSharpnessSelectedHE;

} Settings;

#endif  // SETTINGS_H
