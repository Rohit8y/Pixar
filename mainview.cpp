#include "mainview.h"

#include <math.h>

#include <QLoggingCategory>
#include <QOpenGLVersionFunctionsFactory>

/**
 * @brief MainView::MainView
 * @param Parent
 */
MainView::MainView(QWidget* Parent) : QOpenGLWidget(Parent), scale(1.0f) {}

/**
 * @brief MainView::~MainView Deconstructs the main view.
 */
MainView::~MainView() {
    debugLogger.stopLogging();
    makeCurrent();
}

/**
 * @brief MainView::initializeGL Initializes the opengl functions and settings,
 * initialises the renderers and sets up the debugger.
 */
void MainView::initializeGL() {
    initializeOpenGLFunctions();
    qDebug() << ":: OpenGL initialized";

    connect(&debugLogger, SIGNAL(messageLogged(QOpenGLDebugMessage)), this,
            SLOT(onMessageLogged(QOpenGLDebugMessage)), Qt::DirectConnection);

    if (debugLogger.initialize()) {
        QLoggingCategory::setFilterRules(
        "qt.*=false\n"
        "qt.text.font.*=false");
        qDebug() << ":: Logging initialized";
        debugLogger.startLogging(QOpenGLDebugLogger::SynchronousLogging);
        debugLogger.enableMessages();
    }

    QString glVersion;
    glVersion = reinterpret_cast<const char*>(glGetString(GL_VERSION));
    qDebug() << ":: Using OpenGL" << qPrintable(glVersion);

    makeCurrent();
    // Enable depth buffer
    glEnable(GL_DEPTH_TEST);
    // Default is GL_LESS
    glDepthFunc(GL_LEQUAL);

    // grab the opengl context
    QOpenGLFunctions_4_1_Core* functions =
        QOpenGLVersionFunctionsFactory::get<QOpenGLFunctions_4_1_Core>(
        this->context());

    // initialize renderers here with the current context
    meshRenderer.init(functions, &settings);

    updateMatrices();
}

/**
 * @brief MainView::resizeGL Handles window resizing.
 * @param newWidth The new width of the window in pixels.
 * @param newHeight The new height of the window in pixels.
 */
void MainView::resizeGL(int newWidth, int newHeight) {
    qDebug() << ".. resizeGL";

    settings.dispRatio = float(newWidth) / float(newHeight);

    settings.projectionMatrix.setToIdentity();
    settings.projectionMatrix.perspective(settings.FoV, settings.dispRatio, 0.1f,
                                            40.0f);
    updateMatrices();
}

/**
 * @brief MainView::updateMatrices Updates the matrices used for the model
 * transforms.
 */
void MainView::updateMatrices() {
    settings.modelViewMatrix.setToIdentity();
    settings.modelViewMatrix.translate(QVector3D(0.0, 0.0, -3.0));
    settings.modelViewMatrix.scale(scale);
    settings.modelViewMatrix.rotate(rotationQuaternion);

    settings.normalMatrix = settings.modelViewMatrix.normalMatrix();

    settings.uniformUpdateRequired = true;

    update();
}

/**
 * @brief MainView::updateBuffers Updates the buffers of the renderers.
 * @param mesh The mesh used to update the buffer content with.
 */
void MainView::updateBuffers(Mesh& mesh) {
    mesh.extractAttributes();
    meshRenderer.updateBuffers(mesh);
    update();
}

/**
 * @brief MainView::paintGL Draw call.
 */
void MainView::paintGL() {
    glClearColor(0.0, 0.0, 0.0, 1.0);
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

    if (settings.wireframeMode) {
        glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
    } else {
        glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
    }

    if (settings.modelLoaded) {
        if (settings.showCpuMesh) {
            meshRenderer.draw();
        }

        if (settings.uniformUpdateRequired) {
            settings.uniformUpdateRequired = false;
        }
    }
}

/**
 * @brief MainView::toNormalizedScreenCoordinates Normalizes the mouse
 * coordinates to screen coordinates.
 * @param x X coordinate.
 * @param y Y coordinate.
 * @return A vector containing the normalized x and y screen coordinates.
 */
QVector2D MainView::toNormalizedScreenCoordinates(float x, float y) {
    float xRatio = x / float(width());
    float yRatio = y / float(height());

    // By default, the drawing canvas is the square [-1,1]^2:
    float xScene = (1 - xRatio) * -1 + xRatio * 1;
    // Note that the origin of the canvas is in the top left corner (not the lower
    // left).
    float yScene = yRatio * -1 + (1 - yRatio) * 1;

    return {xScene, yScene};
}

/**
 * @brief MainView::mouseMoveEvent Handles the dragging and rotating of the mesh
 * by looking at the mouse movement.
 * @param event Mouse event.
 */
void MainView::mouseMoveEvent(QMouseEvent* event) {
    if (event->buttons() == Qt::LeftButton) {
        QVector2D sPos = toNormalizedScreenCoordinates(event->position().x(),
                                                       event->position().y());
        QVector3D newVec = QVector3D(sPos.x(), sPos.y(), 0.0);

        // project onto sphere
        float sqrZ = 1.0f - QVector3D::dotProduct(newVec, newVec);
        if (sqrZ > 0) {
            newVec.setZ(sqrt(sqrZ));
        } else {
            newVec.normalize();
        }

        QVector3D v2 = newVec.normalized();
        // reset if we are starting a drag
        if (!dragging) {
            dragging = true;
            oldVec = newVec;
            return;
        }

        // calculate axis and angle
        QVector3D v1 = oldVec.normalized();
        QVector3D N = QVector3D::crossProduct(v1, v2).normalized();
        if (N.length() == 0.0f) {
            oldVec = newVec;
            return;
        }
        float angle = 180.0f / M_PI * acos(QVector3D::dotProduct(v1, v2));
        rotationQuaternion =
                QQuaternion::fromAxisAndAngle(N, angle) * rotationQuaternion;
        updateMatrices();

        // for next iteration
        oldVec = newVec;
    } else {
        // to reset drag
        dragging = false;
        oldVec = QVector3D();
    }
}

/**
 * @brief MainView::mousePressEvent Handles presses by the mouse. Currently does
 * raycasting.
 * @param event Mouse event.
 */
void MainView::mousePressEvent(QMouseEvent* event) {
    setFocus();
    // Only works when vertex selection is checked
    if (event->buttons() == Qt::LeftButton) {

        settings.isEdgeSelected=true;
        // Get screen space coordinates
        int mouse_x = event->position().x();
        int mouse_y = event->position().y();
        GLfloat depth;
        glReadPixels(mouse_x, height() - 1 - mouse_y,1, 1,
                            GL_LEQUAL, GL_FLOAT, &depth);
        // Get NDC
        QVector3D ray_nds = toNormalizedDeviceCoordinates(mouse_x,mouse_y);
        qDebug()<<ray_nds;

        // Clipping
        QVector4D ray_clip = QVector4D(ray_nds.x(),ray_nds.y(), -1.0 , 1.0);

        // Eye
        QVector4D ray_eye = settings.projectionMatrix.inverted() * ray_clip;
        //  float w =   1;
        QVector4D ray_eye_view = QVector4D(ray_eye.x(),ray_eye.y(),-1.0,0.0);

        // World
        QVector4D ray_eye_inverted = settings.modelViewMatrix.inverted() * ray_eye_view;
        QVector3D ray_wor = QVector3D(ray_eye_inverted.x(),ray_eye_inverted.y(),ray_eye_inverted.z());
        ray_wor = ray_wor.normalized();

        // Find closest point in the current mesh
        findClosestHalfEdge(ray_wor,0.5f);

        updateMatrices();
        updateBuffers(settings.meshes[settings.subDivValue]);
        update();
    }
}

/**
 * @brief MainView::wheelEvent Handles zooming of the view.
 * @param event Mouse event.
 */
void MainView::wheelEvent(QWheelEvent* event) {
    // Delta is usually 120
    float phi = 1.0f + (event->angleDelta().y() / 2000.0f);
    scale = fmin(fmax(phi * scale, 0.01f), 100.0f);
    updateMatrices();
}

/**
 * @brief MainView::keyPressEvent Handles keyboard shortcuts. Currently support
 * 'Z' for wireframe mode and 'R' to reset orientation.
 * @param event Mouse event.
 */
void MainView::keyPressEvent(QKeyEvent* event) {
    switch (event->key()) {
        case 'Z':
            settings.wireframeMode = !settings.wireframeMode;
            update();
            break;
        case 'R':
            scale = 1.0f;
            rotationQuaternion = QQuaternion();
            updateMatrices();
            update();
            break;
    }
}

/**
 * @brief MainView::onMessageLogged Helper function for logging messages.
 * @param message The message to log.
 */
void MainView::onMessageLogged(QOpenGLDebugMessage Message) {
    qDebug() << " → Log:" << Message;
}

/**
 * @brief MainView::toNormalizedDeviceCoordinates Transforms the screen corrdinates to
 * 3D normalised device coordinates.
 * @param mouse_x X coordinates of screen space.
 * @param mouse_y Y coordinates of screen space.
 */
QVector3D MainView::toNormalizedDeviceCoordinates(int mouse_x, int mouse_y) {
    QVector3D ray_nds;

    // Scale the range of x and y [-1:1] and reverse the direction of y.
    float x = (2.0f * mouse_x) / width() - 1.0f;
    float y = 1.0f - (2.0f * mouse_y) / height();
    float z = 0.0f;
    ray_nds = QVector3D(x, y, z);
    return ray_nds;
}

QVector3D MainView::extractCameraPos() {
    // Get the 3 basis vector planes at the camera origin and transform them into model space.
    //
    // NOTE: Planes have to be transformed by the inverse transpose of a matrix
    //       Nice reference here: http://www.opengl.org/discussion_boards/showthread.php/159564-Clever-way-to-transform-plane-by-matrix
    //
    //       So for a transform to model space we need to do:
    //            inverse(transpose(inverse(MV)))
    //       This equals : transpose(MV) - see Lemma 5 in http://mathrefresher.blogspot.com.au/2007/06/transpose-of-matrix.html
    //
    // As each plane is simply (1,0,0,0), (0,1,0,0), (0,0,1,0) we can pull the data directly from the transpose matrix.
    //
    QMatrix4x4 modelViewT = settings.modelViewMatrix.transposed();

    // Get plane normals
    QVector3D n1(modelViewT.column(0));
    QVector3D n2(modelViewT.column(1));
    QVector3D n3(modelViewT.column(2));

    // Get plane distances
    float d1(modelViewT.column(0).w());
    float d2(modelViewT.column(1).w());
    float d3(modelViewT.column(2).w());

    // Get the intersection of these 3 planes
    // http://paulbourke.net/geometry/3planes/
    QVector3D n2n3 = QVector3D::crossProduct(n2, n3);
    QVector3D n3n1 = QVector3D::crossProduct(n3, n1);
    QVector3D n1n2 = QVector3D::crossProduct(n1, n2);

    QVector3D top = (n2n3 * d1) + (n3n1 * d2) + (n1n2 * d3);
    float denom = QVector3D::dotProduct(n1, n2n3);

    return top / -denom;
}

void MainView::findClosestHalfEdge(const QVector3D& p, const float maxDist){

    int heIndex= -1;
    float currentDist, minDist = 8;
    float maxCP=100000;
    float finalOriginDistance=10;
    Mesh& currentMesh = settings.meshes[settings.subDivValue];
    QVector<HalfEdge>& heList = currentMesh.getHalfEdges();
    QVector3D cameraPos = extractCameraPos();

    for (int i=0;i<heList.size();i++) {
        QVector3D v1 = heList[i].origin->coords;
        QVector3D v2 = heList[i].next->origin->coords;

        // Find distance from the ray
        float dr1 = v1.distanceToLine(cameraPos,p);
        float dr2 = v2.distanceToLine(cameraPos,p);

        float sumdr = dr1 + dr2;

        float sumcp = v1.distanceToPoint(cameraPos) + v2.distanceToPoint(cameraPos);

        if (sumdr < minDist && sumcp <= maxCP){
            minDist = sumdr;
            maxCP = sumcp;
            heIndex = i;
        }
    }
    qDebug() << "Closest HE points" << heList[heIndex].origin->coords << "and" << heList[heIndex].next->origin->coords;
    QVector<unsigned int> vertexCoords ;
    vertexCoords.append(heList[heIndex].origin->index);
    vertexCoords.append(heList[heIndex].next->origin->index);
    settings.edgeSlected = vertexCoords;
    settings.selectedHE = &heList[heIndex];
    qDebug()<< heList[heIndex].sharpness;
}

