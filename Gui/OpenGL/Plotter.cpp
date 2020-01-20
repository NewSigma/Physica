#include "../Header/Plotter.h"
#include <QGuiApplication>
#include <QScreen>
#include <iostream>

Plotter::Plotter(QWidget* parent, const char* vertexShaderSource, const char* fragmentShaderSource) {
    Q_UNUSED(parent)
    setUpdateBehavior(QOpenGLWidget::NoPartialUpdate);

    this->vertexShaderSource = tr(vertexShaderSource);
    this->fragmentShaderSource = tr(fragmentShaderSource);
    timer = new QBasicTimer();
}

Plotter::~Plotter() {
    delete timer;
    delete shader_program;
}

void Plotter::initializeGL() {
    initializeOpenGLFunctions();
    shader_program = new QOpenGLShaderProgram(this);
    shader_program->create();
    shader_program->addShaderFromSourceFile(QOpenGLShader::Vertex, vertexShaderSource);
    shader_program->addShaderFromSourceFile(QOpenGLShader::Fragment, fragmentShaderSource);
    shader_program->link();

    doInitialize();
}

void Plotter::paintGL() {
    const double retinaScale = devicePixelRatio();

    glClearColor(0.65,0.9,1,1);
    glClear(GL_COLOR_BUFFER_BIT);
    glViewport(0, 0, width() * retinaScale, height() * retinaScale);

    shader_program->bind();

    doPaint();

    shader_program->release();
}

void Plotter::resizeGL(int w, int h) {
    doResize(w, h);
}

void Plotter::reloadTimer(int msec) {
    timer->start(msec, this);
}