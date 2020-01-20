#include "../Header/CubePlayer.h"
#include <QGuiApplication>
#include <QScreen>
#include <cmath>
#include <iostream>

/*
 * Visualize the solution of cube.
 * Copyright (c) 2019 NewSigma@163.com.All rights reserved.
 */
CubePlayer::CubePlayer(QWidget *parent) : Plotter(parent,"../Resources/Shader/cubePlayerVS.glsl","../Resources/Shader/cubePlayerFS.glsl") {
    index = points = 0;
    arrayBuf = new QOpenGLBuffer();
    colorBuf = new QOpenGLBuffer();
    matrix = new QMatrix4x4();
    press_pos = new QVector2D();
    screen_diagonal_length = 1;
    vertices = new GLfloat[0];
    colors = new GLfloat[0];
}

CubePlayer::~CubePlayer() {
    arrayBuf->destroy();
    colorBuf->destroy();
    delete arrayBuf;
    delete colorBuf;
    delete matrix;
    delete press_pos;
    delete[] vertices;
    delete[] colors;
}

void CubePlayer::doInitialize() {
    initVertexAttrib(3);

    setVertexAttrib(1,0,0,1,1,1);
    setVertexAttrib(1,0.5,0,1,1,1);
    setVertexAttrib(0,1,0,1,1,1);

    arrayBuf->create();
    arrayBuf->bind();
    arrayBuf->allocate(vertices, 3 * points * (int)sizeof(float));
    arrayBuf->release();

    colorBuf->create();
    colorBuf->bind();
    colorBuf->allocate(colors, 3 * points * (int)sizeof(float));
    colorBuf->release();

    posAttr = shader_program->attributeLocation("posAttr");
    colAttr = shader_program->attributeLocation("colAttr");
    shader_program->enableAttributeArray(posAttr);
    shader_program->enableAttributeArray(colAttr);
    matrixUniform = shader_program->uniformLocation("matrix");
}

void CubePlayer::doPaint() {
    shader_program->setUniformValue(matrixUniform, *matrix);

    arrayBuf->bind();
    shader_program->setAttributeBuffer(posAttr, GL_FLOAT, 0, 3, 3 * sizeof(float));

    colorBuf->bind();
    shader_program->setAttributeBuffer(colAttr, GL_FLOAT, 0, 3, 3 * sizeof(float));

    glDrawArrays(GL_TRIANGLE_STRIP, 0, points);

    arrayBuf->release();
    colorBuf->release();
}

void CubePlayer::doResize(int w, int h) {
    screen_diagonal_length = std::sqrt(double(width() * width() + height() * height()));
}

void CubePlayer::mousePressEvent(QMouseEvent *e) {
    press_pos->setX((float)e->x());
    press_pos->setY((float)e->y());
}

void CubePlayer::mouseReleaseEvent(QMouseEvent *e) {
    press_pos->setX((float)e->x() - press_pos->x());
    press_pos->setY((float)e->y() - press_pos->y());
    float length = press_pos->length();
    matrix->rotate(float(length / screen_diagonal_length) * 360, press_pos->y(), press_pos->x(), 0);
}

void CubePlayer::timerEvent(QTimerEvent* e) {
    Q_UNUSED(e)
    update();
}

void CubePlayer::initVertexAttrib(int i) {
    delete[] vertices;
    delete[] colors;
    this->points = i;
    vertices = new GLfloat[3 * points];
    colors = new GLfloat[3 * points];
}

void CubePlayer::setVertexAttrib(float x, float y, float z, float r, float g, float b) {
    if(index > points - 1)
        return;
    vertices[3 * index] = x;
    vertices[3 * index + 1] = y;
    vertices[3 * index + 2] = z;
    colors[3 * index] = r;
    colors[3 * index + 1] = g;
    colors[3 * index + 2] = b;
    index += 1;
}