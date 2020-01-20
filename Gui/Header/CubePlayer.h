#ifndef PHYSICA_CUBEPLAYER_H
#define PHYSICA_CUBEPLAYER_H

#include "Plotter.h"
#include <QMouseEvent>
#include <QtGui/QOpenGLBuffer>

class CubePlayer : public Plotter {
public:
    explicit CubePlayer(QWidget *parent);
    ~CubePlayer() override;
protected:
    void doInitialize() override;
    void doPaint() override;
    void doResize(int w, int h) override;
    void timerEvent(QTimerEvent* e) override;
    void mousePressEvent(QMouseEvent *e) override;
    void mouseReleaseEvent(QMouseEvent *e) override;
    void initVertexAttrib(int points);
    void setVertexAttrib(float x, float y, float z, float r, float g, float b);
private:
    QOpenGLBuffer* arrayBuf;
    QOpenGLBuffer* colorBuf;
    GLuint posAttr{};
    GLuint colAttr{};
    GLuint matrixUniform{};
    QMatrix4x4* matrix;
    QVector2D* press_pos;
    double screen_diagonal_length;

    GLfloat* vertices;
    GLfloat* colors;
    int points;
    int index;
};

#endif
