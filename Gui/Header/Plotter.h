#ifndef PHYSICA_C_PLOTER_H
#define PHYSICA_C_PLOTER_H

#include <QtGui/QOpenGLFunctions>
#include <QtWidgets/QOpenGLWidget>
#include <QtGui/QOpenGLShaderProgram>
#include <QtCore/QBasicTimer>

class Plotter : public QOpenGLWidget, protected QOpenGLFunctions
{
Q_OBJECT
public:
    explicit Plotter(QWidget* parent, const char* vertexShaderSource, const char* fragmentShaderSource);
    ~Plotter() override;
    void reloadTimer(int msec);
protected:
    void initializeGL() override;
    void paintGL() override;
    void resizeGL(int w, int h) override;
    QOpenGLShaderProgram* shader_program{};
    QString vertexShaderSource;
    QString fragmentShaderSource;
private:
    virtual void doInitialize() = 0;
    virtual void doPaint() = 0;
    virtual void doResize(int w, int h) = 0;
    QBasicTimer* timer;
};

#endif