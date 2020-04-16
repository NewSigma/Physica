/*
 * Copyright (c) 2019 NewSigma@163.com. All rights reserved.
 */
#include <QApplication>
#include "QTime"
#include "Const.h"
#include "PhysicaMain.h"
#include "Util.h"

const BasicConst* basicConst;
const MathConst* mathConst;

static QtMessageHandler handler;

void init();
void deInit();

int main(int argc, char** argv) {
    init();

    QApplication::setApplicationName("Physica");
    QApplication::setApplicationVersion("0.0.1");
    QApplication::setOrganizationName("NewSigma@163.com");

    QCoreApplication::setAttribute(Qt::AA_UseHighDpiPixmaps);
    QApplication::setAttribute(Qt::AA_EnableHighDpiScaling);

    QApplication app(argc, argv);

    new PhysicaMain();

    int exit_code = QApplication::exec();

    deInit();
    return exit_code;
}

void init() {
    handler = qInstallMessageHandler([](QtMsgType type, const QMessageLogContext &context, const QString &msg) {
        QString prefix{};
        prefix.push_back('[');
        prefix.push_back(QTime::currentTime().toString("hh:mm:ss.zzz"));
        prefix.push_back("] [");
        switch(type) {
            case QtDebugMsg:
                prefix.push_back("Debug");
                break;
            case QtWarningMsg:
                prefix.push_back("Warning");
                break;
            case QtCriticalMsg:
                prefix.push_back("Critical");
                break;
            case QtFatalMsg:
                prefix.push_back("Fatal");
                break;
            case QtInfoMsg:
                prefix.push_back("Info");
                break;
        }
        prefix.push_back("] [");
        prefix.push_back(FILENAME(context.file));
        prefix.push_back(':');
        prefix.push_back(QString::fromStdString(std::to_string(context.line)));
        prefix.push_back("]: ");

        if(handler)
            handler(type, context, prefix + msg);
    });
    basicConst = new BasicConst();
    mathConst = new MathConst();
}

void deInit() {
    delete basicConst;
    delete mathConst;
}