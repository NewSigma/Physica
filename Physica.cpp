/*
 * Copyright (c) 2019 NewSigma@163.com. All rights reserved.
 */
#include "Physica.h"
#include <QApplication>
#include "QTime"
#include "Const.h"
#include "PhysicaMain.h"
#include "ConstTest.h"
#include "Numerical.h"

const BasicConst* basicConst;
const MathConst* mathConst;

static QtMessageHandler handler;

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

#include <iostream>
#include <Core/Header/CalcBasic.h>
#include "Numerical.h"
void init() {
    handler = qInstallMessageHandler([](QtMsgType type, const QMessageLogContext &context, const QString &msg) {
        QString prefix{};
        prefix.push_back('[');
        prefix.push_back(QTime::currentTime().toString("hh:mm:ss"));
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
        auto name = getNameFromFile(context.file);
        prefix.push_back(name);
        delete[] name;
        prefix.push_back(':');
        prefix.push_back(QString::fromStdString(std::to_string(context.line)));
        prefix.push_back("]: ");

        if(handler)
            handler(type, context, prefix + msg);
    });
    basicConst = new BasicConst();
    mathConst = new MathConst();
}
//e.g turn /home/user/Physica/Physica.cpp into Physica
char* getNameFromFile(const char* file) {
    int length = strlen(file);
    int temp = length;
    while(--temp)
        if(file[temp] == '/')
            break;
    length = length - temp - 5;
    auto result = new char[length + 1];
    result[length] = '\0';
    memcpy(result, file + temp + 1, length);
    return result;
}

void deInit() {
    delete basicConst;
    delete mathConst;
}