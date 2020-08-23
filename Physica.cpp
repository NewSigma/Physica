/*
 * Copyright 2019 WeiBo He.
 *
 * This file is part of Physica.

 * Physica is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * Physica is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with Physica.  If not, see <https://www.gnu.org/licenses/>.
 */
#include <QApplication>
#include <QTime>
#include "Physica/Physica.h"
#include "Physica/Core/MultiPrecition/Scalar.h"
#include "Physica/Gui/PhysicaMain.h"
#include "Physica/Core/Utils/DebugUtil.h"

using namespace Physica::Core;
static QtMessageHandler handler;

int main(int argc, char** argv) {
    initPhysica();

    QApplication::setApplicationName("Physica");
    QApplication::setApplicationVersion("0.0.1");
    QApplication::setOrganizationName("NewSigma@163.com");

    QCoreApplication::setAttribute(Qt::AA_UseHighDpiPixmaps);
    QApplication::setAttribute(Qt::AA_EnableHighDpiScaling);

    QApplication app(argc, argv);

    new Physica::Gui::PhysicaMain();

    int exit_code = QApplication::exec();

    deInitPhysica();
    return exit_code;
}

void initPhysica() {
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
    BasicConst::init();
    MathConst::init();
}

void deInitPhysica() {
    BasicConst::deInit();
    MathConst::deInit();
}