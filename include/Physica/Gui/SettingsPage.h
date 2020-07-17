#ifndef PHYSICA_SETTINGSPAGE_H
#define PHYSICA_SETTINGSPAGE_H

#include <QtWidgets/QVBoxLayout>
#include <QtWidgets/QLabel>
#include <QtWidgets/QComboBox>
#include "QWidget"

namespace Physica::Gui {
    class SettingsPage : public QWidget {
    public:
        SettingsPage(QWidget* parent);
    private:
        QVBoxLayout* default_layout;
        QLabel* themeLabel;
        QComboBox* themeBox;
    };
}

#endif
