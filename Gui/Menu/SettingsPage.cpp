/*
 * Copyright (c) 2019 NewSigma@163.com. All rights reserved.
 */

#include "../Header/SettingsPage.h"

SettingsPage::SettingsPage(QWidget* parent) : QWidget(parent) {
    default_layout = new QVBoxLayout(this);

    themeLabel = new QLabel("Theme", this);
    default_layout->addWidget(themeLabel);

    themeBox = new QComboBox(this);
    default_layout->addWidget(themeBox);
}