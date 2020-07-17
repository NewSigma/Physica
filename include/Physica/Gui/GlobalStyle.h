/*
 * Copyright (c) 2019 NewSigma@163.com. All rights reserved.
 */
#ifndef PHYSICA_GLOBALSTYLE_H
#define PHYSICA_GLOBALSTYLE_H

#include <QtWidgets/QProxyStyle>

namespace Physica::Gui {
    class GlobalStyle : public QProxyStyle
    {
    public:
        GlobalStyle();
        int styleHint(StyleHint hint, const QStyleOption* option = nullptr, const QWidget* widget = nullptr, QStyleHintReturn* returnData = nullptr) const override;
    };
}

#endif
