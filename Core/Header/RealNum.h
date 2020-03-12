#ifndef _Physica_C_ReakNumber_H
#define _Physica_C_ReakNumber_H

#include "ComplexNum.h"
/*
 * The following two classes have the same statistics. The differences between them are their method to
 * handle the statistics. RealNum will not consider the accuracy because it thinks the statistics is
 * accurate or the accuracy will not change. However, RealNumberA will calculate the accuracy.
 *
 * They can be convented to each other safely.
 */
class RealNum : public ComplexNum {
public:
    RealNum(Numerical* n);
    RealNum(const RealNum& n);
    RealNum(const RealNum* n);
    virtual ~RealNum();

    virtual AbstractNum* operator+ (const AbstractNum& n) const;
    virtual AbstractNum* operator- (const AbstractNum& n) const;
    virtual AbstractNum* operator* (const AbstractNum& n) const;
    virtual AbstractNum* operator/ (const AbstractNum& n) const;
    virtual AbstractNum* operator^ (const AbstractNum& n) const;
    virtual AbstractNum* operator- () const;
    virtual bool operator== (const AbstractNum& n) const;
    bool operator> (const RealNum& n) const;
    bool operator< (const RealNum& n) const;
    bool operator>= (const RealNum& n) const;
    bool operator<= (const RealNum& n) const;
};

#endif