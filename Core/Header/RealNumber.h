#ifndef _Physica_C_ReakNumber_H
#define _Physica_C_ReakNumber_H

#include <ostream>
#include "Number.h"
/*
 * The following two classes have the same statistics. The differences between them are their method to
 * handle the statistics. RealNumber will not consider the accuracy because it thinks the statistics is
 * accurate or the accuracy will not change. However, RealNumberA will calculate the accuracy.
 *
 * They can be convented to each other safely.
 */
class RealNumber : private Number {
public:
    //Store effective digits.
    unsigned char* byte;
    //Length of byte.
    int length;
    //Number = (x0 +- a * 10^(1-length)) * 10^power
    int power;
    //True if RealNumber > 0 and false if RealNumber < 0. (RealNumber != 0)
    bool sign;
    //Accuracy
    unsigned char a = 0;

    RealNumber();
    explicit RealNumber(double d, unsigned char acc = 0);
    explicit RealNumber(std::wstring s, unsigned char acc = 0);
    RealNumber(const RealNumber& n);
    RealNumber(unsigned char* byte, int length, int power, bool sign = true, unsigned char acc = 0);
    explicit RealNumber(const RealNumber* n);
    ~RealNumber();

    std::string toString() const;
    friend std::ostream& operator<<(std::ostream& os, const RealNumber& n);
    void operator<<(RealNumber& n);
    explicit operator double() const;
    RealNumber& operator= (const RealNumber& n);
    RealNumber* getAccuracy() const;
    RealNumber* getMaximum() const;
    RealNumber* getMinimum() const;
    bool applyError(const RealNumber* copy_error);
    bool isZero() const { return byte[0] == 0; }
    bool isPositive() const { return !isZero() && sign; }
    bool isNegative() const { return !isZero() && !sign; }
    bool isInteger() const { return length == power + 1; }
};
//////////////////////////////////Operators///////////////////////////////////////
RealNumber* operator+ (const RealNumber& n1, const RealNumber& n2);
RealNumber* operator- (const RealNumber& n1, const RealNumber& n2);
RealNumber* operator* (const RealNumber& n1, const RealNumber& n2);
RealNumber* operator/ (const RealNumber& n1, const RealNumber& n2);
void operator+= (RealNumber& n1, const RealNumber& n2);
void operator-= (RealNumber& n1, const RealNumber& n2);
void operator*= (RealNumber& n1, const RealNumber& n2);
void operator/= (RealNumber& n1, const RealNumber& n2);
bool operator> (const RealNumber& n1, const RealNumber& n2);
bool operator< (const RealNumber& n1, const RealNumber& n2);
bool operator>= (const RealNumber& n1, const RealNumber& n2);
bool operator<= (const RealNumber& n1, const RealNumber& n2);
bool operator== (const RealNumber& n1, const RealNumber& n2);
bool operator!= (const RealNumber& n1, const RealNumber& n2);
RealNumber* operator^ (const RealNumber& n1, const RealNumber& n2);
void operator^= (RealNumber& n1, const RealNumber& n2);
RealNumber* operator- (const RealNumber& n);
////////////////////////////////Helper functions/////////////////////////////////////
RealNumber* randomRealNumber();
RealNumber* randomRealNumber(RealNumber* lowerBound, RealNumber* upperBound);
//////////////////////////////Process functions////////////////////////////////////////
RealNumber* add (const RealNumber& n1, const RealNumber& n2);
RealNumber* subtract (const RealNumber& n1, const RealNumber& n2);
RealNumber* multiply (const RealNumber& n1, const RealNumber& n2);
RealNumber* divide (const RealNumber& n1, const RealNumber& n2);
bool cutLength(RealNumber* n);
void cutZero(RealNumber* n);

#endif