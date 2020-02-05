#ifndef _Physica_C_ReakNumber_H
#define _Physica_C_ReakNumber_H

#include <ostream>
/*
 * The following two classes have the same statistics. The differences between them are their method to
 * handle the statistics. RealNumber will not consider the accuracy because it thinks the statistics is
 * accurate or the accuracy will not change. However, RealNumberA will calculate the accuracy.
 *
 * They can be convented to each other safely.
 */
class RealNumber {
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
    explicit RealNumber(double d);
    RealNumber(const RealNumber& n);
    RealNumber(unsigned char* byte, int length, int power, bool sign = true);
    explicit RealNumber(const RealNumber* n);
    ~RealNumber();
    bool isZero() const { return byte[0] == 0; }
    bool isPositive() const { return !isZero() && sign; }
    bool isNegative() const { return !isZero() && !sign; }

    friend std::ostream& operator<<(std::ostream& os, const RealNumber& n);
    explicit operator double() const;
    RealNumber& operator= (const RealNumber& n);
};
/*
 * Real number with accuracy.
 */
class RealNumberA : public RealNumber {
public:
    RealNumberA();
    explicit RealNumberA(double d, unsigned char acc = 0);
    RealNumberA(const RealNumberA& n);
    RealNumberA(unsigned char* byte, int length, int power, bool sign = true, unsigned char acc = 0);
    explicit RealNumberA(const RealNumber* n, unsigned char a = 0);
    friend std::ostream& operator<<(std::ostream& os, const RealNumberA& n);
    RealNumberA& operator= (const RealNumberA& n);
    //Helper functions
    RealNumber* getAccuracy() const;
    RealNumber* getMaximum() const;
    RealNumber* getMinimum() const;
    bool applyError(const RealNumber* copy_error);
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

RealNumberA* operator+ (const RealNumberA& n1, const RealNumberA& n2);
RealNumberA* operator- (const RealNumberA& n1, const RealNumberA& n2);
RealNumberA* operator* (const RealNumberA& n1, const RealNumberA& n2);
RealNumberA* operator/ (const RealNumberA& n1, const RealNumberA& n2);
void operator+= (RealNumberA& n1, const RealNumberA& n2);
void operator-= (RealNumberA& n1, const RealNumberA& n2);
void operator*= (RealNumberA& n1, const RealNumberA& n2);
void operator/= (RealNumberA& n1, const RealNumberA& n2);
////////////////////////////////Helper functions/////////////////////////////////////
RealNumber* randomRealNumber();
RealNumber* randomRealNumber(RealNumber* lowerBound, RealNumber* upperBound);
inline bool isInteger(const RealNumber* n) { return n->length == n->power + 1; };
//////////////////////////////Process functions////////////////////////////////////////
RealNumber* add (const RealNumber* n1, const RealNumber* n2);
RealNumber* subtract (const RealNumber* n1, const RealNumber* n2);
RealNumber* multiply (const RealNumber* n1, const RealNumber* n2);
RealNumber* divide (const RealNumber* n1, const RealNumber* n2);
bool cutArray(RealNumber* n);

#endif