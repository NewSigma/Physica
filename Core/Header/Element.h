#ifndef _Physica_C_Element_H
#define _Physica_C_Element_H

#include "Geography.h"

class Element
{
public:
	Point3D* p1;
	Point3D* p2;
	Point3D* p3;
	Point3D* p5;

	Element(Point3D* p1, Point3D* p2, Point3D* p3, Point3D* p5);
};

class Element2D
{
public:
	Point2D* p1;
	Point2D* p2;
	Point2D* p3;

	Element2D(Point2D* p1, Point2D* p2, Point2D* p3);
};

Point2D* get_P4_inEle2D(Element2D* element);
Point3D* get_P4_inEle3D(Element* element);
Point3D* get_P6_inEle3D(Element* element);
Point3D* get_P7_inEle3D(Element* element);
Point3D* get_P8_inEle3D(Element* element);

#endif