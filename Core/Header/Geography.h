#ifndef _Physica_C_Geography_H
#define _Physica_C_Geography_H

class Point2D
{
public:
	double x;
	double y;

	Point2D(double x, double y);
};

class Point3D
{
public:
	double x;
	double y;
	double z;

	Point3D(double x, double y, double z);
};

class Vector2D
{
public:
	double x;
	double y;

	Vector2D(double x, double y);
};

class Vector3D
{
public:
	double x;
	double y;
	double z;

	Vector3D(double x, double y, double z);
};

Point2D* Point_3Dto2D(Point3D* p);
Point3D* Point_2Dto3D(Point2D* p);
Vector2D* getVectorByPoints(Point2D* start, Point2D* end);
Vector3D* getVectorByPoints3D(Point3D* start, Point3D* end);
void toUnitVector(Vector2D* v);

#endif