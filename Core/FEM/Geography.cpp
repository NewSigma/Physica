#include "../Header/Geography.h"
#include <math.h>

Point2D::Point2D(double x, double y)
{
	(*this).x = x;
	(*this).y = y;
}

Point3D::Point3D(double x, double y, double z)
{
	(*this).x = x;
	(*this).y = y;
	(*this).z = z;
}

Vector2D::Vector2D(double x, double y)
{
	(*this).x = x;
	(*this).y = y;
}

Vector3D::Vector3D(double x, double y, double z)
{
	(*this).x = x;
	(*this).y = y;
	(*this).z = z;
}

Point2D* Point_3Dto2D(Point3D* p)
{
	return new Point2D(p->x, p->y);
}

Point3D* Point_2Dto3D(Point2D* p)
{
	return new Point3D(p->x, p->y, 0);
}

Vector2D* getVectorByPoints(Point2D* start, Point2D* end)
{
	return new Vector2D(end->x - start->x, end->y - start->y);
}

Vector3D* getVectorByPoints3D(Point3D* start, Point3D* end)
{
	return new Vector3D(end->x - start->x, end->y - start->y, end->z - start->z);
}

void toUnitVector(Vector2D* v)
{
	double norm = pow(v->x * v->x + v->y * v->y, 0.5);
	v->x = v->x / norm;
	v->y = v->y / norm;
}

void toUnitVector(Vector3D* v)
{
	double norm = pow(v->x * v->x + v->y * v->y + v->z * v->z, 0.5);
	v->x = v->x / norm;
	v->y = v->y / norm;
	v->z = v->z / norm;
}