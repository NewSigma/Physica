#include "../../Header/Element.h"

/*
p1, p2, p3, p4 consist a skew system where p1 is the origin, as shown in figure below (layer 2 is higher than layer 1):
	layer 1		p3-----p4				layer 2		p7-----p8
					  |          |								  |          |
					p1-----p2								p5-----p6
*/
Element::Element(Point3D* p1, Point3D* p2, Point3D* p3, Point3D* p5)
{
	(*this).p1 = p1;
	(*this).p2 = p2;
	(*this).p3 = p3;
	(*this).p5 = p5;
}

/*
p1, p2, p3 consist a skew system where p1 is the origin, as shown in figure below:
p3-----p4
  |          |
p1-----p2
*/
Element2D::Element2D(Point2D* p1, Point2D* p2, Point2D* p3)
{
	(*this).p1 = p1;
	(*this).p2 = p2;
	(*this).p3 = p3;
}

Point2D* get_P4_inEle2D(Element2D* element)
{
	return new Point2D(element->p1->x + element->p2->x + element->p3->x,
		element->p1->y + element->p2->y + element->p3->y);
}

Point3D* get_P4_inEle3D(Element* element)
{
	return new Point3D(element->p1->x + element->p2->x + element->p3->x,
		element->p1->y + element->p2->y + element->p3->y,
		element->p1->z + element->p2->z + element->p3->z);
}

Point3D* get_P6_inEle3D(Element* element)
{
	return new Point3D(element->p1->x + element->p2->x + element->p5->x,
		element->p1->y + element->p2->y + element->p5->y,
		element->p1->z + element->p2->z + element->p5->z);
}

Point3D* get_P7_inEle3D(Element* element)
{
	return new Point3D(element->p1->x + element->p5->x + element->p3->x,
		element->p1->y + element->p5->y + element->p3->y,
		element->p1->z + element->p5->z + element->p3->z);
}

Point3D* get_P8_inEle3D(Element* element)
{
	return new Point3D(element->p1->x + element->p2->x + element->p3->x + element->p5->x,
		element->p1->y + element->p2->y + element->p3->y + element->p5->y,
		element->p1->z + element->p2->z + element->p3->z + element->p5->z);
}