#ifndef OBJECT_H
#define OBJECT_H

#define PRECISION 1e-10
#define PRECISION2 1e-16
#include <vector>
#include <cmath>
#include "vec3.h"
using namespace std;

enum Type {HOLE, ENTITY};
enum BoundaryType {REFLECT, TRANSMIT, PERIODIC};
Vec3 randomEmitSpace();
Vec3 randomEmitSurface(int normalAxis);
Vec3 randomEmitSurface(Vec3 normal);
Vec3 randomPos(Vec3 boundaries[2], int axis, double value);
Vec3 randomPos(Vec3 boundaries[2]);
double randomNumber();


class Object {
public:
	double transparency;
	Type objectType;
	double specularity;
	Vec3 boundaries[2];
	Object(double alpha, Type t, double s): transparency(alpha), objectType(t), specularity(s) {}
	Vec3 reflection(const Vec3& rayDirection, Vec3& point, bool& scatterWithBoundary);
	virtual bool isWithinBoundaries(Vec3 point, double tolerance) = 0;
	virtual double intersectionDistance(const Vec3 &rayOrigin, const Vec3 &rayDirection) = 0;
	virtual Vec3 getNormal(Vec3& point) = 0;
	virtual bool periodicTransmit(Vec3& point) = 0;
	virtual Vec3 randomPos() = 0;
	virtual Vec3 randomPos(int axis, double value) = 0;
};

void getOutmostBoundaries(vector<Object*> objects, Vec3* boundaries);

class Sphere : public Object {
public:
	Vec3 center;
	double radius;
	BoundaryType bcondition;
	Sphere(Vec3 center, double radius, double transparency, Type objectType, double specularity);
	double intersectionDistance(const Vec3 &rayOrigin, const Vec3 &rayDirection);
	Vec3 getNormal(Vec3& point);
	bool periodicTransmit(Vec3& point){
		return false;
	}
	Vec3 randomPos(){
		double r = randomNumber() * radius;
		return randomEmitSpace() * r + center;
	};
	Vec3 randomPos(int axis, double value){
		Vec3 pos;
		double z = value - center[axis];
		// Assume axis = 2
		double r = randomNumber() * sqrt(radius * radius - z * z);
		double theta = randomNumber() * 2 * M_PI;
		pos[axis] = z;
		pos[(axis+1) % 3] = r * cos(theta);
		pos[(axis+2) % 3] = r * sin(theta);
		return pos + center;
	}
	bool isWithinBoundaries(Vec3 point, double tolerance){
		Vec3 d = point - center;
		if (d.length() <= radius + tolerance)
			return true;
		else
			return false;
	}
};

class Cylinder : public Object {
public:
	Vec3 center, bottomCenter, topCenter;
	double radius, height;
	Cylinder(Vec3 center, double radius, double height, double transparency, Type objectType, double specularity);
	double intersectionDistance(const Vec3 &rayOrigin, const Vec3 &rayDirection);
	bool periodicTransmit(Vec3& point);
	Vec3 getNormal(Vec3& point);
	BoundaryType bcondition; // boundary conditions in z direction
	Vec3 randomPos(){
		double z = (randomNumber() - 0.5) * height;
		double r = randomNumber() * radius;
		double theta = randomNumber() * 2 * M_PI;
		return Vec3(r * cos(theta), r * sin(theta), z) + center;
	}
	Vec3 randomPos(int axis, double value){
		if (axis == 2)
		{
			double z = value;
			double r = randomNumber() * radius;
			double theta = randomNumber() * 2 * M_PI;
			return Vec3(r * cos(theta), r * sin(theta), z) + center;
		}
		else
		{
			// only consider axis=0, otherwise swap x and y
			double x = value - center[axis];
			double y = (2 * randomNumber()  - 1) * sqrt(radius * radius - x * x);
			double z = (randomNumber() - 0.5) * height;
			if (axis == 1) swap(x, y);
			return Vec3(x, y, z) + center;
		}
	}
	bool isWithinBoundaries(Vec3 point, double tolerance){
		Vec3 d = point - center;
		double z = d.z;
		d.z = 0;
		if (d.length() > radius + tolerance)
			return false;
		if (fabs(z) > height / 2 + tolerance)
			return false;
		return true;
	}
};

class Cube : public Object{
public:
	double length, width, height;
	Vec3 center;

	Cube(Vec3 center, double length, double width, double height, double transparency, Type objectType, double specularity);
	double intersectionDistance(const Vec3 &rayOrigin, const Vec3 &rayDirection);
	Vec3 getNormal(Vec3& point);
	bool periodicTransmit(Vec3& point);
	BoundaryType bconditions[3]; // boundary conditions in x, y & z directions

	Vec3 randomPos(){
		Vec3 position=center;
		Vec3 size(length, width, height);
		for (int i = 0; i < 3; i++)
			position[i] += (randomNumber() - 0.5) * size[i];
		return position;
	}
	Vec3 randomPos(int axis, double value){
		Vec3 position=center;
		Vec3 size(length, width, height);
		for (int i = 0; i < 3; i++)
			if (axis != i)
				position[i] += (randomNumber() - 0.5) * size[i];
			else
				position[i] = value;
		return position;
	}
	bool isWithinBoundaries(Vec3 point, double tolerance){
		Vec3 d = point - center;
		if (fabs(d[0]) > length / 2 + tolerance)
			return false;
		if (fabs(d[1]) > width / 2 + tolerance)
			return false;
		if (fabs(d[2]) > height / 2 + tolerance)
			return false;
		return true;
	}
};

#endif
