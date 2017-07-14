#include <iostream>
#include "object.h"
#include <vector>
#include <cmath>
#include <cstdlib>
#include<ctime>
#include "matrix.h"
using namespace std;

void swap(double& i, double& j){
	double tmp = i;
	i = j;
	j = tmp;
}

double randomNumber(){

	return (double) rand() / RAND_MAX;
}

Vec3 randomEmitSpace(){
	double theta = randomNumber();
	double v = randomNumber();
	double cosphi = 1 - 2 * v;
	double sinphi = sqrt(1 - cosphi * cosphi);
	theta *= 2 * M_PI;
	Vec3 dir;
	dir[0] = cosphi;
	dir[1] = cos(theta) * sinphi;
	dir[2] = sin(theta) * sinphi;
	return dir.normalize();
}

Vec3 randomPos(Vec3 boundaries[2], int axis, double value){
	Vec3 size = boundaries[1] - boundaries[0];
	Vec3 position = boundaries[0];
	for (int i = 0; i < 3; i++)
		if (axis != i)
			position[i] += randomNumber() * size[i];
		else
			position[i] = value;
	return position;
}

Vec3 randomPos(Vec3 boundaries[2]){
	Vec3 size = boundaries[1] - boundaries[0];
	Vec3 position = boundaries[0];
	for (int i = 0; i < 3; i++)
		position[i] += randomNumber() * size[i];
	return position;
}

Vec3 randomEmitSurface(int normalAxis){ // normalAxis: 0, 1, 2
	double theta = randomNumber();
	double v = randomNumber();
	double sinphi = sqrt(v);
	double cosphi = sqrt(1-v);
	theta *= 2 * M_PI;
	Vec3 dir;
	dir[normalAxis] = cosphi;
	dir[(normalAxis+1)%3] = cos(theta) * sinphi;
	dir[(normalAxis+2)%3] = sin(theta) * sinphi;
	return dir.normalize();
}

void getOutmostBoundaries(vector<Object*> objects, Vec3* boundaries)
{
	boundaries[0] = INFINITY; boundaries[1] = -INFINITY;
	for (int i= 0; i < objects.size(); i++)
	{
		for (int j = 0; j < 3; j++){
			if (boundaries[0][j] > objects[i]->boundaries[0][j])
				boundaries[0][j] = objects[i]->boundaries[0][j];
			if (boundaries[1][j] < objects[i]->boundaries[1][j])
				boundaries[1][j] = objects[i]->boundaries[1][j];
		}//find the outmost boundaries
	}
}

Vec3 randomEmitSurface(Vec3 normal){
	int i;
	for (i = 0; i < 3; i++){
		if (fabs(fabs(normal[i]) - 1) < PRECISION){
			return randomEmitSurface(i);
		}
	} // For the case the normal vector is parallel to an axis
	// The following has reference:
	//https://math.stackexchange.com/questions/180418/calculate-rotation-matrix-to-align-vector-a-to-vector-b-in-3d
	int argmin = 0;
	double normalMin = INFINITY;
	for (i = 0; i < 3; i++){
		if (fabs(normal[i]) < normalMin){
			normalMin = fabs(normal[i]);
			argmin = i;
		}
	}
	// avoid collinear situation between normal and u;
	Vec3 u(0, 0, 0);
	u[argmin] = 1;
	Vec3 emitz = randomEmitSurface(argmin);
	Matrix U = rotationFromTwoVectors(u, normal);
	return U.dot(emitz);
}
// Object
/*********************************************************************************************************************/

Vec3 Object::reflection(const Vec3& rayDirection, Vec3& point, bool& isScatterWithBoundary){
	Vec3 normal = getNormal(point);
	Vec3 newDir;
	double randomNum;
	bool isSpecular;
	if (specularity == 0 )
		isSpecular = false;
	else if (specularity == 1)
		isSpecular = true;
	else
	{
		randomNum = (double) (rand() / RAND_MAX);
		if (randomNum > specularity)
			isSpecular = true;
		else
			isSpecular = false;
	}
	isScatterWithBoundary = (!isSpecular);
	if (isSpecular)
	{
		newDir = rayDirection.project(normal) * (-2) + rayDirection;
		return newDir.normalize();
	}
	else
	{
		newDir = randomEmitSurface(normal);
		if (newDir.dot(normal) < 0)
			newDir = newDir * -1;
		return newDir;
	}

}

/*********************************************************************************************************************/
// Sphere Object
/*********************************************************************************************************************/
Sphere::Sphere(Vec3 center, double radius, double transparency, Type objectType, double specularity)
	:Object(transparency, objectType, specularity) {
	this->center = center;
	this->radius = radius;
	this->boundaries[0] = Vec3(center-radius);
	this->boundaries[1] = Vec3(center+radius);
}

/*
	Sphere intersection
		- solving equation for a point on ray which is at distance r from center
*/
double Sphere::intersectionDistance(const Vec3 &ro, const Vec3 &rd) {
	/*Transformation done to ray origin and ray direction for world  -> object coordinates*/
	Vec3 rayOrigin = ro;
	Vec3 rayDirection = rd;
	double Ird = rayDirection.length();
	rayDirection.normalize();

	vector<double> points;
	double a, b, c;
	a = rayDirection.length2();
	b = 2 * rayDirection.dot(rayOrigin - center);
	c = (rayOrigin - center).length2() - radius*radius;
	double discriminant = b*b - 4*a*c;
	if (discriminant < 0) return -1;
	else {
		discriminant = sqrt(discriminant);
		double t1 = ((-1 * b) + discriminant) / (2 * a);
		double t2 = ((-1 * b) - discriminant) / (2 * a);
		points.push_back(t1);
		points.push_back(t2);
	}
	double minT = INFINITY;
	bool flag = false;
	for(int i=0;i<points.size();i++){
		if(points[i]>=PRECISION2 && minT > points[i]) {
			minT = points[i];
			flag = true;
		}
	}

	if (flag)
		return minT/Ird;
	else
		return -1;
}
Vec3 Sphere::getNormal(Vec3& point){
	/*for normal first convert p to object space than find normal and then normal converted to world space and normalized*/
	return (center - point).normalize();
}

/*********************************************************************************************************************/

// Cylinder Object
/*********************************************************************************************************************/
Cylinder::Cylinder(Vec3 center, double radius, double height, double transparency, Type objectType, double specularity)
	:Object(transparency, objectType, specularity) {
	this->center = center;
	this->radius = radius;
	this->height = height;
	this->bottomCenter = center; this->bottomCenter[2] -= height / 2;
	this->topCenter = center; this->topCenter[2] += height / 2;
	this->boundaries[0][0] = center[0] - radius;
	this->boundaries[0][1] = center[1] - radius;
	this->boundaries[0][2] = center[2] - height / 2;
	this->boundaries[1][0] = center[0] + radius;
	this->boundaries[1][1] = center[1] + radius;
	this->boundaries[1][2] = center[2] + height / 2;
}

// source : http://mrl.nyu.edu/~dzorin/rend05/lecture2.pdf
/*
	Cylinder Intersection 
		- First solving a point on ray at distance r from a infinite line which is axis of cylinder
		- then cut cylinder through two planes and check if the intersection is still there
		- if not then check if intersection on plane is within range of cylinder
*/
double Cylinder::intersectionDistance(const Vec3 &ro, const Vec3 &rd) {
	Vec3 rayOrigin = ro;
	Vec3 rayDirection = rd;
	Vec3 center1 = center, center2 = center;
	center1[2] -= height / 2; center2[2] += height / 2;
	double a = rayDirection[0] * rayDirection[0] + rayDirection[1] * rayDirection[1];
	double b = 2 * (rayOrigin[0] - center[0]) * rayDirection[0] + 2 * (rayOrigin[1] - center[1]) * rayDirection[1];
	double c = (rayOrigin[0] - center[0]) * (rayOrigin[0] - center[0]) + (rayOrigin[1] - center[1]) * (rayOrigin[1] - center[1]) - radius*radius;
	double discriminant = b*b - 4*a*c;
	vector<double> points;
	if (discriminant < 0)
		return -1;

	discriminant = sqrt(discriminant);
	double t1 = ((-1 * b) + discriminant) / (2 * a);
	double t2 = ((-1 * b) - discriminant) / (2 * a);
	double intersectz;
	if(t1 >= PRECISION2){
		intersectz = rayOrigin[2] + rayDirection[2] * t1;
		if( intersectz>(center1[2]) && intersectz<(center2[2]))
			points.push_back(t1);
	}
	if(t2 >= PRECISION2)
	{
		intersectz = rayOrigin[2] + rayDirection[2] * t2;
		if( intersectz>(center1[2]) && intersectz<(center2[2]))
			points.push_back(t2);
	}

	double denom = rayDirection[2];
	if (denom > PRECISION || denom < -PRECISION) { // not parallel with bottom plane
		Vec3 co = center1 - rayOrigin;
		double t3 = co[2] / denom; // distance to the bottom plane
		if(t3 > PRECISION2 && (rayDirection * t3 - co).length2() <= radius*radius)
			points.push_back(t3);
		Vec3 co2 = center2 - rayOrigin;
		double t4 = co2[2] / denom;
		if(t4 > PRECISION2 && (rayDirection * t4 - co2).length2() <= radius*radius)
			points.push_back(t4);
	}
	double minT = INFINITY;
	bool flag = false;
	for(int i=0;i<points.size(); i++){
		if(minT > points[i] && points[i]>=0) {
			minT = points[i];
			flag = true;
		}
	}
	if(flag)
		return minT;
	else
		return -1;
}
Vec3 Cylinder::getNormal(Vec3& p){
	if(fabs(p[2] - bottomCenter[2]) < 1e-6) // bottom
	{
		p[2] = bottomCenter[2]; // Avoid accumulation error
		if (objectType == ENTITY)
			return Vec3(0, 0, 1);
		else
			return Vec3(0, 0, -1);
	}
	else if(fabs(p[2] - topCenter[2]) < 1e-6) // top
	{
		p[2] = topCenter[2];
		if (objectType == ENTITY)
			return Vec3(0, 0, -1);
		else
			return Vec3(0, 0, 1);
	}

	Vec3 normal(p[0], p[1], 0);
	double correction = radius / normal.length();
	p[0] *= correction;// Avoid accumulation error
	p[1] *= correction; // Avoid accumulation error
	return normal.normalize() * -1;
}

bool Cylinder::periodicTransmit(Vec3& point)
{
	if (bcondition != PERIODIC)
		return false;
	if(fabs(point[2] - bottomCenter[2]) < PRECISION) // bottom
	{
		point[2] = topCenter[2];
		return true;
	}
	else if (fabs(point[2] - topCenter[2]) < PRECISION)
	{
		point[2] = bottomCenter[2];
		return true;
	}
	return false;
}
/*********************************************************************************************************************/

/*********************************************************************************************************************/
Cube::Cube(Vec3 center, double length, double width, double height, double transparency, Type objectType, double specularity)
	:Object(transparency, objectType, specularity){
	this->center = center;
	this->length = length;
	this->width = width;
	this->height = height;
	Vec3 boundaryMin(center[0]-length/2, center[1]-width/2, center[2]-height/2);
	Vec3 boundaryMax(center[0]+length/2, center[1]+width/2, center[2]+height/2);
	this->boundaries[0] = boundaryMin;
	this->boundaries[1] = boundaryMax;
}

double Cube::intersectionDistance(const Vec3& ro, const Vec3 & rd)
{
	Vec3 tmin, tmax;
	tmin = (boundaries[0] - ro) / rd;
	tmax = (boundaries[1] - ro) / rd;
	for (int i = 0; i < 3; i++)
		if (tmin[i] > tmax[i])
			swap(tmin[i], tmax[i]);
	double minDistance = tmin.max();
	if (minDistance > tmax.min())
		return -1;
	if (minDistance < PRECISION)
		return tmax.min();
	else
		return minDistance;
}

Vec3 Cube::getNormal(Vec3& p){
	Vec3 normal;
	int i;
	for (i = 0; i < 3; i++){
		if (fabs(p[i] - boundaries[0][i]) < PRECISION){
			normal[i] = 1;
			p[i] = boundaries[0][i]; // avoid accumulation error
			break;
			}
		else if (fabs(p[i] - boundaries[1][i]) < PRECISION){
			normal[i] = -1;
			p[i] = boundaries[1][i];
			break;
		}
	}
	if (objectType == HOLE)
		normal[i] = -normal[i];
	return normal;
}

bool Cube::periodicTransmit(Vec3& p){
	int i;
	for (i = 0; i < 3; i++){
		if (bconditions[i] != PERIODIC)
			continue;
		if (fabs(p[i] - boundaries[0][i]) < PRECISION){
			p[i] = boundaries[1][i];
			return true;
		}
		else if (fabs(p[i] - boundaries[1][i]) < PRECISION){
			p[i] = boundaries[0][i];
			return true;
		}
	}
	return false;
}
/*********************************************************************************************************************/


