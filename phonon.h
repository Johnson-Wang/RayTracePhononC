#ifndef PHONON_H
#define PHONON_H
#include<cmath>
#include<cassert>
#include <vector>
#include "vec3.h"
#include "object.h"
using namespace std;

class EmitorCamera{
public:
	EmitorCamera(): axis(0), value(0), direction(1){}
	EmitorCamera(int a, double v, int d): axis(a), value(v), direction(d){}
	int axis;
	double value;
	int direction;//+1: see along + direction; -1: along - direction; 0: random (no direction is specified)
};

class Phonon {
public:
	Phonon(double mfp) {
		this->mfp = mfp;
	}
	void initPositionAndDirection(Vec3 pos, Vec3 dir);
	void initPositionAndDirection(vector<Object*> obj);
	void initPositionAndDirection(vector<Object*> objects, int axis, double value, int direction);
	void initPositionAndDirection(vector<Object*> objects, EmitorCamera* emitor);
	double scatterDistance(vector<Object*> objects, EmitorCamera* light, EmitorCamera* camera);
	double scatterDistance(vector<Object*> objects);
	double scatterDistanceWithBoundary(vector<Object*> objects, int& objectIndex);
	double scatterDistanceWithPhonon();
	bool cameraCapture(EmitorCamera* light, EmitorCamera* camera);
	bool lightReturn(EmitorCamera* light, EmitorCamera* camera);
	Vec3 pos0, dir0; // Initial position and direction
	Vec3 pos, dir; // Real-time position and direction
	Vec3 distPeriod; // distance vector (periodic boundaries)
	double mfp;
};

void Phonon::initPositionAndDirection(Vec3 pos, Vec3 dir){
	this->pos = pos;
	this->dir = dir.normalize();
	this->pos0 = pos;
	this->dir0 = dir;
	this->distPeriod = 0;
}

void Phonon::initPositionAndDirection(vector<Object*> objects){
	Vec3 boundaries[2];
	getOutmostBoundaries(objects, boundaries);
	Vec3 posInit;
	Vec3 dirInit;
	bool isWithin, isEntity, isAccept;
	do{
		isAccept=true;
		posInit = randomPos(boundaries);
		for (int i = 0; i < objects.size(); i++)
		{
			isWithin = (objects[i]->isWithinBoundaries(posInit, -PRECISION)); // strictly within
			isEntity = (objects[i]->objectType != ENTITY);
			if (isWithin ^ isEntity){
				isAccept = false;
				break;
			}
		}
	}while(!isAccept);

	dirInit = randomEmitSpace();
	initPositionAndDirection(posInit, dirInit);
} // Space

void Phonon::initPositionAndDirection(vector<Object*> objects, int axis, double value, int direction){
	Vec3 boundaries[2];
	getOutmostBoundaries(objects, boundaries);

	Vec3 posInit, posInitDisrupt;
	Vec3 dirInit;
	bool isWithin, isEntity, isAccept;

	int count=0;
	do{
		isAccept=false;
		posInit = randomPos(boundaries, axis, value);
		posInitDisrupt = posInit;
		posInitDisrupt[axis] += PRECISION;
		for (int i = 0; i < objects.size(); i++)
		{
			isWithin = (objects[i]->isWithinBoundaries(posInitDisrupt, 0)); // make sure phonons only emit to entity
			isEntity = (objects[i]->objectType == ENTITY);
			if (isWithin && isEntity)
				isAccept = true; // make sure find at least one entity
			else if (isWithin && !isEntity) // within a hole
			{
				isAccept = false;
				break;
			}
		}
		if (count > 100)
			cout << "Error";
	}while(!isAccept);

	dirInit = randomEmitSurface(axis);
	if (direction != 0)
		if (dirInit[axis] * direction < 0)
			dirInit = dirInit * -1;
	initPositionAndDirection(posInit, dirInit);
} //Surface


void Phonon::initPositionAndDirection(vector<Object*> objects, EmitorCamera* emitor){
	initPositionAndDirection(objects, emitor->axis, emitor->value, emitor->direction);
} //Surface

double Phonon::scatterDistance(vector<Object*> objects){
	return scatterDistance(objects, NULL, NULL);
}

double Phonon::scatterDistance(vector<Object*> objects, EmitorCamera* light, EmitorCamera* camera){
	double distanceBoundary, dist = 0, distancePhonon;
	bool scatterWithBoundary = false, scatterPhonon=false;
	int scatterObjectIndex=-1;
	Vec3 posTmp, dirTmp;
	Object* obj;
	distancePhonon = scatterDistanceWithPhonon();
	do
	{
		distanceBoundary = scatterDistanceWithBoundary(objects, scatterObjectIndex);
		if (dist + distanceBoundary > distancePhonon)
		{
			scatterPhonon = true;
			break;
		}
		dist += distanceBoundary;
		obj = objects[scatterObjectIndex];
		pos += dir * distanceBoundary;
		posTmp = pos;
		if (obj->periodicTransmit(pos)) //passing through a periodic boundary
			distPeriod += (posTmp - pos);
		else
			dir = obj->reflection(dir, pos, scatterWithBoundary);
		if (cameraCapture(light, camera) || lightReturn(light, camera))
			break;
	} while(!scatterWithBoundary);
	if (scatterPhonon)
	{
		pos += dir * (distancePhonon - dist);
		dir = randomEmitSpace();
		dist = distancePhonon;
	}
	return dist;
}

double Phonon:: scatterDistanceWithBoundary(vector<Object*> objects, int& objectIndex){
	double minDistance=INFINITY, tmpDistance;
	for (int i = 0; i < objects.size(); i++)
	{
		tmpDistance = objects[i]->intersectionDistance(pos, dir);
		if (tmpDistance > 0 && tmpDistance < minDistance)
		{
			minDistance = tmpDistance;
			objectIndex = i;
		}
	}
	return minDistance;
}

double Phonon::scatterDistanceWithPhonon(){
	double randNum = randomNumber();
	return -log(randNum) * mfp;
}

bool Phonon::cameraCapture(EmitorCamera* light, EmitorCamera* camera){
	if (camera == NULL || light == NULL)
		return false;
	Vec3 distVec = (pos - pos0) + distPeriod;
	assert(light->axis == camera->axis);
	double distanceThreshold = camera->value - light->value;
	bool isCapture = false;
	if (camera->direction < 0){
		if (distVec[camera->axis] > distanceThreshold - PRECISION)
			isCapture = true;
	}
	else if (camera->direction > 0){
		if (distVec[camera->axis] < distanceThreshold + PRECISION)
			isCapture = true;
	}
	return isCapture;
}

bool Phonon::lightReturn(EmitorCamera* light, EmitorCamera* camera){
	if (camera == NULL || light == NULL)
		return false;
	Vec3 distVec = pos - pos0 + distPeriod;
	bool isCapture = false;
	if (camera->direction > 0)
		if (distVec[light->axis] > -PRECISION)
			isCapture = true;
	if (camera->direction < 0)
		if (distVec[light->axis] < PRECISION)
			isCapture = true;
	return isCapture;
}


#endif
