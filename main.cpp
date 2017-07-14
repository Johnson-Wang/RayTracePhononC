#include <iostream>
#include <cmath>
#include <vector>
#include <fstream>
#include <string>
#include <sstream>
#include "vec3.h"
#include "object.h"
#include "matrix.h"
#include "phonon.h"
#include <algorithm>
#include <cstdlib>
#include "phonon.h"
#include <omp.h>
using namespace std;

//#define DEBUG
#ifndef DEBUG

bool mfpCompare(vector<double>a, vector<double> b){
	return (a[0] < b[0]);
}

bool read_line(ifstream& ifs, stringstream& ss){
	string line;
	bool success = getline(ifs, line);
//	line.erase(find(line.begin(), line.end(), '#'), line.end()); // erase anything after #
	ss.str(line);
	return success;
}

vector< vector<double> > generateRandomMFPs(vector< vector<double> > cumkappa, int nseg){
	double totalK, kappa;
	int i, j;
	vector < vector<double> > mfpk;
	mfpk.assign(nseg, vector<double>(2,0));
	totalK = cumkappa.back()[1];
	for (i = 0; i < nseg; i++){
		kappa = randomNumber() * totalK;
		mfpk[i][1] = kappa;
		int argmin=0, left, right;
		double min=INFINITY;
		for (j = 0; j < cumkappa.size(); j++){
			if (fabs(cumkappa[j][1] - kappa) < min)
			{
				min = fabs(cumkappa[j][1] - kappa);
				argmin = j;
			}
		}
		if (cumkappa[argmin][1] > kappa) {
			while(cumkappa[argmin-1][1] > kappa){
				argmin--;
			}
			right = argmin;
			left = argmin - 1;
		}
		else{
			while(cumkappa[argmin+1][1] < kappa){
				argmin++;
			}
			left = argmin;
			right = argmin+1;
		}
		// find the closest two nodes
		mfpk[i][0] = pow(10,
				(kappa - cumkappa[left][1]) / (cumkappa[right][1] - cumkappa[left][1])
				*(log10(cumkappa[right][0]) - log10(cumkappa[left][0]))  + log10(cumkappa[left][0]));
//		if (isinf(mfpk[i][0]))
//			cout << "Something Wrong";
	}
	sort(mfpk.begin(), mfpk.end(), mfpCompare);
	return mfpk;
}

int main(int argc, char** argv){
	if(argc!=2){
		cout<<"Usage: "<<argv[0]<<"+ Description file!"<<endl;
		return -1;
	}
	/*Opening and parsing of scene description file*/
	ifstream ifs;
	ifs.open(argv[1], ios::in);
	if(!ifs.is_open()) {
		cout<<"Error opening file"<<endl;
		return -1;
	}
	srand((unsigned)time(0));
	vector<Object*> objects;
	EmitorCamera light;
	EmitorCamera camera;
	int i, j;
	bool isLight=false, isCamera=false, isMFPFile=false, isMFP=false;
	bool isObject=false, isNPhonon=false, isNSeg=false;
	double tolerance = 1e-8;
	float minCameraDistance = 1;
	int mfpseg;
	int nPhonon;
	string mfpfile, line;
	double kfactor, kappa, mfp, totalK, cameraMax=INFINITY;
	Phonon* phon;

	while(ifs.good()){
		string type;
		stringstream ss;
		if (!read_line(ifs, ss)) break;
		ss>>type;
		if (type[0] == '#');
		else if(type=="LIGHT"){
			isLight = true;
			ss>>light.axis>>light.value>>light.direction;
		}
		else if(type=="CAMERA"){
			isCamera = true;
			ss>>camera.axis>>camera.value>>camera.direction;
		}

		else if(type == "OBJECTS") {
			isObject = true;
			while(ifs.good()){
				stringstream ss;
				read_line(ifs, ss);
				ss>>type;
				if (type == "BEGIN");
				else if (type == "END") break;
				else if(type == "SPHERE"){
					Vec3 center;
					double radius;
					string sur;
					Type surface;
					string specularity;
					double sp;
					read_line(ifs, ss);
					ss>>center.x>>center.y>>center.z;
					read_line(ifs, ss);
					ss>>radius;
					read_line(ifs, ss);
					ss >> sur;
					if (sur[0] == 'E')
						surface = ENTITY;
					else if (sur[0] == 'H')
						surface = HOLE;
					else{
						cerr<<"Surface type set incorrectly!";
						exit(1);
					}
					read_line(ifs, ss);
					ss>>specularity;
					try{
						sp = atof(specularity.c_str());
					}
					catch(...){
						if (specularity == "DIFFUSE")
							sp = 0;
						else if (specularity == "SPECULAR")
							sp = 1;
					}
					Sphere* sphere = new Sphere(center, radius, 0, surface, sp);
					string bcondition;
					read_line(ifs, ss);
					ss>> bcondition;
					if (bcondition[0] == 'P')
						sphere->bcondition = PERIODIC;
					else if (bcondition[0] == 'R')
						sphere->bcondition = REFLECT;
					else if (bcondition[0] == 'T')
						sphere->bcondition = TRANSMIT;
					else{
						cerr<<"Boundary type set incorrectly!"<<endl;
						exit(1);
					}
					objects.push_back(sphere);
				}
				else if(type == "CYLINDER"){
					Vec3 center;
					double radius;
					double height;
					Type surface;
					string specularity;
					double sp;
					read_line(ifs, ss);
					ss>>center.x>>center.y>>center.z;
					read_line(ifs, ss);
					ss>>radius>>height;
					string sur;
					read_line(ifs, ss);
					ss>>sur;
					if (sur[0] == 'E')
						surface = ENTITY;
					else if (sur[0] == 'H')
						surface = HOLE;
					else{
						cerr<<"Surface type set incorrectly!";
						exit(1);
					}
					read_line(ifs, ss);
					ss>>specularity;
					try{
						sp = atof(specularity.c_str());
					}
					catch(...){
						if (specularity == "DIFFUSE")
							sp = 0;
						else if (specularity == "SPECULAR")
							sp = 1;
					}
					Cylinder* cylinder = new Cylinder(center, radius, height, 0, surface, sp);
					string bcondition;
					read_line(ifs, ss);
					ss >> bcondition;
					if (bcondition[0] == 'P')
						cylinder->bcondition = PERIODIC;
					else if (bcondition[0] == 'R')
						cylinder->bcondition = REFLECT;
					else if (bcondition[0] == 'T')
						cylinder->bcondition = TRANSMIT;
					else{
						cerr<<"Boundary type set incorrectly!"<<endl;
						exit(1);
					}
					objects.push_back(cylinder);
				}
				else if(type == "CUBE"){
					Vec3 center;
					double length, width, height;
					Type surface;
					string specularity;
					double sp;
					read_line(ifs, ss);
					ss>>center.x>>center.y>>center.z;
					read_line(ifs, ss);
					ss>>length>>width>>height;
					string sur;
					read_line(ifs, ss);
					ss>>sur;
					if (sur[0] == 'E')
						surface = ENTITY;
					else if (sur[0] == 'H')
						surface = HOLE;
					else{
						cerr<<"Surface type set incorrectly!";
						exit(1);
					}
					read_line(ifs, ss);
					ss>>specularity;
					try{
						sp = atof(specularity.c_str());
					}
					catch(...){
						if (specularity == "IFFUSE")
							sp = 0;
						else if (specularity == "SPECULAR")
							sp = 1;
						else
						{
							cout << "Specularity set incorrectly!"<<endl;
							exit(1);
						}
					}
					Cube* cube = new Cube(center, length, width, height, 0, surface, sp);
					vector<string> bconditions;
					bconditions.assign(3, "");
					read_line(ifs, ss);
					ss>> bconditions[0] >> bconditions[1] >> bconditions[2];
					for (i = 0; i < 3; i++){
						if (bconditions[i][0] == 'P')
							cube->bconditions[i] = PERIODIC;
						else if (bconditions[i][0] == 'R')
							cube->bconditions[i] = REFLECT;
						else if (bconditions[i][0] == 'T')
							cube->bconditions[i] = TRANSMIT;
						else{
							cerr<<"Boundary type set incorrectly!"<<endl;
							exit(1);
						}
					}
					objects.push_back(cube);
				}
				else {
					cout << "Wrong objects settings;"<<endl;
					break;
				}
			}
		}
		else if (type == "MFPFILE")
		{
			isMFPFile = true;
			ss>>mfpfile;
		}
		else if (type == "MFP"){
			isMFP = true;
			ss >> mfp;
		}
		else if (type == "MFPSEG")
		{
			isNSeg = true;
			ss>>mfpseg;
		}
		else if (type == "NPHONON")
		{
			isNPhonon = true;

			ss>>nPhonon;
		}
		else if (type == "TOLERANCE")
		{
			ss>>tolerance;
		}

		else if (type == "CAMERARATIO") //MINIMUM CAMERA DISTANCE (in regards to mfp)
		{
			ss>>minCameraDistance;
		}
		else if (type == "CAMERAMAX")
		{
			ss>>cameraMax;
		}

		else if(type == "CUT"){
			break;
		}
		else {
			cerr << "type " << type << " not recognizable!"<<endl;
			return 1;
		}
	}
	ifs.close();
	if (!(isLight && isCamera&&isObject&&isNPhonon))
	{
		cout<<"Setup file error!"<<endl;
		exit(1);
	}
	vector< vector<double> > mfpks;
	if (isMFPFile){ // mfp sampling mode
		ifs.open(mfpfile.c_str(), ios::in);
		if(!ifs.is_open()) {
			cerr<<"Error opening file"<<endl;
			return -1;
		}
		vector< vector<double> > cumkappa;
		while(ifs.good()){
			string line;
			if (!getline(ifs, line)) break;
			if (line == "") continue;
			stringstream ss(line);
			vector<double> mfp_kappa;
			mfp_kappa.assign(2, 0);
			ss>>mfp_kappa[0]>>mfp_kappa[1];
			cumkappa.push_back(mfp_kappa);
		}
		ifs.close();
		totalK = cumkappa.back()[1];
		mfpks = generateRandomMFPs(cumkappa, mfpseg);
	}
	else if(isMFP){ // single MFP mode
		vector<double> mfp_kappa;
		mfp_kappa.assign(2, 1);
		mfp_kappa[0] = mfp;
		mfpks.push_back(mfp_kappa);
		mfpseg = 1;
		totalK = 1;
	}
	else{
		cerr<<"MFP mode has to be set! (MFP or MFPFILE)";
		return 1;
	}

	kfactor = 0;
	for (i = 0; i < mfpks.size(); i++){
		mfp = mfpks[i][0];
		kappa = mfpks[i][1];
		// Insertion

		int collect = 0;
		double cameratolight = camera.value - light.value;
		if (fabs(cameratolight) < mfp * minCameraDistance){
			int ratio = (int) (mfp*minCameraDistance / cameratolight);
			if (ratio > 1 && fabs(cameratolight) <= cameraMax )
			{
				camera.value = light.value + cameratolight * ratio;
				cout << "Increased the camera-light distance to: " << (camera.value - light.value) <<endl;
			}

		}
#pragma omp parallel reduction(+ : collect) private(phon, j)
		{
			phon = new Phonon(mfp);
#pragma omp for
			for (j = 0; j < nPhonon; j++)
			{
				phon->initPositionAndDirection(objects, &light);
				while(true){
					phon->scatterDistance(objects, &light, &camera);
					if (phon->cameraCapture(&light, &camera))
					{
						collect += 1;
						break;
					}
					if (phon->lightReturn(&light, &camera))
						break;
				}
			}
			delete phon;
		}
		double mfp_eff = 0.75 * (double) collect *  (camera.value - light.value) / nPhonon;
		cout <<"("<< i+1 <<"/"<<mfpks.size() << ")"<<"MFP: "<<mfp;
		cout<< ";  MFP(eff):"<<mfp_eff<<endl;
		kfactor += mfp_eff / mfp;
	}
	kfactor /= mfpseg;
	cout << "Final K: "<< totalK * kfactor <<"(W/m-K)" <<endl;
	// Read mfp-k accumulation file
	for (i = 0; i < objects.size(); i++)
		delete objects[i];
}
#endif

#ifdef DEBUG
int main(int argc, char** argv){
	srand((unsigned)time(0));
	Vec3 center(0,0,0);

	Cube box(center, 100, 100, 100, 0, ENTITY, 1);
	box.bconditions[2] = PERIODIC;
	box.bconditions[0] = PERIODIC;
	box.bconditions[1] = PERIODIC;

//	Cylinder box(center, 0.5, 1, 0, ENTITY, 0);
//	box.bcondition = PERIODIC;
	Phonon phon(1);
	vector <Object*> objects;
	objects.push_back(&box);
	long int i, N=100000;
	double distance;
	Vec3 distVector;
	double collect=0;
	EmitorCamera* emitor = new EmitorCamera(2, box.boundaries[0][2], 1);
	EmitorCamera* camera = new EmitorCamera(2, box.boundaries[1][2], -1);
	phon.light = emitor;
	phon.camera = camera;
#define RAY_TRACING
#ifdef RAY_TRACING
	int nmaxCol=0;
	vector<Object*> objs;
	objs.push_back(&box);
	for (i = 0; i < N; i++)
	{
		phon.initPositionAndDirection(objs, emitor);
		int ncol=0;
		while(true){
			distance = phon.scatterDistance(objects);
			if (phon.cameraCapture())
			{
				collect += 1;
				break;
			}
			if (phon.lightReturn())
			{
				break;
			}
			ncol++;
		}
		if (ncol>nmaxCol) nmaxCol = ncol;

	}
	cout<< 0.75 * (double) collect *  (camera->value - emitor->value) / N<<endl;
	cout <<"Maximum number of boundary collisions: "<<nmaxCol<<endl;
	delete emitor;
	delete camera;
}
#else
double dist = 0, dz = 0;

for (i = 0; i < N; i++)
{
	phon.initPositionAndDirection(&box);
	dz = phon.dir[2];
	distance = phon.scatterDistance(objects);
	dist += distance * dz * dz;
}
cout<<3 * dist / N<<endl;
}
#endif
#endif
