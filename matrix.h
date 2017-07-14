#ifndef MATRIX_H
#define MATRIX_H

#include <iostream>
#include <vector>
#include "vec3.h"
using namespace std;
enum Affine {ROTATEX,ROTATEY,ROTATEZ,SCALE,TRANSLATE,SHEAR};

class Matrix {
public:
	vector<vector<double> > mat;
	Matrix();
	Matrix(const Matrix& t);
	Matrix(const double* V);
	Matrix(Vec3 u, Vec3 v, Vec3 w);
	void swapRows(int i, int j);
	void divideRow(int i, double temp);
	int getRow(int i);
	void subRow(int i, int j, double temp);
	Matrix dot(const Matrix& v) const;
	Vec3 dot(Vec3& v) const;
	Matrix& operator = (const Matrix& v);
	double det() const;
	Matrix inverse() const;
	void print();
	Matrix transpose();
};

//The rotation matrix that maps from the vector 'from' to vector 'to'
Matrix rotationFromTwoVectors(Vec3 from, Vec3 to);  // from and to are unit vectors

#endif
