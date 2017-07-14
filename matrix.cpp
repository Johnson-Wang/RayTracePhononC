#include <iostream>
#include <vector>
#include <cmath>
#include <algorithm>
#include "matrix.h"
using namespace std;

Matrix::Matrix() {
	mat.assign(3,vector<double>(3,0));
	for (int i=0;i<3;i++) mat[i][i] = 1;
}

Matrix::Matrix(const Matrix& t) {
	mat.assign(3,vector<double>(3,0));
	for (int i=0;i<3;i++) mat[i][i] = 1;
	for (int i=0;i<3;i++) {
		for (int j=0;j<3;j++)
			mat[i][j] = t.mat[i][j];
	}
}

Matrix::Matrix(Vec3 u, Vec3 v, Vec3 w) {
	mat.assign(3,vector<double>(3,0));
	mat[0][0] = u[0]; mat[0][1] = u[1]; mat[0][2] = u[2];
	mat[1][0] = v[0]; mat[1][1] = v[1]; mat[1][2] = v[2];
	mat[2][0] = w[0]; mat[2][1] = w[1]; mat[2][2] = w[2];
}

Matrix::Matrix(const double V[9]){
	mat.assign(3,vector<double>(3,0));
	for(int i=0; i<9; i++){
		mat[i/3][i%3] =  V[i];
	}
}

Matrix Matrix::dot(const Matrix& v) const {
	Matrix temp;
	for (int i=0;i<3;i++) {
		for (int j=0;j<3;j++) {
			if (i==j) temp.mat[i][j]--;
			for (int k=0;k<3;k++) {
				temp.mat[i][j] += mat[i][k]*v.mat[k][j];
			}
		}
	}
	return temp;
}

double Matrix::det() const{
	return mat[0][0] * (mat[1][1] * mat[2][2] - mat[1][2] * mat[2][1])
	    + mat[0][1] * (mat[1][2] * mat[2][0] - mat[1][0] * mat[2][2])
	    + mat[0][2] * (mat[1][0] * mat[2][1] - mat[1][1] * mat[2][0]);
}

Matrix Matrix::inverse() const{
	  Matrix inv;
	  double d = det();
	  inv.mat[0][0] = (mat[1][1] * mat[2][2] - mat[1][2] * mat[2][1]) / d;
	  inv.mat[1][0] = (mat[1][2] * mat[2][0] - mat[1][0] * mat[2][2]) / d;
	  inv.mat[2][0] = (mat[1][0] * mat[2][1] - mat[1][1] * mat[2][0]) / d;
	  inv.mat[0][1] = (mat[2][1] * mat[0][2] - mat[2][2] * mat[0][1]) / d;
	  inv.mat[1][1] = (mat[2][2] * mat[0][0] - mat[2][0] * mat[0][2]) / d;
	  inv.mat[2][1] = (mat[2][0] * mat[0][1] - mat[2][1] * mat[0][0]) / d;
	  inv.mat[0][2] = (mat[0][1] * mat[1][2] - mat[0][2] * mat[1][1]) / d;
	  inv.mat[1][2] = (mat[0][2] * mat[1][0] - mat[0][0] * mat[1][2]) / d;
	  inv.mat[2][2] = (mat[0][0] * mat[1][1] - mat[0][1] * mat[1][0]) / d;
	  return inv;
}

Vec3 Matrix::dot(Vec3& v) const{
	Vec3 v2;
	int i,j;
	for (i = 0; i < 3; i++){
		for (j = 0; j < 3; j++){
			v2[i] += mat[i][j] * v[j];
		}
	}
	return v2;
}

Matrix& Matrix:: operator = (const Matrix& v){
	for (int i = 0; i < 3; i++)
		for (int j = 0; j < 3; j++)
			mat[i][j] = v.mat[i][j];
	return *this;
}

Matrix Matrix::transpose() {
	Matrix temp;
	for (int i=0;i<3;i++) {
		for (int j=0;j<3;j++) {
			temp.mat[i][j] = mat[j][i];
		}
	}
	return temp;
}

void Matrix::swapRows(int i, int j) {
	for (int k=0;k<3;k++)
		swap(mat[i][k],mat[j][k]);
}

void Matrix::divideRow(int i, double temp) {
	for (int j=0;j<3;j++)
		mat[i][j]/=temp;
}

void Matrix::subRow(int i, int j, double temp) {
	for (int k=0;k<3;k++)
		mat[i][k] -= mat[j][k]*temp;
}

int Matrix::getRow(int i) {
	for (int j=i;j<3;j++) {
		if (mat[j][i]!=0)
			return j;
	}
	return NULL;
}

void Matrix::print() {
	for (int i=0;i<3;i++) {
		for (int j=0;j<3;j++) {
			cout<<mat[i][j]<<" ";
		}
		cout<<endl;
	}
}

Matrix rotationFromTwoVectors(Vec3 from, Vec3 to){ // from and to are both unit vectors
	double fromdotto = from.dot(to);
	Vec3 v = (to - from * fromdotto).normalize();
	Vec3 w = to.cross(from);
	Matrix F(from, v, w);
	double g[9] = {fromdotto, -w.length(), 0, w.length(), fromdotto, 0, 0, 0, 1};
	Matrix G(g); // Rotation matrix in the coordinate formed by u, v, w
	return F.transpose().dot(G.dot(F)); // Coordinate transformation matrix
}


