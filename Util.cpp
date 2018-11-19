#include "Util.h"
#include <cmath>
#include <list>

void zeroMat(double mat[3][3])
{
	for (int i = 0; i < 3; i++)
	{
		for (int j = 0; j < 3; j++)
		{
			mat[i][j] = 0;
		}
	}
}

void multiplyMat2Vec(double vec[3], double mat[3][3], double result[3])
{
	for (int i = 0; i < 3; i++) 
	{
		for (int j = 0; j < 3; j++)
		{
			result[i] += mat[i][j] * vec[j];
		}
	}
}

double angleBetweenVectors(double ovec1[3], double nvec1[3], double ovec2[3], double nvec2[3])
{
	double dotResult = 0.0;
	double tvec1 = 0.0;
	double tvec2 = 0.0;
	
	// compute dot product
	for (int i = 0; i < 3; i++)
	{
		double dvec1 = nvec1[i] - ovec1[i];
		double dvec2 = nvec2[i] - ovec2[i];
		dotResult += dvec1 * dvec2;		
		tvec1 += dvec1*dvec1;
		tvec2 += dvec2*dvec2;
	}
	double multipliedLengths = sqrt(tvec1) * sqrt(tvec2);
	
	return acos(dotResult/multipliedLengths);
}

