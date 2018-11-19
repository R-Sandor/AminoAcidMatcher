#ifndef UTIL_H
#define UTIL_H

#include <list>

extern void zeroMat(double mat[3][3]);

extern void multiplyMat2Vec(double vec[3], double mat[3][3], double result[3]);

extern double angleBetweenVectors(double ovec1[3], double nvec1[3], double ovec2[3], double nvec2[3]);

template <class T> inline 
T& getListElement(std::list<T> & searchList, int index)
{
	auto  iterator = searchList.begin();
	for (int i = 0; i < index; i++) iterator++;
	return *iterator;
};

#endif
