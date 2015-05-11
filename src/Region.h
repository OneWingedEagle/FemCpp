/*
 * Region.h
 *
 *  Created on: 2015/01/23
 *      Author: Hassan
 */
#include <vector>
#include "Vect.h"



#ifndef REGION_H_
#define REGION_H_

class Region{
public:
	bool hasJ;
	Region();

	Region(int dim);

	void setJ0(Vect J);

	void setMu(Vect  Mu);

	void setSigma(Vect sigma);

	void setHasJ(bool b);

	void setIsNonlin(bool b);

	void setIsConduct(bool b);

	Vect getNu();

	Vect getSigma();

	bool getHasJ();

	bool getIsNonlin();

	bool getIsConduct();

	bool getIsRotor();

private:
	int dim,BHnumber;
	string material;
	string name;
	Vect mur,sigma,nu;
	Vect J,M;
	static const double  mu0=4e-7*3.14159265;
	int firstElement,lastElement;
	bool isConductor;
	double Jz;


};

#endif /* REGION_H_ */
