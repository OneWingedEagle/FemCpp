/*
 * Element.h
 *
 *  Created on: 2015/01/23
 *      Author: Hassan
 */

#include <vector>
#include "Vect.h"

#ifndef ELEMENT_H_
#define ELEMENT_H_

class Element{
public:
	Element();
	Element(char elemCode);
	void setJ0(Vect J);
	void setNu(Vect  Nu);
	void setSigma(Vect sigma);
	void setHasJ(bool b);
	void setIsNonlin(bool b);
	void setIsConduct(bool b);

	void setHasM(bool b);
	void setIsRotor(bool b);
	void setB(Vect    B);
	void setJe(Vect    Je);
	void setJ(Vect J);

	void setRegion(int nr);

	void setEdgeNumb(int ne[]);

	void setVertNumb(int ne[]);




	Vect getNu();

	Vect getSigma();

	bool getHasJ();

	bool getIsNonlin();

	bool getIsConduct();

	bool getIsRotor();

	Vect getJ0();

	Vect getB();

	Vect getJe();

	int getRegion();

	vector<int> getVertNumb();

	vector<int> getEdgeNumb();


private:
	char dim,elCode;
	vector<int> vertNumb, edgeNumb, edgeXYNumb;
	int nRegion;
	Vect nu,sigma;
	Vect B,J,Je;
	bool hasJ,isConductor;
};


#endif /* ELEMENT_H_ */
