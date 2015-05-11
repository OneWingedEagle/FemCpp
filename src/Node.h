/*
 * Node.h
 *
 *  Created on: 2015/01/23
 *      Author: Hassan
 */

#include <vector>
#include "Vect.h"

#ifndef NODE_H_
#define NODE_H_

class Node {

public:
	Node();
	Node(char dim);

	void setPhi(double phi);

	void setCoord(Vect coord);

	void setCoord(double u,int i);

	void setIsPhiVar(bool b);

	void setIsPhiKnown(bool b);

	void setPBC(int nPBC);

	bool hasPBC();

	int getnPBC();


	void setMap(int map);

	bool getIsAknown();
	bool getIsPhiKnown();
	bool getIsPhiVar();
	double getPhi();
	Vect getCoord();
	double getCoord(int i);
	int getMap();

private:
	Vect F;
	Vect coord;
	double phi;
	bool onBound[];
	bool Aknown,phiKnown,phiVar,hasF,hasFms;
	bool common,sPBC,aPBC,inUse,rotor;
	char dim;
	int map;

};

#endif /* NODE_H_ */
