/*
 * Edge.h
 *
 *  Created on: 2015/01/23
 *      Author: Hassan
 */

#ifndef EDGE_H_
#define EDGE_H_

class Edge {

private:
	char direction;
	double A,Ap,length;
	int endNodeNumber[2];
	bool edgeKnown,hasJ,sPBC,aPBC,common;
	int map;

public:
	Edge(int n1,int n2);

	void setKnownA(double A);
	 void setSolvedAL(double A);
	 void saveAp();
	 void setLength(double length);
	 void setDirection(char i);
	 double getDiffA();
	 void setA(double A);

	 double getA();

	 void setPBC(int nPBC);

	 bool hasPBC();

	 int getnPBC();


};

#endif /* EDGE_H_ */
