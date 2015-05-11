/*
 * Edge.cpp
 *
 *  Created on: 2015/01/23
 *      Author: Hassan
 */

#include "Edge.h"

Edge::Edge(int n1,int n2)
	{

	if(n1<n2){
	endNodeNumber[0]=n1;
	endNodeNumber[1]=n2;
	}
	else{
		endNodeNumber[0]=n2;
		endNodeNumber[1]=n1;
	}

	}

/*	void Edge::setKnownA(double Ax){
		edgeKnown=true;
		A=Ax;

	}

	 void Edge::setSolvedAL(double Ax){
		 A=Ax;

	}

	 void Edge::saveAp(){
		Ap=A;
	}

	 void Edge::setLength(double l){

		length=l;
	}

	 void Edge::setDirection(char i){

		direction=i;
	}

	 double Edge::getDiffA(){
		return A-Ap;
	}

	 void Edge::setA(double Ax) {

		 A=Ax;

	}

	 double Edge::getA() {

		return A;

	}

	 void Edge::setPBC(int nPBC){
		sPBC=(nPBC==1);
		aPBC=(nPBC==-1);

	}

	 bool Edge::hasPBC(){
		return (sPBC || aPBC);

	}

	 int Edge::getnPBC(){

		if(aPBC) return -1;
		 return 1;

		}*/


