/*
 * Element.cpp
 *
 *  Created on: 2015/01/23
 *      Author: Hassan
 */

#include "Element.h"


	Element::Element(){};

	Element::Element(char elemCode){
		elCode=elemCode;

		if(elCode=0){

			vertNumb.assign (3,0);
			edgeNumb.assign (3,0);

			dim = 2;
		}
		else if(elCode=1){
			vertNumb.assign (3,0);
			edgeNumb.assign (3,0);
			dim = 2;
		}
		else if(elCode=2){
			vertNumb.assign (4,0);
			edgeNumb.assign (6,0);
			dim = 3;
		}
		else if(elCode=3){
			vertNumb.assign (6,0);
			edgeNumb.assign (9,0);
			dim = 3;
		}
		else if(elCode=4){
			vertNumb.assign (8,0);
			edgeNumb.assign (12,0);
			dim = 3;
		}
	}

	void Element::setJ0(Vect J){
		this->J=J.copy();
	}
	void Element::setNu(Vect  Nu){

		this->nu=Nu.copy();
	}
	void Element::setSigma(Vect sigma){
		this->sigma=sigma.copy();
	}
	void Element::setHasJ(bool b){
		hasJ=b;
	}

	void Element::setIsConduct(bool b){
		isConductor=b;
	}


	void Element::setB(Vect    B){
		this->B=B.copy();
	}
	void Element::setJe(Vect    Je){
		this->Je=Je.copy();;
	}
	void Element::setJ(Vect J){
		this->J=J.copy();;
	}

	void Element::setRegion(int nr){
		nRegion=nr;
	}
	void Element::setEdgeNumb(int ne[]){
		for(int i=0;i<(int)edgeNumb.size();i++){
			edgeNumb[i]=ne[i];
		}
	}
	void Element::setVertNumb(int ne[]){
		for(int i=0;i<(int)vertNumb.size();i++){
			vertNumb[i]=ne[i];
		}
	}




	Vect Element::getNu(){
		return nu;
	}
	Vect Element::getSigma(){;
	return sigma;
	}

	bool Element::getHasJ(){
		return hasJ;
	}

	bool Element::getIsConduct(){
		return isConductor;
	}


	Vect Element::getJ0(){
		return J;
	}
	Vect Element::getB(){
		return B;
	}
	Vect Element::getJe(){
		return Je;
	}

	int Element::getRegion(){
		return nRegion;
	}
	vector<int> Element::getVertNumb(){
		return vertNumb;
	}

	vector<int> Element::getEdgeNumb(){
		return edgeNumb;
	}



