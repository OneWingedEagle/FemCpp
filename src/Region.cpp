/*
 * Region.cpp
 *
 *  Created on: 2015/01/23
 *      Author: Hassan
 */

#include "Region.h"


 Region::Region(){

	};

 Region::Region(int dim){
		this->dim=dim;
	};

	void  Region::setJ0(Vect J){
		this->J=J.copy();
	}
	void  Region::setMu(Vect  Mu){
		this->mur=Mu.copy();

		this->nu=(Mu*mu0).inv();
	}
	void Region::setSigma(Vect sigma){
		this->sigma=sigma.copy();
	}
	void Region::setHasJ(bool b){
		hasJ=b;
	}

	void Region::setIsConduct(bool b){
		isConductor=b;
	}

	Vect Region::getNu(){
		return nu;
	}
	Vect Region::getSigma(){;
	return sigma;
	}

	bool Region::getHasJ(){
		return hasJ;
	}

	bool Region::getIsConduct(){
		return isConductor;
	}


