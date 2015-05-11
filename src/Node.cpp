/*
 * Node.cpp
 *
 *  Created on: 2015/01/23
 *      Author: Hassan
 */

#include "Node.h"

	Node::Node(){};

	Node::Node(char dim){
		this->dim=dim;
	}
	void Node::setPhi(double phi){
		this->phi=phi;

	}
	void Node::setCoord(Vect coord){
		this->coord=coord.copy();
	}

	void Node::setCoord(double u,int i){
		coord.set(u,i);
	}
	void Node::setIsPhiVar(bool b){
		phiVar=b;
	}
	void Node::setIsPhiKnown(bool b){
		phiKnown=b;
	}

	void Node::setPBC(int nPBC){
		this->sPBC=(nPBC==1);
		this->aPBC=(nPBC==-1);

	}

	bool Node::hasPBC(){
		return (this->sPBC || this->aPBC);

	}

	int Node::getnPBC(){
		if(this->aPBC) return -1;
		 return 1;

		}


	void Node::setMap(int map){
		this->map=map;
	}
	bool Node::getIsAknown(){
		return this->Aknown;
	};

	bool Node::getIsPhiKnown(){
		return phiKnown;
	};

	bool Node::getIsPhiVar(){
		return phiVar;
	}

	double Node::getPhi(){
		return phi;
	}
	Vect Node::getCoord(){
		return coord;
	}
	double Node::getCoord(int i){
		return coord.get(i);
	}
	int Node::getMap(){
		return map;
	}

