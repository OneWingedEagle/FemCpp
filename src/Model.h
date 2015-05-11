/*
 * Model.h
 *
 *  Created on: 2015/01/23
 *      Author: Hassan
 */

#include <vector>
#include "Vect.h"
#include "Region.h"
#include "Element.h"
#include "Node.h"
#include "Edge.h"

#include "BoundarySet.h"
#include "MagMatrix.h"
#include "SpMat.h"
#include "SpMatSolver.h"
#include "FEMSolver.h"


#ifndef MODEL_H_
#define MODEL_H_

class Model {


private:

	 MagMatrix magMat;
	 SpMatSolver solver;
	 FEMsolver femSolver;
	 BoundarySet bcSet;


	int dim,iterMax,nonLinIterMax,cpb,nRotReg;
	double cpm,Rg;
	int numberOfRegions, numberOfNodes,numberOfElements,numberOfEdges,nXYedges;
	int nBlocks,nBH,nLam;
	int nElVert,nBoundary;
	int nElEdge;
	double spaceBoundary[];
	int BCtype[];
	Vect* diricB;
	Region* region;
	Node* node;
	Element* element;
	Edge* edge;

	vector <int> phiUnknownIndex,unknownPhiNumber,A_unknownIndex,
	unknownAnumber,nodeVarIndex,varNodeNumber;
	vector<double> knownEdgeValue[];

	double scaleFactor;
	int numberOfUnknownEdges,numberOfMappedEdges,numberOfKnownEdges,numberOfVarNodes;
	int numberOfKnownPhis,numberOfUnknowns,analysisMode;

	bool hasJ,hasM,nonLin,hasPBC,rotIndex;
	int nEdEd,nNodNod,nEdNod,nNodEd,nEdgeHasJ;

	char elCode;
	double nu0;
	double rotAng,rotSubAng,rotStep,meshAngStep,freq,dt,errMax;
	int nTsteps,nBegin,nEnd,nInc,nRotorElements,tagx;
	int coordCode,timeIntegMode,eddyTimeIntegMode;


	Vect xp,Ci,Cj;
	string meshFilePath,dataFilePath,fluxFilePath,resultFolder;
	double H2[],H3[];
	double C[],Cj2d[];

	 void setElementB(int i);

	 Vect getElementB(int i, Vect* rotNe);


public:

	Model();
	Model(string bun);
	Model(int nRegions, int nElements,int nNodes, char elCode);
	string elType;
	SpMat Hs,Ls,Ss;
	Vect b,HpAp,HkAk;
	bool circuit,stranded;
	int PBCpair[];


	 void setFemCalc();

	 void loadMesh(string bunFilePath);

	 void loadData(string dataFilePath);

	 double edgeLength(int i);

	 void setHasJ();


	 void setFemSolver();


	 void setMagMat(int iter);

	 void writeMesh(string bunFilePath);

	 void writeData(string dataFilePathOut) ;


	 Node* elementNodes(int i);

	 Edge* elementEdges(int i);

	 void setSolvedAL(Vect x);


	 void setElType(string type);

	 void setMagBC();

	 void setEdge();

	 void setBounds();


	void setElementJ0(int i);

	 void setJ0();


	 void setFreq(double f);
	 void setDt(double dt);

	 void setnTsteps(int N);


	 void setJe();

	 void setElementJe(int i);


	 Vect getElementA(int ie);

	 double getElementPhi(int ie);


	 Vect getElementdA(Node *vertexNode,double *Ae);

	 void setElementsParam();


	 void scaleKnownEdgeAL(double c);


	 void saveAp();



	 Vect solveMagLin(int n);


	 void setNodePhi(Vect x);

	 void setSolution(Vect x);

	 Vect getUnknownA();

	 Vect getUnknownAp();

	 int *getRegNodes(int ir);

	 void setB();

	 Vect getBAt(Vect P);

	 double *getApAnAt(Vect P);

	 void setJ(double t);


	 void writeB(string file);


	 Vect getElementCenter(int ie);


};

#endif /* MODEL_H_ */
