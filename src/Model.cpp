/*
 * Model.cpp
 *
 *  Created on: 2015/01/23
 *      Author: Hassan
 */

#include "Model.h"
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


	Model::Model(){};

	Model::Model(string bun){
		Model();

		loadMesh(bun);
	}

	Model::Model(int nRegions, int nElements,int nNodes, char elCode){

		this->numberOfRegions=nRegions;
		this->numberOfElements=nElements;

		this->numberOfNodes=nNodes;

		//setElType(elType);

		region=new Region[nRegions];
		for(int i=1;i<=this->numberOfRegions;i++)
			region[i]=Region(dim);

		element=new Element[this->numberOfElements+1];
		for(int i=1;i<=this->numberOfElements;i++)
			element[i]=Element(elCode);

		node=new Node[this->numberOfNodes+1];
		for(int i=1;i<=this->numberOfNodes;i++)
			node[i]=Node(dim);
	}



/*	 void Model::loadData(string dataFilePath)
	{
		loader.loadData(this,dataFilePath);
		this->elementMatrix=ElementMatrix(this);
		setMagMech();

	}*/


	 void Model::loadMesh(string bun){
		//loader.loadMesh(this);

	}

/*
	 double Model::edgeLength(int i){
		double length=node[edge[i].endNodeNumber[1]].getCoord().sub(node[edge[i].endNodeNumber[0]].getCoord()).norm();
		return length;
	}



	 void Model::setHasJ(){
		for(int i=1;i<=numberOfRegions;i++)
			if(region[i].getHasJ()){
				hasJ=true;
				break;
			}
	}



	 void Model::femSolver(){
		this->femSolver=new FEMsolver(this);

	}


	 void Model::setMagMat(int iter){
		magMat.setMagMat(this,iter);

	}


	 void Model::writeMesh(String bunFilePath){
		writer.writeMesh(this,bunFilePath);
	}


	 Node* Model::elementNodes(int i){
		Node* elementNode=new Node[nElVert];
		int* vertNumb=element[i].getVertNumb();
		for(int j=0;j<nElVert;j++){

			elementNode[j]=node[vertNumb[j]];
		}

		return elementNode;
	}

	 Edge* Model::elementEdges(int i){
		Edge* elementEdge=new Edge[nElEdge];
		int* edgeNumb=element[i].getEdgeNumb();
		for(int j=0;j<nElEdge;j++)
			elementEdge[j]=edge[edgeNumb[j]];
		return elementEdge;
	}



	 void Model::setSolvedAL(Vect x){

		for(int i=1;i<=numberOfUnknownEdges;i++){


			edge[unknownEdgeNumber[i]].setSolvedAL(x.el[i-1]);

		}


		if(this->hasPBC || this->hasTwoNodeNumb){

			for(int i=1;i<=numberOfEdges;i++){
				if(edge[i].map>0 && !edge[i].edgeKnown)
				{
					if(edge[i].aPBC){

						edge[i].setSolvedAL(this->cpb*edge[edge[i].map].A);

					}


					else{

						edge[i].setSolvedAL(edge[edge[i].map].A);


					}
				}

			}
		}




	}


	 void Model::setElementB(int i){


			Node* vertexNode=elementNodes(i);
			Vect zero=new Vect(3);
			Mat jac=femCalc.jacobian(vertexNode,zero);
			Vect B;
			Vect* rotNe=femCalc.rotNe(jac,zero);
			B=getElementB(i,rotNe);

		// B=new Vect(0,0,getElementA(i).el[1]);

			element[i].setB(B);



	}




	 void Model::setElType(String type){
		elType=type;
		if(type.equals("triangle") ){
			elCode=0;
			nElVert=3;
			nElEdge=3;
			this->dim=2;
		}
		else if(type.equals("quadrangle") ){
			elCode=1;
			nElVert=4;
			nElEdge=4;
			this->dim=2;
		}
		else if(type.equals("tetrahedron") ){
			elCode=2;
			nElVert=4;
			nElEdge=6;
			dim=3;
		}
		else if(type.equals("prism") ){
			elCode=3;
			nElVert=6;
			nElEdge=9;
			dim=3;
		}
		else if(type.equals("hexahedron") ){
			elCode=4;
			nElVert=8;
			nElEdge=12;
			dim=3;
		}

		else if(type.equals("pyramid") ){
			elCode=5;
			nElVert=5;
			nElEdge=8;
			dim=3;
		}
		nBoundary=2*dim;
	}




	 Vect Model::getElementB(int i, Vect* rotNe){

		Edge* edge=elementEdges(i);
		Vect B=new Vect(dim);
		for(int j=0;j<nElEdge;j++)		{

			B=B.add(rotNe[j].times(edge[j].A));
		}

		return B;

	}



	 double Model::getElementVolume(int i){

		double vol=0;
		Node* vertexNode=elementNodes(i);
		Mat jac;
		double detJac,ws;
		Vect localCo=new Vect(dim);
		ws=8;
		jac=femCalc.jacobian(vertexNode,localCo);
		detJac=abs(jac.determinant());

		vol=detJac*ws;

		return vol;

	}




	 void Model::setMagBC(){
		this->bcSet.setMagBC(this);
	}


	 void Model::setEdge(){

		EdgeSet edgeSet=new EdgeSet();
		edgeSet.setEdge(this);

	}


	 void Model::setBounds(){
		if(this->coordCode==1)
		this->bcSet.setSliceBounds(this);
		else if(coordCode==0){

			double spb[2*this->dim];


			for(int i=0;i<this->dim;i++){
				spb[2*i]=1e10;
				spb[2*i+1]=-1e10;
			}


			for(int i=1;i<=this->numberOfNodes;i++){

				if(this->node[i].getCoord(0)<spb[0]) spb[0]=this->node[i].getCoord(0);
				else if(this->node[i].getCoord(0)>spb[1]) spb[1]=this->node[i].getCoord(0);

				if(this->node[i].getCoord(1)<spb[2]) spb[2]=this->node[i].getCoord(1);
				else if(this->node[i].getCoord(1)>spb[3]) spb[3]=this->node[i].getCoord(1);
				if(this->dim==3){
					if(this->node[i].getCoord(2)<spb[4]) spb[4]=this->node[i].getCoord(2);
					else if(this->node[i].getCoord(2)>spb[5]) spb[5]=this->node[i].getCoord(2);
				}


				}

			this->spaceBoundary=spb;

		}
	}


	 void Model::setNodeOnBound(){
		this->bcSet.setNodeOnBound(this);
	}


	 void Model::setElementJ0(int i){

		if(this->elCode!=4) return;

			Node* vertexNode=elementNodes(i);
			Vect zero=new Vect(3);
			Mat jac=femCalc.jacobian(vertexNode,zero);
			Vect* rotNe=femCalc.rotNe(jac,zero);

			Edge* edge=elementEdges(i);
			Vect J=new Vect(dim);
			for(int j=0;j<nElEdge;j++)		{
				J=J.add(rotNe[j].times(edge[j].T));
			}

			element[i].setJ(J);



	}



 void Model::setFreq(double f){
		this->freq=f;

	}
	 void Model::setDt(double dt){

		this->dt=dt;
	}

	 void Model::setnTsteps(int N){

		this->nTsteps=N;
	}





	 void Model::setJe(){
		double Jn2=0,Jmax2=0,Jmin2=0;
		for(int i=1;i<=numberOfElements;i++){

			if(!element[i].isConductor()) continue;

			setElementJe(i);

			Vect Je=element[i].getJe();

			if(element[i].isConductor()){
				Jn2=Je.dot(Je);
				if(Jn2>Jmax2)
					Jmax2=Jn2;
				if(Jn2<Jmin2)
					Jmin2=Jn2;

			}

		}

		Jmax=sqrt(Jmax2);
		Jmin=sqrt(Jmin2);
	}

	 void Model::setElementJe(int i){

		Node* vertex=elementNodes(i);
		Edge* edge=elementEdges(i);
		double dAe[nElEdge];
		Vect dA=new Vect(dim);
		for(int j=0;j<nElEdge;j++)
			dAe[j]=edge[j].getDiffA();

		dA=getElementdA(vertex,dAe);

		double rdt=1.0/dt;

		if(analysisMode==1){
			element[i].setJe(dA.times(element[i].getSigma()).times(-rdt));
		}
		else if(analysisMode==2){
			double nodePhi[nElVert];
			Vect gradPhi(dim);

			for(int j=0;j<nElVert;j++)
				nodePhi[j]=vertex[j].getPhi();
			gradPhi=femCalc.gradPhi(vertex,nodePhi);
			element[i].setJe((dA.times(rdt).add(gradPhi)).times(element[i].getSigma()).times(-1));
		}

	}


	 Vect Model::getElementA(int ie){

		Vect A=new Vect(dim);
		Vect zero=new Vect(dim);
		Node* vertexNode=this->elementNodes(ie);
		Edge* edge=this->elementEdges(ie);
		Mat jac=femCalc.jacobian(vertexNode,zero);
		Vect* Ne=femCalc.Ne(jac,zero);

		for(int j=0;j<nElEdge;j++)	{
			A= A.add(Ne[j].times(edge[j].A));
		}
		return  A;

	}

	 double Model::getElementPhi(int ie){

		double phi=0;
		Vect zero=new Vect(dim);
		Node* vertexNode=this->elementNodes(ie);
		double* N=femCalc.N(zero);

		for(int j=0;j<nElVert;j++)	{
			phi+=N[j]*vertexNode[j].getPhi();
		}
		return  phi;

	}




	 Vect Model::getElementdA(Node* vertexNode,double* Ae){

		Vect dA=Vect(dim);
		Vect zero=Vect(dim);
		Mat jac=femCalc.jacobian(vertexNode,zero);
		Vect* Ne=femCalc.Ne(jac,zero);

		for(int j=0;j<nElEdge;j++)	{
			dA= dA.add(Ne[j].times(Ae[j]));
		}
		return  dA;

	}





	 void Model::setElementsParam(){



		for(int ir=1;ir<=numberOfRegions;ir++){

			boolean regCond=region[ir].isConductor;

			for(int i=region[ir].getFirstEl();i<=region[ir].getLastEl();i++){

				element[i].setRo(region[ir].getRo());
				element[i].setNu(region[ir].getNu());
				element[i].setSigma(region[ir].getSigma());

				if(regCond)
					element[i].setSigma(region[ir].getSigma());

				element[i].setRegion(ir);

				element[i].setYng(region[ir].getYng());

				element[i].setPois(region[ir].getPois());

				element[i].setShear(region[ir].getShear());


			}
		}



		if(elType.equals("triangle")) elCode=0;
		else if(elType.equals("quadrangle")) elCode=1;
		else if(elType.equals("tetrahedron")) elCode=2;
		else if(elType.equals("prism")) elCode=3;
		else if(elType.equals("hexahedron")) elCode=4;

	}





	 void Model::scaleKnownEdgeAL(double c){
		for(int i=1;i<=numberOfKnownEdges;i++)
			if(knownEdgeValue[i]!=0)
				edge[knownEdgeNumber[i]].setSolvedAL(knownEdgeValue[i]*c);
	}


	 void Model::saveAp(){
		for(int i=1;i<=this->numberOfEdges;i++){
			edge[i].saveAp();


		}

	}



	 Vect Model::solveMagLin(int n){

		return femSolver.solveMagLin(this, n);
	}



	 void Model::setNodePhi(Vect x){
		int nodeNumber;

		for(int i=1;i<=numberOfVarNodes;i++){
			nodeNumber=varNodeNumber[i];
			if(nodeNumber>0){
				util.pr(x.el[i+numberOfUnknownEdges-1]);

				node[varNodeNumber[i]].setPhi(x.el[i+numberOfUnknownEdges-1]);
			}
		}

	}



	 void Model::setSolution(Vect x){

		setSolvedAL(x);

		if(analysisMode==2)
			setNodePhi(x);

	}

	 Vect Model::getUnknownA(){
		Vect x(this->numberOfUnknownEdges);
		for(int i=0;i<x.length;i++){
			x.el[i]=edge[unknownEdgeNumber[i+1]].A;
		}

		return x;
	}

		 Vect Model::getUnknownAp(){
		Vect x(this->numberOfUnknownEdges);
		for(int i=0;i<x.length;i++)
			x.el[i]=edge[unknownEdgeNumber[i+1]].Ap;

		return x;
	}


	 void Model::setB(){

		double Bn2,Bmax2=0,Bmin2=0;

		for(int i=1;i<=numberOfElements;i++){
			setElementB(i);

			Bn2=element[i].getB().dot(element[i].getB());
			if(Bn2>Bmax2)
				Bmax2=Bn2;
			if(Bn2<Bmin2)
				Bmin2=Bn2;}

		Bmax=sqrt(Bmax2);
		Bmin=sqrt(Bmin2);



	}




	 void Model::setJ(double t){

	for(int ir=1;ir<=numberOfRegions;ir++){

			if(!region[ir].circuit){

				if(region[ir].hasJ  )
					setJfJ(ir,t);
			else if(this->region[ir].terminalVoltage0!=0){
				setJfV(ir,t);
			}

		}
	}


	if(this->threePhaseRegs[0]!=0 && this->threePhaseRegs[1]!=0 && this->threePhaseRegs[2]!=0)
		setJ3phaseStrandedCircuit(t);
	else if(this->circuit)
		setJStrandedCircuit(t);

	}

	 void Model::setJfJ(int ir,double t){


		Vect J=new Vect();

		J=this->region[ir].getJ().times(cos(this->region[ir].omega*t+this->region[ir].phase0+this->region[ir].beta));



		for(int i=region[ir].getFirstEl();i<=region[ir].getLastEl();i++){


			this->element[i].setJ(J);


		}
	}


	 void Model::writeB(String file){
		writer.writeB(this,file);

	}


*/










