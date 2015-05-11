/*
 * MagMatrix.h
 *
 *  Created on: 2015/01/23
 *      Author: Hassan
 */

#ifndef MAGMATRIX_H_
#define MAGMATRIX_H_

class MagMatrix {

	private:
 int nEdEd,nNodEd,nEdNod;
 int nNodNod,numberOfRegions;
 int nElVert,nElEdge,dim;
 //Calculator calc;

public:


	  MagMatrix(){}

	/*   MagMatrix(Model model){


		this.nElVert=model.nElVert;
		this.nElEdge=model.nElEdge;
		this.dim=model.dim;

		this.nEdEd=model.nEdEd;
		this.nNodEd=model.nNodEd;
		this.nNodNod=model.nNodNod;
		this.nEdNod=model.nEdNod;
		this.numberOfRegions=model.numberOfRegions;

		this.calc=new Calculator(model);
	}


	 void setMagMat(Model model, int step){

		if(model.analysisMode>0 && model.circuit && model.eddyTimeIntegMode<=-2){

			setMagMatEdgeCircuit(model,step);
			return;
		}


		boolean fillSs=(model.analysisMode>0 && model.Ss==null);


		double eps=1e-6,cPB=model.cpb;
		boolean nonLinear,eddy;
		double[][] H1=new double[this.nElEdge][this.nElEdge];
		Vect J=new Vect();
		double[] Cj=new double[this.nElEdge];
		int m,columnEdgeNumb,columnIndex,matrixRow=0, rowEdgeNumb,nBH=0,nLam=0,ext=6;
		boolean hasJ,hasM;

		int[] nz=new int[model.numberOfUnknowns];
		int[] nzs=new int[model.numberOfUnknowns];
		model.Hs=new SpMat(model.numberOfUnknowns);

		if(fillSs){
		model.Ss=new SpMat(model.numberOfUnknowns);
		for(int i=0;i<model.numberOfUnknowns;i++){

			model.Ss.row[i]=new SpVect(model.numberOfUnknowns,this.nEdEd);
		}
		}

		for(int i=0;i<model.numberOfUnknownEdges;i++){

			model.Hs.row[i]=new SpVect(model.numberOfUnknowns,this.nEdEd);
		}




		if(model.eddyTimeIntegMode==1){
		if(model.b!=null)
		model.bT=model.b.deepCopy();
		else
			model.bT=null;
		}

		model.b=new Vect(model.numberOfUnknowns);
		model.HkAk= new Vect(model.numberOfUnknowns);



		for(int ir=1;ir<=this.numberOfRegions;ir++){

			nBH=model.region[ir].BHnumber;
			nLam=model.region[ir].lamBNumber;
			for(int i=model.region[ir].getFirstEl();i<=model.region[ir].getLastEl();i++){

				hasJ=model.element[i].hasJ();
				if(hasJ){
					J=model.element[i].getJ();

				}



				hasM=model.element[i].hasM();
				nonLinear=model.element[i].isNonlin();

				if(fillSs && step==0 &&  model.element[i].isConductor())
					eddy=true;
				else
					eddy=false;

				H1=this.calc.He(model,nBH,nLam,i,nonLinear,eddy,hasJ,hasM);


				if(hasJ ){

					for(int j=0;j<this.nElEdge;j++){

						if(model.dim==2)
						Cj[j]=J.el[2]*model.Cj2d[j];
						else
							Cj[j]=J.dot(model.Cj[j]);
					}
				}

				int[] edgeNumb=model.element[i].getEdgeNumb();

				for(int j=0;j<this.nElEdge;j++){
					rowEdgeNumb=edgeNumb[j];

					if(model.edge[rowEdgeNumb].edgeKnown ) continue;


					matrixRow=model.edgeUnknownIndex[rowEdgeNumb]-1;


					//===========  right-hand side
					if( hasM   ){
						model.b.el[matrixRow]+=model.C[j];
					}

					if(hasJ ){

						model.b.el[matrixRow]+=Cj[j];

					}

					for(int k=0;k<this.nElEdge;k++){

						columnEdgeNumb=edgeNumb[k];



						//===== Periodic BC ================

						if(model.edge[columnEdgeNumb].aPBC ||model.edge[rowEdgeNumb].aPBC){

							cPB=model.edge[columnEdgeNumb].getnPBC()*model.edge[rowEdgeNumb].getnPBC();

							H1[j][k]*=cPB;
							if(nonLinear){
								model.H2[j][k]*=cPB;
							}
							if(eddy)
								model.H3[j][k]*=cPB;

						}

						//===========================

						//===========  right-hand side 1

						double Ak=0;
						if(!model.edge[columnEdgeNumb].hasPBC()) Ak= model.edge[columnEdgeNumb].A;
						else {

							Ak= model.edge[model.edge[columnEdgeNumb].map].A;
						}


						if(model.edge[columnEdgeNumb].edgeKnown || model.nonLin ){

							model.HkAk.el[matrixRow]+=H1[j][k]*Ak;


								if(model.edge[columnEdgeNumb].edgeKnown )
								continue;
						}



						//=======================

						if(nonLinear){

							H1[j][k]+= model.H2[j][k];
						}


						columnIndex=model.edgeUnknownIndex[columnEdgeNumb]-1;
						if(columnIndex>matrixRow) continue;
						m=util.search(model.Hs.row[matrixRow].index,nz[matrixRow]-1,columnIndex);
						if(m<0)
						{

							if(abs(H1[j][k])>eps  ){

								model.Hs.row[matrixRow].index[nz[matrixRow]]=columnIndex;

								model.Hs.row[matrixRow].el[nz[matrixRow]]=H1[j][k];


								if(fillSs && model.H3!=null && abs(model.H3[j][k])>eps ){

									model.Ss.row[matrixRow].index[nz[matrixRow]]=columnIndex;

									model.Ss.row[matrixRow].el[nz[matrixRow]]=model.H3[j][k];
									nzs[matrixRow]++;

									if(nzs[matrixRow]==model.Ss.row[matrixRow].nzLength-1){
										model.Ss.row[matrixRow].extend(ext);
									}

								}

								nz[matrixRow]++;

								//===========================
								if(nz[matrixRow]==model.Hs.row[matrixRow].nzLength-1){
									model.Hs.row[matrixRow].extend(ext);
								}
								//===========================
							}
						}
						else{

							model.Hs.row[matrixRow].addToNz(H1[j][k],m);
							if(fillSs)
								model.Ss.row[matrixRow].addToNz(model.H3[j][k],m);
						}


					}
				}
			}
		}


		for(int i=0;i<model.numberOfUnknownEdges;i++){
			model.Hs.row[i].sortAndTrim(nz[i]);
		}


		if(fillSs){

		for(int i=0;i<model.numberOfUnknownEdges;i++){
			model.Ss.row[i].sortAndTrim(nzs[i]);

			//model.Ss.showcl();
		}


		if(model.eddyTimeIntegMode==0 || model.eddyTimeIntegMode==1){

		model.Ss.times(1.0/model.dt);

		}


		}


		if( model.analysisMode==2 && this.dim==3){
			SpMat T;
			if(step==0){
				T=new SpMat(model.numberOfVarNodes, model.numberOfUnknownEdges,0);

			}
			else
			{
				T=getPs(model);
				Vect x=T.smul(model.getUnknownA());
				for( int i=0;i<x.length;i++){
					model.b.el[model.numberOfUnknownEdges+i]+=x.el[i];
				}
			}
			for(int i=0;i<model.numberOfVarNodes;i++)
				model.Hs.row[i+model.numberOfUnknownEdges]=T.row[i].deepCopy();

			T=getQs(model);

			for(int i=0;i<model.numberOfVarNodes;i++){
				model.Hs.row[i+model.numberOfUnknownEdges]=model.Hs.row[i+model.numberOfUnknownEdges].augh(T.row[i]);
			}

		}

*/


};

#endif /* MAGMATRIX_H_ */
