package fem;


import math.Mat;
import math.SpMat;
import math.SpMatSolver;
import math.SpVect;
import math.Vect;
import math.util;

/**
 * TODO Put here a description of what this class does.
 *
 * @author Hassan.
 *         Created Aug 20, 2012.
 */
public class SeepMatrix {

	public SeepMatrix(){}

	public SeepMatrix(Model model){


		model.femCalc=new Calculator(model);

	}

	public void setSeepMat(Model model){
		setSeepMat(model,false);
	}

	public void setSeepMat(Model model,boolean massNeeded){

		System.out.println(" Structural analysis... ");

		System.out.println();
		System.out.println(" Number of unknown non-fixed points : "+model.numberOfUnknownT);
		System.out.println();

		System.out.println(" Calculating stiffness matrix ...");
		Mat Se=new Mat(model.nElVert,model.nElVert);
		Mat Me=new Mat(model.nElVert,model.nElVert);
		int m,row,column,rowNodeNumber,colNodeNumber,ext=4;
		int[] nz=new int[model.numberOfUnknownT];
		model.Hs=new SpMat(model.numberOfUnknownT);
		model.Ms=new SpMat(model.numberOfUnknownT);
		
		for(int i=0;i<model.numberOfUnknownT;i++){
			model.Hs.row[i]=new SpVect(model.numberOfUnknownT,model.nNodNod);
		
		if(massNeeded)
			model.Ms.row[i]=new SpVect(model.numberOfUnknownT,model.nNodNod);
		}

		
		model.RHS=new Vect(model.numberOfUnknownT);

		for(int i=1;i<=model.numberOfElements;i++){

			if(!model.element[i].isConductor()) {continue;}

			int[] vertNumb=model.element[i].getVertNumb();

			if(massNeeded){

				Se=model.femCalc.Se(model,i,Me);


			}
			else{
				Se=model.femCalc.Se(model,i);

			}

			for(int j=0;j<model.nElVert;j++){

				rowNodeNumber=vertNumb[j];

				if(model.node[rowNodeNumber].is_T_known()) continue;

				row=model.T_unknownIndex[rowNodeNumber]-1;

				for(int k=0;k<model.nElVert;k++){

					colNodeNumber=vertNumb[k];


					//=========================

					if(model.node[colNodeNumber].is_T_known() ){
						double T;
							T=model.node[colNodeNumber].T;									
			
				
							model.RHS.el[row]-=Se.el[j][k]*T;
							
										continue;
					}

				

					column=model.T_unknownIndex[colNodeNumber]-1;

					if(column>row) continue;

					m=util.search(model.Hs.row[row].index,nz[row]-1,column);

					if(m<0)
					{	
						model.Hs.row[row].index[nz[row]]=column;															
						model.Hs.row[row].el[nz[row]]=Se.el[j][k];
						if(massNeeded){
							model.Ms.row[row].index[nz[row]]=column;	
							model.Ms.row[row].el[nz[row]]=Me.el[j][k];
						}

						nz[row]++;

						if(nz[row]==model.Hs.row[row].nzLength-1){
							model.Hs.row[row].extend(ext);
							if(massNeeded)
								model.Hs.row[row].extend(ext);
						}

					}

					else{

						model.Hs.row[row].addToNz(Se.el[j][k],m);
						if(massNeeded)
							model.Ms.row[row].addToNz(Me.el[j][k],m);


					}
				}

			}
		
		}
	
		


	model.Hs.sortAndTrim(nz);



	if(massNeeded){

		model.Ms.sortAndTrim(nz);	
		}

}



	public Vect getbUt(Model model,int mode){

		Vect bUt=new Vect(model.bU.length);
		int[][] index=new int[1+model.numberOfNodes][model.dim];
		Mat R=new Mat();
		if(model.dim==2)
			R=util.rotMat2D(model.cpm);		
		else{
			Mat R2D=util.rotMat2D(model.cpm);
			R=new Mat(model.dim,model.dim);
			for(int m=0;m<2;m++)
				for(int n=0;n<2;n++)
					R.el[m][n]=R2D.el[m][n];

			R.el[2][2]=1;
		}

		int ix=0,nn;
		Vect F;

		if(mode>0){
			for(int i=1;i<=model.numberOfNodes;i++){

				if(!model.node[i].isDeformable() || 
						model.node[i].is_U_known() || 			
						model.node[i].getMap()>0) continue;

				for(int p=0;p<model.dim;p++)
					if(!model.node[i].is_U_known(p)){

						index[i][p]=ix++;
					}
			}


			for(int i=1;i<=model.numberOfNodes;i++){

				if(!model.node[i].isDeformable() || model.node[i].is_U_known()) continue;
				nn=i;
				F=model.node[i].getNodalVect(mode);

				if(F==null) continue;

				if(model.node[i].hasPBC())
				{
					nn=model.node[i].getMap();
					F=R.mul(F);
				}

				if(!model.node[i].hasPBC()){
					for(int p=0;p<model.dim;p++)
						if(!model.node[i].is_U_known(p)){
							bUt.el[index[nn][p]]+=F.el[p];
						}

				}
			}

		}

		return bUt;
	}


	public Vect getHead( Model model,SpMatSolver solver,int mode){

		solver.terminate(false);
		System.out.println(" Calculating deformation....");


		Vect bU1=model.bU.add(this.getbUt(model, mode));


		if(model.Ci==null)
			model.Ci=model.Ks.scale(bU1);
		else
			bU1.timesVoid(model.Ci);
		Vect u;
		boolean firstu=false;
		if(model.xp==null){
			model.xp=new Vect(bU1.length);
			firstu=true;
		}

		if(model.Ls==null)
			model.Ls=model.Ks.ichol();

		if(firstu)
			u=solver.ICCG(model.Ks,model.Ls, bU1,1e-5,5000);
		else{
			//	u=solver.ICCG(model.Ks,model.Ls, bU1,2e-3,5000,model.xp);
			u=model.solver.err0ICCG(model.Ks,model.Ls, bU1,1e-6,5000,model.xp);	

		}

		model.xp=u.deepCopy();

		u.timesVoid(model.Ci);

		return u;
	}



}


