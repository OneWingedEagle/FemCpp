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
public class AcoustMatrix {

	public AcoustMatrix(){}

	public AcoustMatrix(Model model){


	}

	public void setSeepMat(Model model){
		setAcoustMat(model);
	}

	public void setAcoustMat(Model model){

		System.out.println(" Accoustic analysis... ");

		System.out.println();
		System.out.println();

		System.out.println(" Calculating stiffness matrix ...");
		Mat He=new Mat(model.nElVert,model.nElVert);
		Mat Se=new Mat(model.nElVert,model.nElVert);
		int m,row,column,rowNodeNumber,colNodeNumber,ext=4;
		int[] nz=new int[model.numberOfUnknownT];
		model.Hs=new SpMat(model.numberOfUnknownT);
		model.Ss=new SpMat(model.numberOfUnknownT);
		
		for(int i=0;i<model.numberOfUnknownT;i++){
			model.Hs.row[i]=new SpVect(model.numberOfUnknownT,model.nNodNod);
	
			model.Ss.row[i]=new SpVect(model.numberOfUnknownT,model.nNodNod);
		}

		
		model.RHS=new Vect(model.numberOfUnknownT);

		for(int i=1;i<=model.numberOfElements;i++){

			if(!model.element[i].isConductor()) {continue;}

			int[] vertNumb=model.element[i].getVertNumb();

				Se=model.femCalc.NiNjQ(model,i);
				He=model.femCalc.gradNi_gradNjQ(model,i);


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
						model.Hs.row[row].el[nz[row]]=He.el[j][k];
							model.Ss.row[row].index[nz[row]]=column;	
							model.Ss.row[row].el[nz[row]]=Se.el[j][k];
						

						nz[row]++;

						if(nz[row]==model.Hs.row[row].nzLength-1){
							model.Hs.row[row].extend(ext);
								model.Hs.row[row].extend(ext);
						}

					}

					else{

						model.Hs.row[row].addToNz(He.el[j][k],m);
							model.Ms.row[row].addToNz(Se.el[j][k],m);


					}
				}

			}
		
		}

	model.Hs.sortAndTrim(nz);

	model.Ss.sortAndTrim(nz);	
		

}




}


