package fem;


import math.Mat;
import math.SpBlockMat;
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
public class PCMatrix {

	private Calculator calc;

	public PCMatrix(){}

	public PCMatrix(Model model){


		this.calc=new Calculator(model);

	}

	public void setMats(Model model,int kH){

		System.out.println(" Photonic analysis... ");

		System.out.println();
		System.out.println(" Number of unknown non-fixed points : "+model.numberOfUnknownT);
		System.out.println();

		System.out.println(" Calculating stiffness matrix ...");
		Mat Se=new Mat(model.nElVert,model.nElVert);
		Mat Te=new Mat(model.nElVert,model.nElVert);
		
		Mat Pe=new Mat(model.nElVert,model.nElVert);
		
		int m,row,column,rowNodeNumber,colNodeNumber,ext=4;
		int[] nz=new int[model.numberOfUnknownT];
		SpMat Ss=new SpMat(model.numberOfUnknownT);
		SpMat Ms=new SpMat(model.numberOfUnknownT);
		
		for(int i=0;i<model.numberOfUnknownT;i++){
			Ss.row[i]=new SpVect(model.numberOfUnknownT,model.nNodNod);
	
			Ms.row[i]=new SpVect(model.numberOfUnknownT,model.nNodNod);
		}

		

		for(int i=1;i<=model.numberOfElements;i++){

			if(!model.element[i].isThermal()) {continue;}

			int[] vertNumb=model.element[i].getVertNumb();

			//	Se=this.calc.PCe(model,i,Te,Pe,kH);
		

			for(int j=0;j<model.nElVert;j++){

				rowNodeNumber=vertNumb[j];

				if(model.node[rowNodeNumber].is_T_known()) continue;

				row=model.T_unknownIndex[rowNodeNumber]-1;

				for(int k=0;k<model.nElVert;k++){

					colNodeNumber=vertNumb[k];

					column=model.T_unknownIndex[colNodeNumber]-1;


					if(column>row) continue;

					m=util.search(Ss.row[row].index,nz[row]-1,column);

					if(m<0)
					{	
						
						Ss.row[row].index[nz[row]]=column;															
						Ss.row[row].el[nz[row]]=Se.el[j][k];
							Ms.row[row].index[nz[row]]=column;	
							Ms.row[row].el[nz[row]]=Te.el[j][k];
						

						nz[row]++;

						if(nz[row]==Ss.row[row].nzLength-1){
							Ss.row[row].extend(ext);
							Ms.row[row].extend(ext);
						}

					}

					else{

						Ss.row[row].addToNz(Se.el[j][k],m);
						
						Ms.row[row].addToNz(Te.el[j][k],m);


					}
				}

			}
		
		}
	
		


	Ss.sortAndTrim(nz);


			model.Ss=Ss.deepCopy();
			Ms.sortAndTrim(nz);	
			model.Ms=Ms.deepCopy();
		}



	}



