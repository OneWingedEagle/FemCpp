package fem;

import static java.lang.Math.PI;

import java.awt.Color;
import java.util.Arrays;

import math.SpMat;
import math.Vect;
import math.util;

public class Geo {
	
	public Vect xp;
	public int[] nn=new int[1];

	
	
	public Geo(){}
	
	public static void main(String[] args){
		new Main();
	}
	
	
	public void runGeo(Model model, Main main){



		String Trq = System.getProperty("user.dir") + "//torque.txt";
		model.resultFolder=System.getProperty("user.dir") + "//geoResults";
		String folder=model.resultFolder;
		String tmpFile =model.resultFolder+"//appen.txt";
		main.prepare();

		int[] nn1=new int[model.numberOfNodes];
		int ix=0;
		for(int i=1;i<=model.numberOfNodes;i++) 
			if(util.getAng(model.node[i].getCoord())<PI/2)
			nn1[ix++]=i;
		
		 nn=Arrays.copyOf(nn1,ix);

		 model.setSeepBC();

			//model.rdt=2;
			
			//==========================


			int nT=50;
	
			
			double[] hh=new double[nT];


			if(nT>1000)
				model.writeFiles=false;
			else
				model.writeFiles=true;

			boolean writeFiles=model.writeFiles;

			main.gui.lbX[0].setText("rot ang. ");
			main.gui.lbX[1].setText("Trq. ");




			model.solver.terminate(false);


		//	if(model.Hs==null)
			 model.setSeepMat(true);
			 
			 double rdt=1.0/model.dt;

			SpMat Hs=model.Hs.addNew(model.Ms.timesNew(rdt));


			for(int i=0;i<nT;i++) {

				System.out.println(" Calculating seepage....time step "+i);

				main.gui.tfX[1].setText(Integer.toString(i)+"/"+(nT));

				Vect bU1;
				if(model.hp==null)
				 bU1=model.b.deepCopy();
				else 
				 bU1=model.b.add(model.Ms.smul(model.hp.times(rdt)));
			
				
			if(model.Ci==null)
					model.Ci=Hs.scale(bU1);
			else
					bU1.timesVoid(model.Ci);
				
				Vect h;

			
				if(model.Ls==null)
					model.Ls=Hs.ichol();

				if(model.xp==null)
					h=model.solver.ICCG(Hs,model.Ls, bU1,1e-5,5000);
				else{
					h=model.solver.err0ICCG(Hs,model.Ls, bU1,1e-6,5000,model.xp);	

				}

				model.xp=h.deepCopy();

				h.timesVoid(model.Ci);

				
				model.hp=h.deepCopy();


				model.setHead(h);
				//model.setVelocity();
				
				//hh[i]=model.node[16].T;
				hh[i]=model.node[model.element[16].getVertNumb(0)].T;
				
				if(writeFiles){

				
			

					String ff=folder+"\\bun"+i+".txt";
					if(i==0)
					model.writeMesh(ff);
					


					ff=folder+"\\Tmp"+i+".txt";
					
				/*	if(i%10==0)*/{
						
						model.appendScalarField(tmpFile,0,nn,i);
				//	model.writeNodalScalar(ff);

				/*	ff=folder+"\\veloc"+i+".txt";
					model.writeB(ff);*/
					}

					}
			}
			
			util.plot(hh);

		String logFilePath = model.resultFolder+ "\\log.txt";
		main.gui.writeLog(logFilePath);
		main.gui.Run.setBackground(Color.green);

	}

}


