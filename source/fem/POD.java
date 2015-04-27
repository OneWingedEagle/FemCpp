package fem;


import static java.lang.Math.PI;
import math.Complex;
import math.DFT;
import math.Eigen;
import math.Mat;
import math.MatSolver;
import math.SpMat;
import math.SpMatSolver;

import math.Vect;
import math.util;

/**
 * TODO Put here a description of what this class does.
 *
 * @author Hassan.
 *         Created Aug 20, 2012.
 */
public class POD {


	public POD(){}

	
	
	public void setVibPOD(Model model, Main main){


		util.pr("mech integ mode: "+model.timeIntegMode);

		util.pr("dt:"+model.dt);
		
		//=======================================
		int[] reg8nodes=new int[1];
		model.mapnr=new int[1];
		if(model.dim==3 && model.transfer2DTo3D){
			String m2dBun= System.getProperty("user.dir") + "\\mot4th2DFine.txt";

			model.m2d=new Model(m2dBun);

			model.mapnr=new int[model.m2d.numberOfNodes+1];
			reg8nodes=model.m2d.getRegNodes(8);
			for(int u=0;u<reg8nodes.length;u++)
				model.mapnr[reg8nodes[u]]=u;
			model.forceLamin=new Vect[1][model.m2d.region[8].getNumbElements()];

			model.m2d.coordCode=1;
		}
		//=======================================


		int N=model.nTsteps;

		Vect T=new Vect(N);
		//	model.writeMesh( model.resultFolder+"\\bun.txt");

		model.setMechBC();

	//	String solutionfile=System.getProperty("user.dir") +"\\eigsReact10\\eigVects.txt";
		//String solutionfile=System.getProperty("user.dir") +"\\eigs40Half\\eigVects.txt";
		String solutionfile=System.getProperty("user.dir") +"\\solutions100motj10.txt";
	//	solutionfile=System.getProperty("user.dir") +"\\solutions.txt";

		int ix=0;

		int mode=1;
		int L=model.numberOfUnknownUcomp;

		int D=100;

		util.pr(" Using "+D+" basis.");



		Vect UU=new Vect(model.nTsteps);



		Mat	W=new Mat(model.loader.loadArrays(L, D, solutionfile));

		W.normalizeColumns();

		model.setStiffMat(true);
		
		
		Mat C=W.transp().mul(W).times(1.0/D);
	
		
		 Eigen eg=new Eigen(C);
		
		 
		 Mat Q=eg.V;
	
		 Mat Phi=W.mul(Q);

		 Mat PhiT=Phi.transp();
		 

		 
		Mat Kr=PhiT.mul(model.Ks.smul(Phi));
		Mat Mr=PhiT.mul(model.Ms.smul(Phi));

		double a1=model.rayAlpha;
		double a2=model.rayBeta;


		Mat Cr=Mr.times(a1).add(Kr.times(a2));


	
			double dt=model.dt;
			double beta=.25;
			double gama=.5;
			double b1=1./beta/Math.pow(dt,2);

			double b2=-1./beta/dt;

			double b3=1-.5/beta;
			double b4=gama*dt*b1;
			double b5=1+gama*dt*b2;
			double b6=dt*(1+gama*b3-gama);

			 Mat Kr2=Kr.add(Mr.times(b1).add(Cr.times(b4)));

			Mat Kr1=Kr.deepCopy();
			
			Kr1.lu();
			Kr2.lu();


		Vect u;


		MatSolver ms=new MatSolver();
		
	

		for(int i=	model.nBegin;i<=	model.nEnd;	i+=model.nInc){

						
			String file=model.forceFile[ix];

			main.gui.tfX[1].setText(Integer.toString(i)+"/"+(model.nEnd));
			
			if(model.dim==3 && model.transfer2DTo3D){


						String[] files=new String[model.forceLamin.length];
						for(int k=0;k<files.length;k++){

							int m=i%1800;
							//files[k]=System.getProperty("user.dir") + "\\forcesMotMSz"+k+"\\force"+m+".txt";
							files[k]=System.getProperty("user.dir") + "\\forcesMotMSzerostress\\force"+m+".txt";

							model.m2d.loadNodalField(files[k], 1);
							for(int j=0;j<reg8nodes.length;j++){
								model.forceLamin[k][j]=model.m2d.node[reg8nodes[j]].F;
							}
						}

							model.transfer2DTo3D( file,false,true);
							
							


						}
			else{
			model.loadNodalField(file,1);
			}
			
			
			for(int n=1;n<=0*model.numberOfNodes;n++){
				
				Vect v=model.node[n].getCoord();
				
				if( v.el[0]>-.4999){
					int m=(ix%180);
					model.node[n].F=new Vect(0,0,-1);//.times(Math.cos(2*PI*m/180)+0*Math.cos(0*PI*m/180));
				}
					}
			
			Vect bU1=PhiT.mul(model.bU.add(model.getbUt(mode)));

			
			Vect bp=new Vect();
			Vect ur;

			if(ix<2 ){


				Vect br=bU1.deepCopy();
				
				ur=ms.solvelu(Kr1, br);
	
				
				if(ix==1)
					model.ud=ur.sub(model.up).times(1.0/dt);	
				
				model.up=ur.deepCopy();
				
				model.ud=new Vect(bU1.length);
				model.udd=new Vect(bU1.length);
			}

			else{
			bp=Mr.mul(model.up.times(b1).add(model.ud.times(-b2)).add(model.udd.times(-b3)))
			.add(Cr.mul(model.up.times(b4).add(model.ud.times(-b5)).add(model.udd.times(-b6))));

		
			
			

			Vect br=bU1.add(bp);
			
			ur=ms.solvelu(Kr2, br);

	

				Vect ud1=model.ud.deepCopy();
				Vect udd1=model.udd.deepCopy();
				Vect up1;
				

				if(model.up==null) up1=new Vect(bU1.length);
				else up1=model.up.deepCopy();
			
					model.ud=ur.sub(up1).times(b4).add(ud1.times(b5)).add(udd1.times(b6));	
					model.udd=ur.sub(up1).times(b1).add(ud1.times(b2)).add(udd1.times(b3));
	

				model.up=ur.deepCopy();
			}

				//u=ur.deepCopy();
					
		u=Phi.mul(ur);


				model.setU(u);
				
	

				int nx0=2352;
				int cmp=0;
				nx0=1558; cmp=1; // reactor
				//nx0=15;
				nx0=24697; cmp=0;// motor half
				

				int nx=Math.min(nx0,model.numberOfNodes);

				UU.el[ix]=model.node[nx].getU(cmp);


				ix++;


				if(model.saveForce)
					model.writeNodalField( model.resultFolder+"\\force"+i+".txt",1);

				if(model.saveDisp)
					model.writeNodalField( model.resultFolder+"\\disp"+i+".txt",-1);

				if(model.solver.terminate) break;



			}
		
		
		UU.timesVoid(1e9);
		util.plot(UU);
		
		UU.show();
		
		}
	
	
	
	}


















