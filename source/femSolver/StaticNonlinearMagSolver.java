package femSolver;

import java.text.DecimalFormat;

import fem.Model;
import math.SpMat;
import math.Vect;


public class StaticNonlinearMagSolver{
	int stepNumb;
	boolean usePrev=false;
	int totalNonlinIter;
	boolean relaxedNR=true;

	public StaticNonlinearMagSolver(){	}


	
	
	public Vect solve(Model model,Vect x,boolean echo,int step){


		DecimalFormat dfB=new DecimalFormat("0.0000");
		DecimalFormat dfe=new DecimalFormat("0.0E00");


		int iterMax=model.iterMax;
		int nonLinIterMax=model.nonLinIterMax;

		Vect Ci;
		SpMat L;
		Vect dA=new Vect(model.numberOfUnknowns);

		Vect[] B1,B2;

		model.solver.terminate(false);

		double errNR=1,resNR=0,resNR0=1;
		double errFlux=1;

		int nonLinIter=0;
		


		for( nonLinIter=0; errNR>model.errNRmax && nonLinIter<nonLinIterMax;nonLinIter++)

			
		{
			totalNonlinIter++;

			if(echo)
			{
				System.out.println();
				System.out.println();
				System.out.println("    nonlinear Iteration: "+nonLinIter+
						"     Bmax: "+dfB.format(model.Bmax)+ "     error: "+dfe.format(errNR));
				System.out.println(" Flux    error: "+dfe.format(errFlux));
				System.out.println();
			}

			SpMat Ks=new SpMat();

			Vect b=new Vect();


			model.setMagMat();
			


				Ks=model.Hs.deepCopy();

				b=model.RHS.sub(model.HkAk);
		

			if(nonLinIter==0){
				resNR0=b.norm();
				model.solver.resRef=resNR0;
			}
			
			resNR=b.norm();
			
		
			errNR=resNR/resNR0;



			Ci=Ks.scale(b);
			L=Ks.ichol();

				dA=model.solver.ICCG(Ks,L, b,model.errCGmax,iterMax);


			if(model.solver.terminate) break;


			dA.timesVoid(Ci);
		
			x=x.add(dA);
			

			
			
			
			B1=model.getAllB();

			model.setSolution(x);	

			model.setB();	
			
			B2=model.getAllB();
			
			errFlux=model.getDiffMax(B1,B2);

			
	

		}
		

	
		
		System.out.println();
		System.out.println("-----------------------------------------------------");
		System.out.println(" Nonlinear Iteration: "+nonLinIter+
				"     Bmax: "+dfB.format(model.Bmax)+ "     error: "+dfe.format(errNR));
		System.out.println("Flux    error: "+dfe.format(errFlux));
		System.out.println("-----------------------------------------------------");
		System.out.println();
		
		System.out.println();
		System.out.println("=======================================================");
		System.out.println("Total number of ICCG iterations: "+model.solver.totalIter);
		System.out.println("Total number of Nonlinear iterations: "+totalNonlinIter);
		
		
		System.out.println();
	
		System.out.println(" Final NR residual: "+resNR);
		
		System.out.println("=======================================================");
		System.out.println();
		


		return x;

	}





}
