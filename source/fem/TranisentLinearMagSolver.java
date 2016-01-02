package fem;

import math.Complex;
import math.SpMat;
import math.SpMatComp;
import math.SpVect;
import math.Vect;
import math.VectComp;
import math.util;

public class TranisentLinearMagSolver {

int stepNumb;
	

	public TranisentLinearMagSolver(){	}

	public Vect solve(Model model, int step ){

		this.stepNumb=step;

		
		SpMat L=new SpMat();

		Vect x=new Vect(model.numberOfUnknowns);

		model.solver.terminate(false);

		model.magMat.setRHS(model);

	if(step==0)
		model.setMagMat();
				

		//=== known values go to right hand side 
	
	
		model.RHS=model.RHS.sub(model.HkAk);
		



		if(model.eddyTimeIntegMode==0) {
				

				
				
				if(step==0)
				model.Hs.addSmaller(model.Ss);

				Vect v1=model.getUnknownA();
				
				if(model.analysisMode==2){
					Vect v=model.Ts.amul(v1);
					
					if(model.RHS!=null)
					for( int i=0;i<v.length;i++){
						model.RHS.el[model.numberOfUnknownEdges+i]+=v.el[i];
					}
					
					v1=v1.aug(new Vect(model.numberOfUnknowns-v1.length));
				}
				

		
	
				model.RHS=model.RHS.add(model.Ss.smul(v1));



	
			SpMat Ks=model.Hs.deepCopy();

			Vect Ci=Ks.scale(model.RHS);


			L=Ks.ichol();

			if(model.RHS.abs().max()>1e-8){

				if(model.xp==null){
					x=model.solver.ICCG(Ks,L, model.RHS,model.errCGmax,model.iterMax);
				//	x=model.solver.err0ICCG(model.Hs,L, model.RHS,1e-3*model.errCGmax,model.iterMax);	

				}
				else{
					x=model.solver.err0ICCG(Ks,L, model.RHS,1e-2*model.errCGmax,model.iterMax,model.xp);	
					//x=model.solver.ICCG(Ks,L, model.RHS,model.errCGmax,model.iterMax,model.xp);	

				}
			}

			else
				x=new Vect(x.length);

			model.xp=x.deepCopy();
	

			x.timesVoid(Ci);

			model.setSolution(x);	

			model.setB();	

			System.out.println("Bmax ( linear analysis): "+model.Bmax);

			return x;


		}



		//=========

		else if( model.eddyTimeIntegMode==1){

//crank linear

			SpMat  Ks=model.Hs.deepCopy();

			Ks.addSmaller(model.Ss);

			Vect bb1=model.RHS.deepCopy();

			model.Ci=Ks.scale(bb1);

			L=Ks.ichol();

			if(this.stepNumb<1){


				if(bb1.norm()>1e-8){
					x=model.solver.ICCG(Ks,L, bb1,model.errCGmax,model.iterMax);
				}
				else
					x=new Vect(bb1.length);

				x.timesVoid(model.Ci);	
				model.up=x.deepCopy();


				model.setSolution(x);	

				model.setB();	

				return x;

			}


			Ks=model.Hs.timesNew(.5).addNew(model.Ss);

			Vect bp=model.Ss.smul(model.up).add(model.Hs.smul(model.up.times(-.5)));

			Vect bU1=new Vect();
			if(model.bT==null)
				bU1=model.RHS.add(bp);
			else{
				bU1=model.RHS.add(model.bT).times(.5).add(bp);

			}


			model.Ci=Ks.scale(bU1);



			L=Ks.ichol();

			if(model.xp==null){
				x=model.solver.ICCG(Ks,L, bU1,model.errCGmax,model.iterMax);
			}
			else{
				x=model.solver.err0ICCG(Ks,L,bU1,1e-1*model.errCGmax,model.iterMax,model.xp);	

			}



			model.xp=x.deepCopy();

			x.timesVoid(model.Ci);

			model.up=x.deepCopy();

			model.setSolution(x);	

			model.setB();	

			return x;


		}



		return null;

	}



}
