package femSolver;


import fem.Model;
import math.SpMat;
import math.SpMatComp;
import math.Vect;
import math.VectComp;

public class ACMagSolver {


	public ACMagSolver(){	}

	public Vect solve(Model model, int step ){


		SpMat L=new SpMat();

		Vect x=new Vect(model.numberOfUnknowns);

		model.solver.terminate(false);
		

		model.magMat.setRHS(model);

		

		model.setMagMat();

	

		
		model.freq=50;
	
						double  w=2*Math.PI*model.freq;

						SpMatComp Ks=new SpMatComp(model.Hs,model.Ss.timesNew(w));

						Ks.setSymHerm(1);

						
						VectComp  v=new VectComp(model.RHS);
						int m=v.length;
							model.Ci=Ks.scale(v);

				

							SpMatComp Ls=Ks.ichol(1.05);
							Ls.setSymHerm(0);
		
						
							VectComp xc;

							if(v.norm()>1e-8){
								xc=model.solver.COICCG(Ks,Ls,v,model.errCGmax,model.iterMax,new VectComp(m),1,true);
							}
							else
								xc=new VectComp(m);
							
							xc.timesVoid(model.Ci);	
	

						Vect vr=new Vect(m);
						for(int i=0;i<m;i++){
							vr.el[i]=xc.el[i].re;
						}

						model.setSolution(vr);	

						model.setB();

						return vr;

					
			
		}


}
