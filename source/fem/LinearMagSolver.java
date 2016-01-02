package fem;

import math.Complex;
import math.SpMat;
import math.SpMatComp;
import math.SpVect;
import math.Vect;
import math.VectComp;
import math.util;

public class LinearMagSolver extends FEMsolver{


	public LinearMagSolver(){	}

	public Vect solve(Model model, int step ){

		this.stpNumb=step;

		
		SpMat L=new SpMat();

		Vect x=new Vect(model.numberOfUnknowns);

		model.solver.terminate(false);

		model.magMat.setRHS(model);

	if(step==0)
		model.setMagMat();
				

		//=== known values go to right hand side 
	
	
		model.RHS=model.RHS.sub(model.HkAk);
		

		if(model.analysisMode==0) {

			Vect Ci=model.Hs.scale(model.RHS);

			L=model.Hs.ichol();

			if(model.RHS.abs().max()>1e-8){

				if(!usePrev || model.xp==null){
					x=model.solver.ICCG(model.Hs,L, model.RHS,model.errCGmax,model.iterMax);
				}
				else{
					x=model.solver.ICCG(model.Hs,L, model.RHS,model.errCGmax,model.iterMax,model.xp);
					//x=model.solver.err0ICCG(model.Hs,L, model.RHS,1e-2*model.errCGmax,model.iterMax,model.xp);	

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



		if(model.eddyTimeIntegMode==0) {


			if( model.analysisMode>0){
				

				
				
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



			}

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

			if(this.stpNumb<1){


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
		else if(model.eddyTimeIntegMode==2){

			double dt=model.dt;

			SpMat  Ks=model.Hs.deepCopy();

			Vect bb1=model.RHS.deepCopy();

			double beta=.25;
			double gama=.5;
			double b1=1./beta/Math.pow(dt,2);

			double b2=-1./beta/dt;

			double b3=1-.5/beta;
			double b4=gama*dt*b1;
			double b5=1+gama*dt*b2;
			double b6=dt*(1+gama*b3-gama);




			if(this.stpNumb<1){

				model.ud=new Vect(Ks.nRow);
				model.udd=new Vect(Ks.nRow);

				Ks=model.Hs.addNew(model.Ss.timesNew(b4));
				model.Ci=Ks.scale(bb1);
				L=Ks.ichol();


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




			Ks=model.Hs.addNew(model.Ss.timesNew(b4));

			Vect bp=model.Ss.smul(model.up.times(b4).add(model.ud.times(-b5)).add(model.udd.times(-b6)));

			Vect bU1=model.RHS.add(bp);


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


			Vect ud1=model.ud.deepCopy();
			Vect udd1=model.udd.deepCopy();
			Vect up1;

			if(model.up==null) up1=new Vect(bU1.length);
			else up1=model.up.deepCopy();

			model.ud=x.sub(up1).times(b4).add(ud1.times(b5)).add(udd1.times(b6));	
			model.udd=x.sub(up1).times(b1).add(ud1.times(b2)).add(udd1.times(b3));


			model.up=x.deepCopy();

			model.setSolution(x);	

			model.setB();	

			return x;

		}

		else if(model.eddyTimeIntegMode==-1){
//time harmonic iccog

						double  w=2*Math.PI*model.freq;

						SpMatComp Ks=new SpMatComp(model.Hs,model.Ss.timesNew(w));

						Ks.setSymHerm(1);

						
						VectComp  v=new VectComp(model.RHS);
						int m=v.length;
							model.Ci=Ks.scale(v);

						//	Ks.show();

							SpMatComp Ls=Ks.ichol(1.0);
							Ls.setSymHerm(0);
						//	Ls.show();
							//  problem with precond
							
						/*	Ls=Ls.timesNew(0);
							Ls.addToDiag(new VectComp(new Vect().ones(m)));*/
						
							VectComp xc;

							if(v.norm()>1e-8){
								xc=model.solver.COICCG(Ks,Ls,v,model.errCGmax,model.iterMax,new VectComp(m),1,true);
								//xc=model.solver.COCG(Ks,Ls,v,model.errCGmax,model.iterMax,new VectComp(m),1,true);
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

		else if(model.eddyTimeIntegMode==-2 || model.eddyTimeIntegMode==3){

			//if(dg) return circuitDG(model, step);
			
			SpMat Hs=getCircuitHs(model,step);
			

		
			int nNeut=model.nNeutral;
			int nUnCur=model.numberOfUnknownCurrents;


			Vect v1=model.getUnknownA();
		

			Vect v2=new Vect(model.numberOfUnknowns);

			for(int j=0;j<v1.length;j++)
				v2.el[j]=v1.el[j];

			
		model.RHS=model.RHS.add(model.Ss.smul(v2));


			for(int i=0;i<nUnCur;i++){
				int nr=model.unCurRegNumb[i];

				
				double vp=model.region[nr].terminalVoltagep;
				if(	this.stpNumb==0)
					vp=model.region[nr].terminalVoltage;
		

				if(model.eddyTimeIntegMode==-3){

					double cf=this.theta;

	
					
								if(model.HpAp!=null){
					model.RHS=model.RHS.sub(model.HpAp.times(1-cf));
						}
				
			model.RHS.el[model.Hs.nRow-nUnCur+i-nNeut]=cf*(model.lastRows[i].dot(v2)
				+(((1-cf)*model.region[nr].getWireRes()*model.region[nr].current
				+cf*model.vNeutral
				-(cf*model.region[nr].terminalVoltage+(1-cf)*vp))*model.dt
				-this.coilInduct*model.region[nr].current)/model.height);
		
				
			model.RHS=model.RHS.sub(model.lastRows[i].times((1-cf)*model.region[nr].current).vectForm());
			
				}
				else{					
		
				
					model.RHS.el[model.Hs.nRow-nUnCur+i-nNeut]=model.lastRows[i].dot(v2)
					-model.region[nr].terminalVoltage*model.dt/model.height
					-this.coilInduct*model.region[nr].current/model.height;
					
					
				}

			}

			
			

			Vect Ci=Hs.scale(model.RHS);

			L=Hs.ichol();


			if(model.RHS.abs().max()>1e-8){

				if(!usePrev  || model.xp==null){
					x=model.solver.ICCG(Hs,L, model.RHS,model.errCGmax,model.iterMax);
				}
				else{
					x=model.solver.err0ICCG(Hs,L, model.RHS,1e-3*model.errCGmax,model.iterMax,model.xp);	

				}
		
			}

			else
				x=new Vect(x.length);

			model.xp=x.deepCopy();


			x.timesVoid(Ci);

			model.up=x.deepCopy();
			
			model.setSolution(x);	

			model.setB();	

			System.out.println("Bmax ( linear analysis): "+model.Bmax);


				v1=model.getUnknownA();

				v2=new Vect(model.numberOfUnknowns);

				for(int j=0;j<v1.length;j++)
					v2.el[j]=v1.el[j];

				Vect v3=model.getUnknownAp();

				Vect v4=new Vect(model.numberOfUnknowns);

				for(int j=0;j<v3.length;j++)
					v4.el[j]=v3.el[j];

				Vect dv=v4.sub(v2);


				for(int k=0;k<nUnCur;k++){
					int nr=model.unCurRegNumb[k];

					model.region[nr].inducedVoltage=model.lastRows[k].dot(dv)/model.dt*model.height;
					

					model.region[nr].currentp=model.region[nr].current;
					model.region[nr].current=x.el[x.length-nUnCur-nNeut+k];

				}



				if(nNeut>0)
					model.vNeutral=x.el[x.length-nNeut];


			return x;


		}
		else if(model.eddyTimeIntegMode==-4){
			
			double w=2*Math.PI*model.freq;
			
		//	Complex jw=new Complex(0,w);
			SpMatComp Ks=new SpMatComp(model.Hs,model.Ss.timesNew(w));
			
			int Nun=Ks.nRow;
			VectComp b=new VectComp(Nun);
			for(int i=0;i<Nun;i++){
				b.el[i]=new Complex(model.RHS.el[i],0);
			}
			
			
			
			Ks.setSymHerm(1);
			
			model.Ci=Ks.scale(b);
			
			SpMatComp Ls=Ks.ichol(1.1);
			Ls.setSymHerm(0);

			VectComp xc=new VectComp();
		
		

				xc=model.solver.COICCG(Ks,Ls,b,model.errCGmax,model.iterMax,new VectComp(Nun),1,true);
				//xc=model.solver.COCG(Ks,Ls,b,model.errCGmax,model.iterMax,new VectComp(Nn),1,true);
	
			
			xc.timesVoid(model.Ci);	


			for(int i=0;i<Nun;i++){
				x.el[i]=xc.el[i].re;
			
				}
			
			model.setSolution(x);	

			model.setB();	

			System.out.println("Bmax ( linear analysis): "+model.Bmax);
			
	
			return x;
			
	
		
		}


		return null;

	}

/*	public SpMat getCircuitHs(Model model,int step){


		SpMat Hs=model.Hs.deepCopy();

		int nUnCur=model.numberOfUnknownCurrents;

		int nNeut=model.nNeutral;
		
	double cf=1;
	
		if(model.eddyTimeIntegMode==-3) cf=this.theta;
		
		
		int nr=Hs.nRow-nUnCur-nNeut;
		

		for(int i=0;i<nr;i++)
		 Hs.row[i]=Hs.row[i].times(cf);
		
		

		for(int i=0;i<nUnCur;i++){


			int jr=Hs.nRow-nUnCur-nNeut+i;

			Hs.row[jr]=model.lastRows[i].times(cf);
			int nz=Hs.row[jr].nzLength;
			Hs.row[jr].extend(1);

	
			double R=model.region[model.unCurRegNumb[i]].getWireRes();

			Hs.row[jr].el[nz]=-cf*(cf*model.dt*R+this.coilInduct)/model.height;

			Hs.row[jr].index[nz]=jr;
		}

		int jr=Hs.nRow-1;

		
		if(nNeut>0){
			Hs.row[jr]=new SpVect(Hs.nRow,1+nUnCur);		

			for(int i=0;i<nUnCur;i++){

				Hs.row[jr].el[i]=-cf*model.dt/model.height;

				Hs.row[jr].index[i]=jr-nUnCur+i;
			}

			Hs.row[jr].el[nUnCur]=cf*model.dt/model.Rg/model.height;

			Hs.row[jr].index[nUnCur]=jr;

		}

		if(model.analysisMode>0){
			Hs.addSmaller(model.Ss);

		}
		
		return Hs;

	}

*/

}
