package fem;

import java.text.DecimalFormat;

import math.Complex;
import math.Mat;
import math.MatSolver;
import math.SpMat;
import math.SpMatComp;
import math.SpVect;
import math.SpVectComp;
import math.Vect;
import math.VectComp;
import math.util;

public class FEMsolver {


	double coilInduct=0e-2;
	boolean dg=false;
	///  important  axisym has prob with nonlinear even noncoupled 
	//   crank nonlinear problem

	double theta=.5;
	
	boolean usePrev=false;
	int stpNumb=0;
	
	public FEMsolver(){	}

	public FEMsolver(Model model)	{}



	public Vect solveMagLin(Model model, int step){

		this.stpNumb=step;
		
		if(!usePrev)
			step=0;
		
		SpMat L=new SpMat();

		Vect x=new Vect(model.numberOfUnknowns);

		model.solver.terminate(false);

		model.setMagMat(step);

		//=== known values go to right hand side 

	
		model.b=model.b.sub(model.HkAk);

		if(model.analysisMode==0) {

			Vect Ci=model.Hs.scale(model.b);

			L=model.Hs.ichol(1.0);
			

			if(model.b.abs().max()>1e-8){

				if(!usePrev || model.xp==null){
					x=model.solver.ICCG(model.Hs,L, model.b,model.errMax,model.iterMax);
				}
				else{
					x=model.solver.err0ICCG(model.Hs,L, model.b,1e-1*model.errMax,model.iterMax,model.xp);	

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


			if(model.analysisMode>0){
				model.Hs.addSmaller(model.Ss);
				

				Vect v1=model.getUnknownAp();


				model.b=model.b.add(model.Ss.smul(v1));



			}


			Vect Ci=model.Hs.scale(model.b);


			L=model.Hs.ichol();



			if(model.b.abs().max()>1e-8){

				if(model.xp==null){
					x=model.solver.ICCG(model.Hs,L, model.b,model.errMax,model.iterMax);
				}
				else{
					x=model.solver.err0ICCG(model.Hs,L, model.b,1e-1*model.errMax,model.iterMax,model.xp);	
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

			Vect bb1=model.b.deepCopy();

			model.Ci=Ks.scale(bb1);

			L=Ks.ichol();

			if(this.stpNumb<1){


				if(bb1.norm()>1e-8){
					x=model.solver.ICCG(Ks,L, bb1,model.errMax,model.iterMax);
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
				bU1=model.b.add(bp);
			else{
				bU1=model.b.add(model.bT).times(.5).add(bp);

			}


			model.Ci=Ks.scale(bU1);



			L=Ks.ichol();

			if(model.xp==null){
				x=model.solver.ICCG(Ks,L, bU1,model.errMax,model.iterMax);
			}
			else{
				x=model.solver.err0ICCG(Ks,L,bU1,1e-1*model.errMax,model.iterMax,model.xp);	

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

			Vect bb1=model.b.deepCopy();

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
					x=model.solver.ICCG(Ks,L, bb1,model.errMax,model.iterMax);
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

			Vect bU1=model.b.add(bp);


			model.Ci=Ks.scale(bU1);

			L=Ks.ichol();

			if(model.xp==null){
				x=model.solver.ICCG(Ks,L, bU1,model.errMax,model.iterMax);
			}
			else{
				x=model.solver.err0ICCG(Ks,L,bU1,1e-1*model.errMax,model.iterMax,model.xp);	

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

						
						VectComp  v=new VectComp(model.b);
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
								xc=model.solver.COICCG(Ks,Ls,v,model.errMax,model.iterMax,new VectComp(m),1,true);
								//xc=model.solver.COCG(Ks,Ls,v,model.errMax,model.iterMax,new VectComp(m),1,true);
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
			
			// Hs.row[Hs.nRow-1].show();

		
			int nNeut=model.nNeutral;
			int nUnCur=model.numberOfUnknownCurrents;


			Vect v1=model.getUnknownA();
		

			Vect v2=new Vect(model.numberOfUnknowns);

			for(int j=0;j<v1.length;j++)
				v2.el[j]=v1.el[j];

			
		model.b=model.b.add(model.Ss.smul(v2));


			for(int i=0;i<nUnCur;i++){
				int nr=model.unCurRegNumb[i];

				
				double vp=model.region[nr].terminalVoltagep;
				if(	this.stpNumb==0)
					vp=model.region[nr].terminalVoltage;
		

				if(model.eddyTimeIntegMode==-3){

					double cf=this.theta;

	
					
								if(model.HpAp!=null){
					model.b=model.b.sub(model.HpAp.times(1-cf));
						}
				
			model.b.el[model.Hs.nRow-nUnCur+i-nNeut]=cf*(model.lastRows[i].dot(v2)
				+(((1-cf)*model.region[nr].getWireRes()*model.region[nr].current
				+cf*model.vNeutral
				-(cf*model.region[nr].terminalVoltage+(1-cf)*vp))*model.dt
				-this.coilInduct*model.region[nr].current)/model.height);
		
				
			model.b=model.b.sub(model.lastRows[i].times((1-cf)*model.region[nr].current).vectForm());
			
				}
				else{					
		
				
					model.b.el[model.Hs.nRow-nUnCur+i-nNeut]=model.lastRows[i].dot(v2)
					-model.region[nr].terminalVoltage*model.dt/model.height
					-this.coilInduct*model.region[nr].current/model.height;
					
					
				}

			}

			
			

			Vect Ci=Hs.scale(model.b);

			L=Hs.ichol();


			if(model.b.abs().max()>1e-8){

				if(!usePrev  || model.xp==null){
					x=model.solver.ICCG(Hs,L, model.b,model.errMax,model.iterMax);
				}
				else{
					x=model.solver.err0ICCG(Hs,L, model.b,1e-3*model.errMax,model.iterMax,model.xp);	

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
				b.el[i]=new Complex(model.b.el[i],0);
			}
			
			
			
			Ks.setSymHerm(1);
			
			model.Ci=Ks.scale(b);
			
			SpMatComp Ls=Ks.ichol(1.1);
			Ls.setSymHerm(0);

			VectComp xc=new VectComp();
		
		

				xc=model.solver.COICCG(Ks,Ls,b,model.errMax,model.iterMax,new VectComp(Nun),1,true);
				//xc=model.solver.COCG(Ks,Ls,b,model.errMax,model.iterMax,new VectComp(Nn),1,true);
	
			
			xc.timesVoid(model.Ci);	


			for(int i=0;i<Nun;i++){
				x.el[i]=xc.el[i].re;
			
				}
			
			model.setSolution(x);	

			model.setB();	

			System.out.println("Bmax ( linear analysis): "+model.Bmax);
			
	
			return x;
			
	
		
		}
			
			/*			
						
//time harmonic old
			double  w=2*Math.PI*model.freq;

			SpMat Ks=model.Hs.deepCopy();
	

			Ks.augh(new SpMat(Ks.nRow,Ks.nRow,0));

			SpMat Rs=model.Ss.timesNew(-w).reflect(model.nEdEd);

			Rs.augh(model.Hs.timesNew(-1));
//Rs.plot();
			Ks=Ks.augv(Rs);


			int m=model.b.length;

			Vect v=model.b.aug(new Vect(m));

			boolean lu=false;

			if(lu){

				Mat A=Ks.matForm(true);

				A.lu();
				x=new MatSolver().solvelu(A,v);

			}

			else{

				model.Ci=Ks.scale(v);

				L=Ks.ichol();

				if(v.norm()>1e-8){
					x=model.solver.ICCG(Ks,L,v,model.errMax,model.iterMax);

					//x=model.solver.gaussIter(Ks, v, model.errMax,2*model.iterMax);
				}
				else
					x=new Vect(2*m);

				x.timesVoid(model.Ci);	

			}

			Vect vr=new Vect(m);
			Vect vm=new Vect(m);
			for(int i=0;i<m;i++){
				vr.el[i]=x.el[i];
				vm.el[i]=x.el[i+m];
			}

			//vm.show();

			
			model.setSolution(vr);	

			model.setB();

			return vr;
}
*/

		return null;

	}

	public Vect solveNonLinear(Model model,Vect x,boolean echo,int step){


		DecimalFormat dfB=new DecimalFormat("0.0000");
		DecimalFormat dfe=new DecimalFormat("0.0E00");


		int iterMax=model.iterMax;
		int nonLinIterMax=model.nonLinIterMax;
		double errMax=1e-2;
		boolean fluxErr=true;
		Vect Ci;
		SpMat L;
		Vect dA=new Vect(model.numberOfUnknowns);

		Vect[] B1,B2;

		model.solver.terminate(false);

		double err=1,err0=1;


		int nonLinIter=0;


		for( nonLinIter=0; err>errMax && nonLinIter<nonLinIterMax;nonLinIter++)

			
		{
			

			if(echo)
			{
				System.out.println();
				System.out.println();
				System.out.println("    nonlinear Iteration: "+nonLinIter+"     Bmax: "+dfB.format(model.Bmax)+ "     error: "+dfe.format(err));
				System.out.println();
			}

			SpMat Ks=new SpMat();

			Vect b=new Vect();


			model.setMagMat(nonLinIter);


			Vect dv=new Vect();

			if(model.analysisMode>0)
			{


				if(model.eddyTimeIntegMode==0/* || step==0*/)
				{


					Ks=model.Hs.addNew(model.Ss);
					Vect vp=model.getUnknownAp();

					Vect vi=model.getUnknownA();

					dv=vp.sub(vi);
					
				

					b=model.b.sub(model.HkAk).add(model.Ss.smul(dv));
				}

				else if(model.eddyTimeIntegMode==1){

//crank
					Vect vp=model.getUnknownAp();
					Vect vi=model.getUnknownA();
			

					b=model.b.add(model.bT).sub(model.HkAk.add(model.HpAp).times(0.5)).add(model.Ss.smul(vp.sub(vi)));

					Ks=model.Hs.timesNew(.5).addNew(model.Ss);

				}

				else  if(model.eddyTimeIntegMode<=-2){


					 Ks=getCircuitHs(model,step);


					int nNeut=model.nNeutral;
					int nUnCur=model.numberOfUnknownCurrents;


					Vect vk=model.getUnknownA();

					Vect vk2=new Vect(model.numberOfUnknowns);

					for(int j=0;j<vk.length;j++)
						vk2.el[j]=vk.el[j];

					
					Vect vp=model.getUnknownAp();

					Vect vp2=new Vect(model.numberOfUnknowns);

					for(int j=0;j<vp.length;j++)
						vp2.el[j]=vp.el[j];

					dv=vp2.sub(vk2);

					model.b=model.b.add(model.Ss.smul(dv));


					
					for(int i=0;i<nUnCur;i++){
						int nr=model.unCurRegNumb[i];

						double vprev=model.region[nr].terminalVoltagep;
						if(	this.stpNumb==0)
							vprev=model.region[nr].terminalVoltage;
						
				
						double ip=model.region[nr].currentp;

		
						if(model.eddyTimeIntegMode==-3){
														
							double cf=this.theta;
							
							if(model.HpAp!=null){
								model.b=model.b.sub(model.HpAp.times(1-cf));
									}
			
					
						model.b.el[model.Hs.nRow-nUnCur+i-nNeut]=cf*(model.lastRows[i].dot(vp2)
						+(((1-cf)*model.region[nr].getWireRes()*ip
								+(1-cf)*model.vNeutral
								-(cf*model.region[nr].terminalVoltage+(1-cf)*vprev)))*model.dt			
						-this.coilInduct*ip)/model.height;
											

						model.b=model.b.sub(model.lastRows[i].times((1-cf)*model.region[nr].currentp).vectForm());

						}
						else{
				
							model.b.el[model.Hs.nRow-nUnCur+i-nNeut]=model.lastRows[i].dot(vp2)
							+((model.vNeutral-model.region[nr].terminalVoltage)*model.dt
							-this.coilInduct*ip)/model.height;

						}

					}

					SpMat Q=new SpMat(model.b.length,model.b.length);
					for(int i=0;i<model.b.length;i++){
						if(i<vk.length)
							Q.row[i]=null;
						else
							Q.row[i]=Ks.row[i].deepCopy();
					}

			
					Vect v3=new Vect(model.numberOfUnknowns);
					for(int j=0;j<vk.length;j++){
						v3.el[j]=vk.el[j];
					}

					for(int i=0;i<nUnCur;i++){
						int nr=model.unCurRegNumb[i];
						v3.el[vk.length+i]=model.region[nr].current;

					}
					if(model.nNeutral>0)
					v3.el[vk.length+nUnCur]=model.vNeutral;
					
			
			
					model.b=model.b.sub(Q.smul(v3));


					if(model.eddyTimeIntegMode==-3){
				
						b=model.b.sub(model.HkAk.times((this.theta)));

					}
					else{

						b=model.b.sub(model.HkAk);
						
					}
			

				}
			}


			else
			{
				Ks=model.Hs.deepCopy();

				b=model.b.sub(model.HkAk);
			}




			Ci=Ks.scale(b);
			L=Ks.ichol();

			if(b.abs().max()>1e-6)
				dA=model.solver.err0ICCG(Ks,L,b,1e-3*model.errMax,iterMax);	

			if(model.solver.terminate) break;


			dA.timesVoid(Ci);
			
			

			int nUnCur=model.numberOfUnknownCurrents;
			int nNeut=model.nNeutral;

			// very important 
		/*	for(int k=0;k<nUnCur+nNeut;k++){
			dA.el[dA.length-1-k]=dA.el[dA.length-1-k];
			}*/
			
		//	dA.show();

			x=x.add(dA);
			
			//x=x.times(.5).add(x.add(dA).times(.5));

			model.up=x.deepCopy();


			if(model.eddyTimeIntegMode<=-2){

				
				Vect vk=model.getUnknownA();

				Vect vk2=new Vect(model.numberOfUnknowns);

				for(int j=0;j<vk.length;j++)
					vk2.el[j]=vk.el[j];

				Vect vp=model.getUnknownAp();

				Vect vp2=new Vect(model.numberOfUnknowns);

				for(int j=0;j<vp.length;j++)
					vp2.el[j]=vp.el[j];

				dv=vp2.sub(vk2);

				if(nNeut>0)
					model.vNeutral=x.el[x.length-nNeut];

				for(int k=0;k<nUnCur;k++){
					int nr=model.unCurRegNumb[k];
					model.region[nr].inducedVoltage=model.lastRows[k].dot(dv)/model.dt*model.height;

					model.region[nr].current=x.el[x.length-nUnCur-nNeut+k];

				}
			}





			if(fluxErr){
				B1=model.getAllB();

				model.setSolution(x);	

				model.setB();	

				B2=model.getAllB();
			
				err=model.getDiffMax(B1,B2);//+Math.abs(dA.el[dA.length-1]);
			}
			else{

				if(nonLinIter==1){
					err0=model.b.norm();
					if(err0<errMax) err0=1;
				}
				err=errMax/1e-6*dA.norm()/err0;
				model.setSolution(x);	
				model.setB();	

			}


		}
		
		int nUnCur=model.numberOfUnknownCurrents;


		for(int k=0;k<nUnCur;k++){
			int nr=model.unCurRegNumb[k];

	 model.region[nr].currentp=model.region[nr].current;
		}


		model.HpAp=model.HkAk.deepCopy();

		System.out.println();

		System.out.println("    nonlinear Iteration: "+nonLinIter+"     Bmax: "+dfB.format(model.Bmax)+ "     error: "+dfe.format(err));

		System.out.println();

		return x;

	}

	public void solveCoupled(Model model, Vect x){
		double t1=System.currentTimeMillis();
		boolean twoloops=false;
		model.setStiffMat();
		model.solver.terminate(false);
		Vect Ci;
		SpMat L;
		Vect dA=new Vect(model.numberOfUnknownEdges);

		Vect[] B1 = null,B2=null;

		double err=1,errMax=1e-3;
		int cascadeIter=0;

		for(cascadeIter=0 ;err>errMax && cascadeIter<100;cascadeIter++)

		{

			if(model.solver.terminate) break;
			System.out.println();

			System.out.println(" >>>>>>  coscade Iter: "+cascadeIter+ "  error: "+err+"  Bmax: "+model.Bmax);
			System.out.println();


			if(!twoloops){

				model.setMagMat(cascadeIter);


				Ci=model.Hs.scale(model.b);

				L=model.Hs.ichol();

				double t3=System.currentTimeMillis();

				dA=model.solver.err0ICCG(model.Hs,L,model.b,1e-8,model.iterMax);	

				double t4=System.currentTimeMillis();

				System.out.println();
				System.out.println(" solving time: "+(t4-t3)/1000+ " seconds ");


				if(model.solver.terminate ) break;

				dA.timesVoid(Ci);

				x=x.add(dA);

			}


			if(twoloops){

				B1=model.getAllB();
				x=solveNonLinear(model,x,true,cascadeIter);
				B2=model.getAllB();
				err=model.getDiffMax(B1,B2);
			}
			else{
				B1=model.getAllB();
				model.setSolution(x);	
				model.setB();	
				B2=model.getAllB();
			}

			err=model.getDiffMax(B1,B2);

			model.setReluctForce();

			model.setMSForce();
			model.setDeformation();

		}


		System.out.println("    cascade Iteration: "+cascadeIter+ "  error: "+err);

		double t2=System.currentTimeMillis();
		util.pr("[][][][][][][][][][]>>>> cascade time(sec): "+(t2-t1)/1000);
	}

	public void solveMech(Model model){

		boolean nonlinU=false;
		boolean def=true;

		if(def){


			Vect u=model.getDeformation();

			if(nonlinU)
				u=solveMechNonLinear(model,u,true);

			model.setU(u);
		//	model.setStress();

		}
	}


	public Vect solveMechNonLinear(Model model,Vect x,boolean echo){

		DecimalFormat dfU=new DecimalFormat("0.0000");
		DecimalFormat dfe=new DecimalFormat("0.0E00");

		int[] nn=new int[1000];
		int ix=0;
		for(int n=1;n<=model.numberOfNodes;n++)
			if(model.node[n].getCoord(0)>.49) nn[ix++]=n;

		int iterMax=model.iterMax;
		int nonLinIterMax=model.nonLinIterMax;
		double errMax=1e-3;
		Vect Ci;
		SpMat L;
		Vect du=new Vect(model.numberOfUnknownUcomp);

		nonLinIterMax=90;

		model.solver.terminate(false);

		double err=1,err0=1;


		int nonLinIter=1;

		for( nonLinIter=1; err>errMax && nonLinIter<=nonLinIterMax;nonLinIter++)

		{

			if(echo)
			{
				System.out.println();
				System.out.println();
				System.out.println("    nonlinear Iteration: "+nonLinIter+"     Bmax: "+dfU.format(model.uMax*1e9)+ "     error: "+dfe.format(err));
				System.out.println();
			}

			//	model.writeMesh(System.getProperty("user.dir") +"//results3D//ssforce"+nonLinIter+".txt");
			model.setStress();
			model.setStressForce();

			model.deformMesh();

			//model.writeNodalField(System.getProperty("user.dir") +"//results3D//ssforce"+nonLinIter+".txt",2);
			model.setStiffMat();

			SpMat Ks=model.Ks.deepCopy();

			System.out.println(" Calculating deformation....");

			double a=nonLinIter;
			if(nonLinIter>50) a=50;



			/*for(int i=0;i<ix;i++){
				int n=nn[i];
				model.node[n].F=new Vect(model.node[n].getCoord(2),0,-model.node[n].getCoord(0)).times(40*a);
			}*/


			Vect bU1=model.bU.add(this.getbUt(model, 1).add(this.getbUt(model, 2)));


			Ci=Ks.scale(bU1);
			L=Ks.ichol();

			double t3=System.currentTimeMillis();
			if(bU1.abs().max()>1e-6)
				du=model.solver.ICCG(Ks,L,bU1,1e-8,iterMax);



			double t4=System.currentTimeMillis();
			System.out.println();
			System.out.println(" solving time: "+(t4-t3)/1000+ " seconds ");


			if(model.solver.terminate) break;

			du.timesVoid(Ci);

			err=du.abs().max();

			//	x=x.add(du);


			model.setU(du);	


		}

		System.out.println();

		System.out.println("    nonlinear Iteration: "+nonLinIter+"     Bmax: "+dfU.format(model.uMax)+ "     error: "+dfe.format(err));

		System.out.println();

		return x;

	}

	public void solveSeepage(Model model){



		model.solver.terminate(false);

		System.out.println(" Calculating seepage....");


		Vect bU1=model.b.deepCopy()/*.add(this.getbQt(model))*/;

		SpMat Hs=model.Hs.deepCopy();

		if(model.Ci==null)
			model.Ci=Hs.scale(bU1);
		else
			bU1.timesVoid(model.Ci);
		Vect h;


		boolean firstu=false;
		if(model.xp==null){
			model.xp=new Vect(bU1.length);
			firstu=true;
		}

		if(model.Ls==null)
			model.Ls=Hs.ichol();

		if(firstu)
			h=model.solver.ICCG(Hs,model.Ls, bU1,1e-5,5000);
		else{
			h=model.solver.err0ICCG(Hs,model.Ls, bU1,1e-6,5000,model.xp);	

		}

		model.xp=h.deepCopy();

		h.timesVoid(model.Ci);


		for(int i=1;i<=model.numberOfNodes;i++){
			if(model.T_unknownIndex[i]!=0)
			{
				model.node[i].T=h.el[model.T_unknownIndex[i]-1];

			}
		}

		model.setVelocity();



	}



	public Vect getbUt(Model model,int mode){

		int dim=model.dim;

		Vect bUt=new Vect(model.bU.length);
		int[][] index=new int[1+model.numberOfNodes][model.dim];
		Mat R=new Mat();
		if(model.dim==2)
			R=util.rotMat2D(model.cpm);		
		else{
			Mat R2D=util.rotMat2D(model.cpm);
			R=new Mat(dim,dim);
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

				for(int p=0;p<dim;p++)
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
					for(int p=0;p<dim;p++)
						if(!model.node[i].is_U_known(p)){
							bUt.el[index[nn][p]]+=F.el[p];
						}

				}
			}

		}

		return bUt;
	}

	public Vect getTmp(Model model){

		Vect bU1=model.bT.deepCopy();

		SpMat Ks=model.Ss.deepCopy();

		Vect Ci=Ks.scale(bU1);
		SpMat L=Ks.ichol();

		Vect T=new Vect();

		if(bU1.abs().max()>1e-6)
			T=model.solver.ICCG(Ks,L,bU1,1e-8,2000);
		else 
			T=new Vect(bU1.length);

		T.timesVoid(Ci);


		return T;
	}

	public SpMat getCircuitHs(Model model,int step){


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
	



}
