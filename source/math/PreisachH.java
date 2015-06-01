package math;

import java.util.Random;

import math.Mat;
import math.Vect;
import math.util;

public class PreisachH {

	public Random r;
	public int M,nphi,kk;
	public long seed;
	public double mean,width,Hs,Bs,DB;
	public boolean[] on;
	public double[] phi;
	public double[][][] a;
	public double[][] K;
	public double  eps=1e-4;


	public PreisachH(){}


	public PreisachH(int M, double mean,double width,double Hmax,double Bs, long seed){

		this.M=M;
		this.mean=mean;
		this.width=width;
		this.seed=seed;
		this.Hs=Hmax;
		this.Bs=Bs;
		this.DB=Bs/M;

		 nphi=18;

		r=new Random(3564656);

		K=new double[M][nphi];
		a=new double[M][2][nphi];
		on=new boolean[M];
		phi=new double[18];
		
		for(int k=0;k<nphi;k++){
			phi[k]=k*10;
			double phirad=phi[k]*Math.PI/180;

		for(int j=0;j<M;j++){

			K[j][k]=1;

			double am=mean*r.nextGaussian()*(1+Math.sin(phirad));
			double d=width*r.nextDouble()*(1+Math.sin(phirad));
			a[j][0][k]=am-d/2;

			a[j][1][k]=am+d/2;


		}
		}


	}

	public static void main2(String[] args)
	{

		double Bs=1.5;
		double Hs=1600;

		int nInit=1;
		int nMajor=1;
		int nSymLoops=1;
		int nDescending=0;
		int nAscending=0;
		int nTot=nInit+nMajor+nSymLoops+nDescending+nAscending;

		int Mp=1000;
		double mean=200;
		double width=500;

		PreisachH ps=new PreisachH(Mp,mean,width,Hs,Bs,3564656);

		util.pr(ps.getRes());
	//	ps.demagnetize();
		Mat BH=ps.getCurve(new Vect().linspace(0, ps.Hs,100));
		//util.plotBunch(BH.el);
		//BH=ps.demagnetize(Hs);
		//util.plotBunch(BH.el);	
		//BH=ps.magnetize();
	//	util.pr(ps.getRes());
		//BH=ps.symMajorFull(100);
	//	BH=ps.symDesc(600,100);
		//BH=ps.symFull(600,100);
	//	BH=ps.symMajorAsc(100);
		/*util.pr(ps.getRes());
		util.plotBunch(BH.el);*/
		/*util.pr(ps.getRes());

		util.plotBunch(BH.el);*/


	}


	public Mat demagnetize(){
		return demagnetize(-Hs,20);
	}
	public Mat demagnetize(double Hin){
		return demagnetize(Hin,20);
	}

	public Mat demagnetize(int nCycles){
		return demagnetize(-Hs,nCycles);
	}

	public Mat demagnetize(double Hin,int nCycles){

		int L=400;
		double t=0;
		double dt=1./L;


		//double cc=2.0/nCycles;

		Vect B=new Vect(nCycles*(L+L/4));
		Vect H=new Vect(nCycles*(L+L/4));

		int  ix=0;
		for(int i=0;i<nCycles;i++)
			for(int j=0;j<L;j++){

				double x=(1-t/nCycles)*Math.cos(2*Math.PI*t);
				t+=dt;
				H.el[ix]=x*Hin;
				ix++;
			}


		Mat BH=this.getCurve(H);

		return BH;

	}

	public Mat magnetize()
	{
		return magnetize(20,this.Hs);	
	}

	public Mat magnetize(double Hmax)
	{
		return magnetize(20,Hmax);	
	}
	public Mat magnetize(int nCycles,double Hmax){

		int L=400;
		double t=0;
		double dt=1./L;

		Vect H=new Vect(nCycles*L+L/4);


		int ix=0;

		for(int i=0;i<nCycles;i++)
			for(int j=0;j<L;j++){

				double x=t/nCycles*Math.sin(2*Math.PI*t);;

				t+=dt;
				H.el[ix]=x*Hmax;
				ix++;
			}




		for(int j=0;j<L/4;j++){

			double x=t/nCycles*Math.sin(2*Math.PI*t);
			t+=dt;
			H.el[ix]=x*Hmax;
			ix++;
		}


		Mat BH=this.getCurve(H);
		return BH;

	}
	

	public Mat getCurve(Vect H){



		int L=H.length;
		Vect B=new Vect(L);
		
		


		for(int i=0;i<L;i++){
			if(i==0){
				B.el[i]=this.getRes();
				H.el[i]=H.el[0];

				continue;
			}

			double dB=0;
			for(int j=0;j<M;j++)
			{

				if(H.el[i]>H.el[i-1] && H.el[i]>a[j][1][kk]){
					if(!on[j]){
						dB+=this.DB*K[j][kk];
						on[j]=true;
					}
				}
				else if(H.el[i]<H.el[i-1] && H.el[i]<a[j][0][kk] ){

					if(on[j]){
						dB-=this.DB*K[j][kk];
						on[j]=false;
					}
				}

			}

			B.el[i]=B.el[i-1]+2*dB;
		}


		Mat BH=new Mat(L,2);
		for(int i=0;i<L;i++){
			BH.el[i][0]=H.el[i];
			BH.el[i][1]=B.el[i];
		}

		return BH;

	}

	public double getRes(){
		double Br=0;

		for(int j=0;j<M;j++){

			if(on[j]) 
				Br+=this.DB;
			else	
				Br+=-this.DB;


		}

		return Br;
	}

	public  Mat initLoop(int L){
		this.demagnetize();

		Vect seqH=new Vect().linspace(0, this.Hs, L);

		Mat BH1=this.getCurve(seqH);

		Mat BH=this.distill(BH1);

		BH.el[BH.nRow-1][1]=this.Bs;


		return BH;
	}

	public  Mat symMajorDesc(int L){

		this.demagnetize();
		this.magnetize(Hs);

		Vect seqH=new Vect().linspace(this.Hs, -this.Hs, L);

		Mat BH1=this.getCurve(seqH);


		Mat BH=this.distill(BH1);



		return BH;
	}


	public  Mat symMajorAsc(int L){

		this.demagnetize();

		this.magnetize(-Hs);
		Vect seqH=new Vect().linspace(-this.Hs, this.Hs, L);

		Mat BH1=this.getCurve(seqH);


		Mat BH=this.distill(BH1);

		return BH;
	}
	
	public  Mat symMajorFull(int L){
		
		Mat BH1=this.symMajorAsc(L);
		Mat BH2=this.symMajorDesc(L);
		int L1=BH1.nRow;
		int L2=BH2.nRow;
		Mat BH3=new Mat(L1+L2+1,2);
		for(int i=0;i<BH1.nRow;i++){
			BH3.el[i]=BH1.el[i];
		}
		for(int i=0;i<BH2.nRow;i++){
			BH3.el[i+L1]=BH2.el[i];
		}

		BH3.el[L1+L2]=BH1.el[0];
		
		return BH3;
		
	}
	
	public  Mat symDesc(double Hpeak,int L){

		this.demagnetize();

		this.magnetize(Hpeak);

		Vect seqH=new Vect().linspace(Hpeak, -Hpeak, L);
	

		Mat BH1=this.getCurve(seqH);


		Mat BH=this.distill(BH1);



		return BH;
	}


	public  Mat symAsc(double Hpeak,int L){

		this.demagnetize();

		this.magnetize(-Hpeak);

		Vect seqH=new Vect().linspace(-Hpeak, Hpeak, L);
	

		Mat BH1=this.getCurve(seqH);


		Mat BH=this.distill(BH1);



		return BH;
	}
	
	public  Mat symFull(double Bpeak,int L){
		
		Mat BH1=this.symAsc(Bpeak,L);
		Mat BH2=this.symDesc(Bpeak,L);
		int L1=BH1.nRow;
		int L2=BH2.nRow;
		Mat BH3=new Mat(L1+L2+1,2);
		for(int i=0;i<BH1.nRow;i++){
			BH3.el[i]=BH1.el[i];
		}
		for(int i=0;i<BH2.nRow;i++){
			BH3.el[i+L1]=BH2.el[i];
		}

		BH3.el[L1+L2]=BH1.el[0];
		
		return BH3;
		
	}
	
	public  Mat revDesc(double Hpeak,int L){

		this.demagnetize();

		this.magnetize(Hpeak);

		Vect seqH=new Vect().linspace(Hpeak, -Hs, L);
	

		Mat BH1=this.getCurve(seqH);


		Mat BH=this.distill(BH1);



		return BH;
	}
	
	public  Mat revAsc(double Hpeak,int L){
		this.demagnetize();
	
		this.magnetize(Hpeak);

		Vect seqH=new Vect().linspace(Hpeak, Hs, L);
	

		Mat BH1=this.getCurve(seqH);


		Mat BH=this.distill(BH1);



		return BH;
	}

	public Mat distill(Mat BH1){

		int Leff=1;

		int L=BH1.nRow;
		boolean[] skip=new boolean[L];
		for(int i=1;i<BH1.nRow;i++){

			if(i<BH1.nRow/2 && Math.abs(BH1.el[i][1]-BH1.el[i-1][1])<eps){
				skip[i-1]=true;
				continue;
			}
			if( Math.abs(BH1.el[i][1]-BH1.el[i-1][1])<eps){
				skip[i]=true;
				continue;
			}


			Leff++;
		}


		Mat BH=new Mat(Leff,2);

		int ix=0;
		for(int i=0;i<BH1.nRow;i++){

			if(!skip[i]){
				BH.el[ix][0]=BH1.el[i][0];
				BH.el[ix][1]=BH1.el[i][1];
				ix++;
			}


		}




		return BH;
	}

}
