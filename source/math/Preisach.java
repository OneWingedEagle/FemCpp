package math;

import java.util.Random;

import math.Mat;
import math.Vect;
import math.util;

public class Preisach {

	public Random r;
	public int M,nphi,kk;
	public long seed;
	public double mean,width,Hs,Bs,DB;
	public boolean[] on;
	public double[] phi;
	public double[][][] a;
	public double[][] K;
	public double  eps=1e-3;


	public Preisach(){}


	public Preisach(int M, double mean,double width,double Hmax,double Bs, long seed){

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

		Preisach ps=new Preisach(Mp,mean,width,Hs,Bs,3564656);

		util.pr(ps.getRes());
	
		int ix=0;
		
		Mat[] BH=new Mat[10];
		
		ps.demagnetize(Hs,20);
		util.pr(ps.getRes());
		//Mat BH=ps.magnetizeUpTo(Bs,500);
		//util.pr(ps.getRes());
		/*Mat BH=ps.magnetizeDownTo(-.75,200);
		util.pr(ps.getRes());
		ps.demagnetize(-Hs,20);
		util.pr(ps.getRes());
		//Mat BH=ps.getCurve(new Vect().linspace(0, ps.Hs,100));
	//	BH.*/
	//	util.plotBunch(BH.el);
		//BH=ps.demagnetize(Hs);
		//util.plotBunch(BH.el);	
		//BH=ps.magnetize();
	//	util.pr(ps.getRes());
		//BH[ix++]=ps.symMajorFull(100);
	/*	BH[ix++]=ps.symDesc(.7,500);
		 BH[ix++]=ps.symAsc(.7,500);
		 BH[ix++]=ps.symFull(.7,500);*/
		 BH[ix++]=ps.revAsc(.7,500);
	//	BH=ps.symMajorAsc(100);
		/*util.pr(ps.getRes());
		util.plotBunch(BH.el);*/
		/*util.pr(ps.getRes());
*/
		util.plotBunch(BH,ix);


	}


	public Mat demagnetize(){
		return demagnetize(Hs,20);
	}
	public Mat demagnetize(double Hm){
		return demagnetize(Hm,20);
	}

	public Mat demagnetize(int nCycles){
		return demagnetize(Hs,nCycles);
	}

	
	
	public Mat demagnetize(double Hm,int nCycles){

		int L=200;
		double t=0;
		double dt=1./L;

		//double cc=2.0/nCycles;

		Vect H=new Vect(nCycles*L);

		int  ix=0;
		for(int i=0;i<nCycles;i++)
			for(int j=0;j<L;j++){

					double x=(1-(t+dt)/nCycles)*Math.sin(2*Math.PI*t);
					//double x=Math.exp(-cc*t)*Math.sin(2*Math.PI*t);
				t+=dt;
				H.el[ix]=x*Hm;
				ix++;
			}


		Mat BH=this.getCurve(H);

		return BH;

	}

	
public Mat magnetizeUpTo(double Bpeak,int L){
		

		Vect H=new Vect().linspace(0, Hs, L);
		Vect B=new Vect(L);

		int ix=0;
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
			
			ix++;
			if(B.el[i]>Bpeak || Math.abs(B.el[i]-Bpeak)<=eps) {
				B.el[i]=Bpeak;
				break;
			}
			
			
		}


		Mat BH1=new Mat(ix+1,2);
		for(int i=0;i<=ix;i++){
			BH1.el[i][0]=H.el[i];
			BH1.el[i][1]=B.el[i];
		}
		 Mat BH=this.distill(BH1);

		return BH;

	
	}
	
public Mat magnetizeDownTo(double Bpeak,int L){
		

		Vect H=new Vect().linspace(0, -Hs, L);
		Vect B=new Vect(L);

		int ix=0;
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
			
			ix++;
			if(B.el[i]<Bpeak || Math.abs(B.el[i]-Bpeak)<=eps) {
				B.el[i]=Bpeak;
				break;
			}
			
			
		}


		Mat BH1=new Mat(ix+1,2);
		for(int i=0;i<=ix;i++){
			BH1.el[i][0]=H.el[i];
			BH1.el[i][1]=B.el[i];
		}
		 Mat BH=this.distill(BH1);

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

	public  Mat initial(int L){
		
		this.demagnetize();
		return 	this.magnetizeUpTo(Bs,L);
	}
	
	public  Mat initLoop(double Bpeak, int L){
		
		return 	this.magnetizeDownTo(Bpeak,L);
	}
	
	public  Mat initLoopUptoBpeak(double Bpeak,int L){
		this.demagnetize();

		Vect seqH=new Vect().linspace(0, this.Hs, L);

		Mat BH1=this.getCurve(seqH);

		Mat BH=this.distill(BH1);

		BH.el[BH.nRow-1][1]=this.Bs;


		return BH;
	}
	

	public  Mat symMajorDesc(int L){
		return symDesc(Bs,L);
	}


	public  Mat symMajorAsc(int L){
		
		return symAsc(-Bs,L);
	}
	
	public  Mat symMajorFull(int L){
		
		return symFull(Bs,L);
	}
	
	public  Mat symDesc(double Bpeak,int L){

		this.demagnetize();

		Mat BH1=this.magnetizeUpTo(Bpeak,L);

		int L1=BH1.nRow;
		double H1=BH1.el[L1-1][0];
		
		Vect seqH=new Vect().linspace(H1, -Hs, L);
	

		Mat BH2=this.getCurve(seqH);
		
		Mat BH3=this.distill(BH2);
		
		int ix=0;
	
		for(int i=0;i<BH3.nRow;i++){
			if(BH3.el[i][1]>=-Bpeak)
				ix++;
		}

		
		Mat BH4=new Mat(ix,2);
		
		for(int i=0;i<ix;i++){
			BH4.el[i]=BH3.el[i];
	
		}
		

		return BH4;
	}


	public  Mat symAsc(double Bpeak,int L){

		this.demagnetize();

		Mat BH1=this.magnetizeDownTo(-Bpeak,L);
		
		int L1=BH1.nRow;
		double H1=BH1.el[L1-1][0];
		
		Vect seqH=new Vect().linspace(H1, Hs, L);
	

		Mat BH2=this.getCurve(seqH);
		
		Mat BH3=this.distill(BH2);
		
		int ix=0;
	
		for(int i=0;i<BH3.nRow;i++){
			if(BH3.el[i][1]<=Bpeak)
				ix++;
		}

		
		Mat BH4=new Mat(ix,2);
		
		for(int i=0;i<ix;i++){
			BH4.el[i]=BH3.el[i];
	
		}
		

		return BH4;

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
	
	public  Mat revDesc(double Bpeak,int L){
		
	this.demagnetize();
	this.magnetizeDownTo(-Bs,L);
	
	Vect H=new Vect().linspace(-Hs, Hs, L);
	
	Vect B=new Vect(L);

	int ix=0;
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
		
		ix++;
		if(B.el[i]>Bpeak || Math.abs(B.el[i]-Bpeak)<=eps) {
			B.el[i]=Bpeak;
			break;
		}
		
		
	}

	double H1=H.el[ix-1];

	H=new Vect().linspace(H1, -Hs, L);

	Mat BH1=this.getCurve(H);

	 Mat BH=this.distill(BH1);

	
	return BH;
	
	}
	
	public  Mat revAsc(double Bpeak,int L){
		
		this.demagnetize();
		this.magnetizeUpTo(Bs,L);
		
		Vect H=new Vect().linspace(Hs, -Hs, L);
	
		Vect B=new Vect(L);

		int ix=0;
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
			
			ix++;
			if(B.el[i]<Bpeak || Math.abs(B.el[i]-Bpeak)<=eps) {
				B.el[i]=Bpeak;
				break;
			}
			
			
		}
		


		double H1=H.el[ix-1];

		H=new Vect().linspace(H1, Hs, L);

		Mat BH1=this.getCurve(H);

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
