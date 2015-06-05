package PlayModel;

import java.util.Random;
import static java.lang.Math.*;
import math.Mat;
import math.Vect;
import math.util;

public class Preisach2D {

	public Random r;
	public int M,nphi,nh,kk,dim;
	public long seed;
	public double cfm,cfw,mean,width,Hs,Bs,DB2D,projCoef,dphiRad;
	public boolean[][] on;
	public double[] phi;
	public double[][][] a;
	public double[][] K;
	public double  epsdB=1e-3;
	//public double  epsdBdH=1e-5;


	public Preisach2D(){}


	public Preisach2D(int M, double mean,double width,double Hmax,double Bs, long seed){

		this.M=M;
		this.mean=mean;
		this.width=width;
		this.seed=seed;
		this.Hs=Hmax;
		this.Bs=Bs;

		dim=2;

		cfm=0;
		cfw=0;

		nh=9;

		nphi=2*nh+1;

		if(nphi==1)  this.dphiRad=0;
		else
			this.dphiRad=Math.PI/(nphi-1);


		double sum=0;
		for(int i=-nh;i<=nh;i++){
			sum+=cos(i*dphiRad);
		}


		sum/=nphi;


		projCoef=sum;


		DB2D=Bs/M/nphi/sum;;


		r=new Random(3564656);

		K=new double[M][nphi];
		a=new double[M][2][nphi];
		on=new boolean[M][nphi];
		phi=new double[nphi];

		double dphiDeg=180.0/(nphi-1);

		for(int k=-nh;k<=nh;k++){
			int kp=k+nh;

			phi[kp]=k*dphiDeg;
			double phirad=k*this.dphiRad;

			for(int j=0;j<M;j++){

				K[j][kp]=1-.0*sin(phirad+0*PI/12);

				double am=mean*r.nextGaussian()*(1+cfm*abs(sin(phirad+0*PI/12)));
				double d=width*abs(r.nextGaussian())*(1+cfw*abs(sin(phirad+0*PI/12)));
				a[j][0][kp]=am-d/2;

				a[j][1][kp]=am+d/2;
		



			}
		}


	}

	public Preisach2D deepCopy(){

		//Random r=new Random();
		long ss=this.seed;
		Preisach2D pr=new Preisach2D(this.M, this.mean,this.width,this.Hs,this.Bs,ss);
		pr.projCoef=this.projCoef;
		pr.dphiRad=this.dphiRad;

		return pr;
	}

	public static void main2(String[] args)
	{

		double Bs=1.8;
		double Hs=1000;

		int nInit=1;
		int nMajor=1;
		int nSymLoops=1;
		int nDescending=0;
		int nAscending=0;
		int nTot=nInit+nMajor+nSymLoops+nDescending+nAscending;

		int Mp=5000;
		double mean=.2*Hs;
		double width=.3*Hs;

		Preisach2D ps=new Preisach2D(Mp,mean,width,Hs,Bs,3564656);

		ps.getRes().hshow();

		ps.demagnetize(10);

		int steps=360;
		int nc=2;
		double t=0;
		double dt=1.0/steps;
		int L=nc*steps;
		Mat H=new Mat(L,ps.dim);
		int jx=0;
		for (int n = 0; n <nc; ++n)
			for (int i = 0; i <steps; ++i){

				double Hm=(1-exp(-4*t/nc))*Hs/4;

				Hm=Hs/2;
				H.el[jx][0] = Hm*sin(2 * PI*t);
				H.el[jx][1] = Hm*cos(2 * PI*t);
				t+=dt;
				jx++;

			}
		
		
	//	Mat BH=ps.magnetizeUpTo(1.8,L);
		
		Mat BH=ps.initial(L);
		
	Vect Hr=new Vect(BH.nRow);
	Vect Br=new Vect(BH.nRow);
	
		for(int i=0;i<BH.nRow;i++){
			Hr.el[i]=new Vect(BH.el[i][0],0).norm();
			Br.el[i]=new Vect(BH.el[i][1],BH.el[i][2]).norm();
		}
		
		util.plot(Hr,Br);
		
/*Mat[] BH2=new Mat[2];
BH2[0]=new Mat(BH.nRow,2);
BH2[1]=new Mat(BH.nRow,2);

		Vect H1=BH.getColVect(0);
		Vect B1=BH.getColVect(3);

		BH2[0].setCol(H1.times(400), 0);
		BH2[0].setCol(B1.times(400), 1);
		
		BH2[1].setCol(BH.getColVect(0), 0);
		BH2[1].setCol(BH.getColVect(1), 1);
		
		
	
		util.plot(BH.getColVect(2),BH.getColVect(3));*/


	}
	
	public  Mat initial(int L){


		return 	initCurveUptoBpeak(this.Bs,L);
	}


	public  Mat initCurveUptoBpeak(double Bpeak,int L){
		
		this.demagnetize();

		Vect seqH=new Vect().linspace(0, this.Hs, L);

		Mat BH1=this.getCurveAlt(seqH);

		Mat BH=this.distill(BH1);


		return BH;
	}
	
/*	public Mat initial(int L){
		
		this.demagnetize();
	Mat BH1=this.getCurveAlt(new Vect().linspace(0,Hs,L));
	
	
	return BH1;
	
	}*/


	public Mat getCurveAlt(Vect H){



		int L=H.length;
		Vect[][] B=new Vect[L][nphi];
		for(int i=0;i<L;i++)
			for(int k=0;k<nphi;k++)
				B[i][k]=new Vect(2);


		for(int kh=-nh;kh<=nh;kh++){


			int k=kh+nh;

			Vect Halt=H.times(cos(kh*dphiRad));


			Vect er=new Vect(cos(kh*dphiRad),sin(kh*dphiRad));

			for(int i=0;i<L;i++){

				if(i==0){

					B[i][k]=B[i][k].add(this.getRes(k));
					//B[i][k]=B[i][k].add(er.times(0));
					continue;
				}

				Vect dB=new Vect(2);
				for(int j=0;j<M;j++)
				{

					if(Halt.el[i]>Halt.el[i-1] && Halt.el[i]>a[j][1][k]){
						if(!on[j][k]){
							dB=dB.add(er.times(this.DB2D*K[j][k]));
							on[j][k]=true;
						}
					}
					else if(Halt.el[i]<=Halt.el[i-1] && Halt.el[i]<a[j][0][k] ){

						if(on[j][k]){
							dB=dB.sub(er.times(this.DB2D*K[j][k]));
							on[j][k]=false;
						}
					}

				}

				B[i][k]=B[i-1][k].add(dB.times(2));
			}
		}

		Vect[] Bsum=new Vect[L];

		for(int i=0;i<L;i++){
			Bsum[i]=new Vect(2);
			for(int j=0;j<nphi;j++)
				Bsum[i]=Bsum[i].add(B[i][j]);

		}


		Mat BH=new Mat(L,3);
		for(int i=0;i<L;i++){
			BH.el[i][0]=H.el[i];

			BH.el[i][1]=Bsum[i].el[0];
			BH.el[i][2]=Bsum[i].el[1];
		}

	//	BH.transp().show();
		return BH;
		
	

	}
	
	
	public Mat getCurveAltxx(Vect H){



		int L=H.length;
		Vect[][] B=new Vect[L][nphi];
		for(int i=0;i<L;i++)
			for(int k=0;k<nphi;k++)
				B[i][k]=new Vect(2);


		for(int kh=-nh;kh<=nh;kh++){


			int k=kh+nh;

			Vect Halt=H.times(cos(kh*dphiRad));


			Vect er=new Vect(cos(kh*dphiRad),sin(kh*dphiRad));

			for(int i=0;i<L;i++){

				if(i==0){

					B[i][k]=B[i][k].add(this.getRes(k));
					//B[i][k]=B[i][k].add(er.times(0));
					continue;
				}

				Vect dB=new Vect(2);
				for(int j=0;j<M;j++)
				{

					if(Halt.el[i]>Halt.el[i-1] && Halt.el[i]>a[j][1][k]){
						if(!on[j][k]){
							dB=dB.add(er.times(this.DB2D*K[j][k]));
							on[j][k]=true;
						}
					}
					else if(Halt.el[i]<=Halt.el[i-1] && Halt.el[i]<a[j][0][k] ){

						if(on[j][k]){
							dB=dB.sub(er.times(this.DB2D*K[j][k]));
							on[j][k]=false;
						}
					}

				}

				B[i][k]=B[i-1][k].add(dB.times(2));
			}
		}

		Vect[] Bsum=new Vect[L];

		for(int i=0;i<L;i++){
			Bsum[i]=new Vect(2);
			for(int j=0;j<nphi;j++)
				Bsum[i]=Bsum[i].add(B[i][j]);

		}


		Mat BH=new Mat(L,3);
		for(int i=0;i<L;i++){
			BH.el[i][0]=H.el[i];

			BH.el[i][1]=Bsum[i].el[0];
			BH.el[i][2]=Bsum[i].el[1];
		}
BH.transp().show();

		return BH;

	}
	

	public Mat getLocus(Mat H){

		int L=H.nRow;
		Vect[][] B=new Vect[L][nphi];
		for(int i=0;i<L;i++)
			for(int k=0;k<nphi;k++)
				B[i][k]=new Vect(2);

		Vect Hr=new Vect(L);


		for(int kh=-nh;kh<=nh;kh++){
			//	for(int kh=0;kh<=0;kh++){


			int k=kh+nh;




			Vect er=new Vect(cos(kh*dphiRad),sin(kh*dphiRad));

			for(int i=0;i<L;i++)
				Hr.el[i]=new Vect(H.el[i]).dot(er);


			for(int i=0;i<L;i++){


				if(i==0){
					//B[i][k]=er.times(this.getRes());
					B[i][k]=B[i][k].add(this.getRes(k));

					continue;
				}

				Vect dB=new Vect(2);
				for(int j=0;j<M;j++)
				{

					if(Hr.el[i]>Hr.el[i-1] && Hr.el[i]>a[j][1][k]){
						if(!on[j][k]){
							dB=dB.add(er.times(this.DB2D*K[j][k]));
							on[j][k]=true;
						}
					}
					else if(Hr.el[i]<=Hr.el[i-1] && Hr.el[i]<a[j][0][k] ){

						if(on[j][k]){
							dB=dB.sub(er.times(this.DB2D*K[j][k]));
							on[j][k]=false;
						}
					}

				}

				B[i][k]=B[i-1][k].add(dB.times(2));
			}
		}

		Vect[] Bsum=new Vect[L];

		for(int i=0;i<L;i++){
			Bsum[i]=new Vect(2);
			for(int j=0;j<nphi;j++)
				Bsum[i]=Bsum[i].add(B[i][j]);
		}


		Mat BH=new Mat(L,4);
		for(int i=0;i<L;i++){
			BH.el[i][0]=H.el[i][0];
			BH.el[i][1]=H.el[i][1];
			BH.el[i][2]=Bsum[i].el[0];
			BH.el[i][3]=Bsum[i].el[1];
		}

		//BH=BH.times(projCoef);


		return BH;

	}


	public  double getH(Mat BH,double B){

		if(B<=0)
			return 0;
		
		int i=BH.nRow-1;
		if(B>=BH.el[i][1])
			return BH.el[i][0]+(B-BH.el[i][1])*(BH.el[i][0]-BH.el[i-1][0])/(BH.el[i][1]-BH.el[i-1][1]);;
			//
			//return BH.el[BH.nRow-1][0];
		
		int j=0;

		while(BH.el[j+1][1]<B){j++;}
		
		double cc=(BH.el[j+1][0]-BH.el[j][0])/(BH.el[j+1][1]-BH.el[j][1]);

		//double H= BH.el[j][0];
		double H=BH.el[j][0]+(B-BH.el[j][1])*cc;

		
		return H;
				
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

		Mat BH=this.getCurveAlt(H);

		return BH;

	}


	public Mat magnetizeUpTo(double Bpeak,int L){

		Vect H=new Vect().linspace(0, Hs, L);

			
		Vect[][] B=new Vect[L][nphi];
		for(int i=0;i<L;i++)
			for(int k=0;k<nphi;k++)
				B[i][k]=new Vect(2);

		Vect Hr=new Vect(L);
		

		for(int kh=-nh;kh<=nh;kh++){


			int k=kh+nh;

			Vect er=new Vect(cos(kh*dphiRad),sin(kh*dphiRad));

		
				Hr=H.times(cos(kh*dphiRad));

			for(int i=0;i<L;i++){


				if(i==0){
					//B[i][k]=er.times(this.getRes());
					B[i][k]=B[i][k].add(this.getRes(k));

					continue;
				}

				Vect dB=new Vect(2);
				for(int j=0;j<M;j++)
				{

					if(Hr.el[i]>Hr.el[i-1] && Hr.el[i]>a[j][1][k]){
						if(!on[j][k]){
							dB=dB.add(er.times(this.DB2D*K[j][k]));
							on[j][k]=true;
						}
					}
					else if(Hr.el[i]<=Hr.el[i-1] && Hr.el[i]<a[j][0][k] ){

						if(on[j][k]){
							dB=dB.sub(er.times(this.DB2D*K[j][k]));
							on[j][k]=false;
						}
					}

				}

				B[i][k]=B[i-1][k].add(dB.times(2));
				

			}
		}
		

		Vect[] Bsum=new Vect[L];

		for(int i=0;i<L;i++){
			Bsum[i]=new Vect(2);
			for(int j=0;j<nphi;j++)
				Bsum[i]=Bsum[i].add(B[i][j]);
		}
		
		int ix=0;
		for(int i=0;i<L;i++){
			if(Bsum[i].el[0]>Bpeak) break;
			ix++;
		}

		Mat BH1=new Mat(ix,3);
		for(int i=0;i<ix;i++){
			BH1.el[i][0]=H.el[i];
			BH1.el[i][1]=Bsum[i].el[0];
			BH1.el[i][2]=Bsum[i].el[1];
		}

		 Mat BH=this.distill(BH1);


		 
		return BH;

	}
	/*	
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
					if(!on[j][kk]){
						dB+=this.DB*K[j][kk];
						on[j][kk]=true;
					}
				}
				else if(H.el[i]<H.el[i-1] && H.el[i]<a[j][0][kk] ){

					if(on[j][kk]){
						dB-=this.DB*K[j][kk];
						on[j][kk]=false;
					}
				}

			}

			B.el[i]=B.el[i-1]+2*dB;

			ix++;
			if(B.el[i]<Bpeak || Math.abs(B.el[i]-Bpeak)<=epsdB) {
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


	}*/




	/*
public Mat getCurveX(Vect H){

	Mat H1=new Mat(H.length,2);
	H1.setCol(H, 0);

	Mat BH=getLocus(H1);

	Mat BH1=new Mat(H.length,2);

	BH1.setCol(H, 0);
	BH1.setCol(BH.getColVect(2), 1);

	return BH1;

}*/



	/*public Mat getLocus2D(Mat H){


		int dim=H.nCol;

		int n1=this.nphi;


		Preisach[] pr=new Preisach[n1];

		int L=H.nRow;

		Vect Hn=new Vect(L);
		Vect[] Bn=new Vect[n1];

		Mat BH=new Mat(L,4);


			Vect er=new Vect(dim);

			for(int k=0;k<n1;k++){

				int kp=k;
				kp=0;
			//	kp=0;

				pr[kp]=this.deepCopy();

				pr[kp].kk=kp;
				pr[kp].kk=kp;

				pr[kp].demagnetize();

				double phirad= phi[kp]*Math.PI/180;

				er.el[0]=Math.cos(phirad);
				er.el[1]=Math.sin(phirad);

				for(int i=0;i<L;i++)	
					Hn.el[i]=new Vect(H.el[i]).dot(er);



				//Mat BH1=pr[kp].getCurveAlt(Hn);
				Mat BH1=this.getCurveAlt(Hn);

				Bn[kp]=BH1.getColVect(1);




				for(int i=0;i<L;i++){

					BH.el[i][2]=BH.el[i][2]+Bn[kp].el[i]*er.el[0];
					BH.el[i][3]=BH.el[i][3]+Bn[kp].el[i]*er.el[1];
				}
			}

			util.pr(projCoef);

			for(int i=0;i<L;i++){
				BH.el[i][2]*=1.0/n1*projCoef;
				BH.el[i][3]*=1.0/n1*projCoef;
				BH.el[i][0]=H.el[i][0];
				BH.el[i][1]=H.el[i][1];
			}


		return BH;

	}*/

	public Vect getRes(){

		Vect Br=new Vect(2);


		for(int kh=-nh;kh<=nh;kh++){

			int k=kh+nh;


			Vect er=new Vect(cos(kh*dphiRad),sin(kh*dphiRad));

			for(int j=0;j<M;j++)
			{

				if(on[j][k])
					Br=Br.add(er.times(this.DB2D*K[j][k]));
				else
					Br=Br.add(er.times(-this.DB2D*K[j][k]));

			}

		}


		return Br;
	}

	public Vect getRes(int k){

		int kh=k-nh;
		Vect Br=new Vect(2);

		Vect er=new Vect(cos(kh*dphiRad),sin(kh*dphiRad));

		for(int j=0;j<M;j++)
		{

			if(on[j][k])
				Br=Br.add(er.times(this.DB2D*K[j][k]));
			else
				Br=Br.add(er.times(-this.DB2D*K[j][k]));

		}




		return Br;
	}


	/*	public Vect getResV(){
		int ix=0;

		Vect Bres= new Vect(dim);
		Vect eph=new Vect(dim);

		for(int k=0;k<this.nphi;k++){

			double phirad=PI*phi[k]/180;

			eph.el[0]=Math.cos(phirad);
			eph.el[1]=Math.sin(phirad);

		for(int j=0;j<M;j++){

			if(on[j][k]) {
				Bres=Bres.add(eph.times(this.DBv));

			}
			else{
				ix++;
				Bres=Bres.sub(eph.times(this.DBv));
			}


		}

	}

		util.pr(ix*this.DBv);

		return Bres;
	}*/

	/*
	public  Mat initial(int L){


		return 	initCurveUptoBpeak(this.Bs,L);
	}


	public  Mat initCurveUptoBpeak(double Bpeak,int L){
		this.demagnetize();

		Vect seqH=new Vect().linspace(0, this.Hs, L);

		Mat BH1=this.getCurveAlt(seqH);

		Mat BH=this.distill(BH1);


		return BH;
	}

	}*/
	
	
		public  Mat symMajorDesc(int L){
			return symDesc(Bs,L);
		}

		public  Mat symDesc(double Bpeak,int L){
			

		this.demagnetize();

		Mat BH1=this.magnetizeUpTo(Bpeak,L);

		int L1=BH1.nRow;
		double H1=BH1.el[L1-1][0];
	
		this.demagnetize();
		Vect seqH=new Vect().linspace(0, H1, L);
		Mat BH2=this.getCurveAlt(seqH);
		
		seqH=new Vect().linspace(H1, -Hs, L);
		
		 BH2=this.getCurveAlt(seqH);
		
		Mat BH3=this.distill(BH2);
		


		int ix=0;

		for(int i=0;i<BH3.nRow;i++){
			if(BH3.el[i][1]>=-Bpeak)
				ix++;
		}


		Mat BH4=new Mat(ix,3);

		for(int i=0;i<ix;i++){
			BH4.el[i][0]=BH3.el[i][0];
			BH4.el[i][1]=BH3.el[i][1];
			BH4.el[i][2]=BH3.el[i][2];
		}


		return BH4;
	}


	/*public  Mat symAsc(double Bpeak,int L){

		this.demagnetize();

		Mat BH1=this.magnetizeDownTo(-Bpeak,L);

		int L1=BH1.nRow;
		double H1=BH1.el[L1-1][0];

		Vect seqH=new Vect().linspace(H1, Hs, L);


		Mat BH2=this.getCurveAlt(seqH);

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
				if(!on[kk][j]){
					dB+=this.DB*K[j][kk];
					on[j][kk]=true;
				}
			}
			else if(H.el[i]<H.el[i-1] && H.el[i]<a[j][0][kk] ){

				if(on[j][kk]){
					dB-=this.DB*K[j][kk];
					on[j][kk]=false;
				}
			}

		}

		B.el[i]=B.el[i-1]+2*dB;

		ix++;
		if(B.el[i]>Bpeak || Math.abs(B.el[i]-Bpeak)<=epsdB) {
			B.el[i]=Bpeak;
			break;
		}


	}

	double H1=H.el[ix-1];

	H=new Vect().linspace(H1, -Hs, L);

	Mat BH1=this.getCurveAlt(H);

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
					if(!on[j][kk]){
						dB+=this.DB*K[j][kk];
						on[j][kk]=true;
					}
				}
				else if(H.el[i]<H.el[i-1] && H.el[i]<a[j][0][kk] ){

					if(on[j][kk]){
						dB-=this.DB*K[j][kk];
						on[j][kk]=false;
					}
				}

			}

			B.el[i]=B.el[i-1]+2*dB;

			ix++;
			if(B.el[i]<Bpeak || Math.abs(B.el[i]-Bpeak)<=epsdB) {
				B.el[i]=Bpeak;
				break;
			}


		}



		double H1=H.el[ix-1];

		H=new Vect().linspace(H1, Hs, L);

		Mat BH1=this.getCurveAlt(H);

		 Mat BH=this.distill(BH1);

		return BH;

		}
*/
	public Mat distill(Mat BH1){

		int Leff=1;

		int L=BH1.nRow;
		boolean[] skip=new boolean[L];
		for(int i=1;i<BH1.nRow;i++){


			
			if( Math.abs(BH1.el[i][1]-BH1.el[i-1][1])/*/Math.abs(BH1.el[i][0]-BH1.el[i-1][0])*/<epsdB)
			//	if( Math.abs(BH1.el[i][2]-BH1.el[i-1][2])/Math.abs(BH1.el[i][0]-BH1.el[i-1][0])<epsdBdH)
				{
				skip[i]=true;

				continue;
			}


			Leff++;
		}

		Mat BH=new Mat(Leff,3);

		int ix=0;
		for(int i=0;i<BH1.nRow;i++){

			if(!skip[i]){
				BH.el[ix][0]=BH1.el[i][0];
				BH.el[ix][1]=BH1.el[i][1];
				BH.el[ix][2]=BH1.el[i][2];
				ix++;

			}


		}




		return BH;
	}
	
}
