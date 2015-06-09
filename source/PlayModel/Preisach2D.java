package PlayModel;

import java.util.Random;

import static java.lang.Math.*;
import math.Mat;
import math.Vect;
import math.util;

public class Preisach2D {

	public Random r;
	public int M,nphi,nh,dim;
	public long seed;
	public double cfm,cfw,mean,width,Hs,Bs,DB2D,projCoef,dphiRad;
	public boolean[][] on;
	public double[] phi;
	public double[][][] a;
	public double[][] K;
	public double  epsdB=1e-2;
	public double  epsdBdH=1e-5;


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

				K[j][kp]=1-.0*sin(2*phirad);

				double am=mean*r.nextGaussian()*(1+cfm*abs(sin(phirad)));
				double d=width*abs(r.nextGaussian())*(1+cfw*abs(sin(phirad)));

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

		//ps.getRes().hshow();


		Mat R=ps.getLocusHRotation(.4);
	
		R.show();
		
		



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
	
	
	public Mat getLocusBRotation(double Hm,int Lc,int Nc){
		
		int L=Nc*Lc;
		Mat Hp=new Mat(L,2);
		for(int j=0;j<L;j++){
			Hp.el[j][0]=Hm*Math.cos(4*j*Math.PI/L);
			Hp.el[j][1]=Hm*Math.sin(4*j*Math.PI/L);
		}
		
		return getLocus(Hp);
		
		
		
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



		return BH;

	}
	
	public Mat getLocusHRotation(double Br){

		this.demagnetize();
		
		int L=10000;
		Vect[] H=new Vect[L];
		
		Vect[][] B=new Vect[L][nphi];
		for(int i=0;i<L;i++){
			H[i]=new Vect(2);
			for(int k=0;k<nphi;k++)
				B[i][k]=new Vect(2);
		}
		
		double dphi=20*PI/L;
		

		int i=0;
		double kr=0;
		double Hrp,Hr=0;

		while(i<L){
			
		
			if(i>0){
				
		double Hp=400*(Math.PI/2+Math.atan(.005*kr));
				//kr=4000;
				
			H[i]=new Vect(cos(i*dphi),sin(i*dphi)).times(Hp);
			}
	//	H[i].hshow();
		for(int kh=-nh;kh<=nh;kh++){
	

			int k=kh+nh;
				
			Vect er=new Vect(cos(kh*dphiRad),sin(kh*dphiRad));

		
			Hrp=Hr;
				Hr=H[i].dot(er);


				if(i==0){

					B[i][k]=B[i][k].add(this.getRes(k));

					continue;
				}

				Vect dB=new Vect(2);
				for(int j=0;j<M;j++)
				{

					if(Hr>Hrp && Hr>a[j][1][k]){
						if(!on[j][k]){
							dB=dB.add(er.times(this.DB2D*K[j][k]));
							on[j][k]=true;
						}
					}
					else if(Hr<=Hrp && Hr<a[j][0][k] ){

						if(on[j][k]){
							dB=dB.sub(er.times(this.DB2D*K[j][k]));
							on[j][k]=false;
						}
					}

				}

				B[i][k]=B[i-1][k].add(dB.times(2));
			
			//}
		}
		
		Vect B1=new Vect(2);
		for(int k=0;k<nphi;k++)
			B1=B1.add(B[i][k]);
		double Bri=B1.norm();///nphi;

		if(Bri>Br) kr--;
		else kr++;
		
		util.pr(Bri);
		//util.pr(kr);
		i++;
		
		}

		Vect[] Bsum=new Vect[L];

		for(int j=0;j<L;j++){
			Bsum[j]=new Vect(2);
			for(int k=0;k<nphi;k++)
				Bsum[j]=Bsum[j].add(B[j][k]);
			
		
		}


		Mat BH=new Mat(L,4);
		for(int j=0;j<L;j++){
			BH.el[j][0]=H[j].el[0];
			BH.el[j][1]=H[j].el[1];
			BH.el[j][2]=Bsum[j].el[0];
			BH.el[j][3]=Bsum[j].el[1];
		}



		return BH;

	}



	public Mat demagnetize(){
		return demagnetize(Hs,200,20);
	}
	public Mat demagnetize(double Hm){
		return demagnetize(Hm,200,20);
	}

	public Mat demagnetize(int nCycles){
		return demagnetize(Hs,200,nCycles);
	}



	public Mat demagnetize(double Hm,int L,int nCycles){


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


	public Mat distill(Mat BH1){

		int Leff=1;

		int L=BH1.nRow;
		boolean[] skip=new boolean[L];
		for(int i=1;i<BH1.nRow;i++){


			
			if( Math.abs(BH1.el[i][1]-BH1.el[i-1][1])/Math.abs(BH1.el[i][0]-BH1.el[i-1][0])<epsdB)
				if( Math.abs(BH1.el[i][2]-BH1.el[i-1][2])/Math.abs(BH1.el[i][0]-BH1.el[i-1][0])<epsdBdH)
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
