package PlayModel;

import java.util.Random;

import static java.lang.Math.*;
import math.Mat;
import math.Vect;
import math.util;

public class Preisach2D {

	public Random r;
	public int M,nphi,nh,dim,kRotated;
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
		cfw=2;

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
		



/*		Mat R=ps.getLocusHRotation(.4);
	
		R.show();
		*/
		



	}
	
	public  Mat initial(double Bpeak,int L){
	
		this.demagnetize();
		Mat BH1=this.magnetizeUpTo(Bpeak,L);

		
		return BH1;
	}

/*	public  Mat initial(double Hm,int L, boolean distill){

	
		this.demagnetize();

		Vect seqH=new Vect().linspace(0, Hm, L);

		Mat BH1=this.getCurveAlt(seqH);
		util.plot(BH1.getColVect(0),BH1.getColVect(1));
		
		if(!distill) return BH1;

		Mat BH=this.distill(BH1);


		return BH;
	}
*/

	public Mat getCurveAlt(Vect H){

		int L=H.length;
		Vect[][] B=new Vect[L][nphi];
		for(int i=0;i<L;i++)
			for(int k=0;k<nphi;k++)
				B[i][k]=new Vect(2);


		for(int kh=-nh;kh<=nh;kh++){


			int k=kh+nh;
			
			int kr=(k+kRotated)%nphi;

			Vect Halt=H.times(cos(kh*dphiRad));


			Vect er=new Vect(cos(kh*dphiRad),sin(kh*dphiRad));

			for(int i=0;i<L;i++){

				if(i==0){

					B[i][k]=B[i][k].add(this.getRes(kr));
					continue;
				}

				Vect dB=new Vect(2);
				for(int j=0;j<M;j++)
				{

					if(Halt.el[i]>Halt.el[i-1] && Halt.el[i]>a[j][1][kr]){
						if(!on[j][kr]){
							dB=dB.add(er.times(this.DB2D*K[j][kr]));
							on[j][kr]=true;
						}
					}
					else if(Halt.el[i]<=Halt.el[i-1] && Halt.el[i]<a[j][0][kr] ){

						if(on[j][kr]){
							dB=dB.sub(er.times(this.DB2D*K[j][kr]));
							on[j][kr]=false;
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

		

		for(int k=0;k<this.nphi;k++)
			for(int j=0;j<this.M;j++)
				on[j][k]=false;
			
		double t=0;
		double dt=1./L;

		//double cc=2.0/nCycles;

		Vect H=new Vect(nCycles*L);

		int  ix=0;
		for(int i=0;i<nCycles;i++)
			for(int j=0;j<L;j++){

				double x=(1-(t+dt)/nCycles)*Math.sin(2*Math.PI*t-PI/2);
				//double x=Math.exp(-cc*t)*Math.sin(2*Math.PI*t);
				t+=dt;
				H.el[ix]=x*Hm;
				ix++;
			}
		
		Mat Hp=new Mat(nCycles*L,2);

		  ix=0;
		  t=0;
		for(int i=0;i<nCycles;i++)
			for(int j=0;j<L;j++){

				double x=(1-(t+dt)/nCycles)*Math.sin(2*Math.PI*t-PI/2);
				double y=(1-(t+dt)/nCycles)*Math.cos(2*Math.PI*t-PI/2);
				//double x=Math.exp(-cc*t)*Math.sin(2*Math.PI*t);
				t+=dt;
				Hp.el[ix][0]=x*Hm;
				Hp.el[ix][1]=y*Hm;
				ix++;
			}


		Mat BH=this.getLocus(Hp);

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
		
		Vect seqH=new Vect().linspace(0, H1, L).aug(new Vect().linspace(H1, -Hs, L));
		BH1=this.getCurveAlt(seqH);

	
		
	int jx=0;
		
		while(jx<BH1.nRow && BH1.el[jx+1][0]>=BH1.el[jx][0]){jx++;}


		
		Mat BH2=new Mat(BH1.nRow-jx,3);
	

		for(int i=0;i<BH2.nRow;i++)
			BH2.el[i]=BH1.el[i+jx];

		
	//util.plot(BH2.getColVect(0),BH2.getColVect(1));

		
		Mat BH3=this.distill(BH2);

		int ix=0;
		
		while(ix<BH3.nRow && BH3.el[ix][1]>=-Bpeak){ix++;}

		Mat BH4=new Mat(ix,3);

		for(int i=0;i<ix;i++){
			BH4.el[i][0]=BH3.el[i][0];
			BH4.el[i][1]=BH3.el[i][1];
			BH4.el[i][2]=BH3.el[i][2];
		}

		//util.plot(BH4.getColVect(0),BH4.getColVect(1));	 
	
		return BH4;
	}
		


		public Mat distill(Mat BH1){

			int Leff=1;
			boolean col3=(BH1.nCol==3);

			int L=BH1.nRow;
			boolean[] skip=new boolean[L];
			for(int i=1;i<BH1.nRow;i++){

				if( Math.abs(BH1.el[i][1]-BH1.el[i-1][1])/Math.abs(BH1.el[i][0]-BH1.el[i-1][0])<epsdB)
					if( Math.abs(BH1.el[i][1]-BH1.el[i-1][1])/Math.abs(BH1.el[i][0]-BH1.el[i-1][0])<epsdBdH)
					{
					skip[i]=true;

					continue;
				}


				Leff++;
			}

			int ncol=2;
			if(col3) ncol=3;
			Mat BH=new Mat(Leff,ncol);

			int ix=0;
			for(int i=0;i<BH1.nRow;i++){

				if(!skip[i]){
					BH.el[ix][0]=BH1.el[i][0];
					BH.el[ix][1]=BH1.el[i][1];
					if(col3)
					BH.el[ix][2]=BH1.el[i][2];
					ix++;

				}


			}




			return BH;
		}


	
}
