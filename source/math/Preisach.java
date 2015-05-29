package math;

import java.util.Random;

import math.Mat;
import math.Vect;
import math.util;

public class Preisach {
	
	public Random r;
	public int M;
	public long seed;
	public double mean,width,Hs,Bs,DB;
	public boolean[] on;
	public double[][] a;
	public double  eps=3e-2;

	


	public Preisach(){}
	

	public Preisach(int M, double mean,double width,double Hmax,double Bs, long seed){
	
	this.M=M;
	this.mean=mean;
	this.width=width;
	this.seed=seed;
	this.Hs=Hmax;
	this.Bs=Bs;
	this.DB=2*Bs/M;

	
	 r=new Random(3564656);
	
	a=new double[M][2];
	on=new boolean[M];
	
	for(int j=0;j<M;j++){
			double am=mean*r.nextGaussian();
		double d=width*r.nextDouble();
		a[j][0]=am-d/2;
		
		a[j][1]=am+d/2;

	
	}
	
	demagnetize(10);
	
	}
	
	public void demagnetize(int nCycles){
		demagnetize(nCycles,-Hs,-Bs);
	}
	
	public void demagnetize(int nCycles,double Hin,double Bin){
		
		int L=400;
		double t=0;
		double dt=1./L;
	
		
		double cc=4.0/nCycles;
		
		double Bs1=2*Bs/M;
		
		Vect B=new Vect(nCycles*(L+L/4));
		Vect H=new Vect(nCycles*(L+L/4));
		
		int  ix=0;
		for(int i=0;i<nCycles;i++)
			for(int j=0;j<L;j++){
			
			double x=Math.exp(-cc*t)*Math.cos(2*Math.PI*t);
			t+=dt;
			H.el[ix]=x*Hin;
			ix++;
		}
	
	/*		for(int j=0;j<L/4;j++){
			
			double x=-1e-3*Math.cos(2*Math.PI*t);
			t+=dt;
			H.el[ix]=x*Hin;
			ix++;
		}
		*/
		

	
		
		
		int n=-1;
		
		for(int i=0;i<H.length;i++){

			n++;

			if(n==0) 
				{
				B.el[n]=Bin;
				continue;
				}
			
			for(int j=0;j<M;j++)
			{

				if(H.el[n]>H.el[n-1] && H.el[n]>a[j][1]){
					if(!on[j]){
					B.el[n]+=Bs1;
					on[j]=true;
					}
					}
				else if(H.el[n]<H.el[n-1] && H.el[n]<a[j][0] ){
				
					if(on[j]){
						B.el[n]-=Bs1;
					on[j]=false;
					}
		
					}
		
					
				
			}
			if(n>0) 
			B.el[n]+=B.el[n-1];

			//if(i<10) util.pr(x+"  "+B.el[n]);
			
	}
	
	//	util.plot(H,B);
		
	}
	
public void magnetize(int nCycles,double Bin,double Hmax){
		
		int L=400;
		double t=0;
		double dt=1./L;

		
		double cc=4.0/nCycles;
	
		Vect B=new Vect(nCycles*(L+L/4));
		Vect H=new Vect(nCycles*(L+L/4));
		
		
		int ix=0;
		
		for(int i=0;i<nCycles;i++)
			for(int j=0;j<L;j++){
			
				double x=(1-Math.exp(-cc*t))*Math.sin(2*Math.PI*t);;
			t+=dt;
			H.el[ix]=x*Hmax;
			ix++;
		}
	
			
		
		for(int j=0;j<L/4;j++){
		
			double x=(1-Math.exp(-cc*t))*Math.sin(2*Math.PI*t);
			t+=dt;
			H.el[ix]=x*Hmax;
			ix++;
		}
		
		int n=-1;
		
		for(int i=0;i<H.length;i++){


			n++;
	
	
			if(n==0) 
				{
				B.el[n]=Bin;
				continue;
				}
			
			for(int j=0;j<M;j++)
			{

				if(H.el[n]>H.el[n-1] && H.el[n]>a[j][1]){
					if(!on[j]){
					B.el[n]+=this.DB;
					on[j]=true;
					}
					}
				else if(H.el[n]<H.el[n-1] && H.el[n]<a[j][0] ){
				
					if(on[j]){
						B.el[n]-=this.DB;
					on[j]=false;
					}
		
					}
		
					
				
			}
			if(n>0) 
			B.el[n]+=B.el[n-1];

			//if(i<10) util.pr(x+"  "+B.el[n]);
			
	}
	
		util.plot(H,B);
		
	}

public Mat getCurve(Vect H, double Bin){
	

	int L=H.length;
	Vect B=new Vect(L);
	double Bs1=2*Bs/M;
	
	for(int i=0;i<L;i++){

	
		if(i==0){
			B.el[i]=Bin;
			continue;
			}
		
		if(i==0){
			if(H.el[i]>=Hs)
			B.el[i]=Bs;
			continue;
			}
		
		for(int j=0;j<M;j++)
		{

			if(H.el[i]>H.el[i-1] && H.el[i]>a[j][1]){
				if(!on[j]){
				B.el[i]+=Bs1;
				on[j]=true;
				}
				}
			else if(H.el[i]<H.el[i-1] && H.el[i]<a[j][0] ){
			
				if(on[j]){
					B.el[i]-=Bs1;
				on[j]=false;
				}
	
				}
	
				
			
		}

		B.el[i]+=B.el[i-1];

		
}

	
	Mat BH=new Mat(L,2);
	BH.setCol(H, 0);
	BH.setCol(B, 1);
	


	
	return BH;
	
}

}
