package fem;

import math.Vect;
import math.util;

public class PlayModel {

	
	public PlayModel()
	{	}

	public static void main(String[] args)
	{
		int L=1000;
		int M=1;
		Vect B=new Vect(2*L);
		Vect H=new Vect(2*L);
		Vect[] pk=new Vect[M];
		Vect[] sk=new Vect[M];
		double pkp=0;
		double skp=0;
		double Bp=0;
		double Xs=4;
		double pkps=0;

		double[] zeta=new double[M];
		double[] ita=new double[M];
		
		for(int k=0;k<M;k++){
			zeta[k]=1*(.5+k*1.0/M);
			ita[k]=1*(.8+k*1.0/M);
			pk[k]=new Vect(2*L);
			sk[k]=new Vect(2*L);
			}
		
		for(int i=-L;i<L;i++){
			int n=i+L;
			//double x=i*4.0/L;
			double x=10*Math.sin(2*Math.PI*n*1.0/L);
			//double x=n*.2;
			H.el[n]=x;
			//double pkp=-1000;
			for(int k=0;k<M;k++){
			if(n>0) pkp=pk[k].el[n-1];
			pkps=pks(pkp,Xs,zeta[k]);
			pk[k].el[n]=pk(pkps,x,zeta[k]);
			
			B.el[n]+=func(pk[k].el[n]);
		}
		}
		util.plot(H,B);
		
		 B=new Vect(2*L);
		 H=new Vect(2*L);
		
		for(int i=-L;i<-L;i++){
			int n=i+L;
			//double x=i*4.0/L;
			double x=Math.sin(2*Math.PI*i*1.0/L)+0*Math.sin(12*Math.PI*i*1.0/L);
			
			B.el[n]=x;
			
			if(n>0){
				Bp=B.el[n-1];
			
			}
			for(int k=0;k<M;k++){
			if(n>0){
				skp=sk[k].el[n-1];
			
			}
			sk[k].el[n]=sk(skp,B.el[n],Bp,ita[k]);
			
			H.el[n]=func2(sk[k].el[n]);
		}
		}
	
		//util.plot(H,B);
		//util.plot(H,sk[0]);
		

	}
	
	public static double func(double x)
	{
		
	double y=Math.cbrt(x);

		return y;

	}
	
	public static double func2(double x)
	{
		
		double y=Math.pow(x, 3);

	
	return y;
		

	}
	
	public static double pk(double pkp,double H,double zeta)
	{
		
	double y=0;
	
	y=Math.max(Math.min(pkp,H+zeta),H-zeta);
	
	return y;

	}
	
	public static double pks(double pkp,double H,double zeta)
	{
		
	double y=0;
	
	y=Math.max(Math.min(pkp,H-zeta),-H+zeta);
	
	return y;

	}
	
	public static double sk(double skp,double B,double B0,double ita)
	{
		
	double y=0;
	
	y=Math.max(Math.min(B-B0+skp,ita),-ita);
	
	return y;
		

	}
}
