package fem;

import java.util.Random;

import math.Vect;
import math.util;

public class PlayModel {

	
	public PlayModel()
	{	}

	public static void main(String[] args)
	{
		int L=100;
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
		
		for(int i=-L;i<-L;i++){
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
		//util.plot(H,B);
		
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
	
	

	reversal();

	}
	
	public static void curves()
	{

		
		int M=1000;
		int L=1000;
		
		Vect B=new Vect(2*L+1);
		Vect H=new Vect(2*L+1);
		Random r=new Random(3564656);
		
		double[][] a=new double[M][2];
		boolean[] on=new boolean[M];
		
		for(int j=0;j<M;j++){
			/*a[j][0]=30*r.nextGaussian()-100;//*(.5-r.nextGaussian());
			a[j][1]=30*r.nextGaussian()+100;//*(.5+.2*r.nextGaussian());
*/			double am=500*r.nextGaussian();
			double d=500*r.nextDouble();
			a[j][0]=am-d/2;
			a[j][1]=am+d/2;
		}
		
	/*	 B=new Vect(2*L+1);
		 H=new Vect(2*L+1);
		*/
		 double Bs=2;
		double Bs1=2*Bs/M;
		
		
		for(int i=-L;i<=L;i++){
			
			int n=i+L;
			double x=-(.1+.001*n)*Math.cos(20*Math.PI*i*1.0/L);
			//double x=-Math.cos(2*Math.PI*i*1.0/L)-.4*Math.cos(8*Math.PI*i*1.0/L)+.2*Math.cos(12*Math.PI*i*1.0/L);
			
			x*=1000;
/*			if(n<=500)
			x=n*.5;
			else if(n<=700)
				x=250-(n-500)*.5;
			else x=150+(n-700)*.5;*/
			
			H.el[n]=x;
	
			if(n==0) 
				{
				B.el[n]=-Bs;
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
			}

		util.plot(H);
	//util.show(a);
		util.plot(H,B);

	}
	
	public static void reversal()
	{

		
		int M=1000;
		int L=8000;
		
		Vect B=new Vect(L);
		Vect H=new Vect(L);
		Random r=new Random(3564656);
		
		double[][] a=new double[M][2];
		boolean[] on=new boolean[M];
		
		for(int j=0;j<M;j++){
			/*a[j][0]=30*r.nextGaussian()-100;//*(.5-r.nextGaussian());
			a[j][1]=30*r.nextGaussian()+100;//*(.5+.2*r.nextGaussian());
*/			double am=500*r.nextGaussian();
			double d=500*r.nextDouble();
			a[j][0]=am-d/2;
			a[j][1]=am+d/2;
		}
		
	/*	 B=new Vect(2*L+1);
		 H=new Vect(2*L+1);
		*/
		 double Bs=2;
		double Bs1=2*Bs/M;
		
		double Hs=1000;
	
		int P=10;
		
	
		//
		double dH= 16*Hs/L;
		
		int cp=1;
		H.el[0]=-Hs;
		B.el[0]=-Bs;
		

		
		double[] Ha=new double[P];
		for(int p=0;p<P;p++)
			Ha[p]=-Hs+Hs*2.0*p/P;
		
		int p=0;
		for(int i=1;i<L;i++){
	
		if(i>1){
		if(cp==1 && H.el[i-1]>Ha[p] )	{cp=-1; p++;}
			
		else if(cp==-1 && H.el[i-1]<=-Hs )	cp=1;
		//double x=-1000*Math.cos(4*Math.PI*i*1.0/L);
		}
		
		H.el[i]=H.el[i-1]+dH*cp;

				
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

		util.plot(H);
	//util.show(a);
		util.plot(H,B);

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
