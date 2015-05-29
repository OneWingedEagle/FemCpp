package math;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.Random;

import math.Mat;
import math.Vect;
import math.util;

public class PlayModel {

	
	public PlayModel()
	{	}

	public static void main(String[] args)
	{

		
		double Bs=1.5;
		double Hs=1600;
	
		int nInit=1;
		int nMajor=1;
		int nSymLoops=1;
		int nDescending=0;
		int nAscending=0;
		int nTot=nInit+nMajor+nSymLoops+nDescending+nAscending;
		
		int Mp=2000;
		double mean=200;
		double width=200;
		
		Preisach ps=new Preisach(Mp,mean,width,Hs,Bs,3564656);
	/*	
		Mat[] BHs=new Mat[nTot];
		
		BHs[0]=initLoop(ps,Bs,100,1.1*Hs);
		ps.demagnetize(50);
		if(nMajor==1)
		BHs[1]=symMajorHalf(ps,200);
		ps.demagnetize(50);
		if(nTot>=2)
		for(int i=2;i<nTot;i++)
		BHs[i]=symLoopHalf(ps,i*.9*Bs/(nTot+1),100,Hs);
		util.plotBunch(BHs);
		//util.plot(BHs[2].el);
		
		String hysdatafile="C:\\Works\\HVID\\hys_dataH";


		
		try{
			PrintWriter pwBun = new PrintWriter(new BufferedWriter(new FileWriter(hysdatafile)));		
	
			pwBun.println(1+"\t"+1);
			pwBun.println("*Bs*Hs*");
			pwBun.println(Bs+"\t"+Hs);

			pwBun.println("* 初磁化曲線数 * メジャーループ数 * 対称ループ数 * 下降曲線数 * 上昇曲線数 *");
			pwBun.println(nInit+"\t"+nMajor+"\t"+nSymLoops+"\t"+nDescending+"\t"+nAscending);
			
			for(int i=0;i<nTot;i++){
				pwBun.println("*xxx");
				pwBun.println(BHs[i].nRow);
				for(int j=0;j<BHs[i].nRow;j++)
					pwBun.println(BHs[i].el[j][0]+"\t"+BHs[i].el[j][1]);
			}
			
			pwBun.println("* ----- 回転ヒステリシス損");
			pwBun.println("* B数 *");
			pwBun.println("0");
			pwBun.println("* B * 損失");
			pwBun.println("* ----- 異方性");
			pwBun.println("* B数 * 角度数 *");
			pwBun.println(0+"\t"+0);
	
			
			pwBun.close();
		}
		catch(IOException e){}
	;*/
	//reversal();
	//curves();
	
/*		Preisach ps=new Preisach(1000,200,500,1000,2,3564656);
		ps.magnetize(10,0);
		ps.demagnetize(20);
		ps.magnetize(10,0)*/;
	/*	ps.demagnetize(10);
		ps.magnetize(10);*/

		
		//ps.magnetize(10);
		
/*		int L=100;
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
	*/
	



	}
	
	public static Mat symLoop(Preisach ps,double Bpeak,int LinitMax, double Hmax){
				

		int M=ps.M;
		
		double[] DM=new double[M];
		
		
		for(int j=0;j<M;j++){

			DM[j]=1;
		}
		

	
		

/*		Vect seqInitial0=new Vect().linspace(0, Math.sqrt(Hs), LinitMax);
		Vect seqInitial1=seqInitial0.times(seqInitial0);*/
		
		Vect seqInitial1=new Vect().linspace(0, Hmax, LinitMax);

		
		Vect H=new Vect(LinitMax);
		Vect B=new Vect(LinitMax);
		
		int Linit=0;
		
		for(int i=0;i<seqInitial1.length;i++){
			
			H.el[i]=seqInitial1.el[i];
	
			if(i==0) 
				{
				B.el[i]=0;
				continue;
				}
			
			for(int j=0;j<M;j++)
			{

				if(H.el[i]>H.el[i-1] && H.el[i]>ps.a[j][1]){
					if(!ps.on[j]){
						B.el[i]+=ps.DB*DM[j];
						ps.on[j]=true;
					}
					}
				else if(H.el[i]>H.el[i-1]  && H.el[i]<ps.a[j][0] ){
				
					if(ps.on[j]){
						B.el[i]-=ps.DB*DM[j];;
						ps.on[j]=false;
					}
		
					}

				
			}
			if(	Linit>0) 
				B.el[Linit]+=B.el[Linit-1];
			
		
		
			
			if(B.el[Linit]>=ps.Bs-ps.eps){
				H.el[Linit]=ps.Hs;
				B.el[Linit]=ps.Bs;
				util.pr(i);
		
				Linit++;
				break;
			}
			
			Linit++;

			}


		
		Vect seqInitial=new Vect(Linit);
		for(int i=0;i<seqInitial.length;i++){
			seqInitial.el[i]=seqInitial1.el[i];
				}

	
		Vect seqHdown=new Vect(2*Linit-2);

		for(int i=0;i<seqHdown.length;i++){
			if(i<Linit-1)
			seqHdown.el[i]=seqInitial.el[Linit-i-2];
			else
			seqHdown.el[i]=-seqInitial.el[i-Linit+2];
		}

		Vect seqHup=new Vect(seqHdown.length-1);
		for(int i=0;i<seqHup.length;i++)
			seqHup.el[i]=seqHdown.el[seqHdown.length-2-i];

		Vect seqH=seqInitial.aug(seqHdown).aug(seqHup);
		
	

	
		ps.demagnetize(50);

		Mat BH1=ps.getCurve(seqH, 0);
		
		Mat BH=new Mat(BH1.nRow-Linit+2,2);
		
		for(int i=0;i<BH.nRow-1;i++){
			BH.el[i][0]=BH1.el[i+Linit-1][0];
			BH.el[i][1]=BH1.el[i+Linit-1][1];
		}

		BH.el[BH.nRow-1][0]=BH1.el[Linit-1][0];
		BH.el[BH.nRow-1][1]=BH1.el[Linit-1][1];


		
		return BH;
	}
	
	public static Mat symLoopHalf(Preisach ps,double Bpeak,int LinitMax,double Hmax){
		
		int M=ps.M;
		
		double[] DM=new double[M];
		
		
		for(int j=0;j<M;j++){

			DM[j]=1;
		}
		

		
		/*		Vect seqInitial0=new Vect().linspace(0, Math.sqrt(Hs), LinitMax);
		Vect seqInitial1=seqInitial0.times(seqInitial0);*/
		
		Vect seqInitial1=new Vect().linspace(0, Hmax, LinitMax);

	
		
		Vect H=new Vect(LinitMax);
		Vect B=new Vect(LinitMax);
		
		int Linit=0;
		
		for(int i=0;i<seqInitial1.length;i++){
			
			H.el[i]=seqInitial1.el[i];
	
			if(i==0) 
				{
				B.el[i]=0;
				continue;
				}
			
			for(int j=0;j<M;j++)
			{

				if(H.el[i]>H.el[i-1] && H.el[i]>ps.a[j][1]){
					if(!ps.on[j]){
						B.el[i]+=ps.DB*DM[j];
						ps.on[j]=true;
					}
					}
				else if(H.el[i]>H.el[i-1]  && H.el[i]<ps.a[j][0] ){
				
					if(ps.on[j]){
						B.el[i]-=ps.DB*DM[j];;
						ps.on[j]=false;
					}
		
					}

				
			}
			
			if(	Linit>0) 
				B.el[Linit]+=B.el[Linit-1];
			
		
		
			
			if(B.el[Linit]>=Bpeak-ps.eps){
				B.el[Linit]=Bpeak;
				Linit++;
				break;
			}
			
			Linit++;

			}


		

		
		Vect seqInitial=new Vect(Linit);
		for(int i=0;i<seqInitial.length;i++){
			seqInitial.el[i]=seqInitial1.el[i];
				}

	
	
		Vect seqHdown=new Vect(2*Linit-2);

		for(int i=0;i<seqHdown.length;i++){
			if(i<Linit-1)
			seqHdown.el[i]=seqInitial.el[Linit-i-2];
			else
			seqHdown.el[i]=-seqInitial.el[i-Linit+2];
		}

	

		Vect seqH=seqInitial.aug(seqHdown);
		
	

	
		ps.demagnetize(50);

		Mat BH1=ps.getCurve(seqH, 0);

		
		Mat BH=new Mat(BH1.nRow-Linit+1,2);
		
		for(int i=0;i<BH.nRow;i++){
			BH.el[i][0]=BH1.el[i+Linit-1][0];
			BH.el[i][1]=BH1.el[i+Linit-1][1];
		}


		
		return BH;
	}
public static Mat symMajorHalf(Preisach ps,int LinitMax){
		
		int M=ps.M;
		
		double[] DM=new double[M];
		
		
		for(int j=0;j<M;j++){

			DM[j]=1;
		}
		

		
		/*		Vect seqInitial0=new Vect().linspace(0, Math.sqrt(Hs), LinitMax);
		Vect seqInitial1=seqInitial0.times(seqInitial0);*/
		
		ps.demagnetize(50);
	//	ps.magnetize(50,0);
	
		
	//	Vect seqH=new Vect().linspace(ps.Hs, -ps.Hs, LinitMax);

		Vect seqH=new Vect().linspace(0, ps.Hs/2, 100);
	
		Mat BH1=ps.getCurve(seqH, ps.Bs);

	
/*		Mat BH=new Mat(Linit,2);
		
		
		for(int i=0;i<BH.nRow;i++){
			BH.el[i][0]=H.el[i];
			BH.el[i][1]=B.el[i];
		}
*/
BH1.show();
		
		return BH1;
	}
	
	public static Mat initLoop(Preisach ps,double Bpeak,int LinitMax, double Hmax){
		
		int M=ps.M;
		
		double[] DM=new double[M];
		
		
		for(int j=0;j<M;j++){

			DM[j]=1;
		}
		

	
		

		/*		Vect seqInitial0=new Vect().linspace(0, Math.sqrt(Hs), LinitMax);
		Vect seqInitial1=seqInitial0.times(seqInitial0);*/
		
		Vect seqInitial1=new Vect().linspace(0, Hmax, LinitMax);
	//	seqInitial1.show();
	
		Vect H=new Vect(LinitMax);
		Vect B=new Vect(LinitMax);
		
		int Linit=0;
		
		for(int i=0;i<seqInitial1.length;i++){
			
			H.el[i]=seqInitial1.el[i];
	
			if(i==0) 
				{
				B.el[i]=0;
				continue;
				}
			
			for(int j=0;j<M;j++)
			{

				if(H.el[i]>H.el[i-1] && H.el[i]>ps.a[j][1]){
					if(!ps.on[j]){
						B.el[i]+=ps.DB*DM[j];
						ps.on[j]=true;
					}
					}
				else if(H.el[i]>H.el[i-1]  && H.el[i]<ps.a[j][0] ){
				
					if(ps.on[j]){
						B.el[i]-=ps.DB*DM[j];;
						ps.on[j]=false;
					}
		
					}

				
			}
			
			if(	Linit>0) 
				B.el[Linit]+=B.el[Linit-1];
			
		
		
			
			if(B.el[Linit]>=ps.Bs-ps.eps ){
		
				B.el[Linit]=ps.Bs;
				Linit++;
				break;
			}
			
			Linit++;
	

			}

		
		Mat BH=new Mat(Linit,2);
		
		for(int i=0;i<BH.nRow;i++){
			BH.el[i][0]=H.el[i];
			BH.el[i][1]=B.el[i];
		}

		
		return BH;
	}
	
	public static void curves()
	{

		int M=1000;
		double mean=200;
		double width=500;
		double Hs=1000;
		
		Preisach ps=new Preisach(M,mean,width,Hs,2,3564656);
		ps.demagnetize(10);
		
		
		int L=10000;
		
		Vect B=new Vect(2*L+1);
		Vect H=new Vect(2*L+1);
		
		double[] DM=new double[M];

		
		for(int j=0;j<M;j++){

			
			DM[j]=1;
		}
		

		 double Bs=2;
		double Bs1=2*Bs/M;


			
		 B=new Vect(2*L+1);
		 H=new Vect(2*L+1);
		
		for(int i=-L;i<=L;i++){
			
			int n=i+L;
			//double x=(.1+.001*n)*Math.sin(20*Math.PI*i*1.0/L);
			double x=Math.sin(2*Math.PI*i*1.0/L)-.0*Math.cos(8*Math.PI*i*1.0/L)+.0*Math.cos(12*Math.PI*i*1.0/L);
		
			x*=500;
			
/*			if(n<=500)
			x=n*.5;
			else if(n<=700)
				x=250-(n-500)*.5;
			else x=150+(n-700)*.5;*/
			
			
			H.el[n]=x;
	
			if(n==0) 
				{
				B.el[n]=0;
				continue;
				}
			
			for(int j=0;j<M;j++)
			{

				if(H.el[n]>H.el[n-1] && H.el[n]>ps.a[j][1]){
					if(!ps.on[j]){
					B.el[n]+=Bs1*DM[j];
					ps.on[j]=true;
					}
					}
				else if(H.el[n]<H.el[n-1] && H.el[n]<ps.a[j][0] ){
				
					if(ps.on[j]){
						B.el[n]-=Bs1*DM[j];;
						ps.on[j]=false;
					}
		
					}
		
					
				
			}
			if(n>0) 
			B.el[n]+=B.el[n-1];

			//if(i<10) util.pr(x+"  "+B.el[n]);
			}

	
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
		int Q=10;
		
		double[][] Ev=new double[P][Q];
		
	
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
