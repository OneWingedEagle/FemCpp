package math;

import io.HystDataLoader;

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
	
	public Mat shapeFunc;
	public int nSym, nDesc,nAsc, nInit, nMajor,nTotCurves;
	public boolean anisot;
	public double Bs,Hs;
	public Mat[] BHraw,BH;
	public double[] zk,pk;
	public int dim=2,nHyst;
	

	
	public PlayModel()
	{	}

	public static void main(String[] args)
	{
		
		PlayModel pm=new PlayModel();
		
		String file=System.getProperty("user.dir") + "\\hys_dataH.txt";

	
		
		pm.loadData(file);
	//	pm.distillData();
		
		pm.simulateData();
/*		for(int i=0;i<pm.BH.length;i++)
			pm.BH[i].transp().show();*/
	//	pm.createData();
	/*	
		Mat[] BH1=new Mat[2*pm.BH.length];
		
		for(int i=0;i<pm.BH.length;i++){
			BH1[i]=pm.BHraw[i];
			BH1[i+pm.BH.length]=pm.BH[i];
	
		}
		
		util.plotBunch(BH1);*/
	}
	
	public void loadData(String file){
		
		HystDataLoader hysLoader=new HystDataLoader();
		
		hysLoader.loadData(this,file);
		
		double db=Bs/(this.nSym);
	
		this.nHyst=2*this.nSym+1;
		
		pk=new double[this.nSym+1];
		zk=new double[this.nHyst];
		
		for(int i=0;i<pk.length;i++){
			pk[i]=i*db;
		}
		
		for(int i=0;i<zk.length;i++){
			zk[i]=i*db;
		}

		
		this.shapeFunc=new Mat(nHyst,nHyst);
		
		for(int i=0;i<shapeFunc.nRow;i++){
			for(int j=0;j<=i;j++){
				shapeFunc.el[i][j]=BH[i].el[j][0]-BH[i].el[j-2][0];

			}
		}
		
		/*for(int i=0;i<shapeFunc.nRow;i++){
			double cc=.1*(i+1);
			for(int j=0;j<shapeFunc.nCol;j++){
				shapeFunc.el[i][j]=Math.pow(j*.1, cc);

			}
		}*/


	}
	

	
	public void createData(){

		
		double Bs=1.5;
		double Hs=1600;

		int nInit=1;
		int nMajor=1;
		int nSymLoops=10;
		int nDescending=0;
		int nAscending=0;
		int nTot=nInit+nMajor+nSymLoops+nDescending+nAscending;

		int Mp=3000;
		double mean=200;
		double width=500;

		Preisach ps=new Preisach(Mp,mean,width,Hs,Bs,3564656);
		
		Mat[] BHs=new Mat[nTot];
		
		int ix=0;
		
		for(int i=0;i<nInit;i++){

		BHs[ix++]=ps.initial(500);
		//BHs[ix++]=ps.symMajorFull(200);
		
		}
		
		for(int i=0;i<nMajor;i++){
		//BHs[ix++]=ps.symMajorFull(500);
		BHs[ix++]=ps.symMajorDesc(500);

		}

		
	
		for(int i=0;i<nSymLoops;i++){
			double Bp=Bs-(i+1)*Bs/(nSymLoops+1);
			
		//BHs[ix++]=ps.symFull(Bp,500);
		BHs[ix++]=ps.symDesc(Bp,500);
		//util.pr(ps.getRes());
	
		}
	
			
		for(int i=0;i<nDescending;i++){
		
			double Bp=Bs*(1-2.0*(i+1)/(nDescending+1));
			util.pr(Bp);
			BHs[ix++]=ps.revDesc(Bp,200);
		}
		
		for(int i=0;i<nAscending;i++){
		
			double Bp=-Bs*(1-2.0*(i+1)/(nAscending+1));
				util.pr(Bp);
			BHs[ix++]=ps.revAsc(Bp,1000);
		}
		
	//	util.plotBunch(BHs);
		//util.plot(BHs[2].el);
		
		String hysdatafile=System.getProperty("user.dir") + "\\hys_dataH.txt";


		
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




	
	}
	
	public double shapeFuncAt(int kz, double p){
		

			
		double f;

/*		
		if(kz>zk.length-1)
			kz=zk.length-1;


		if(p<=this.pk[0] )
			return shapeFunc.el[kz][0];

		if(p>=this.pk[pk.length-1] )
			return shapeFunc.el[kz][pk.length-1];
*/
		
		int jp=getj(this.pk,p);

		
		f=shapeFunc.el[kz][jp]*(p-this.pk[jp])+shapeFunc.el[kz][jp+1]*(this.pk[jp+1]-p);
	
		f*=1.0/(this.pk[jp+1]-this.pk[jp]);
		
	
		f=Math.cbrt(p);
		
		return f;
		
	}
	
	public int getj(double[] array, double p){
		
		int j=0;
		if(array.length==1) return j;
		while(array[j+1]<p){j++;}
		return j;
	}
	
	public double shapeFuncAt(double z, double p){
		
		double f=0;
		int jz=-1;
		while(this.zk[jz+1]>z){jz++;}
		
		int jp=-1;
		while(this.pk[jp+1]>z){jp++;}
		
		
		double fQ11=shapeFunc(jz,jp);
		double fQ12=shapeFunc(jz+1,jp);
		double fQ21=shapeFunc(jz,jp+1);
		double fQ22=shapeFunc(jz+1,jp+1);
		
		f+=fQ11*(this.pk[jp+1]-p)*(this.zk[jz+1]-z);
		
		f+=fQ21*(p-this.pk[jp])*(this.zk[jz+1]-z);
		
		f+=fQ12*(this.pk[jp+1]-p)*(z-this.zk[jz]);
		
		f+=fQ22*(p-this.pk[jp])*(z-this.zk[jz]);
		
		f*=1.0/(this.pk[jp+1]-this.pk[jp])*(this.zk[jz+1]-this.zk[jz]);
		

		return f;
		
	}
	
private double shapeFunc(int jz, int jp) {
		// TODO Auto-generated method stub
		return 0;
	}

public void simulateData(){


		//ps.magnetize(10);
		
		int L=100;
		int M=this.nHyst;
		Vect B=new Vect(2*L);
		Vect H=new Vect(2*L);
		Vect[] pk=new Vect[M];
		Vect[] sk=new Vect[M];
		double pkp=0;
		double skp=0;
		double Bp=0;
		double Xs=1.5;
		double pkps=0;

	//	this=new double[M];
		//double[] ita=new double[M];
		
		for(int k=0;k<M;k++){
			pk[k]=new Vect(2*L+1);
			//zeta[k]=1*k;
		}
	
		
		for(int i=-L;i<L;i++){
			int n=i+L;
	
			double x=Math.exp(-.01*n)*Math.sin(2*Math.PI*n*1.0/L);
			//double x=n*.2;
			//x=this.pk[i];
			
			B.el[n]=x;
			//double pkp=-1000;
			for(int k=0;k<this.nHyst;k++){
				
			if(n>0) pkp=pk[k].el[n-1];

			pk[k].el[n]=hystron1D(x,this.zk[k],pkp);
		//	pk[k].el[n]=pk(x,this.zk[k],pkp);
		
			H.el[n]+=this.shapeFuncAt(k,pk[k].el[n]);
			//H.el[n]+=this.func(pk[k].el[n]);
	
		}
		}


		util.plot(B,H);
	



	
	}
	
	public void distillData(){
		
		BH=new Mat[this.nTotCurves];
		
		Vect divBinit=new Vect().linspace(0, Bs, this.nTotCurves);
		Vect divB=new Vect().linspace(Bs, -Bs, 2*this.nTotCurves);
		//Vect[] divB=new Vect[this.nTotCurves];
		
		int L0=divB.length;
		
	/*	for(int i=0;i<divB.length;i++){
			divB=new Vect(2*this.nTotCurves+1-2*i);
		}
		
		=new Vect().linspace(Bs, -Bs, 2*this.nTotCurves);*/
		int L=0;
		
		int ix=0;
		
		for(int i=0;i<this.nInit;i++){
			
			L=divBinit.length;
			
			BH[ix]=new Mat(L,2);
		
			for(int j=0;j<L;j++){

					BH[ix].el[j]=this.gethbAsc(ix,divBinit.el[j]);
				
		
			}
			
			ix++;
		}
		

		for(int i=0;i<this.nMajor;i++){
			
			L=divB.length;
			
			BH[ix]=new Mat(L,2);
		
			for(int j=0;j<L;j++){
			

				BH[ix].el[j]=this.gethbDesc(ix,divB.el[j]);
				//util.hshow(BH[ix].el[j]);
			}
			
			ix++;
		}
		

		for(int i=0;i<this.nSym;i++){
			
			L=divB.length-2*i;
			
			BH[ix]=new Mat(L,2);
		
			for(int j=i;j<L0-i;j++){

				BH[ix].el[j-i]=this.gethbDesc(ix,divB.el[j]);
			
			}
			
			ix++;
		}
	
		
	
	
	}
	
	public double[] gethbDesc(int i,double B ){
		double[] hb=new double[2];
			int L=this.BHraw[i].nRow;
	
			
			if(B>=this.BHraw[i].el[0][1]){
				hb=BHraw[i].el[0]; 
			
			}
			
			else if(B<=this.BHraw[i].el[L-1][1]){
				hb=BHraw[i].el[L-1]; 
			}
			else{	
			int j=-1;
			while(this.BHraw[i].el[j+1][1]>B){j++;}
			double m1=(B-this.BHraw[i].el[j][1])/(this.BHraw[i].el[j+1][1]-this.BHraw[i].el[j][1]);
			double dH=this.BHraw[i].el[j+1][0]-this.BHraw[i].el[j][0];
			
			double H=m1*dH+this.BHraw[i].el[j][0];
			hb[0]=H; 
			hb[1]=B; 
			
			
			}
		
			return hb;
	}
	
	public double[] gethbAsc(int i,double B ){
		double[] hb=new double[2];
			int L=this.BHraw[i].nRow;
	
			if(B<=this.BHraw[i].el[0][1]){
				hb=BHraw[i].el[0]; 
			
			}
			
			else if(B>=this.BHraw[i].el[L-1][1]){
				hb=BHraw[i].el[L-1]; 
			}
			else{	
			int j=-1;
			while(this.BHraw[i].el[j+1][1]<B){j++;}
			double m1=(B-this.BHraw[i].el[j][1])/(this.BHraw[i].el[j+1][1]-this.BHraw[i].el[j][1]);
			double dH=this.BHraw[i].el[j+1][0]-this.BHraw[i].el[j][0];
			
			double H=m1*dH+this.BHraw[i].el[j][0];
			hb[0]=H; 
			hb[1]=B; 
			
			
			}
		
			return hb;
	}
	

	public static double func(double x)
	{
		util.pr(x);
	double y=Math.cbrt(x);

		return y;

	}
	
	public static double func2(double x)
	{
		
		double y=Math.pow(x, 3);

	
	return y;
		

	}
	
	public static double pk(double H,double zeta,double pkp)
	{
		
	double y=0;
	
	y=Math.max(Math.min(pkp,H+zeta),H-zeta);
	
	return y;

	}
	
	public static double hystron1D(double B,double zeta,double pp)
	{
		
	double y=0;
	
	if(zeta==0) y=B;
	else
	y=B-(B-pp)/Math.max(Math.abs(B-pp)/zeta,1);
	
	return y;

	}
	
	public static double sk(double skp,double B,double B0,double ita)
	{
		
	double y=0;
	
	y=Math.max(Math.min(B-B0+skp,ita),-ita);
	
	return y;
		

	}
}
