package PlayModel;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.text.DecimalFormat;
import java.util.Random;

import math.Mat;
import math.Vect;
import math.util;

public class PlayModel2D {

	public Mat shapeFunc;
	public int nSym, nDesc,nAsc, nInit, nMajor,nTotCurves;
	public boolean anisot;
	public double Bs,Hs;
	public Mat[] BHraw,BH;
	public double[] zk,pk;
	public int dim=2,nHyst;



	public PlayModel2D()
	{	}

	public static void main(String[] args)
	{

		PlayModel2D pm=new PlayModel2D();

		String file=System.getProperty("user.dir") + "\\hys_dataH.txt";
		//String file="C:\\Works\\HVID\\hys_dataH";

		//pm.createData(file);

		pm.rotation();
		
		
		try {
			Thread.sleep(200);
		} catch (InterruptedException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}

		//pm.loadData(file);


	}
	
	public void rotation(){


		double Bs=1.8;
		double Hs=1800;

		int nInit=0;
		int nMajor=0;
		int nSymLoops=0;
		int nDescending=0;
		int nAscending=0;
		int nAni=1;
		int Lani=18;
		
		int nTot=nInit+nMajor+nSymLoops+nDescending+nAscending;

		int Mp=2000;
		

		double Hs0=1200;
		double mean=.2*Hs0;
		double width=.3*Hs0;

		Preisach2D ps=new Preisach2D(Mp,mean,width,Hs,Bs,3564656);
		
	//	Mat R=ps.getLocusHRotation(1.3);
		
		Mat R=ps.getLocusBRotation(1000,100,4);
		
		//R.show();
		
		//Mat BxBy=this.getHBij(R, 0);
		
		util.plot(R.getColVect(2),R.getColVect(3));
		//util.plot(R.getColVect(0),R.getColVect(1));
	}


	public void loadData(String file){/*

		HystDataLoader hysLoader=new HystDataLoader();

		hysLoader.loadData(this,file);

		this.distillBHData();



		this.nHyst=2*(this.nMajor+this.nSym);
		double db=Bs/this.nHyst;

		pk=new double[this.nHyst];
		zk=new double[this.nHyst];

		for(int i=0;i<pk.length;i++){
			pk[i]=i*db;
		}

		for(int i=0;i<zk.length;i++){
			zk[i]=i*db;
		}

		doIdentification();





	*/}



	public void createData(String file){


		double Bs=1.8;
		double Hs=1800;

		int nInit=0;
		int nMajor=0;
		int nSymLoops=0;
		int nDescending=0;
		int nAscending=0;
		int nAni=1;
		int Lani=18;
		
		int nTot=nInit+nMajor+nSymLoops+nDescending+nAscending;

		int Mp=2000;
		

		double Hs0=1200;
		double mean=.2*Hs0;
		double width=.3*Hs0;

		Preisach2D ps=new Preisach2D(Mp,mean,width,Hs,Bs,3564656);

		int LL=1000;
	
		double Bseff=.98*Bs;

		Mat[] BHs=new Mat[nTot];
		Mat[] BHsr=new Mat[nTot];

		int ix=0;

		for(int i=0;i<nInit;i++){

			Mat BHtemp=ps.initial(LL);
				
			int L=200;
			Vect B0=new Vect().linspace(0,.1,4);
			Vect B1=new Vect().linspace(.125,Bseff,L-4);
			Vect B=B0.aug(B1);
			
			BHs[ix]=new Mat(L,3);
			for(int j=0;j<L;j++){
				BHs[ix].el[j][0]=getH(BHtemp,B.el[j]);
				BHs[ix].el[j][1]=B.el[j];

			}
		

		BHsr[ix]=this.getHrBr(BHs[ix]);
		


			ix++;

		}
		
	
		


		for(int i=0;i<nMajor;i++){

		Mat BHtemp=ps.symMajorDesc(LL);
		
		int L=400;
		Vect B=new Vect().linspace(Bseff,-Bseff,L);
		BHs[ix]=new Mat(L,3);
		for(int j=0;j<L;j++){
			BHs[ix].el[j][0]=getH(BHtemp,B.el[j]);
			BHs[ix].el[j][1]=B.el[j];

		}
		

		BHsr[ix]=this.getHBij(BHs[ix],0);
			
			
			ix++;

		}



		for(int i=0;i<nSymLoops;i++){
			
			double Bp=Bseff-(i+1)*Bseff/(nSymLoops+1);	

			Mat BHtemp=ps.symDesc(Bp,LL);
			
			int L=400;
			Vect B=new Vect().linspace(Bp,-Bp,L);
			BHs[ix]=new Mat(L,3);
			for(int j=0;j<L;j++){
				BHs[ix].el[j][0]=getH(BHtemp,B.el[j]);
				BHs[ix].el[j][1]=B.el[j];

			}
			

			BHsr[ix]=this.getHBij(BHs[ix],0);
				
				
				ix++;

			

		}


		for(int i=0;i<nDescending;i++){
			
		}

		for(int i=0;i<nAscending;i++){
			
		}

		if(nTot>0)
		util.plotBunch(BHsr);
			

		
		Vect B=new Vect().linspace(0, Bseff, Lani);
	
		Mat[] BHani=new Mat[nAni];
		Mat[] BHaniT=new Mat[nAni];
		
		
		for(int i=0;i<nAni;i++){

			double phiRad=i*10*Math.PI/180;
			phiRad=0*Math.PI/2;

			Vect Hr=new Vect().linspace(0,1.2*Hs,LL);
			
			int Lp=2000;
			Mat Hp=new Mat(Lp,2);
			for(int j=0;j<Lp;j++){
				Hp.el[j][0]=.8*Hs*Math.cos(4*j*Math.PI/Lp)*Math.cos(phiRad);
				Hp.el[j][1]=.8*Hs*Math.sin(4*j*Math.PI/Lp)*Math.cos(phiRad);
			}
			
						
			Mat H=new Mat(LL,2);
			for(int j=0;j<LL;j++){
				H.el[j][0]=Hr.el[j]*Math.cos(phiRad);
				H.el[j][1]=Hr.el[j]*Math.sin(phiRad);
			}
			
			
			for(int k=0;k<ps.M;k++)
				for(int j=0;j<ps.nphi;j++){
					ps.on[k][j]=false;
				}
			
			ps.demagnetize(5*Hs,100,20);
		/*	util.pr(M.nCol);

			util.plot(getHBij(M,0).el);
*/

		//	ps.getRes().hshow();
			
			Mat BH2=ps.getLocus(Hp);
	
			Mat BH3=this.getHrBr(BH2);
			
			//util.plot(getHBij(BH2,3).el);
			util.plot(BH2.getColVect(2),BH2.getColVect(3));
			
			BHaniT[i]=BH3.deepCopy();
				
			BHaniT[i]=new Mat(Lani,2);

			BHani[i]=new Mat(Lani,2);
		
		
			for(int j=0;j<Lani;j++){
				double Ht=this.getH(BH3, B.el[j]);
				BHaniT[i].el[j][0]=Ht;
				BHaniT[i].el[j][1]=B.el[j];
				BHani[i].el[j][0]=Ht*Math.cos(phiRad);
				BHani[i].el[j][1]=Ht*Math.sin(phiRad);
		
			}

		}
		
		//util.plot(BHaniT[0].el);

	//util.plotBunch(BHaniT);
	//util.plotBunch(BHs);

	boolean write =false;

if(write){


	DecimalFormat dfB=new DecimalFormat("#.00");
	DecimalFormat dfH=new DecimalFormat("#.0");


		try{
			PrintWriter pwBun = new PrintWriter(new BufferedWriter(new FileWriter(file)));		

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
			pwBun.println(Lani+"\t"+nAni);
			pwBun.println("* B * H ･････ *　磁化容易軸");
			
			for(int i=0;i<Lani;i++){
				pwBun.print(dfB.format(B.el[i])+"\t");
				for(int j=0;j<nAni;j++){
					pwBun.print(dfH.format(BHani[j].el[i][0])+"\t");
				}
				pwBun.println();
			}
			pwBun.println("* B * H ･････ *　磁化困難軸");
			for(int i=0;i<Lani;i++){
				pwBun.print(dfB.format(B.el[i])+"\t");
				for(int j=0;j<nAni;j++){
					pwBun.print(dfH.format(BHani[j].el[i][1])+"\t");
				}
				pwBun.println();
			}
				

			util.pr("Simulated hysteresis data was written to "+file+".");

			pwBun.close();
		}
		catch(IOException e){}

}



	}
	

	public void doIdentification(){
		this.shapeFunc=new Mat(nHyst+1,nHyst+1);

		int j = nHyst;



		for (int i=0; i<=nHyst; ++i) {
			this.shapeFunc.el[i][j - i] =this.BH[1].el[i][0];//.getMajorLoop()->H[i];

		}




		for (int k=0; k<this.nSym; ++k) {


			j--;
			for (int i=0; i<j-k-1; ++i) {

				this.shapeFunc.el[i][j-i] = this.BH[k+1].el[i][0];
			}
			if ((j-k-1) <= 0) {
				break;
			}

			this.shapeFunc.el[j-k-1][k+1] = this.BH[k+1].el[j-k-1][0];
		}

		this.shapeFunc.el[0][this.nSym+1] = 0.;

		j = nHyst;

		for (int k=0; k<this.nSym; ++k) {
			this.shapeFunc.el[0][j] -= this.shapeFunc.el[1][j-1];
			for (int i=1; i<this.nHyst-2*k-1; ++i) {
				this.shapeFunc.el[i][j-i]
						= this.shapeFunc.el[i][j-i] +this.shapeFunc.el[i][j-i-1]
								- this.shapeFunc.el[i-1][j-i] - this.shapeFunc.el[i+1][j-i-1];
			}
			this.shapeFunc.el[this.nHyst-2*k-1][j-this.nHyst+2*k+1]
					= 2. * (this.shapeFunc.el[this.nHyst-2*k-1][j-this.nHyst+2*k+1]
							- this.shapeFunc.el[this.nHyst-2*k-2][j-this.nHyst+2*k+1]);
			j--;
		}


		for (int i=0; i<this.nHyst-1; ++i) {
			for (int k=0; k<this.nSym+0.5-0.5*i; ++k) {
				this.shapeFunc.el[i][k+1] = this.shapeFunc.el[i][this.nHyst-i-k];
			}
		}


		for (int i=0; i<this.nHyst; ++i) {
			for (int k=2; k<this.nHyst+1-i; ++k) {
				this.shapeFunc.el[i][k] += this.shapeFunc.el[i][k-1];
			}
			this.shapeFunc.el[i][0] = 0.;
		}


		double center;
		for (int i=0; i<this.nHyst; ++i) {
			center = 0.5 * this.shapeFunc.el[i][this.nHyst-i];
			for (int k=0; k<this.nHyst+1-i; ++k) {
				this.shapeFunc.el[i][k] -= center;
			}
		}


		Vect pp=new Vect().linspace(-Bs, Bs, nHyst+1);
	//p.show();
		shapeFunc.show();

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
	
	public  double getH(Mat BH,double B){
		return getH(BH,B,0);
	}
	
	public  double getH(Mat BH,double B, int mode){

		int ih, ib;
		if(mode==0){
			ih=0;ib=1;
		}
		else{
			ih=0;ib=2;
		}
		
		boolean ascending=true;
		if (BH.el[1][ih]<BH.el[0][ih]) ascending=false;
		
		int i1=0;
		int i2=BH.nRow-1;
		
		if(ascending){

		if(B<=BH.el[i1][ib])
			return BH.el[i1][ih];
		if(B>=BH.el[i2][ib])
			return BH.el[i2][ih];
			/*return BH.el[i][ih]+(B-BH.el[i][ib])*(BH.el[i+1][ih]-BH.el[i][ih])/(BH.el[i+1][ib]-BH.el[i][ib]);
		i=BH.nRow-1;
		if(B>=BH.el[i][ib])
			return BH.el[i][ih]+(B-BH.el[i][ib])*(BH.el[i][ih]-BH.el[i-1][ih])/(BH.el[i][ib]-BH.el[i-1][ib]);
		*/
		}
		else{

			if(B<=BH.el[i2][ib])
				return BH.el[i2][ih];
			if(B>=BH.el[i1][ib])
				return BH.el[i1][ih];
			}

		int j=0;

		if(ascending)
			while(BH.el[j+1][ib]<B){j++;}
		else
			while(BH.el[j+1][ib]>B){j++;}
		
		double cc=(BH.el[j+1][ih]-BH.el[j][ih])/(BH.el[j+1][ib]-BH.el[j][ib]);

		double H=BH.el[j][ih]+(B-BH.el[j][ib])*cc;

		
		return H;
				
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

	public void distillBHData(){

		BH=new Mat[this.nTotCurves];

		int nDB=2*(this.nSym+nMajor)+1;

		Vect divBinit=new Vect().linspace(0, Bs, (this.nSym+nMajor)+1);
		Vect divB=new Vect().linspace(Bs, -Bs,nDB);



		int L0=divB.length;


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
	
	public Mat getHBij(Mat BH,int k){
		
		Mat BHr=new Mat(BH.nRow,2);
		
		if(BH.nCol==4){
			
			if(k==0)
			for(int i=0;i<BHr.nRow;i++){
			BHr.el[i][0]=BH.el[i][0];
			
			BHr.el[i][1]=BH.el[i][2];
			}
		else if(k==1)
			for(int i=0;i<BHr.nRow;i++){
			BHr.el[i][0]=BH.el[i][0];
			
			BHr.el[i][1]=BH.el[i][3];
			}
		else if(k==2)
			for(int i=0;i<BHr.nRow;i++){
			BHr.el[i][0]=BH.el[i][1];
		
			BHr.el[i][1]=BH.el[i][2];
		}
		else if(k==3)
			for(int i=0;i<BHr.nRow;i++){
			BHr.el[i][0]=BH.el[i][1];
		
			BHr.el[i][1]=BH.el[i][3];
		}
		else {
			throw new IllegalArgumentException("Matrix has 4 columns. No more than 4 combination possible.");
		}
}
	
		else if(BH.nCol==3){
			
			if(k==0)
				for(int i=0;i<BHr.nRow;i++){
			BHr.el[i][0]=BH.el[i][0];
	
			BHr.el[i][1]=BH.el[i][1];
			}
		else if(k==1)
			for(int i=0;i<BHr.nRow;i++){
			BHr.el[i][0]=BH.el[i][0];
			
			BHr.el[i][1]=BH.el[i][2];
			}
		else {
			throw new IllegalArgumentException("Matrix has 3 columns. No more than 2 combination possible.");
		}
}
		
		return BHr;

	}
	
	public Mat getHrBr(Mat BH){
		
		Mat BHr=new Mat(BH.nRow,2);
		
		for(int i=0;i<BHr.nRow;i++)
		{
			if(BH.nCol==4)
			{
			BHr.el[i][0]=new Vect(BH.el[i][0],BH.el[i][1]).norm();
			
			BHr.el[i][1]=new Vect(BH.el[i][2],BH.el[i][3]).norm();
			}
			else if(BH.nCol==3){
				
				BHr.el[i][0]=BH.el[i][0];
				
				BHr.el[i][1]=new Vect(BH.el[i][1],BH.el[i][2]).norm();
			}
		}

		
		return BHr;

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
