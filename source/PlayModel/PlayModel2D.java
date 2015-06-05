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

		//String file=System.getProperty("user.dir") + "\\hys_dataH.txt";
		String file="C:\\Works\\HVID\\hys_dataH";

		pm.createData(file);

		try {
			Thread.sleep(200);
		} catch (InterruptedException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}

		//pm.loadData(file);


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


		double Bs=1.6;
		double Hs=4500;

		int nInit=1;
		int nMajor=1;
		int nSymLoops=3;
		int nDescending=0;
		int nAscending=0;
		int nTot=nInit+nMajor+nSymLoops+nDescending+nAscending;

		int Mp=2000;
		

		double Hs0=1200;
		double mean=.2*Hs0;
		double width=.3*Hs0;

		Preisach2D ps=new Preisach2D(Mp,mean,width,Hs,Bs,3564656);

		ps.getRes().hshow();

		//ps.demagnetize(10);
		
int LL=2000;
		int nAni=0;
		int Lani=0;

		Mat[] BHs=new Mat[nTot];
		Mat[] BHsr=new Mat[nTot];

		int ix=0;

		for(int i=0;i<nInit;i++){

			BHs[ix]=ps.initial(LL);

		
			
			double b0=1e10;
			for(int j=0;j<BHs[ix].nRow;j++)
				if(BHs[ix].el[j][1]<b0)
					b0=BHs[ix].el[j][1];
	
			
		for(int j=0;j<BHs[ix].nRow;j++)
				BHs[ix].el[j][1]+=-b0;
		
		BHsr[ix]=this.getHrBr(BHs[ix]);
		
		/*	for(int j=0;j<BHs[ix].nRow;j++)
				BHs[ix].el[j][1]=Math.abs(BHs[ix].el[j][1]);*/
			
			//BHs[ix++]=ps.symMajorFull(200);
			
			ix++;

		}
		
	
		


		for(int i=0;i<nMajor;i++){
			//BHs[ix++]=ps.symMajorFull(500);
			
		//	ps.getRes().hshow();
		BHs[ix]=ps.symMajorDesc(LL);
	/*		BHsr[ix]=new Mat(BHs[ix].nRow,2);
			BHsr[ix].setCol(BHs[ix].getColVect(0), 0);
			BHsr[ix].setCol(BHs[ix].getColVect(1), 1);*/
	
			BHsr[ix]=this.getHBij(BHs[ix],0);
			//BHsr[ix].transp().show();
			
			ix++;

		}



		for(int i=0;i<nSymLoops;i++){
			double Bp=Bs-(i+1)*Bs/(nSymLoops+1);
			//BHs[ix++]=ps.symFull(Bp,500);
			BHs[ix]=ps.symDesc(Bp,LL);
			//util.pr(ps.getRes());
			BHsr[ix]=this.getHBij(BHs[ix],0);
			//BHs
		//	BHsr[ix]=this.getHBij(BHs[ix],0);BHsr[ix].transp().show();
			
			ix++;

		}


		for(int i=0;i<nDescending;i++){

			double Bp=Bs*(1-2.0*(i+1)/(nDescending+1));

		}

		for(int i=0;i<nAscending;i++){

			double Bp=-Bs*(1-2.0*(i+1)/(nAscending+1));
	
		}

			
			util.plotBunch(BHsr);
			

		Hs=800;
		

		
		Vect B=new Vect(Lani);
		
		for(int i=0;i<Lani;i++)
			B.el[i]=i*.1;
		
		Mat[] BHani=new Mat[nAni];
		Mat[] BHaniT=new Mat[1000];
		for(int i=0;i<nAni;i++){
			ps.kk=i;


			ps.demagnetize();
			//util.pr(ps.getRes());
			Mat BH2=ps.initial(1000);
			
			BHaniT[i]=new Mat(BH2.nRow,2);
			BHaniT[i].setCol(BH2.getColVect(0), 0);
			BHaniT[i].setCol(BH2.getColVect(1), 1);

			BHani[i]=new Mat(Lani,2);
		
		
			for(int j=0;j<Lani;j++){
				double Ht=this.getH(BH2, B.el[j]);
				BHaniT[i].el[j][0]=Ht;
				BHaniT[i].el[j][1]=B.el[j];
				BHani[i].el[j][0]=Ht*Math.cos(i*Math.PI/18);
				BHani[i].el[j][1]=Ht*Math.sin(i*Math.PI/18);
			//	BHani[i].el[j][1]=B.el[j];
			}

			
		}


//	util.plotBunch(BHaniT,5);
	//util.plotBunch(BHs);





	DecimalFormat dfB=new DecimalFormat("#.00");
	DecimalFormat dfH=new DecimalFormat("#.0");


		try{
			PrintWriter pwBun = new PrintWriter(new BufferedWriter(new FileWriter(file)));		

			pwBun.println(1+"\t"+1);
			pwBun.println("*Bs*Hs*");
			pwBun.println(Bs+"\t"+Hs);

			pwBun.println("* �������Ȑ��� * ���W���[���[�v�� * �Ώ̃��[�v�� * ���~�Ȑ��� * �㏸�Ȑ��� *");
			pwBun.println(nInit+"\t"+nMajor+"\t"+nSymLoops+"\t"+nDescending+"\t"+nAscending);

			for(int i=0;i<nTot;i++){
				pwBun.println("*xxx");
				pwBun.println(BHs[i].nRow);
				for(int j=0;j<BHs[i].nRow;j++)
					pwBun.println(BHs[i].el[j][0]+"\t"+BHs[i].el[j][1]);
			}

			pwBun.println("* ----- ��]�q�X�e���V�X��");
			pwBun.println("* B�� *");
			pwBun.println("0");
			pwBun.println("* B * ����");
			pwBun.println("* ----- �ٕ���");
			pwBun.println("* B�� * �p�x�� *");
			pwBun.println(Lani+"\t"+nAni);
			pwBun.println("* B * H ����� *�@�����e�Վ�");
			
			for(int i=0;i<Lani;i++){
				pwBun.print(dfB.format(B.el[i])+"\t");
				for(int j=0;j<nAni;j++){
					pwBun.print(dfH.format(BHani[j].el[i][0])+"\t");
				}
				pwBun.println();
			}
			pwBun.println("* B * H ����� *�@�������");
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
