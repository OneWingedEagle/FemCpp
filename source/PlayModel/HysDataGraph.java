package PlayModel;

import static java.lang.Math.PI;
import static java.lang.Math.acos;
import static java.lang.Math.cos;
import static java.lang.Math.sqrt;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.text.DecimalFormat;
import java.util.Scanner;

import math.Mat;
import math.Vect;
import math.util;
import fem.Model;



public class HysDataGraph {

	Mat[][] BH;
	Mat rotLoss;
	Mat[][] BHAni;
	int nSet;
	double[] Bs, Hs;
	int[] nInitial,nMajor,nSymLoop,nDescending,nAscending,nTotCurves,nAni,Lani;
	String regex="[:; ,\\t]+";


	public static void main(String[] args){

		HysDataGraph pg=new HysDataGraph();

		//pg.loadHysData();
		
		//pg.loadRotHysData();
		pg.loadAngData();
		//pg.loadAngSymData();
		
	}

	public boolean loadHysData(){

		/*String file=util.getFile();
	if(file==null || file.equals("") )return false;*/
		//	String file="C:\\Works\\HVID\\folder1\\data\\A_B入力対称ループhts_data\\hys_data";

		//String file="C:\\Users\\Hassan Ebrahimi\\JavaWorks\\MagFem\\hys_data";
		String file="C:\\Works\\HVID\\hys_data";
	//	String file=System.getProperty("user.dir") + "\\hys_dataH.txt";
		


		try{
			FileReader fr=new FileReader(file);
			BufferedReader br = new BufferedReader(fr);
			String line;
			String s;
			String[] sp;

			
			line=br.readLine();
			
			sp=line.split(regex);
			if(sp.length<3) this.nSet=1;
			else this.nSet=Integer.parseInt(sp[2]);
			
			this.Bs=new double[this.nSet];
			this.Hs=new double[this.nSet];
			this.nInitial=new int[this.nSet];
			this.nMajor=new int[this.nSet];
			this.nSymLoop=new int[this.nSet];
			this.nDescending=new int[this.nSet];
			this.nAscending=new int[this.nSet];
			this.nTotCurves=new int[this.nSet];
			this.nAni=new int[this.nSet];
			this.Lani=new int[this.nSet];

			BH=new Mat[nSet][];
			BHAni=new Mat[nSet][];
	
		
		
			fr=new FileReader(file);
			br = new BufferedReader(fr);

			for(int ia=0;ia<this.nSet;ia++){

				line=br.readLine();
				line=br.readLine();
				line=br.readLine();
				sp=line.split(regex);	

				this.Bs[ia]=Double.parseDouble(sp[0]);
				this.Hs[ia]=Double.parseDouble(sp[1]);
				line=br.readLine();
				line=br.readLine();
				sp=line.split(regex);	
				this.nInitial[ia]=Integer.parseInt(sp[0]);
				this.nMajor[ia]=Integer.parseInt(sp[1]);
				this.nSymLoop[ia]=Integer.parseInt(sp[2]);
				this.nDescending[ia]=Integer.parseInt(sp[3]);
				this.nAscending[ia]=Integer.parseInt(sp[4]);

				this.nTotCurves[ia]=this.nInitial[ia]+this.nMajor[ia]+this.nSymLoop[ia]+this.nDescending[ia]+this.nAscending[ia];
				BH[ia]=new Mat[this.nTotCurves[ia]];

				line=br.readLine();
				line=br.readLine();
				int L=Integer.parseInt(line);

				BH[ia][0]=new Mat(L,2);

				for( int p=0;p<L;p++){
					line=br.readLine();
					sp=line.split(regex);	
					double[] bh=getCSV(line);
					BH[ia][0].el[p][0]=bh[0];
					BH[ia][0].el[p][1]=bh[1];

				}



				line=br.readLine();
				int L1=0;
				for( int ip=1;ip<this.nTotCurves[ia];ip++){
					line=br.readLine();
					if(line.startsWith("*")) {line=br.readLine();};
					sp=line.split(regex);	
					//if(sp.length==1)
					L1=Integer.parseInt(sp[0]);

					BH[ia][ip]=new Mat(L1,2);

					for( int i=0;i<L1;i++){
						line=br.readLine();

						double[] bh=getCSV(line);

						BH[ia][ip].el[i][0]=bh[0];
						BH[ia][ip].el[i][1]=bh[1];
					}


				}

				line=br.readLine();
				line=br.readLine();
				line=br.readLine();

				int nLoss=Integer.parseInt(line);
				rotLoss=new Mat(nLoss,2);
				line=br.readLine();

				for( int i=0;i<nLoss;i++){
					line=br.readLine();
					double[] BL=getCSV(line);
					rotLoss.el[i][0]=BL[0];
					rotLoss.el[i][1]=BL[1];
				}

				line=br.readLine();
				line=br.readLine();
				line=br.readLine();
				double[] dd=this.getCSV(line);
				
				Lani[ia]=(int)dd[0];
				nAni[ia]=(int)dd[1];
			
				line=br.readLine();

				BHAni[ia]=new Mat[nAni[ia]];

				for( int i=0;i<nAni[ia];i++)
					BHAni[ia][i]=new Mat(Lani[ia],3);

				for( int i=0;i<Lani[ia];i++){
					line=br.readLine();
			
					double[] BL=getCSV(line);
					double B=BL[0];
		
					for(int j=0;j<nAni[ia];j++){
						BHAni[ia][j].el[i][0]=B;
						BHAni[ia][j].el[i][1]=BL[j+1];
					}
				}

				line=br.readLine();
				for( int i=0;i<Lani[ia];i++){
					line=br.readLine();
					double[] BL=getCSV(line);
			
					for(int j=0;j<nAni[ia];j++){
						BHAni[ia][j].el[i][2]=BL[j+1];
					}
				}


			}
			
			//util.plotBunch(BH[0]);
			//util.plotBunch(BH[9],3);
			
		createAngleDepData();

			
			boolean ani=false;
			
			if(ani){
				
				Mat[] BHaniT=new Mat[BHAni[0].length];
				
	
				for(int i=0;i<BHAni[0].length;i++){
					
					BHaniT[i]=new Mat(BHAni[0][i].nRow,2);
					
					double ang=i*Math.PI/18;
					
					Vect er=new Vect(Math.cos(ang),Math.sin(ang));

				for(int j=0;j<BHAni[0][i].nRow;j++){
					double Ht=new Vect(BHAni[0][i].el[j][1], BHAni[0][i].el[j][2]).dot(er);
					BHaniT[i].el[j][0]=Ht;
					BHaniT[i].el[j][1]=BHAni[0][i].el[j][0];
					
				}
			
				}
				
			
						
			
		//	util.plotBunch(BHaniT,0);

			}
	
			boolean distill=false;
			if(distill){
			Mat[][] BHdist=distHysData(BH);
			
			
/*			for(int i=0;i<BHdist.length;i++)
				for(int j=0;j<BHdist[i].length;j++)
					for(int k=0;k<BHdist[i][j].nRow;k++)
						BHdist[i][j].el[k][1]*=1+Math.abs(k-BHdist[i][j].nRow/2)*.005;*/
						
			//util.plotBunch(BHdist[0]);
			
		
			if(BHdist.length==1){
				String filed="C:\\Works\\HVID\\hys_dataHAvdist";
			this.writeHystDataAv(BHdist[0], BHAni[0], filed);
			}
			else{
				String filed="C:\\Works\\HVID\\hys_dataHdist";
			this.writeHystData(BHdist, filed);
			}
			}
		
			return true;

		}

		catch(IOException e){System.err.println("Error in loading BH data file.");
		return false;
		}

	}	
	
	
	public void createAngleDepData(){
		
		int nSet=18;
		int nTot=16;
		
		Mat[][] BHs=new Mat[nSet][nTot];
		
		
		for(int ia=0;ia<nSet;ia++)
			for(int i=0;i<BHs[ia].length;i++){
			BHs[ia][i]=this.BH[0][i].deepCopy();
			for(int j=0;j<BHs[ia][i].nRow;j++){
				//BHs[ia][i].el[j][0]*=(1+1*Math.sin(ia*Math.PI/nSet));
			//if(i>0)
			BHs[ia][i].el[j][0]*=1+2.5*Math.sin(ia*Math.PI/nSet);//*(BHs[ia][i].el[j][0]-BHs[ia][i].el[0][0])*(BHs[ia][i].el[j][0]-BHs[ia][i].el[BHs[ia][i].nRow-1][0])/Math.pow(BHs[ia][i].el[0][0]-BHs[ia][i].el[BHs[ia][i].nRow-1][0], 2);
			//BHs[ia][i].el[j][1]*=1+.5*Math.sin(ia*Math.PI/nSet)*(BHs[ia][i].el[j][1]-BHs[ia][i].el[0][1])*(BHs[ia][i].el[j][1]-BHs[ia][i].el[BHs[ia][i].nRow-1][1])/Math.pow(BHs[ia][i].el[0][1]-BHs[ia][i].el[BHs[ia][i].nRow-1][1], 2);
			}
			
			}
		
		util.plotBunch(BHs[0],16);
		util.plotBunch(BHs[9],16);

		
		String file="C:\\Works\\HVID\\hys_dataGen";
		
		writeHystData(BHs,  file);
		
	}
	
	public void loadRotHysData(){


		String file="C:\\Works\\HVID\\hysRotation";

		
	 HystDataLoader loader=new HystDataLoader();
	 

	Mat BHij=new Mat(loader.loadArrays(1024,69,file));
		


		
		Mat[] BB=new Mat[17];
		Mat[] HH=new Mat[17];
		
		for(int i=0;i<BB.length;i++){
			BB[i]=new Mat(BHij.nRow,2);
			HH[i]=new Mat(BHij.nRow,2);
			
			BB[i].setCol(BHij.getColVect(4*i+1), 0);
			BB[i].setCol(BHij.getColVect(4*i+3), 1);
			
			HH[i].setCol(BHij.getColVect(4*i+2), 0);
			HH[i].setCol(BHij.getColVect(4*i+4), 1);
		}
	
		
	util.plotBunch(HH);
		
		//util.plot(HH[9]);

		


	}
	
	public void loadAngSymData(){
		
		int deg=0;
		
		 HystDataLoader loader=new HystDataLoader();
		
		 
		 
		 int nset=7;
	
		 
		 Mat[][] hys=new Mat[nset][];
	
		 
		 for(int ia=0;ia<nset;ia++){
		 String file="C:\\Works\\HVID\\KitaoData\\symmetricData\\sym"+ia*15+".dat";
	;

		 Mat[][] syms=loader.loadDataSym(file);
		 
		 hys[ia]=new Mat[syms[1].length+1];
				 
		 Mat init=new Mat(syms[1].length+1,2);
		 for(int i=1;i<init.nRow;i++){
			 init.el[i][0]=syms[1][i-1].el[0][0];
			 init.el[i][1]=syms[1][i-1].el[0][1];
		 }
		 
		 hys[ia][0]=init;
		 
		 for(int i=1;i<hys[ia].length;i++)
		 hys[ia][hys[ia].length-i]=syms[1][i-1];
		 
		 }
		 
		//util.plotBunch(hys[0]);
		
		PlayModel2D pm=new PlayModel2D();

		
		Mat[][] BHs=new Mat[nset][];

		 int L=50;
		 
		 for(int ia=0;ia<nset;ia++){
			 int nc=hys[ia].length;
			 BHs[ia]=new Mat[nc];
			 
		 for(int i=0;i<nc;i++)
		 {
			 BHs[ia][i]=new Mat(L,2);
			 Vect B=new Vect().linspace(hys[ia][i].el[0][1], hys[ia][i].el[hys[ia][i].nRow-1][1], L);
	
			 for(int j=0;j<L;j++){
				 BHs[ia][i].el[j][1]=B.el[j];
				 BHs[ia][i].el[j][0]=pm.getH(hys[ia][i], B.el[j]);
			 }

		 }
		 }
		
		
		 util.plotBunch(BHs[0]);
		 util.plotBunch(BHs[3]);
		 
		// util.plotBunch(BHs[0]);
		 
		//String fileout="C:\\Works\\HVID\\KitaoData\\symmetricData\\hys_data"+deg;
		String fileout="C:\\Works\\HVID\\KitaoData\\symmetricData\\hys_dataAll";
		pm.writeHystData(BHs, fileout);
		
}
	
	
	
	public void loadAngData(){

		int deg=60;
		
		Vect er=new Vect(Math.cos(deg*180.0/Math.PI),Math.sin(deg*180/Math.PI));

		//String file="C:\\Works\\HVID\\KitaoData\\deg45";

		String file="C:\\Works\\HVID\\KitaoData\\deg"+deg;


	 HystDataLoader loader=new HystDataLoader();
	 

	Mat BHij=new Mat(loader.loadArrays(1024,72,file));
	
		int nc=72/4;
		
		double err=1e-4;
		
		Mat[] BB=new Mat[nc];
		Mat[] HH=new Mat[nc];
		
		Mat[] BH=new Mat[nc];
		Mat[] BHs=new Mat[nc];
		
		for(int i=0;i<BB.length;i++){
			BB[i]=new Mat(BHij.nRow,2);
			HH[i]=new Mat(BHij.nRow,2);
			
			BB[i].setCol(BHij.getColVect(4*i+0), 0);
			BB[i].setCol(BHij.getColVect(4*i+1), 1);
			
			HH[i].setCol(BHij.getColVect(4*i+2), 0);
			HH[i].setCol(BHij.getColVect(4*i+3), 1);
			
			int L=BB[i].nRow;
			
			BH[i]=new Mat(L,2);
			
			for(int j=0;j<L;j++){
				
				BH[i].el[j][0]=new Vect(HH[i].el[j]).dot(er);
				BH[i].el[j][1]=new Vect(BB[i].el[j]).dot(er);
			}
			

			double Bmax=BH[i].getColVect(1).max();
			double Bmin=BH[i].getColVect(1).min();
			
			int n1=0;
			int n2=0;
			int jx=0;

			
			while(BH[i].el[jx][1]<-err+Bmax /*|| BH[i].el[jx+1][0]>=BH[i].el[jx][0]*/){jx++;}
			n1=jx;

			while(BH[i].el[jx][1]>err+Bmin /*|| BH[i].el[jx+1][0]<=BH[i].el[jx][0]*/){jx++;}
			
			n2=jx;
			
			int Ls=n2-n1;
			
			BHs[i]=new Mat(Ls,2);
			
			for(int j=0;j<Ls;j++){
				
				BHs[i].el[j]=BH[i].el[j+n1];
			}
			
			
		}
	

		
	//util.plotBunch(HH,1);
	util.plotBunch(BH,10);
/*	BHs[0].show();	
	BH[0].show();*/
		//util.plot(HH[9]);

		


	}
	

	
	
	Mat[][] distHysData(Mat[][] BH){
		
		Mat[][] distBH=new Mat[BH.length][];
		
		PlayModel2D pm=new PlayModel2D();
		for(int k=0;k<distBH.length;k++){
			
			distBH[k]=new Mat[BH[k].length];
		
			
		for(int i=0;i<distBH[k].length;i++){
			int Lx=50;
			if (i==0) Lx=Lx/2;
			distBH[k][i]=new Mat(Lx,2);
			Vect B=new Vect();
			if(i==0){
				int jx=0;
				while(jx<BH[k][i].nRow-1 &&BH[k][i].el[jx+1][0]>=BH[k][i].el[jx][0]){jx++;};
		
				B=new Vect().linspace(0, BH[k][i].el[jx][1], Lx);
				
				if(B.el[Lx-1]<BH[k][i].el[BH[k][i].nRow-1][1]) B.el[Lx-1]=BH[k][i].el[BH[k][i].nRow-1][1];
			}
			else if(i==1){
				
				int jx1=0;
				while(jx1<BH[k][i].nRow-1&& BH[k][i].el[jx1+1][0]>=BH[k][i].el[jx1][0]){jx1++;};
				
				int r=BH[k][i].nRow;
				int jx2=r-1;
				while(jx2>0 &&BH[k][i].el[jx2-1][0]<=BH[k][i].el[jx2][0]){jx2--;};

				B=new Vect().linspace(BH[k][i].el[jx1][1], BH[k][i].el[jx2][1], Lx);
				
				if(B.el[Lx-1]>-BH[k][0].el[BH[k][0].nRow-1][1]) B.el[Lx-1]=-BH[k][0].el[BH[k][0].nRow-1][1];
				if(B.el[0]<BH[k][0].el[BH[k][0].nRow-1][1]) B.el[0]=BH[k][0].el[BH[k][0].nRow-1][1];
			}
			else{
				B=new Vect().linspace(BH[k][i].el[0][1], BH[k][i].el[BH[k][i].nRow-1][1], Lx);
			}

			for(int j=0;j<Lx;j++)
			{
				double H=pm.getH(BH[k][i], B.el[j]);
				distBH[k][i].el[j][0]=H;
			
				distBH[k][i].el[j][1]=B.el[j];
			}
			
			if(i==0) distBH[k][i].el[0][0]=0;
			
		}
		}
		

		return distBH;

		
	}
	
	public void writeHystDataAv(Mat[] BHs,Mat[] BHani, String file){
		
		int nSet=1;
		
		double Bseff= BHs[0].el[BHs[0].nRow-1][1];
		
		int nAni=BHani.length;
		int Lani=BHani[0].nRow;
		int nInit=1;
		int nMajor=1;
		int nTot=BHs.length;
		int nSymLoops=nTot-2;
		int nDescending=0;
		int nAscending=0;
		
		double Hseff=BHs[0].el[BHs[0].nRow-1][0];

		DecimalFormat dfB=new DecimalFormat("#.00");
		DecimalFormat dfH=new DecimalFormat("#.0");


			try{
				PrintWriter pwBun = new PrintWriter(new BufferedWriter(new FileWriter(file)));		


				pwBun.println(1+"\t"+1+"\t"+nSet+"\t"+0);
				pwBun.println("*Bs*Hs*");
				pwBun.println(Bseff+"\t"+Hseff);

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
					pwBun.print(dfB.format(BHani[0].el[i][0])+"\t");
					for(int j=0;j<nAni;j++){
						pwBun.print(dfH.format(BHani[j].el[i][1])+"\t");
					}
					pwBun.println();
				}
				pwBun.println("* B * H ･････ *　磁化困難軸");
				for(int i=0;i<Lani;i++){
					pwBun.print(dfB.format(BHani[0].el[i][0])+"\t");
					for(int j=0;j<nAni;j++){
						pwBun.print(dfH.format(BHani[j].el[i][1])+"\t");
					}
					pwBun.println();
				}
				BHani[0].show();

				util.pr("Simulated angle-dependent hysteresis data was written to "+file+".");

				pwBun.close();
			}
			catch(IOException e){}
			

	}

	public void writeHystData(Mat[][] BHs, String file){
		
		int nSet=BHs.length;
		
		double Bseff= BHs[0][0].el[BHs[0][0].nRow-1][1];
		
		Mat[] BHani=new Mat[1];
		int nAni=0;
		int Lani=0;
		int nInit=1;
		int nMajor=1;
		int nTot=BHs[0].length;
		int nSymLoops=nTot-2;
		int nDescending=0;
		int nAscending=0;
		
	

		DecimalFormat dfB=new DecimalFormat("#.00");
		DecimalFormat dfH=new DecimalFormat("#.0");


			try{
				PrintWriter pwBun = new PrintWriter(new BufferedWriter(new FileWriter(file)));		

	for(int ia=0;ia<nSet;ia++){

		double Hsefft=BHs[ia][0].el[BHs[ia][0].nRow-1][0];
		double Bsefft=BHs[ia][0].el[BHs[ia][0].nRow-1][1];
				pwBun.println(1+"\t"+1+"\t"+nSet+"\t"+ia*10);
				pwBun.println("*Bs*Hs*");
				pwBun.println(Bsefft+"\t"+Hsefft);

				pwBun.println("* 初磁化曲線数 * メジャーループ数 * 対称ループ数 * 下降曲線数 * 上昇曲線数 *");
				pwBun.println(nInit+"\t"+nMajor+"\t"+nSymLoops+"\t"+nDescending+"\t"+nAscending);

				for(int i=0;i<nTot;i++){
					pwBun.println("*xxx");
					pwBun.println(BHs[ia][i].nRow);
					for(int j=0;j<BHs[ia][i].nRow;j++)
						pwBun.println(BHs[ia][i].el[j][0]+"\t"+BHs[ia][i].el[j][1]);
				}

				pwBun.println("* ----- 回転ヒステリシス損");
				pwBun.println("* B数 *");
				pwBun.println("0");
				pwBun.println("* B * 損失");
				pwBun.println("* ----- 異方性");
				pwBun.println("* B数 * 角度数 *");
				pwBun.println(0+"\t"+0); 		//	pwBun.println(Lani+"\t"+nAni);
				pwBun.println("* B * H ･････ *　磁化容易軸");
				
				for(int i=0;i<0*Lani;i++){
					pwBun.print(dfB.format(BHani[0].el[i][0])+"\t");
					for(int j=0;j<nAni;j++){
						pwBun.print(dfH.format(BHani[j].el[i][0])+"\t");
					}
					pwBun.println();
				}
				pwBun.println("* B * H ･････ *　磁化困難軸");
				for(int i=0;i<0*Lani;i++){
					pwBun.print(dfB.format(BHani[0].el[i][0])+"\t");
					for(int j=0;j<nAni;j++){
						pwBun.print(dfH.format(BHani[j].el[i][1])+"\t");
					}
					pwBun.println();
				}
	/*			pwBun.println();
				pwBun.println("End of hysteresis data set "+ia);
				pwBun.println();*/
	}
					

				util.pr("Simulated angle-dependent hysteresis data was written to "+file+".");

				pwBun.close();
			}
			catch(IOException e){}
			

	}


	private double[] getCSV(String line){

		String[] sp=line.split(regex);	

		int p0=0;
		if(sp[0].equals(""))
		{
			p0=1;
		}
		int L=sp.length-p0;

		double[] v=new double[L];

		for( int p=0;p<L;p++){

			v[p]=Double.parseDouble(sp[p+p0]);
		}

		return v;
	}


}
