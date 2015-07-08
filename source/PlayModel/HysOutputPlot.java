package PlayModel;

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



public class HysOutputPlot {
	
Mat BHs[];

String regex="[:; ,\\t]+";


public static void main(String[] args){

	HysOutputPlot pg=new HysOutputPlot();
	
	pg.loadData();
	
}
	
public  HysOutputPlot(){}

public void loadData(){

	/*String file=util.getFile();
	if(file==null || file.equals("") )return false;*/
//	String file="C:\\Works\\HVID\\folder1\\data\\A_B入力対称ループhts_data\\hys_data";

	//String file="C:\\Users\\Hassan Ebrahimi\\JavaWorks\\MagFem\\hys_data";
	String file="C:\\Works\\HVID\\output";
	
	

		try{
			FileReader fr=new FileReader(file);
			BufferedReader br = new BufferedReader(fr);
			String line;
			String s;
			String[] sp;

int Nmax=20;

Mat[] BB=new Mat[Nmax];
Mat[] HH=new Mat[Nmax];

BHs=new Mat[Nmax];


Mat[] XX=new Mat[2*Nmax];

int numbCurves=0;
for(int i=0;i<Nmax;i++){
	

		line=br.readLine();
		
		if(line==null) break;
		
			line=br.readLine();
			line=br.readLine();

int L=Integer.parseInt(line);

numbCurves++;

Mat bbhh=new Mat(L,4);
	
	for( int p=0;p<L;p++)
	{
		line=br.readLine();
		double[] ar=this.getCSV(line);
		bbhh.el[p]=ar;
	}
		
			
			BB[i]=new Mat(L,2);
			HH[i]=new Mat(L,2);
			
		
		BB[i].setCol(bbhh.getColVect(0).times(1), 0);
		BB[i].setCol(bbhh.getColVect(1).times(1), 1);
			
			HH[i].setCol(bbhh.getColVect(2), 0);
			HH[i].setCol(bbhh.getColVect(3), 1);
	

	
			XX[2*i]=BB[i].times(100);
			XX[2*i+1]=HH[i];
			
			BHs[i]=new Mat(L,2);
			
		//	HH.show();
			
				Vect Hr=new Vect(HH[i].nRow);
			Vect Br=new Vect(BB[i].nRow);
			
		
				
		double ang=Math.atan(BB[i].el[1][1]/BB[i].el[1][0]);
		//util.pr(ang/Math.PI*180);
		//util.pr(ang/Math.PI*180);
		
		Vect er=new Vect(Math.cos(ang),Math.sin(ang));
		for(int j=0;j<Hr.length;j++){
			Hr.el[j]=new Vect(HH[i].el[j][0],HH[i].el[j][1]).dot(er);
			Br.el[j]=new Vect(BB[i].el[j][0],BB[i].el[j][1]).dot(er);
			//util.pr(Hr.el[i]+"\t"+Br.el[i]);
		}

		BHs[i].setCol(Hr,0);
		BHs[i].setCol(Br,1);

		
		line=br.readLine();

}
		
			//BH[1].show();
//util.plotBunch(BHs,1);
//BHs[0].show();
			br.close();
			fr.close();
	
			//util.plotBunch(XX,2);


			int L=HH[0].nRow/2;
			Mat M=new Mat(L,2);
			for(int i=0;i<L;i++){
				M.el[i]=HH[0].el[i+L];
			}
			
			//M.show();
		//	util.plot(M);
			util.plotBunch(BHs,1);
		//	BHs[0].show();
			String file1="C:\\Works\\HVID\\b_times";
			
			if(2>5)
			try{
				PrintWriter pwBun = new PrintWriter(new BufferedWriter(new FileWriter(file1)));		

				int K=3;
				double dB=.1;
				Vect div=	new Vect().linspace(0,1.7,K);
				//util.pr(div.length);
			int nAng=18;
			
			pwBun.println(1);
			pwBun.println((K-2)*nAng);
			
		for(int k=0;k<nAng;k++){
			
			double rad=k*Math.PI/18;
			double cos=Math.cos(rad);
			double sin=Math.sin(rad);
			
				for(int i=2;i<div.length;i++){
					Vect v=new Vect(7);
					v.el[0]=-div.el[i]*sin;
					v.el[1]=-div.el[i]*cos;
					v.el[2]=div.el[i]*sin;
					v.el[3]=div.el[i]*cos;
					v.el[4]=-div.el[i]*sin;
					v.el[5]=-div.el[i]*cos;
					
					v.el[6]=dB;
				
					for(int p=0;p<v.length;p++)
						pwBun.print(v.el[p]+"\t");
					pwBun.println();
				}
		}		
				
					

				util.pr("b_Time was written to "+file1+".");

				pwBun.close();
			}
			catch(IOException e){}
		
			}
		
			catch(IOException e){System.err.println("Error in loading BH data file.");

			}
		
		

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
			pwBun.println(Lani+"\t"+nAni); 		//	pwBun.println(Lani+"\t"+nAni);
			pwBun.println("* B * H ･････ *　磁化容易軸");
			
			for(int i=0;i<Lani;i++){
				pwBun.print(dfB.format(BHani[0].el[i][0])+"\t");
				for(int j=0;j<nAni;j++){
					pwBun.print(dfH.format(BHani[j].el[i][0])+"\t");
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
	
	double Hseff=BHs[0][0].el[BHs[0][0].nRow-1][0];

	DecimalFormat dfB=new DecimalFormat("#.00");
	DecimalFormat dfH=new DecimalFormat("#.0");


		try{
			PrintWriter pwBun = new PrintWriter(new BufferedWriter(new FileWriter(file)));		

for(int ia=0;ia<nSet;ia++){

			pwBun.println(1+"\t"+1+"\t"+nSet+"\t"+ia);
			pwBun.println("*Bs*Hs*");
			pwBun.println(Bseff+"\t"+Hseff);

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
