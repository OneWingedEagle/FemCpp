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
	
Mat BB,HH;

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


			line=br.readLine();
			line=br.readLine();
			line=br.readLine();

int L=Integer.parseInt(line);

	Mat bbhh=new Mat(L,4);
	
	for( int p=0;p<L;p++)
	{
		line=br.readLine();
		double[] ar=this.getCSV(line);
		bbhh.el[p]=ar;
	}
		
			
			BB=new Mat(L,2);
			HH=new Mat(L,2);
			
		
		BB.setCol(bbhh.getColVect(0).times(1), 0);
		BB.setCol(bbhh.getColVect(1).times(1), 1);
			
			HH.setCol(bbhh.getColVect(2), 0);
			HH.setCol(bbhh.getColVect(3), 1);
			
/*			int Lx=90;
			BB.el[Lx][0]*=2;
			BB.el[Lx][1]*=2;
			HH.el[Lx][0]*=2;
			HH.el[Lx][1]*=2;
			*/
					
			Mat[] XX=new Mat[2];
			XX[0]=BB.times(60);
			XX[1]=HH;
			//util.plotBunch(XX);
			
			
				Vect Hr=new Vect(BB.nRow);
			Vect Br=new Vect(BB.nRow);
			
				for(int i=0;i<Hr.length;i++){
					Hr.el[i]=new Vect(HH.el[i][0],HH.el[i][1]).norm();
					Br.el[i]=new Vect(BB.el[i][0],BB.el[i][1]).norm();
				}
				//Hr.show();
			//	Br.show();

		
				
		double ang=Math.atan(XX[0].el[1][1]/XX[0].el[1][0]);
		//util.pr(ang/Math.PI*180);
		
		Vect er=new Vect(Math.cos(ang),Math.sin(ang));
		for(int i=0;i<Hr.length;i++){
			Hr.el[i]=new Vect(HH.el[i][0],HH.el[i][1]).dot(er);
			Br.el[i]=new Vect(BB.el[i][0],BB.el[i][1]).dot(er);
			util.pr(Hr.el[i]+"\t"+Br.el[i]);
		}

		util.plot(Hr,Br);

	
			//BH[1].show();
			
			br.close();
			fr.close();
	
		
		
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
