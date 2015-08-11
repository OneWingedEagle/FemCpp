package PlayModel;

import static java.lang.Math.sqrt;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
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
	
public  HysOutputPlot(){


		}

public void loadData(){


	/*String file=util.getFile();
	if(file==null || file.equals("") )return false;*/
//	String file="C:\\Works\\HVID\\folder1\\data\\A_B“ü—Í‘Î�Ìƒ‹�[ƒvhts_data\\hys_data";

	//String file="C:\\Users\\Hassan Ebrahimi\\JavaWorks\\MagFem\\hys_data";
	//String file="C:\\Works\\HVID\\output";
	String file="C:\\Works\\PlayModel\\output";
	//String file="C:\\Works\\HVIDConv\\output";
	

		try{
			FileReader fr=new FileReader(file);
			BufferedReader br = new BufferedReader(fr);
			String line;
			String s;
			String[] sp;

int Nmax=20;

Mat[] BB=new Mat[Nmax];
Mat[] HH=new Mat[Nmax];

Mat[] BHs1=new Mat[Nmax];


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
	

	
			XX[2*i]=BB[i].times(1);
			XX[2*i+1]=HH[i];
			
			BHs1[i]=new Mat(L,2);
			
		//	HH.show();
			
				Vect Hr=new Vect(HH[i].nRow);
			Vect Br=new Vect(BB[i].nRow);
			
		
				
		double ang=Math.atan(BB[i].el[1][1]/BB[i].el[1][0]);


		//util.pr(ang/Math.PI*180);
		//util.pr(ang/Math.PI*180);
		
		Vect er=new Vect(Math.cos(ang),Math.sin(ang));
		
		//er.hshow();
		for(int j=0;j<Hr.length;j++){
			Hr.el[j]=new Vect(HH[i].el[j][0],HH[i].el[j][1]).dot(er);
			Br.el[j]=new Vect(BB[i].el[j][0],BB[i].el[j][1]).dot(er);
			//util.pr(Hr.el[i]+"\t"+Br.el[i]);
		}
		

		BHs1[i].setCol(Hr,0);
		BHs1[i].setCol(Br,1);

		
		line=br.readLine();

}

			//BH[1].show();
//util.plotBunch(BHs1,numbCurves);

Mat BH=BHs1[0];
Mat HB=new Mat(BH.size());
HB.setCol(BH.getColVect(1), 0);
HB.setCol(BH.getColVect(0), 1);
//HB.show();
//BHs1[0].show();
			br.close();
			fr.close();
			
		util.plotBunch(XX,2);
	

			int L=HH[0].nRow/2;
			Mat M=new Mat(L,2);
			for(int i=0;i<L;i++){
				M.el[i]=HH[0].el[i+L];
			}
			
			M.show();	
			
			Mat[] BHs=new Mat[numbCurves];
			
			for(int i=0;i<numbCurves;i++)
				BHs[i]=BHs1[i];
			//M.show();
		//	util.plot(M);
			//util.plotBunch(BHs);
		//	BHs[0].show();
			String file1="C:\\Works\\HVID\\b_times";

		
			}
		
			catch(IOException e){System.err.println("Error in loading BH data file.");

			}

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
