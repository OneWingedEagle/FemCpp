package math;

import static java.lang.Math.sqrt;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.Scanner;

import fem.Model;



public class PlayModelGraph {
	
Mat[] BH;
int numb;
double Bs, Hs;
int numb2,numb3,numb4,numb5;
String regex="[:; ,\\t]+";


public static void main(String[] args){

	PlayModelGraph pg=new PlayModelGraph();
	
	pg.loadShapeFunc();
}
	
public boolean loadShapeFunc(){

	/*String file=util.getFile();
	if(file==null || file.equals("") )return false;*/

	String file="C:\\Works\\HVID\\shape";

		try{
			FileReader fr=new FileReader(file);
			BufferedReader br = new BufferedReader(fr);
			String line;
			String s;
			String[] sp;


			
			line=br.readLine();
			sp=line.split(regex);	
			
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
			
			this.numb=Integer.parseInt(sp[0]);
			this.Bs=Double.parseDouble(sp[1]);
			this.Hs=Double.parseDouble(sp[2]);
			this.numb2=Integer.parseInt(sp[3]);
			this.numb3=Integer.parseInt(sp[4]);
			this.numb4=Integer.parseInt(sp[5]);
			this.numb5=Integer.parseInt(sp[6]);
			
		
			//util.show(sp);
			
			
			Mat[] BH1=new Mat[200];
			int[] LBH1=new int[200];
			
			int Lx=10000;
			for( int p=0;p<BH1.length;p++)
			BH1[p]=new Mat(Lx,2);
		
			int iloop=0;
			int jx=0;
			while(br.ready()){
				line=br.readLine();
				sp=line.split(regex);	
				if(sp.length==1) break;
				double[] bh=getCSV(line);
				if(jx>0 && bh[0]<BH1[iloop].el[jx-1][0])
				{
					iloop++;
					jx=0;
				}
				BH1[iloop].el[jx][0]=bh[0];
				BH1[iloop].el[jx][1]=bh[1];
				
				LBH1[iloop]++;
				jx++;
			}
			
			this.BH=new Mat[iloop];
			for( int p=0;p<BH.length;p++){
				BH[p]=new Mat(LBH1[p],2);
				for( int i=0;i<BH[p].nRow;i++)
				{
					BH[p].el[i][0]=BH1[p].el[i][0];
					BH[p].el[i][1]=BH1[p].el[i][1];
				}
			}
			
			
			util.plotBunch(BH);
			
		//	BH1[0].show();
			return true;
		
			}
			catch(IOException e){System.err.println("Error in loading BH data file.");
			return false;
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
