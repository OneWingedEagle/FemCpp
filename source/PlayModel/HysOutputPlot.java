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
	

String regex="[:; ,\\t]+";


public static void main(String[] args){

	HysOutputPlot pg=new HysOutputPlot();
	
}
	
public  HysOutputPlot(){

	/*String file=util.getFile();
	if(file==null || file.equals("") )return false;*/
//	String file="C:\\Works\\HVID\\folder1\\data\\A_Bì¸óÕëŒèÃÉãÅ[Évhts_data\\hys_data";

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
		
			
			
			
			Mat[] XX=new Mat[2];
			XX[0]=new Mat(L,2);
			XX[1]=new Mat(L,2);
			
		
			XX[0].setCol(bbhh.getColVect(0).times(1), 0);
			XX[0].setCol(bbhh.getColVect(1).times(1), 1);
			
			XX[1].setCol(bbhh.getColVect(2), 0);
			XX[1].setCol(bbhh.getColVect(3), 1);
			
			
			//util.plot(XX[1].el);
			
			
				Vect Hr=new Vect(XX[0].nRow);
			Vect Br=new Vect(XX[0].nRow);
			
				for(int i=0;i<Hr.length;i++){
					Hr.el[i]=new Vect(XX[1].el[i][0],XX[1].el[i][1]).norm();
					Br.el[i]=new Vect(XX[0].el[i][0],XX[0].el[i][1]).norm();
				}
				//Hr.show();
			//	Br.show();

			//	util.plot(Hr,Br);

			util.plotBunch(XX);
			//BH[1].show();
			
			br.close();
			fr.close();
	
		
		
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
