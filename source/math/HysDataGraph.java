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



public class HysDataGraph {
	
Mat[] BH;
int numb1;
double Bs, Hs;
int numb2,numb3,numb4,numb5;
String regex="[:; ,\\t]+";


public static void main(String[] args){

	HysDataGraph pg=new HysDataGraph();
	
	pg.loadHysData();
}
	
public boolean loadHysData(){

	/*String file=util.getFile();
	if(file==null || file.equals("") )return false;*/
//	String file="C:\\Works\\HVID\\folder1\\data\\A_B“ü—Í‘ÎÌƒ‹[ƒvhts_data\\hys_data";

	String file="C:\\Works\\HVID\\hys_data";

		try{
			FileReader fr=new FileReader(file);
			BufferedReader br = new BufferedReader(fr);
			String line;
			String s;
			String[] sp;


			line=br.readLine();
			line=br.readLine();
			line=br.readLine();
			sp=line.split(regex);	

	
		
			this.Bs=Double.parseDouble(sp[0]);
			this.Hs=Double.parseDouble(sp[1]);
			line=br.readLine();
			line=br.readLine();
			sp=line.split(regex);	
			this.numb1=Integer.parseInt(sp[0]);
			this.numb2=Integer.parseInt(sp[1]);
			this.numb3=Integer.parseInt(sp[2]);
			this.numb4=Integer.parseInt(sp[3]);
			this.numb5=Integer.parseInt(sp[4]);
			
			BH=new Mat[numb4];
			
			line=br.readLine();
			line=br.readLine();
			int L=Integer.parseInt(line);
			BH[0]=new Mat(L,2);
			
			for( int p=0;p<L;p++){
				line=br.readLine();
				sp=line.split(regex);	
				double[] bh=getCSV(line);
				BH[0].el[p][0]=bh[0];
				BH[0].el[p][1]=bh[1];
			
			}
		
		
	;
			line=br.readLine();
			int L1=0;
			for( int ip=1;ip<numb4;ip++){
				line=br.readLine();
				if(line.startsWith("*")) {line=br.readLine();};
				sp=line.split(regex);	
				//if(sp.length==1)
					L1=Integer.parseInt(sp[0]);
						
		
				BH[ip]=new Mat(L1,2);

				for( int i=0;i<L1;i++){
					line=br.readLine();
		
					double[] bh=getCSV(line);

					BH[ip].el[i][0]=bh[0];
					BH[ip].el[i][1]=bh[1];
				}
			
						
			}
		
			util.plotBunch(BH,6);
	
			
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
