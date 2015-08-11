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



public class EMSolOutputReader {

	String regex="[:; ,\\t]+";


	public static void main(String[] args){

		EMSolOutputReader x=new EMSolOutputReader();

		Vect  NRError=x.loadOutput();
		
	}

	public Vect loadOutput(){

		/*String file=util.getFile();
	if(file==null || file.equals("") )return false;*/
		//	String file="C:\\Works\\HVID\\folder1\\data\\A_Bì¸óÕëŒèÃÉãÅ[Évhts_data\\hys_data";

		String file="C:\\Works\\EMSolBuild_C\\EMSolBatch\\Large model_iso\\output";
	//	String file="C:\\Works\\EMSolBuild_C\\EMSolBatch\\Small model\\output";
		//String file="C:\\Works\\HVID\\hys_data";
	//	String file=System.getProperty("user.dir") + "\\hys_dataH.txt";
		
		Vect tempICCGEr=new Vect(100000);
		Vect tempNREr=new Vect(100000);
		Vect tempdB=new Vect(100000);
		try{
			FileReader fr=new FileReader(file);
			BufferedReader br = new BufferedReader(fr);
			String line;
			String s;
			String[] sp;
			int ix=0;
			while((line=br.readLine())!=null){
			while((line=br.readLine())!=null && !line.startsWith("*** ICCG converges ***")){}
			if(line==null) break;
			line=br.readLine();
			sp=line.split(regex);
			tempICCGEr.el[ix]=Double.parseDouble(sp[sp.length-1]);
			line=br.readLine();
					line=br.readLine();
			
			line=br.readLine();
		//	line=br.readLine();
			sp=line.split(regex);
			tempNREr.el[ix]=Double.parseDouble(sp[sp.length-1]);
		
			ix++;
			
			continue;
			
			}
	br.close();
	fr.close();
	
	Vect errs=new Vect(ix);
	Vect NRerr=new Vect(ix);
	for(int i=0;i<ix;i++){
		errs.el[i]=tempICCGEr.el[i];
		NRerr.el[i]=tempNREr.el[i];
	}
	
	//util.plot(dBs);
	util.plot(NRerr);
	//errs.show();
	return errs;
		}

		catch(IOException e){System.err.println("Error in loading output file.");
		return null;
		}

		
	}	
	


}
