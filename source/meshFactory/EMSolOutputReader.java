package meshFactory;

import static java.lang.Math.PI;
import static java.lang.Math.acos;
import static java.lang.Math.cos;
import static java.lang.Math.sqrt;

import java.awt.Image;
import java.io.BufferedInputStream;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.ByteArrayOutputStream;
import java.io.File;
import java.io.FileOutputStream;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStream;
import java.io.PrintWriter;
import java.net.MalformedURLException;
import java.net.URI;
import java.net.URL;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.text.DecimalFormat;
import java.util.Scanner;

import javax.imageio.ImageIO;
import javax.tools.FileObject;

import math.Mat;
import math.Vect;
import math.util;
import fem.Model;

import java.io.File;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.InputStream;
import java.io.OutputStream;


public class EMSolOutputReader {

	private static final int nSourceFields = 0;
	String regex="[:; ,\\t]+";
	private Vect[][] sourceFieldCurrent;


	public static void main(String[] args) throws IOException{

		EMSolOutputReader x=new EMSolOutputReader();
		
		String folder="//192.168.12.103/Users/hassan/Documents/Large Scale Test/FieldSources/results/version 20151110 225548/periodic_source/independent/Igiven/refine 0/REGULARIZATION 1";
		x.compareOutputs(folder);

	}

	public Vect loadOutput(){

		/*String file=util.getFile();
	if(file==null || file.equals("") )return false;*/
		//	String file="C:\\Works\\HVID\\folder1\\data\\A_B入力対称ループhts_data\\hys_data";

		//String file="C:\\Works\\EMSolBuild_C\\EMSolBatch\\Large model_iso\\output";
		String file="C:\\Works\\EMSolBuild_C\\EMSolBatch\\Small model\\output";
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
			sp=line.split(regex);
			tempNREr.el[ix]=Double.parseDouble(sp[sp.length-1]);
			
			line=br.readLine();

			sp=line.split(regex);
			tempdB.el[ix]=Double.parseDouble(sp[4]);
			
			ix++;

			
			}
	br.close();
	fr.close();
	
	Vect errs=new Vect(ix);
	Vect NRerr=new Vect(ix);
	Vect dBs=new Vect(ix);
	for(int i=0;i<ix;i++){
		errs.el[i]=tempICCGEr.el[i];
		NRerr.el[i]=tempNREr.el[i];
		dBs.el[i]=tempdB.el[i];
	}
	
	util.plot(dBs);
	//util.plot(NRerr);
	//errs.show();
	return errs;
		}

		catch(IOException e){System.err.println("Error in loading output file.");
		return null;
		}

		
	}
	
	private void compareOutputs(String folder){
		
/*		//Vect  NRError=x.loadOutput();
		
		
		//String file="//192.168.12.103/Users/hassan/Documents/Large Scale Test/FieldSources/domain/output";
			
			//String file="//192.168.12.103/Users/hassan/Documents/Large Scale Test/FieldSources/results/transient_linear_series_voltage_input/release/refine 0/output";
			String file="//192.168.12.103/Users/hassan/Documents/Large Scale Test/FieldSources/results/transient_linear_series_voltage_input/release/refine 0/output";
			
		//String file="C:\\Works\\Problems and Models\\Large-scale Test\\FieldSources\\transient_linear_series_voltage_input\\release\\refine 1-2coils\\output";
			
			file="//192.168.12.103/Users/hassan/Documents/Large Scale Test/FieldSources/results/version 20151102 111751/transient_linear_series_current_input/release/refine 0/output";
			

			file="//192.168.12.103/Users/hassan/Documents/Large Scale Test/FieldSources/results/version 20151102 111751/parallel/transient_linear_parallel_curren_input/release/refine 0/output";
			//file="\\Hd_dmz\\em\\COMMON\\計画実績\\年度単位\\FY2015\\試験研究\\大規模化Pj\\大規模MPI版評価検証\\verification-program_201510\\result\\reduced_potential_on_periodc_boundary\\release\\output";

			
			file="//192.168.12.103/Users/hassan/Documents/Large Scale Test/FieldSources/domain/output";
			
			//Vect  vt=x.loadSourceVoltage(file,1);
			//util.plot(vt);
			//file="//192.168.12.103/Users/hassan/Desktop/release2/output";
*/
		 folder="//192.168.12.103/Users/hassan/Documents/Large Scale Test/FieldSources/results/version 20151110 225548/essential_para_8_bug";
		 
		int nFiles=3;
		
		String[] tip=new String[nFiles];
		tip[0]="/release/ELCUR_only/output";
		tip[1]="/sequential/ELCUR_only/output";
		tip[2]="/MPI_para_1/ELCUR_only/output";
		
		
		String[] file=new String[nFiles];
		
	/*	file[0]=folder+"/release/3sources_no_phi/output";
		file[1]=folder+"/MPI para 1/3sources_no_phi/output";
		file[2]=folder+"/MPI para 8/3sources_no_phi/output";
		nFiles=3;*/
		for(int i=0;i<nFiles;i++)
		file[i]=folder+tip[i];


			
	//util.pr(file[0]);
			int nPowerSources=0;
			int nFieldSources=0;
			int nFluxes=0;
			int nEnergies=0;
			
			int nMax=100;
			
			for(int i=0;i<nMax;i++){
				Vect v=loadSourceCurrent(file[0],i+1);
			
				if(v==null) break;
				nPowerSources++;
			}
			
			Vect[][] sourceVoltage=new Vect[nFiles][nPowerSources];
			Vect[][] sourceCurrent=new Vect[nFiles][nPowerSources];
			for(int nfile=0;nfile<nFiles;nfile++){
				for(int i=0;i<nPowerSources;i++){
					sourceVoltage[nfile][i]=loadSourceVoltage(file[nfile],i+1);
					sourceCurrent[nfile][i]=loadSourceCurrent(file[nfile],i+1);

			}
			}
			

			
		//--------------------------------------------
			
		
			
			for(int i=0;i<nMax;i++){
				Vect v=loadFieldSourceCurrent(file[0],i+1);
				if(v==null) break;
				nFieldSources++;
			
			}
			
			Vect[][] fieldSourceVoltage=new Vect[nFiles][nFieldSources];
			Vect[][] fieldSourceCurrent=new Vect[nFiles][nFieldSources];
			
			for(int nfile=0;nfile<nFiles;nfile++)
				for(int i=0;i<nFieldSources;i++){
					fieldSourceVoltage[nfile][i]=loadFieldSourceVoltage(file[nfile],i+1);
						fieldSourceCurrent[nfile][i]=loadFieldSourceCurrent(file[nfile],i+1);

			}

			//--------------------------------------------	
				Vect[][] flux=new Vect[nFiles][nMax];
				
				for(int i=0;i<nMax;i++){
					Vect v=loadMagFlux(file[0],i+1);
					if(v==null) break;
					nFluxes++;
					flux[0][i]=v.deepCopy();
				}
				
				for(int nfile=1;nfile<nFiles;nfile++)
					for(int i=0;i<nFluxes;i++){

						flux[nfile][i]=loadMagFlux(file[nfile],i+1);

				}
				
			
		//--------------------------------------------	
				Vect[][] energy=new Vect[nFiles][nMax];
				
				for(int i=0;i<nMax;i++){
					Vect v=loadMagEnergy(file[0],i+1);
					if(v==null) break;
					nEnergies++;
					energy[0][i]=v.deepCopy();
				}
				
				for(int nfile=1;nfile<nFiles;nfile++)
					for(int i=0;i<nEnergies;i++){

						energy[nfile][i]=loadMagEnergy(file[nfile],i+1);

				}

		
				int nT=0;
				
				if(nPowerSources>0)
					nT=sourceCurrent[0][0].length;
				else if(nFluxes>0) nT=flux[0][0].length;
				else if(nEnergies>0) nT=energy[0][0].length;
				
				for(int i=0;i<nPowerSources;i++){
					for(int nfile=0;nfile<nFiles;nfile++){
						if(sourceCurrent[nfile][i]==null) sourceCurrent[nfile][i]=new Vect(nT);
						if(sourceVoltage[nfile][i]==null) sourceVoltage[nfile][i]=new Vect(nT);
					}
				}
				
				for(int i=0;i<nFieldSources;i++){
					for(int nfile=0;nfile<nFiles;nfile++){
						if(fieldSourceCurrent[nfile][i]==null) fieldSourceCurrent[nfile][i]=new Vect(nT);
						if(fieldSourceVoltage[nfile][i]==null) fieldSourceVoltage[nfile][i]=new Vect(nT);
					}
				}
				
				for(int i=0;i<nFluxes;i++){
					for(int nfile=0;nfile<nFiles;nfile++){
						if(flux[nfile][i]==null) flux[nfile][i]=new Vect(nT);
					}
				}
				for(int i=0;i<nEnergies;i++){
					for(int nfile=0;nfile<nFiles;nfile++){
						if(energy[nfile][i]==null) energy[nfile][i]=new Vect(nT);
					}
				}
			
				
			Vect time=	loadTimesSteps(file[0]);

		
			
			Mat tables=new Mat(nT,1+(2*nPowerSources+2*nFieldSources+nFluxes+nEnergies)*(2*nFiles-1));
			
			int ix=0;
			tables.setCol(time, ix++);
			
			
			for(int i=0;i<nPowerSources;i++){
				for(int nfile=0;nfile<nFiles;nfile++){
					tables.setCol(sourceCurrent[nfile][i], ix++);
				}
				
				for(int nfile=1;nfile<nFiles;nfile++){
					tables.setCol(sourceCurrent[nfile][i].sub(sourceCurrent[0][i]).abs().div(sourceCurrent[0][i].abs()).times(100), ix++);
				}
			}
			//--------------------------
			for(int i=0;i<nPowerSources;i++){
				for(int nfile=0;nfile<nFiles;nfile++){
				tables.setCol(sourceVoltage[nfile][i], ix++);
				}
				
				for(int nfile=1;nfile<nFiles;nfile++){
					tables.setCol(sourceVoltage[nfile][i].sub(sourceVoltage[0][i]).abs().div(sourceVoltage[0][i].abs()).times(100), ix++);
				}
			}
			
			//--------------------
			for(int i=0;i<nFieldSources;i++){
				for(int nfile=0;nfile<nFiles;nfile++){
				tables.setCol(fieldSourceCurrent[nfile][i], ix++);
				}
				
				for(int nfile=1;nfile<nFiles;nfile++){
					tables.setCol(fieldSourceCurrent[nfile][i].sub(fieldSourceCurrent[0][i]).abs().div(fieldSourceCurrent[0][i].abs()).times(100), ix++);
				}
			}
			//--------------------------
			//--------------------
			for(int i=0;i<nFieldSources;i++){
				for(int nfile=0;nfile<nFiles;nfile++){
				tables.setCol(fieldSourceVoltage[nfile][i], ix++);
				}
				
				for(int nfile=1;nfile<nFiles;nfile++){
					tables.setCol(fieldSourceVoltage[nfile][i].sub(fieldSourceVoltage[0][i]).abs().div(fieldSourceVoltage[0][i].abs()).times(100), ix++);
				}
			}
			
			for(int i=0;i<nFluxes;i++){
				for(int nfile=0;nfile<nFiles;nfile++){
				tables.setCol(flux[nfile][i].times(1000), ix++);
				}
				
				for(int nfile=1;nfile<nFiles;nfile++){
					tables.setCol(flux[nfile][i].sub(flux[0][i]).abs().div(flux[0][i].abs()).times(100), ix++);
				}
			}
			
			for(int i=0;i<nEnergies;i++){
				for(int nfile=0;nfile<nFiles;nfile++){
				tables.setCol(energy[nfile][i], ix++);
				}
				
				for(int nfile=1;nfile<nFiles;nfile++){
					tables.setCol(energy[nfile][i].sub(energy[0][i]).abs().div(energy[0][i].abs()).times(100), ix++);
				}
			}
			
			String[] titles1=new String[tables.nCol];
			String[] titles2=new String[tables.nCol];
			
			tip[0]="release";
			tip[1]="sequential";
			tip[2]="MPI_para_1";
			
			ix=0;
			titles1[ix]="Time(sec.)";	
			titles2[ix++]="-";
			for (int i=0;i<nPowerSources;i++){
				for(int nfile=0;nfile<nFiles;nfile++){
					titles2[ix]=tip[nfile];	
				titles1[ix++]="PS"+(i+1)+"(A)";	
			}
			for(int nfile=1;nfile<nFiles;nfile++){
				titles2[ix]=tip[nfile];	
				titles1[ix++]="PS"+(i+1)+"(A)err.(%)";	
			}
			}
			
			for (int i=0;i<nPowerSources;i++){
				for(int nfile=0;nfile<nFiles;nfile++){
					titles2[ix]=tip[nfile];	
				titles1[ix++]="PS"+(i+1)+"(V)";	
			}
			for(int nfile=1;nfile<nFiles;nfile++){
				titles2[ix]=tip[nfile];	
				titles1[ix++]="PS"+(i+1)+"V)err.(%)";	
			}
			}
			//============================
			for (int i=0;i<nPowerSources;i++){
				for(int nfile=0;nfile<nFiles;nfile++){
					titles2[ix]=tip[nfile];	
				titles1[ix++]="FS"+(i+1)+"(A)";	
			}
			for(int nfile=1;nfile<nFiles;nfile++){
				titles2[ix]=tip[nfile];	
				titles1[ix++]="FS"+(i+1)+"(A)err.(%)";	
			}
			}
			
			for (int i=0;i<nPowerSources;i++){
				
				for(int nfile=0;nfile<nFiles;nfile++){
					titles2[ix]=tip[nfile];	
				titles1[ix++]="FS"+(i+1)+"(V)";	
			}
			for(int nfile=1;nfile<nFiles;nfile++){
				titles2[ix]=tip[nfile];	
				titles1[ix++]="FS"+(i+1)+"(V)err.(%)";	
			}
			}
			
			//============================
			for (int i=0;i<nFluxes;i++){
				for(int nfile=0;nfile<nFiles;nfile++){
					titles2[ix]=tip[nfile];	
				titles1[ix++]="Flux"+(i+1)+"(mWb)";	
			}
			for(int nfile=1;nfile<nFiles;nfile++){
				titles2[ix]=tip[nfile];	
				titles1[ix++]="Flux"+(i+1)+"Err.(%)";	
			}
			}
			
			for (int i=0;i<nEnergies;i++){
				for(int nfile=0;nfile<nFiles;nfile++){
					titles2[ix]=tip[nfile];	
				titles1[ix++]="Energy(J)";	
			}
			for(int nfile=1;nfile<nFiles;nfile++){
				titles2[ix]=tip[nfile];	
				titles1[ix++]="Energy"+(i+1)+"Err.(%)";	
			}
			}
	
			for (int i=0;i<titles1.length;i++){
			System.out.print("\t"+titles1[i]);
			}
			
			System.out.println();
			
			for (int i=0;i<titles2.length;i++){
				System.out.print("\t"+titles2[i]);
			}
			System.out.println();
		tables.show();
		
	}
	
	public Vect loadCheck(){


		/*String file=util.getFile();
	if(file==null || file.equals("") )return false;*/
		//	String file="C:\\Works\\HVID\\folder1\\data\\A_B入力対称ループhts_data\\hys_data";

		//String file="C:\\Works\\EMSolBuild_C\\EMSolBatch\\Large model_iso\\output";
		//String file="C:\\Works\\EMSolBuild_C\\EMSolBatch\\Large model_Angs\\output";
		String file="C:\\Works\\EMSolBuild_C\\EMSolBatch\\Small model\\check";
		//String file="C:\\Works\\HVID\\hys_data";
	//	String file=System.getProperty("user.dir") + "\\hys_dataH.txt";
		
		Vect tempICCGEr=new Vect(100000);
	
		try{
			FileReader fr=new FileReader(file);
			BufferedReader br = new BufferedReader(fr);
			String line;
			String s;
			String[] sp;
			int ix=0;
			while((line=br.readLine())!=null){
			while((line=br.readLine())!=null && !line.startsWith("*** ICCG iteration")){}
			if(line==null) {break;}
			line=br.readLine();
		
			while(!(line=br.readLine()).startsWith("*")){
		
				sp=line.split(regex);
				//util.show(sp);
				tempICCGEr.el[ix]=Double.parseDouble(sp[1]);
				ix++;
			}
			
			
		
			
			continue;
			
			}
	br.close();
	fr.close();
	
	
	Vect errs=new Vect(ix);
	//Vect NRerr=new Vect(ix);
	//Vect dBs=new Vect(ix);
	for(int i=0;i<ix;i++){
		errs.el[i]=tempICCGEr.el[i];
		//NRerr.el[i]=tempNREr.el[i];
		//dBs.el[i]=tempdB.el[i];
	}
	
	util.plot(errs);
	//util.plot(NRerr);
	//errs.show();
	return errs;
		}

		catch(IOException e){System.err.println("Error in loading output file.");
		return null;
		}

		
	
	}
	
	public Vect loadEnergy(){

	
		String file="C:\\Works\\EMSolBuild_C\\EMSolBatch\\ringCompositAngDep\\output";

		
		Vect tempEnrgy=new Vect(10000);

		try{
			FileReader fr=new FileReader(file);
			BufferedReader br = new BufferedReader(fr);
			String line;
			String s;
			String[] sp;
			int ix=0;
			while((line=br.readLine())!=null){
			while((line=br.readLine())!=null && !line.startsWith("*** Total magnetic energy ")){}
			if(line==null) {break;}
			line=br.readLine();
			line=br.readLine();
		//	line=br.readLine();
			sp=line.split(regex);
			tempEnrgy.el[ix]=Double.parseDouble(sp[2]);
			
			ix++;

			
			}
			
	br.close();
	fr.close();
	
	if(ix==0) return null;
	
	Vect energy=new Vect(ix);

	for(int i=0;i<ix;i++){
		energy.el[i]=tempEnrgy.el[i];

	}
	
	util.plot(energy);
	//util.plot(NRerr);
	energy.times(1e6).show();
	return energy;
		}

		catch(IOException e){System.err.println("Error in loading output file.");
		return null;
		}

		
	}	
	
	public Vect loadSourceCurrent(String file, int sourceID){

	
		
		Vect powerSourceCurrent=new Vect(10000);

		try{
			FileReader fr=new FileReader(file);
			BufferedReader br = new BufferedReader(fr);
			String line;
			String s;
			String[] sp;
			int ix=0;
			while((line=br.readLine())!=null){
			while((line=br.readLine())!=null && !line.startsWith("***  Power Sources")){}
			if(line==null) {break;}

			line=br.readLine();
			for(int i=0;i<sourceID;i++){
			line=br.readLine();
			if(line.startsWith("*"))
			 break;
			}
			//util.pr(line);
		//	line=br.readLine();
			sp=line.split(regex);
			if(sp.length<=2) break;
			//util.pr(sp[2]);
			powerSourceCurrent.el[ix]=Double.parseDouble(sp[2]);

			ix++;

			
			}
			
	br.close();
	fr.close();
	
	if(ix==0) return null;
	
	Vect current=new Vect(ix);

	for(int i=0;i<ix;i++){
		current.el[i]=powerSourceCurrent.el[i];

	}

	//util.plot(current);
	//util.plot(NRerr);

	return current;
		}

		catch(IOException e){System.err.println("Error in loading output file.");
		return null;
		}

		
	}	
	
	
	public Vect loadFieldSourceCurrent(String file,int sourceID){

		
		Vect powerSourceCurrent=new Vect(10000);

		try{
			FileReader fr=new FileReader(file);
			BufferedReader br = new BufferedReader(fr);
			String line;
			String s;
			String[] sp;
			int ix=0;
			while((line=br.readLine())!=null){
			while((line=br.readLine())!=null && !line.startsWith("***  Sources ")){}
			if(line==null) {break;}
			line=br.readLine();
			for(int i=0;i<sourceID;i++){
			line=br.readLine();
			if(line.startsWith("**")) break;
			}
			if(line.startsWith("**")) break;
			sp=line.split(regex);
			if(sp.length<=2) break;
			powerSourceCurrent.el[ix]=Double.parseDouble(sp[2]);
			
			ix++;

			
			}
			
	br.close();
	fr.close();
	
	if(ix==0) return null;
	
	Vect voltage=new Vect(ix);

	for(int i=0;i<ix;i++){
		voltage.el[i]=powerSourceCurrent.el[i];

	}
	
	//util.plot(voltage);
	//util.plot(NRerr);

	return voltage;
		}

		catch(IOException e){System.err.println("Error in loading output file.");
		return null;
		}

		
	}	
	
	
	
	public Vect loadSourceVoltage(String file,int sourceID){

	
		Vect powerSourceCurrent=new Vect(10000);

		try{
			FileReader fr=new FileReader(file);
			BufferedReader br = new BufferedReader(fr);
			String line;
			String s;
			String[] sp;
			int ix=0;
			while((line=br.readLine())!=null){
			while((line=br.readLine())!=null && !line.startsWith("***  Power Sources")){}
			if(line==null) {break;}
			line=br.readLine();
			for(int i=0;i<sourceID;i++){
			line=br.readLine();
			if(line.startsWith("**")) break;
			}
			if(line.startsWith("**")) break;
			sp=line.split(regex);
			if(sp.length<=2) break;
			powerSourceCurrent.el[ix]=Double.parseDouble(sp[3]);
			
			ix++;

			
			}
			
	br.close();
	fr.close();
	
	if(ix==0) return null;
	
	Vect voltage=new Vect(ix);

	for(int i=0;i<ix;i++){
		voltage.el[i]=powerSourceCurrent.el[i];

	}
	
	//util.plot(voltage);
	//util.plot(NRerr);

	return voltage;
		}

		catch(IOException e){System.err.println("Error in loading output file.");
		return null;
		}

		
	}	
	
	public Vect loadFieldSourceVoltage(String file,int sourceID){

		
		Vect powerSourceCurrent=new Vect(10000);

		try{
			FileReader fr=new FileReader(file);
			BufferedReader br = new BufferedReader(fr);
			String line;
			String s;
			String[] sp;
			int ix=0;
			while((line=br.readLine())!=null){
			while((line=br.readLine())!=null && !line.startsWith("***  Sources ")){}
			if(line==null) {break;}
			line=br.readLine();
			for(int i=0;i<sourceID;i++){
			line=br.readLine();
			if(line.startsWith("**")) break;
			}
			if(line.startsWith("**")) break;
			sp=line.split(regex);
			if(sp.length<=3) break;
			powerSourceCurrent.el[ix]=Double.parseDouble(sp[3]);
			
			ix++;

			
			}
			
	br.close();
	fr.close();
	
	if(ix==0) return null;
	
	Vect voltage=new Vect(ix);

	for(int i=0;i<ix;i++){
		voltage.el[i]=powerSourceCurrent.el[i];

	}
	
	//util.plot(voltage);
	//util.plot(NRerr);

	return voltage;
		}

		catch(IOException e){System.err.println("Error in loading output file.");
		return null;
		}

		
	}	
	
	public Vect loadMagFlux(String file,int loopID){

		
		
		Vect magFluxTemp=new Vect(10000);

		try{
			FileReader fr=new FileReader(file);
			BufferedReader br = new BufferedReader(fr);
			String line;
			String s;
			String[] sp;
			int ix=0;
			while((line=br.readLine())!=null){
			while((line=br.readLine())!=null && !line.startsWith("*            Magnetic fluxes of flux loops")){}
			if(line==null) {break;}
			line=br.readLine();
			line=br.readLine();
			line=br.readLine();
			line=br.readLine();
		
			for(int i=0;i<loopID;i++){
			line=br.readLine();
			if(line.startsWith("**")) break;
			}
			if(line.startsWith("**")) break;
			sp=line.split(regex);
			if(sp.length<=2) break;
			magFluxTemp.el[ix]=Double.parseDouble(sp[2]);
			
			ix++;

			
			}
			
	br.close();
	fr.close();
	
	if(ix==0) return null;
	
	Vect magFlux=new Vect(ix);

	for(int i=0;i<ix;i++){
		magFlux.el[i]=magFluxTemp.el[i];

	}
	
	//util.plot(voltage);
	//util.plot(NRerr);

	return magFlux;
		}

		catch(IOException e){System.err.println("Error in loading output file.");
		return null;
		}

		
	}	

	
public Vect loadMagEnergy(String file,int loopID){

		
		
		Vect magEnergyTemp=new Vect(10000);

		try{
			FileReader fr=new FileReader(file);
			BufferedReader br = new BufferedReader(fr);
			String line;
			String s;
			String[] sp;
			int ix=0;
			while((line=br.readLine())!=null){
			while((line=br.readLine())!=null && !line.startsWith("*** Total magnetic energy")){}
			if(line==null) {break;}
			line=br.readLine();
			for(int i=0;i<loopID;i++){
			line=br.readLine();
			if(line.startsWith("**")) break;
			}
			if(line.startsWith("**")) break;
			sp=line.split(regex);
			if(sp.length<=3) break;
			magEnergyTemp.el[ix]=Double.parseDouble(sp[3]);
			
			ix++;

			
			}
			
		
			
	br.close();
	fr.close();
	
	if(ix==0) return null;
	
	Vect magFnergy=new Vect(ix);

	for(int i=0;i<ix;i++){
		magFnergy.el[i]=magEnergyTemp.el[i];

	}
	
	//util.plot(voltage);
	//util.plot(NRerr);

	return magFnergy;
		}

		catch(IOException e){System.err.println("Error in loading output file.");
		return null;
		}

		
	}	

public Vect loadTimesSteps(String file){

	
	
	Vect powerSourceCurrent=new Vect(10000);

	try{
		FileReader fr=new FileReader(file);
		BufferedReader br = new BufferedReader(fr);
		String line;
		String s;
		String[] sp;
		int ix=0;
		while((line=br.readLine())!=null){
		while((line=br.readLine())!=null && !line.startsWith("----- Time step")){}
		if(line==null) {break;}

		//line=br.readLine();

		//util.pr(line);
	//	line=br.readLine();
		sp=line.split(regex);

		powerSourceCurrent.el[ix]=Double.parseDouble(sp[5]);

		ix++;

		
		}
		
br.close();
fr.close();

if(ix==0) return null;

Vect current=new Vect(ix);

for(int i=0;i<ix;i++){
	current.el[i]=powerSourceCurrent.el[i];

}

//util.plot(current);
//util.plot(NRerr);

return current;
	}

	catch(IOException e){System.err.println("Error in loading output file.");
	return null;
	}

	
}	


}
