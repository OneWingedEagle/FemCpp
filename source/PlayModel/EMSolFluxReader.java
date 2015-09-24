package PlayModel;

import static java.lang.Math.PI;
import static java.lang.Math.acos;
import static java.lang.Math.cos;
import static java.lang.Math.sqrt;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
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



public class EMSolFluxReader {

	String regex="[:; ,\\t]+";


	public static void main(String[] args){

		EMSolFluxReader x=new EMSolFluxReader();
		
		//x.fetchB();
		x.fetchJQ();
	
	}
	
	public void fetchB(){
		
		int elNumb=26942;

		
			int nSteps=20;
	
			String fluxFile="C:\\Works\\Problems and Models\\Large-scale Test\\classic\\magnetic";

			Mat BB=extractFlux( fluxFile,3,nSteps,  elNumb);

			BB.show();
			
			util.plot(BB.getColVect(0));

		
		
	}
	
	
	public void fetchBH(){
				
		int elNumb=390;
			int nSteps=21;
			
		String bhfolder="C:\\Works\\Problems and Models\\Large-scale Test";
		Mat BH=	getBHcurve( bhfolder,3,nSteps, elNumb,0);
//		BH.show(); 
		Mat HB=new Mat(BH.size());
		HB.setCol(BH.getColVect(1), 0);
		HB.setCol(BH.getColVect(0), 1);
		HB.show();
		util.plot(BH);
////		util.plot(BH.getColVect(0));
		util.plot(BH.getColVect(1));
	}
	
	public void fetchJQ(){
		
		int elNumb=6513;

			int nSteps=20;
		
			String fluxFile="C:\\Works\\Problems and Models\\Large-scale Test\\classic\\current";

			Mat JJ=extractCurrent( fluxFile,3,nSteps,  elNumb);

			JJ.show();
			
			util.plot(JJ.getColVect(0));

		
		
	}

	public Mat extractFlux(String bbf,int dim, int nSteps, int nelem){

		Vect[] B=new Vect[nSteps];
		
		for(int i=0;i<nSteps;i++)
		B[i]=new Vect(dim);
		
		String regex="[ ,\\t]+";
		try{

			File f=new File(bbf);
			FileReader fr=new FileReader(f);
			BufferedReader br = new BufferedReader(fr);
			String line="";
			String[] sp=new String[15];
	
			for(int i=0;i<nSteps;i++){
						
			while(!line.startsWith("STEP")){
			line=br.readLine();
					}
		
			String ss="";
			while(!ss.equals(Integer.toString(nelem))){
				line=br.readLine();
				sp=line.split(regex);

				ss=sp[1];
				
						}
			sp=line.split(regex);
			for(int j=0;j<dim;j++)
				B[i].el[j]=Double.parseDouble(sp[2+j]);
		
			}

				
		
	br.close();
	fr.close();
	
	Mat result=new Mat(nSteps,dim);
	for(int i=0;i<nSteps;i++)
		result.el[i]=B[i].el;
	
return result;
	
		}


		catch(Exception e){
			System.err.println("error");	e.printStackTrace(); 
			return null;
		}
		

	
	}
	
	
	public Mat getBHcurve(String bhfolder,int dim, int nSteps, int nelem,double angdeg){

		
		String fileB=bhfolder+"\\magnetic";
		String fileH=bhfolder+"\\magnetization";
		
		Mat BH1=new Mat(nSteps,2);
		Mat BH=null;
		Vect B=new Vect(dim);
		Vect H=new Vect(dim);
		Vect er=new Vect(dim);
		er.el[0]=Math.cos(angdeg*Math.PI/180);
		er.el[1]=Math.sin(angdeg*Math.PI/180);
		
			String regex="[ ,\\t]+";
		try{

			
			
			File f=new File(fileB);
			FileReader fr=new FileReader(f);
			BufferedReader br = new BufferedReader(fr);
			String line="";
			String[] sp=new String[15];
	
			int ix=0;
			
			for(int i=0;i<nSteps;i++){
					
			while(!line.startsWith("STEP") ){
				line=br.readLine();
	

				if(line==null) break;
			
					}
			if(line==null) break;

			String ss="";
			while(!ss.equals(Integer.toString(nelem))){
	
				line=br.readLine();

				sp=line.split(regex);

				ss=sp[1];
				
						}
			sp=line.split(regex);
			for(int j=0;j<dim;j++)
				B.el[j]=Double.parseDouble(sp[2+j]);
		
	ix++;

			BH1.el[i][1]=B.dot(er);
			//BH.el[i][1]=B.el[1];
		
		
			}
			
			 f=new File(fileH);
			 fr=new FileReader(f);
			 br = new BufferedReader(fr);
			 line="";
			 sp=new String[15];
	
				for(int i=0;i<nSteps;i++){
					
					while(!line.startsWith("STEP") ){
						line=br.readLine();
				

						if(line==null) break;
					
							}
					if(line==null) break;
				
					String ss="";
					while(!ss.equals(Integer.toString(nelem))){

				line=br.readLine();
				sp=line.split(regex);

				ss=sp[1];
				
						}
			sp=line.split(regex);
			for(int j=0;j<dim;j++)
				H.el[j]=Double.parseDouble(sp[2+j]);
		
	

			BH1.el[i][0]=H.dot(er);
			//BH.el[i][0]=H.el[0];
		
		
			}
	
		
				BH=new Mat(ix,2);
				for(int i=0;i<ix;i++)
					BH.el[i]=BH1.el[i];
			
			//util.plot(BH.getColVect(1));

		//	BH.getColVect(1).show();
			
			br.close();
			fr.close();

	
	
		}


		catch(Exception e){System.err.println("error");	e.printStackTrace(); }

return BH;
	
	
	}
	
	public Mat extractCurrent(String bbf,int dim, int nSteps, int nelem){

		Vect[] B=new Vect[nSteps];
		
		double[] Q=new double[nSteps];
		
		for(int i=0;i<nSteps;i++)
		B[i]=new Vect(dim);
		
		String regex="[ ,\\t]+";
		try{

			File f=new File(bbf);
			FileReader fr=new FileReader(f);
			BufferedReader br = new BufferedReader(fr);
			String line="";
			String[] sp=new String[15];
		
			
			for(int i=0;i<nSteps;i++){
						
			while(!line.startsWith("STEP")){
			line=br.readLine();
	
					}
		
			String ss="";
			while(!ss.equals(Integer.toString(nelem))){
				line=br.readLine();
				sp=line.split(regex);

				ss=sp[1];
				
						}
			sp=line.split(regex);
			for(int j=0;j<dim;j++)
				B[i].el[j]=Double.parseDouble(sp[2+j]);
			
			Q[i]=Double.parseDouble(sp[2+dim]);

		
			}

		
	br.close();
	fr.close();
	
	Mat result=new Mat(nSteps,dim+1);
	for(int i=0;i<nSteps;i++){
		for(int j=0;j<dim;j++)
		result.el[i][j]=B[i].el[j];
		
		result.el[i][dim]=Q[i];
	}
	
return result;
	
		}


		catch(Exception e){
			System.err.println("error");	e.printStackTrace(); 
			return null;
		}
		

	
	}
}
