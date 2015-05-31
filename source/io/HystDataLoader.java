package io;

import java.io.BufferedReader;

import java.io.FileReader;
import java.io.IOException;

import java.util.Arrays;

import math.Mat;
import math.PlayModel;
import math.Vect;
import math.util;
import fem.*;

/**
 * TODO Put here a description of what this class does.
 *
 * @author Hassan.
 *         Created Aug 15, 2012.
 */
public class HystDataLoader {

	private String regex="[:; ,=\\t]+";
	private String regex2="[\\[\\]\\s: ,=\\t]+";

	
public static void main2(String[] args){}
	

	public void loadData(PlayModel pm, String file){

		

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

					
				pm.Bs=Double.parseDouble(sp[0]);
				pm.Hs=Double.parseDouble(sp[1]);
				line=br.readLine();
				
				line=br.readLine();
				int[] numb=this.getCSInt(line);
				
				pm.nInit=numb[0];
				pm.nMajor=numb[1];
				pm.nSym=numb[2];
				pm.nDesc=numb[3];
				pm.nAsc=numb[4];
				
			pm.nTotCurves=pm.nInit+pm.nMajor+pm.nSym+pm.nDesc+pm.nAsc;
				
				
				pm.BHraw=new Mat[pm.nTotCurves];
				
				//line=br.readLine();
	/*			line=br.readLine();
				int L=Integer.parseInt(line);*/
			
				//pm.BH1[0]=new Mat(L,2);
				
/*				for( int p=0;p<L;p++){
					line=br.readLine();
					sp=line.split(regex);	
					double[] bh=getCSV(line);
					BH[0].el[p][0]=bh[0];
					BH[0].el[p][1]=bh[1];
				
				}
			
			
		
				line=br.readLine();
				*/
				int L1=0;
				for( int ip=0;ip<pm.nTotCurves;ip++){
					line=br.readLine();
					if(line.startsWith("*")) {line=br.readLine();};
					L1=Integer.parseInt(line);
			
					pm.BHraw[ip]=new Mat(L1,2);

					for( int i=0;i<L1;i++){
						line=br.readLine();
			
						double[] bh=getCSV(line);

						pm.BHraw[ip].el[i][0]=bh[0];
						pm.BHraw[ip].el[i][1]=bh[1];
					}
				
							
				}

			//	util.plotBunch(pm.BHraw);
	
			
				}
			
				catch(IOException e){System.err.println("Error in loading BH data file.");
			
				}

	}	
	


	private Vect getVectData(String line, int dim){

		String[] sp=line.split(regex2);	
	
		Vect v=new Vect(dim);
	
		int k=0;
		while(sp[k].equals("")){k++;}	

		k++;
			
		for( int p=k;p<sp.length;p++){
		v.el[p-k]=Double.parseDouble(sp[p]);
		}

		return v;
	}



	
	private double[] getTabedData(String line){
		String[] sp=line.split(regex);	
		int L=sp.length;
		double[] v=new double[L];
		for( int p=0;p<L;p++)
			v[p]=Double.parseDouble(sp[p]);

		return v;
	}
	
	private double[] getPair(String line){
		String[] sp=line.split(regex);	
		double[] v=new double[2];
		int k=0;
		if(sp[k].equals("")) k++;

		for( int p=0;p<2;p++)
			v[p]=Double.parseDouble(sp[k+p]);

		return v;
	}
	
	private double getScalarData(String line){
		String[] sp=line.split(regex);	
		return Double.parseDouble(sp[sp.length-1]);
	}

	private int getIntData(String line){
		String[] sp=line.split(regex);	
		return Integer.parseInt(sp[sp.length-1]);
	}

	private boolean getBooleanData(String line){
		boolean b=false;
		String[] sp=line.split(regex);	
		
		if(sp[sp.length-1].startsWith("t"))	
			b=true;
		
		return b;

	}

	private String getStringData(String line){
		String[] sp=line.split(regex);	
		String[] sp2=sp[sp.length-1].split(" ");
		return sp2[sp2.length-1];
	}



public double[] loadArray(){
	String file=util.getFile();
	if(file==null || file.equals("") )  throw new NullPointerException("file not found.");
	return loadArray(file);
}

public double[] loadArray(String arrayPath){

	try{
		FileReader fr=new FileReader(arrayPath);
		BufferedReader br = new BufferedReader(fr);
		String line;
		String s;
		String[] sp;

		int N=100000;
		
		double[] x1=new double[N];
		
		int i=0;
		line=br.readLine();
		while(line!=null){
			if(i>N) break;
			x1[i++]=Double.parseDouble(line);
			line=br.readLine();
			
		}

		double[] x=Arrays.copyOf(x1, i);
		
			return x;
			
	}
	catch(IOException e){
		e.printStackTrace();//System.err.println("Error in loading model file.");
	}


	return null;
}	

public double[][] loadArrays(int n, int m,String arrayPath){

	try{
		FileReader fr=new FileReader(arrayPath);
		BufferedReader br = new BufferedReader(fr);
		String line;
		String s;
		String[] sp;

	
		
		double[][] A=new double[n][m];
		
		for(int i=0;i<n;i++){
			line=br.readLine();
			if(line==null) continue;
			double[] x=getCSV(line);
			for(int j=0;j<m;j++)
				A[i][j]=x[j];
	
			
			
		}

		
			return A;
			
	}
	catch(IOException e){
		e.printStackTrace();//System.err.println("Error in loading model file.");
	}


	return null;
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

private int[] getCSInt(String line){
	String[] sp=line.split(regex);	
	int L=sp.length;
	int[] v=new int[L];
	for( int p=0;p<L;p++)
				v[p]=Integer.parseInt(sp[p]);

	return v;
}


}
