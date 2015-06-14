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
import math.util;
import fem.Model;



public class HysDataGraph {

	Mat[][] BH;
	Mat rotLoss;
	Mat[][] BHAni;
	int nSet;
	double[] Bs, Hs;
	int[] nInitial,nMajor,nSymLoop,nDescending,nAscending,nTotCurves,nAni;
	String regex="[:; ,\\t]+";


	public static void main(String[] args){

		HysDataGraph pg=new HysDataGraph();

		pg.loadHysData();
	}

	public boolean loadHysData(){

		/*String file=util.getFile();
	if(file==null || file.equals("") )return false;*/
		//	String file="C:\\Works\\HVID\\folder1\\data\\A_Bì¸óÕëŒèÃÉãÅ[Évhts_data\\hys_data";

		//String file="C:\\Users\\Hassan Ebrahimi\\JavaWorks\\MagFem\\hys_data";
		//String file="C:\\Works\\HVID\\hys_dataWithLoss";
		String file=System.getProperty("user.dir") + "\\hys_dataH.txt";
		


		try{
			FileReader fr=new FileReader(file);
			BufferedReader br = new BufferedReader(fr);
			String line;
			String s;
			String[] sp;

			
			line=br.readLine();
			
			sp=line.split(regex);

			this.nSet=Integer.parseInt(sp[2]);
			this.Bs=new double[this.nSet];
			this.Hs=new double[this.nSet];
			this.nInitial=new int[this.nSet];
			this.nMajor=new int[this.nSet];
			this.nSymLoop=new int[this.nSet];
			this.nDescending=new int[this.nSet];
			this.nAscending=new int[this.nSet];
			this.nTotCurves=new int[this.nSet];
			this.nAni=new int[this.nSet];

			BH=new Mat[nSet][];
			BHAni=new Mat[nSet][];
		
		
			fr=new FileReader(file);
			br = new BufferedReader(fr);

			for(int ia=0;ia<this.nSet;ia++){

		
				line=br.readLine();
				line=br.readLine();
				line=br.readLine();
				sp=line.split(regex);	

				this.Bs[ia]=Double.parseDouble(sp[0]);
				this.Hs[ia]=Double.parseDouble(sp[1]);
				line=br.readLine();
				line=br.readLine();
				sp=line.split(regex);	
				this.nInitial[ia]=Integer.parseInt(sp[0]);
				this.nMajor[ia]=Integer.parseInt(sp[1]);
				this.nSymLoop[ia]=Integer.parseInt(sp[2]);
				this.nDescending[ia]=Integer.parseInt(sp[3]);
				this.nAscending[ia]=Integer.parseInt(sp[4]);

				this.nTotCurves[ia]=this.nInitial[ia]+this.nMajor[ia]+this.nSymLoop[ia]+this.nDescending[ia]+this.nAscending[ia];
				BH[ia]=new Mat[this.nTotCurves[ia]];

				line=br.readLine();
				line=br.readLine();
				int L=Integer.parseInt(line);

				BH[ia][0]=new Mat(L,2);

				for( int p=0;p<L;p++){
					line=br.readLine();
					sp=line.split(regex);	
					double[] bh=getCSV(line);
					BH[ia][0].el[p][0]=bh[0];
					BH[ia][0].el[p][1]=bh[1];

				}



				line=br.readLine();
				int L1=0;
				for( int ip=1;ip<this.nTotCurves[ia];ip++){
					line=br.readLine();
					if(line.startsWith("*")) {line=br.readLine();};
					sp=line.split(regex);	
					//if(sp.length==1)
					L1=Integer.parseInt(sp[0]);

					BH[ia][ip]=new Mat(L1,2);

					for( int i=0;i<L1;i++){
						line=br.readLine();

						double[] bh=getCSV(line);

						BH[ia][ip].el[i][0]=bh[0];
						BH[ia][ip].el[i][1]=bh[1];
					}


				}

				line=br.readLine();
				line=br.readLine();
				line=br.readLine();

				int nLoss=Integer.parseInt(line);
				rotLoss=new Mat(nLoss,2);
				line=br.readLine();

				for( int i=0;i<nLoss;i++){
					line=br.readLine();
					double[] BL=getCSV(line);
					rotLoss.el[i][0]=BL[0];
					rotLoss.el[i][1]=BL[1];
				}

				line=br.readLine();
				line=br.readLine();
				line=br.readLine();
				double[] dd=this.getCSV(line);

				nAni[ia]=(int)dd[0];
				int Lani=(int)dd[1];
				line=br.readLine();

				BHAni[ia]=new Mat[nAni[ia]];

				for( int i=0;i<nLoss;i++)
					BHAni[ia][i]=new Mat(Lani,3);

				for( int i=0;i<Lani;i++){
					line=br.readLine();
					double[] BL=getCSV(line);
					double B=BL[0];
					for(int j=0;j<nAni[ia];j++){
						BHAni[ia][j].el[i][0]=B;
						BHAni[ia][j].el[i][1]=BL[j+1];
					}
				}

				line=br.readLine();
				for( int i=0;i<Lani;i++){
					line=br.readLine();
					double[] BL=getCSV(line);
					for(int j=0;j<nAni[ia];j++){
						BHAni[ia][j].el[i][2]=BL[j+1];
					}
				}
				
				line=br.readLine();
				line=br.readLine();
				line=br.readLine();

			}
			
		

			util.plotBunch(BH[0]);
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
