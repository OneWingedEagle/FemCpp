package meshFactory;

import static java.lang.Math.PI;
import static java.lang.Math.abs;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStreamReader;
import java.io.PrintWriter;
import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashSet;
import java.util.List;
import java.util.Set;








import fem.BoundarySet;
import fem.Element;
import fem.Model;
import fem.Node;
import fem.Region;
import io.Loader;
import io.Writer;
import math.Mat;
import math.SpBlockMat;
import math.SpMat;
import math.Vect;
import math.util;



public class MeshFormatConverter {

	String regex="[ : ,=\\t]+";
	BoundarySet bs=new BoundarySet();
	Writer writer=new Writer();

	public static void main(String[] args){
		MeshFormatConverter mfc=new MeshFormatConverter();
		//mfc.getPreHexAtlasOrig();
		
		
		//mfc.getPostHexAtlasOrig();
		
		//mfc.getPostHexAtlas();
		mfc.getPostHexaNeu();
		
		
		
		
	}




	public void getNeuMeshQ(int mode){
		String regex="[ ,\\t]+";
		String s=util.getFile();
		try{
			File f=new File(s);
			FileReader fr=new FileReader(f);
			BufferedReader br = new BufferedReader(fr);
			String line="";
			String[] sp=new String[15];

	if(mode==0)
	while(!br.readLine().startsWith("<<< End Solid Transmit <<<")){		}

	for(int i=1;i<100000;i++)
			{

				line=br.readLine();
				sp=line.split(regex);
				if(sp.length==15 && !sp[0].equals("0")) break;
			
				
			}
			
	

			int nnMax=0,nn=1;
			Vect[] coord1=new Vect[1000000];
			int[] map=new int[1000000];
			int nx=0;
			for(int i=1;i<1000000;i++)
			{
				sp=line.split(regex);
				
				line=br.readLine();

				if(sp.length!=15) break;

			
				nn=Integer.parseInt(sp[0]);
				nx++;
			
				map[nn]=nx;
		

				coord1[nx]=new Vect(Double.parseDouble(sp[11]),Double.parseDouble(sp[12]),Double.parseDouble(sp[13]));
				

			}
			
			nnMax=nx;


			int[][] vernumb=new int[10*nnMax+1][4];
			int[] nReg=new int[1000000+1];
			
			for(int i=0;i<nReg.length;i++)
				nReg[i]=-1;
				
			int ix=0;
		
				line=br.readLine();
			//	line=br.readLine();
			

			sp=new String[13];
			for(int i=1;i<=vernumb.length;i++)
			{
				line=br.readLine();

				sp=line.split(regex);
				if(sp.length<5) break;
				ix++;
				nReg[ix]=Integer.parseInt(sp[2]);
				line=br.readLine();
		
				sp=line.split(regex);
				
				line=br.readLine();
				line=br.readLine();
				line=br.readLine();
				line=br.readLine();
				line=br.readLine();

				vernumb[ix][0]=map[Integer.parseInt(sp[0])];
				vernumb[ix][1]=map[Integer.parseInt(sp[1])];
				vernumb[ix][2]=map[Integer.parseInt(sp[2])];
				vernumb[ix][3]=map[Integer.parseInt(sp[3])];

				for(int j=0;j<4;j++){
					if(vernumb[ix][j]==0) {
						if(j>0) vernumb[ix][j]=vernumb[ix][j-1];
						//ix--;
						break;
					}
				}

			}
			
			int nEl=ix;

		
			
			List<Integer> list1=new ArrayList<Integer>();
			for(int i=1;i<nReg.length;i++){
				if(nReg[i]!=-1)
					list1.add(nReg[i]);
			}
		
				Set<Integer> set = new HashSet<Integer>(list1);
				
				ArrayList<Integer> regNumbs = new ArrayList<Integer>(set);
				
				int nRegions=regNumbs.size();
				
				util.pr(nRegions);

				
				int[] regNumber=new int[nRegions+1];
				for(int ir=1;ir<=nRegions;ir++)
					regNumber[ir]=regNumbs.get(ir-1);
				

					
			int[] elOrd=new int[nEl+1];
			int[][]regEnd=new int[nRegions+1][2];
			
			 nx=0;
			for(int ir=1;ir<=nRegions;ir++)
			{
				regEnd[ir][0]=nx+1;
				for(int i=1;i<=nEl;i++)
					if(nReg[i]==regNumber[ir])
						elOrd[++nx]=i;
				regEnd[ir][1]=nx;
				
			}


			int nNodes=nnMax;
			double scaleFactor=1;
			double scx=1;
			DecimalFormat formatter;
			if(scaleFactor==1)
				formatter= new DecimalFormat("0.000000000");
			else 
				formatter= new DecimalFormat("0.000000");
			

			try{


				String fout=System.getProperty("user.dir")+"\\EMSol\\quad.txt";
				PrintWriter pwBun = new PrintWriter(new BufferedWriter(new FileWriter(fout)));		

				pwBun.println("quadrangle");
				pwBun.println("//Number_of_Node");
				pwBun.println(nNodes);

				pwBun.println("//Number_of_Element");
				pwBun.println(nEl);

				pwBun.println("//Number_of_Region");
				pwBun.println(nRegions);
				pwBun.println("//Factor");
				pwBun.println(scx);

				for(int ir=1;ir<=nRegions;ir++)
					for(int i=regEnd[ir][0];i<=regEnd[ir][1];i++){
						for(int j=0;j<4;j++){
							int mpn=vernumb[elOrd[i]][j];
							pwBun.print(mpn+",");
						}
					
						pwBun.println();
					}

			Vect v;
			for(int i=1;i<=nnMax;i++){ 
				if(coord1[i]==null)
					v=new Vect(0,0);
				else
					v=coord1[i].deepCopy();
						for(int j=0;j<2;j++){
							pwBun.print(formatter.format(v.el[j]*scaleFactor)+" ,");
						}
					pwBun.println();	

				}

				for(int ir=1;ir<=nRegions;ir++)
					pwBun.println(regEnd[ir][0]+","+regEnd[ir][1]+","+"region"+ir);
				
				pwBun.close();
				br.close();
				fr.close();
				
				System.out.println();
				System.out.println(" Bun data was written to:");
				System.out.println("    "+fout);
			}
			catch(IOException e){ System.err.println("error");e.printStackTrace(); }

		
		
		}
		
		

		catch(Exception e){System.err.println("error");	e.printStackTrace(); }
	}
	
	
	
	public void getNeuMeshHexa(){
		String regex="[ ,\\t]+";
		String s=util.getFile();
		try{
			File f=new File(s);
			FileReader fr=new FileReader(f);
			BufferedReader br = new BufferedReader(fr);
			String line="";
			String linep="";
			String[] sp=new String[15];
			
			

	while(!br.readLine().startsWith("<<< End Solid Transmit <<<")){		}

	for(int i=1;i<100000;i++)
			{

				line=br.readLine();
				sp=line.split(regex);
				if(sp.length==15 && !sp[0].equals("0")) break;
			
				
			}
			
	

			int nnMax=0,nn=1;
			Vect[] coord1=new Vect[1000000];
			int[] map=new int[1000000];
			int nx=0;
			for(int i=1;i<1000000;i++)
			{
				sp=line.split(regex);
				
				line=br.readLine();

				if(sp.length!=15) break;

			
				nn=Integer.parseInt(sp[0]);
				nx++;
			
				map[nn]=nx;
		

				coord1[nx]=new Vect(Double.parseDouble(sp[11]),Double.parseDouble(sp[12]),Double.parseDouble(sp[13]));
		

			}
			
			nnMax=nx;


			int[][] vernumb=new int[10*nnMax+1][8];
			int[] nReg=new int[1000000+1];
			
			for(int i=0;i<nReg.length;i++)
				nReg[i]=-1;
				
			int ix=0;
		
				line=br.readLine();
			//	line=br.readLine();
	for(int i=0;i<7;i++)
		line=br.readLine();

			sp=new String[13];
			for(int i=1;i<=vernumb.length;i++)
			{
				
				linep=br.readLine();
				line=br.readLine();
	
				
				sp=line.split(regex);
	
				if(sp.length<5) break;
				
				if(sp[6].equals("0")) {
					for(int j=0;j<5;j++)
						line=br.readLine();
					continue;
				}
				

				
				sp=linep.split(regex);
				ix++;
				nReg[ix]=Integer.parseInt(sp[2]);
			
			
				sp=line.split(regex);


				for(int j=0;j<8;j++){

					vernumb[ix][j]=map[Integer.parseInt(sp[j])];
				}
				
			

				for(int j=0;j<8;j++){
					if(vernumb[ix][j]==0) {
						if(j>0) vernumb[ix][j]=vernumb[ix][j-1];
						ix--;
						break;
					}
				}
				
	
				
				for(int j=0;j<5;j++)
					line=br.readLine();

			}
			
			int nEl=ix;

		
			
			List<Integer> list1=new ArrayList<Integer>();
			for(int i=1;i<nReg.length;i++){
				if(nReg[i]!=-1)
					list1.add(nReg[i]);
			}
		
				Set<Integer> set = new HashSet<Integer>(list1);
				
				ArrayList<Integer> regNumbs = new ArrayList<Integer>(set);
				
				int nRegions=regNumbs.size();
				
				util.pr(nRegions);

				
				int[] regNumber=new int[nRegions+1];
				for(int ir=1;ir<=nRegions;ir++)
					regNumber[ir]=regNumbs.get(ir-1);
				

					
			int[] elOrd=new int[nEl+1];
			int[][]regEnd=new int[nRegions+1][2];
			
			 nx=0;
			for(int ir=1;ir<=nRegions;ir++)
			{
				regEnd[ir][0]=nx+1;
				for(int i=1;i<=nEl;i++)
					if(nReg[i]==regNumber[ir])
						elOrd[++nx]=i;
				regEnd[ir][1]=nx;
				
			}


			int nNodes=nnMax;
			double scaleFactor=1;
			double scx=1;
			DecimalFormat formatter;
			if(scaleFactor==1)
				formatter= new DecimalFormat("0.000000000");
			else 
				formatter= new DecimalFormat("0.000000");
			

			try{


				String fout=System.getProperty("user.dir")+"\\EMSol\\Hexa.txt";

				PrintWriter pwBun = new PrintWriter(new BufferedWriter(new FileWriter(fout)));		

				pwBun.println("hexahedron");
				pwBun.println("//Number_of_Node");
				pwBun.println(nNodes);

				pwBun.println("//Number_of_Element");
				pwBun.println(nEl);

				pwBun.println("//Number_of_Region");
				pwBun.println(nRegions);
				pwBun.println("//Factor");
				pwBun.println(scx);

				for(int ir=1;ir<=nRegions;ir++)
					for(int i=regEnd[ir][0];i<=regEnd[ir][1];i++){
						for(int j=0;j<8;j++){
							int mpn=vernumb[elOrd[i]][j];
							pwBun.print(mpn+",");
						}
					
						pwBun.println();
					}

			Vect v;
			for(int i=1;i<=nnMax;i++){ 
				if(coord1[i]==null)
					v=new Vect(0,0,0);
				else
					v=coord1[i].deepCopy();
						for(int j=0;j<3;j++){
							pwBun.print(formatter.format(v.el[j]*scaleFactor)+" ,");
						}
					pwBun.println();	

				}

				for(int ir=1;ir<=nRegions;ir++)
					pwBun.println(regEnd[ir][0]+","+regEnd[ir][1]+","+"region"+ir);
				
				pwBun.close();
				br.close();
				fr.close();
				
				System.out.println();
				System.out.println(" Bun data was written to:");
				System.out.println("    "+fout);
			}
			catch(IOException e){ System.err.println("error");e.printStackTrace(); }

		
		
		}
		
		

		catch(Exception e){System.err.println("error");	e.printStackTrace(); }
	}
	
	
	public void getPostMeshQ(){
		String regex="[ ,\\t]+";
		String s=util.getFile();
		try{
			File f=new File(s);
			FileReader fr=new FileReader(f);
			BufferedReader br = new BufferedReader(fr);
			String line="";
			String[] sp=new String[15];
			
/*	while(!br.readLine().startsWith("<<< End Solid Transmit <<<")){
	
			}*/

	for(int i=1;i<1000;i++)
			{

				line=br.readLine();
				sp=line.split(regex);
				//if(sp.length==15 && !sp[0].equals("0")) break;
				if(sp.length>5) break;
			
				
			}
			
	

			int nnMax=0,nn=1;
			Vect[] coord1=new Vect[1000000];
			//int[] map=new int[1000000];
			int nx=0;
			for(int i=1;i<1000000;i++)
			{
				sp=line.split(regex);
				
				line=br.readLine();

				if(sp.length!=14) break;

			
				nn=Integer.parseInt(sp[0]);
				nx++;
			
			//	map[nn]=nx;
		

				coord1[nn]=new Vect(Double.parseDouble(sp[11]),Double.parseDouble(sp[12]),Double.parseDouble(sp[13]));
				

			}
			
			nnMax=nx;


			int[][] vernumb=new int[10*nnMax+1][4];
			int[] nReg=new int[1000000+1];
			
			for(int i=0;i<nReg.length;i++)
				nReg[i]=-1;
				
			int ix=0;
		
				line=br.readLine();
			//	line=br.readLine();
			

			sp=new String[13];
			for(int i=1;i<=vernumb.length;i++)
			{
				line=br.readLine();

				sp=line.split(regex);
				if(sp.length<5) break;
				ix++;
				nReg[ix]=Integer.parseInt(sp[2]);
				line=br.readLine();
		
				sp=line.split(regex);
				
				line=br.readLine();
				line=br.readLine();
				line=br.readLine();
				line=br.readLine();
				line=br.readLine();

			/*	vernumb[ix][0]=map[Integer.parseInt(sp[0])];
				vernumb[ix][1]=map[Integer.parseInt(sp[1])];
				vernumb[ix][2]=map[Integer.parseInt(sp[2])];
				vernumb[ix][3]=map[Integer.parseInt(sp[3])];*/
				
				vernumb[ix][0]=Integer.parseInt(sp[0]);
				vernumb[ix][1]=Integer.parseInt(sp[1]);
				vernumb[ix][2]=Integer.parseInt(sp[2]);
				vernumb[ix][3]=Integer.parseInt(sp[3]);

				for(int j=0;j<4;j++){
					if(vernumb[ix][j]==0) {
						if(j>0) vernumb[ix][j]=vernumb[ix][j-1];
						//ix--;
						break;
					}
				}

			}
			
			int nEl=ix;

		
			
			List<Integer> list1=new ArrayList<Integer>();
			for(int i=1;i<nReg.length;i++){
				if(nReg[i]!=-1)
					list1.add(nReg[i]);
			}
		
				Set<Integer> set = new HashSet<Integer>(list1);
				
				ArrayList<Integer> regNumbs = new ArrayList<Integer>(set);
				
				int nRegions=regNumbs.size();
				
				util.pr(nRegions);

				
				int[] regNumber=new int[nRegions+1];
				for(int ir=1;ir<=nRegions;ir++)
					regNumber[ir]=regNumbs.get(ir-1);
				

					
			int[] elOrd=new int[nEl+1];
			int[][]regEnd=new int[nRegions+1][2];
			
			 nx=0;
			for(int ir=1;ir<=nRegions;ir++)
			{
				regEnd[ir][0]=nx+1;
				for(int i=1;i<=nEl;i++)
					if(nReg[i]==regNumber[ir])
						elOrd[++nx]=i;
				regEnd[ir][1]=nx;
				
			}


			int nNodes=nnMax;
			double scaleFactor=1;
			double scx=1;
			DecimalFormat formatter;
			if(scaleFactor==1)
				formatter= new DecimalFormat("0.000000000");
			else 
				formatter= new DecimalFormat("0.000000");
			

			try{


				String fout=System.getProperty("user.dir")+"\\quadExtracted.txt";
				PrintWriter pwBun = new PrintWriter(new BufferedWriter(new FileWriter(fout)));		

				pwBun.println("quadrangle");
				pwBun.println("//Number_of_Node");
				pwBun.println(nNodes);

				pwBun.println("//Number_of_Element");
				pwBun.println(nEl);

				pwBun.println("//Number_of_Region");
				pwBun.println(nRegions);
				pwBun.println("//Factor");
				pwBun.println(scx);

				for(int ir=1;ir<=nRegions;ir++)
					for(int i=regEnd[ir][0];i<=regEnd[ir][1];i++){
						for(int j=0;j<4;j++){
							int mpn=vernumb[elOrd[i]][j];
							pwBun.print(mpn+",");
						}
					
						pwBun.println();
					}

			Vect v;
			for(int i=1;i<=nnMax;i++){ 
				if(coord1[i]==null)
					v=new Vect(0,0);
				else
					v=coord1[i].deepCopy();
						for(int j=0;j<2;j++){
							pwBun.print(formatter.format(v.el[j]*scaleFactor)+" ,");
						}
					pwBun.println();	

				}

				for(int ir=1;ir<=nRegions;ir++)
					pwBun.println(regEnd[ir][0]+","+regEnd[ir][1]+","+"region"+ir);
				
				pwBun.close();
				br.close();
				fr.close();
				
				System.out.println();
				System.out.println(" Bun data was written to:");
				System.out.println("    "+fout);
			}
			catch(IOException e){ System.err.println("error");e.printStackTrace(); }

		
		
		}
		
		

		catch(Exception e){System.err.println("error");	e.printStackTrace(); }
	}
	
	
	public void getPostMeshHex(){
		String regex="[ ,\\t]+";
		String s=util.getFile();
		try{
			File f=new File(s);
			FileReader fr=new FileReader(f);
			BufferedReader br = new BufferedReader(fr);
			String line="";
			String linep="";
			String[] sp=new String[15];

			int numbAddedNodes=0;
	//if(mode==0)
	//while(!br.readLine().startsWith("<<< End Solid Transmit <<<")){		}

			for(int i=0;i<7;i++)
				line=br.readLine();
			
			line=br.readLine();

			int nnMax=0,nn=1;
			Vect[] coord1=new Vect[1000000];
			int[] map=new int[1000000];
			int nx=0;
			for(int i=1;i<1000000;i++)
			{
				
				
				
				line=br.readLine();

				sp=line.split(regex);
			
				if(sp.length<14) break;

			
				nn=Integer.parseInt(sp[0]);
			//	nx++;
			nx=nn;
				
			map[nn]=nx;
			
		

				coord1[nx]=new Vect(Double.parseDouble(sp[11]),Double.parseDouble(sp[12]),Double.parseDouble(sp[13]));
		

			}
			
			
			nnMax=0;
			for(int k=0;k<map.length;k++)
				if(map[k]>nnMax) nnMax=map[k];


			int[][] vernumb=new int[10*nnMax+1][8];
			int[] nReg=new int[1000000+1];
			int firstElement=0;
			
			for(int i=0;i<nReg.length;i++)
				nReg[i]=-1;
				
			int ix=0;
		int nEl=0;
	for(int i=0;i<2;i++)
		line=br.readLine();

			sp=new String[13];
			for(int i=1;i<=vernumb.length;i++)
			{
				
				linep=br.readLine();

				line=br.readLine();
				if(line==null) break;
					
				sp=line.split(regex);
	
				if(sp.length<5) break;
				
				if(sp[6].equals("0")) {
					for(int j=0;j<5;j++)
						line=br.readLine();
					continue;
				}
				

				
				sp=linep.split(regex);
		
				ix=Integer.parseInt(sp[0]);
				if(firstElement==0) firstElement=ix;
				
				if(ix>nEl) nEl=ix;

				nReg[ix]=Integer.parseInt(sp[2]);
			
				sp=line.split(regex);


				for(int j=0;j<8;j++){

					vernumb[ix][j]=map[Integer.parseInt(sp[j])];
				}
				
			

				for(int j=0;j<8;j++){
					if(vernumb[ix][j]==0) {

						numbAddedNodes++;
				
						 vernumb[ix][j]=nnMax+numbAddedNodes;//vernumb[ix][j-1];
							if(j>=0){
						 coord1[numbAddedNodes]=coord1[vernumb[ix][j-1]].deepCopy().times(0);
							}
						//ix--;
						//break;
					}
				}
				
				// change order
				
				int[] tmp=new int[8];
				for(int j=0;j<8;j++){
					tmp[j]=vernumb[ix][j];
				
				}
				vernumb[ix][0]=tmp[0];
				vernumb[ix][1]=tmp[3];
				vernumb[ix][2]=tmp[2];
				vernumb[ix][3]=tmp[1];
				
				vernumb[ix][4]=tmp[4];
				vernumb[ix][5]=tmp[7];
				vernumb[ix][6]=tmp[6];
				vernumb[ix][7]=tmp[5];

				//======
				
				for(int j=0;j<5;j++)
					line=br.readLine();

			}
			
		

		
			
			List<Integer> list1=new ArrayList<Integer>();
			for(int i=1;i<=nEl;i++){
				if(nReg[i]==-1){
					nReg[i]=1;
					for(int j=0;j<8;j++)
						vernumb[i][j]=vernumb[firstElement][j];
				}
				else
					nReg[i]+=1;
				
					list1.add(nReg[i]);

			}

				Set<Integer> set = new HashSet<Integer>(list1);
				
				ArrayList<Integer> regNumbs = new ArrayList<Integer>(set);
				
				int nRegions=regNumbs.size();
				
				util.pr(nRegions);

				
				int[] regNumber=new int[nRegions+1];
				for(int ir=1;ir<=nRegions;ir++)
					regNumber[ir]=regNumbs.get(ir-1);

					
			int[] elOrd=new int[nEl+1];
			int[][]regEnd=new int[nRegions+1][2];
			
			 nx=0;
			for(int ir=1;ir<=nRegions;ir++)
			{
	
				regEnd[ir][0]=nx+1;
				for(int i=1;i<=nEl;i++)
					if(nReg[i]==regNumber[ir]){
						elOrd[i]=i;
						//elOrd[++nx]=i;
			
				nx++;
					}
				regEnd[ir][1]=nx;
				
				util.hshow(regEnd[ir]);
			}

			int nNodes=nnMax+numbAddedNodes;
			double scaleFactor=1;
			double scx=1;
			DecimalFormat formatter;
			if(scaleFactor==1)
				formatter= new DecimalFormat("0.000000000");
			else 
				formatter= new DecimalFormat("0.000000");
			

			try{


				String fout=System.getProperty("user.dir")+"\\EMSol\\Hexa.txt";
				PrintWriter pwBun = new PrintWriter(new BufferedWriter(new FileWriter(fout)));		

				pwBun.println("hexahedron");
				pwBun.println("//Number_of_Node");
				pwBun.println(nNodes);

				pwBun.println("//Number_of_Element");
				pwBun.println(nEl);

				pwBun.println("//Number_of_Region");
				pwBun.println(nRegions);
				pwBun.println("//Factor");
				pwBun.println(scx);

				for(int ir=1;ir<=nRegions;ir++)
					for(int i=regEnd[ir][0];i<=regEnd[ir][1];i++){
						for(int j=0;j<8;j++){
							int mpn=vernumb[elOrd[i]][j];
							pwBun.print(mpn+",");
						}
					
						pwBun.println();
					}

			Vect v;
			for(int i=1;i<=nnMax+numbAddedNodes;i++){ 
				if(coord1[i]==null)
					v=new Vect(0,0,0);
				else
					v=coord1[i].deepCopy();
						for(int j=0;j<3;j++){
							pwBun.print(formatter.format(v.el[j]*scaleFactor)+" ,");
						}
					pwBun.println();	

				}

				for(int ir=1;ir<=nRegions;ir++)
					pwBun.println(regEnd[ir][0]+","+regEnd[ir][1]+","+"region"+ir);
				
				pwBun.close();
				br.close();
				fr.close();
				
				System.out.println();
				System.out.println(" Bun data was written to:");
				System.out.println("    "+fout);
			}
			catch(IOException e){ System.err.println("error");e.printStackTrace(); }

		
		
		}
		
		

		catch(Exception e){System.err.println("error");	e.printStackTrace(); }
	}
	
	public void getPreHexAtlas(){

		String s=util.getFile();
		
		String line;
		int max=1000000;
		Vect[] coord1=new Vect[max];
		int[][] vertNumb=new int[max][8];
		int[] nodeMap=new int[max];

		int[] elMap=new int[max];
		int[] regNumb=new int[max];
		
		int[] nRegEls=new int[15];



		int nnx=1;
		int nex=1;
		
		try{
			File f=new File(s);
			FileReader fr=new FileReader(f);
			BufferedReader br = new BufferedReader(fr);

	
			while(true){
				line=br.readLine();
				
				if(line==null) break;
				if(line.startsWith("GRID")){
					nnx=readAtlasNodes(br,nodeMap,coord1,nnx);


				}
				
				
				if(line.startsWith("CONC")){	
		
					nex=readAtlasElements(br,elMap,regNumb,nRegEls,vertNumb,nex);
				}

			}
		
		}
			catch(Exception e){System.err.println("error");	e.printStackTrace(); }

			int nNodes=0;
			for(int i=0;i<nodeMap.length;i++)
				if(nodeMap[i]>0) nNodes++;


			int nEls=0;
		
			int nRegions=0;
			
			for(int i=0;i<nRegEls.length;i++)
				if(nRegEls[i]>0) {
					nRegions++;
					nEls+=nRegEls[i];
				}
			
		
			
			int[] regMap=new int[nRegEls.length];
			int nr=0;
			int[][]regEnd=new int[nRegions+1][2];
			regEnd[0][0]=1;
					
			for(int i=0;i<nRegEls.length;i++){
				
				if(nRegEls[i]>0) {
					nr++;
					regMap[i]=nr;
					regEnd[nr][0]=regEnd[nr-1][1]+1;
				
					regEnd[nr][1]=regEnd[nr][0]+nRegEls[i]-1;
				}
			}
			
			
			Model model=new Model(nRegions,nEls,nNodes,"hexahedron");
				
			for(int ir=1;ir<=nRegions;ir++)
			{
				model.region[ir].setFirstEl(regEnd[ir][0]);;
				model.region[ir].setLastEl(regEnd[ir][1]);;
				model.region[ir].setName("region"+regMap[ir]);
			for(int i=model.region[ir].getFirstEl();i<=model.region[ir].getLastEl();i++){
				for(int j=0;j<8;j++){
					int mpn=nodeMap[vertNumb[i][j]];
					model.element[i].setVertNumb(j,mpn);
					
				}

			}
			
		
			}
			
			for(int i=1;i<=nNodes;i++){
				model.node[i].setCoord(coord1[i]);
				}
		
			model.scaleFactor=1;
			
			for(int i=1;i<=nEls;i++)
			model.element[i].setRegion(regMap[regNumb[i]]);
			
			reRegionGroupEls(model);
			
			String fout=System.getProperty("user.dir")+"\\EMSol\\Hexa.txt";
			
			model.writeMesh(fout);



	}
	
	
	public void getPostHexAtlas(){

		String s=util.getFile();
		
		String line;
		int max=1000000;
		Vect[] coord1=new Vect[max];
		int[][] vertNumb=new int[max][8];
		int[] nodeMap=new int[max];

		int[] elMap=new int[max];
		int[] regNumb=new int[max];
		
		int[] nRegEls=new int[1000];


		int nElMax=0;

		int nnx=1;
		int nex=1;
		
		try{
			File f=new File(s);
			FileReader fr=new FileReader(f);
			BufferedReader br = new BufferedReader(fr);

	
			while(true){
				line=br.readLine();
				
				if(line==null) break;
				if(line.startsWith("GRID")){
					nnx=readAtlasNodes(br,nodeMap,coord1,nnx);


				}
				
				
				if(line.startsWith("CONC")){	
		
					nex=readAtlasElements(br,elMap,regNumb,nRegEls,vertNumb,nex);
				}

			}
		
		}
			catch(Exception e){System.err.println("error");	e.printStackTrace(); }

		
			int nNodes=0;
			for(int i=0;i<nodeMap.length;i++)
				if(nodeMap[i]>0) nNodes++;
			
			for(int i=1;i<elMap.length;i++)
				if(elMap[i]>nElMax) nElMax=elMap[i];


			int nEls=0;
		
			int nRegions=0;
			
			for(int i=0;i<nRegEls.length;i++)
				if(nRegEls[i]>0) {
					nRegions++;
					nEls+=nRegEls[i];
				}
			
		
			
			int[] regMap=new int[nRegEls.length];
			int nr=0;
			int[][]regEnd=new int[nRegions+1][2];
			regEnd[0][0]=1;
					
			for(int i=0;i<nRegEls.length;i++){
				
				if(nRegEls[i]>0) {
					nr++;
					regMap[i]=nr;
					regEnd[nr][0]=regEnd[nr-1][1]+1;
				
					regEnd[nr][1]=regEnd[nr][0]+nRegEls[i]-1;
				}
			}
			
			
			Model model=new Model(nRegions,nEls,nNodes,"hexahedron");
				
			for(int ir=1;ir<=nRegions;ir++)
			{
				model.region[ir].setFirstEl(regEnd[ir][0]);;
				model.region[ir].setLastEl(regEnd[ir][1]);;
				model.region[ir].setName("region"+regMap[ir]);
			for(int i=model.region[ir].getFirstEl();i<=model.region[ir].getLastEl();i++){
				for(int j=0;j<8;j++){
					int mpn=nodeMap[vertNumb[i][j]];
					model.element[i].setVertNumb(j,mpn);
					
				}

			}
			
		
			}
			
			for(int i=1;i<=nNodes;i++){
				model.node[i].setCoord(coord1[i]);
				}
		
			model.scaleFactor=1;
			
			for(int i=1;i<=nEls;i++){
			model.element[i].setRegion(regMap[regNumb[i]]);
			}
			
			int[] elMapReorderd=reRegionGroupEls(model);
			

			
			int[] elMapReorderdRev=new int[elMapReorderd.length];
			for(int i=1;i<elMapReorderd.length;i++){
				elMapReorderdRev[elMapReorderd[i]]=i;
	
			}
	
			
			String fout=System.getProperty("user.dir")+"\\EMSol\\Hexa.txt";

			model.writeMesh(fout);
			
			String map=System.getProperty("user.dir")+"\\EMSol\\map.txt";
			
			try {
				PrintWriter pwBun = new PrintWriter(new BufferedWriter(new FileWriter(map)));
				
				for(int i=1;i<elMap.length;i++)
					if(elMap[i]>0)
						pwBun.println(i+"\t"+elMapReorderdRev[elMap[i]]);
				
				pwBun.close();
			} catch (IOException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}	
			




	}
	
	public void getPostHexaNeu(){


		String file=util.getFile();
		
		String line;
		int max=300000;

		Vect[] coord1=new Vect[max];
		int[][] vertNumb=new int[max][8];
		int[] nodeMap=new int[max];

		int[] elMap=new int[max];
		int[] regNumb=new int[max];
		
		
		int[] nRegEls=new int[1000];



		int nElMax=0;

		int nnx=1;
		int nex=1;
		
		try{
			File f=new File(file);
			FileReader fr=new FileReader(f);
			BufferedReader br = new BufferedReader(fr);

	
	
			while(true){
	
				line=br.readLine();
			
				if(line==null) break;
		
				if(first(line).equals("403")){
		
					nnx=readNeuNodes(br,nodeMap,coord1,nnx);

				}

				
				if(first(line).equals("404")){

					nex=readNeuElements(br,elMap,regNumb,nRegEls,vertNumb,nex);


				}

			}

				
	
		}
			catch(Exception e){System.err.println("error");	e.printStackTrace(); }


			int nNodes=nnx-1;

	
			
			for(int i=1;i<elMap.length;i++)
				if(elMap[i]>nElMax) nElMax=elMap[i];


			int nEls=0;
		
			int nRegions=0;
			
			for(int i=0;i<nRegEls.length;i++)
				if(nRegEls[i]>0) {
					nRegions++;
					nEls+=nRegEls[i];
				}
			
		
			
			int[] regMap=new int[nRegEls.length];
			int nr=0;
			int[][]regEnd=new int[nRegions+1][2];
			regEnd[0][0]=1;
					
			for(int i=0;i<nRegEls.length;i++){
		
				if(nRegEls[i]>0) {
				
					nr++;
					regMap[i]=nr;

					regEnd[nr][0]=regEnd[nr-1][1]+1;
				
					regEnd[nr][1]=regEnd[nr][0]+nRegEls[i]-1;
				}
			}

			int[] regMapRev=new int[regMap.length];
			for(int ir=1;ir<regMap.length;ir++)
				regMapRev[regMap[ir]]=ir;

			
			Model model=new Model(nRegions,nEls,nNodes,"hexahedron");
			
			for(int ir=1;ir<=nRegions;ir++)
			{
				model.region[ir].setFirstEl(regEnd[ir][0]);;
				model.region[ir].setLastEl(regEnd[ir][1]);;
				model.region[ir].setName("region"+regMapRev[ir]);
			for(int i=model.region[ir].getFirstEl();i<=model.region[ir].getLastEl();i++){
				for(int j=0;j<8;j++){
					int mpn=nodeMap[vertNumb[i][j]];
					model.element[i].setVertNumb(j,mpn);
					
				}

			}
			
		
			}
			

			
			for(int i=1;i<=model.numberOfNodes;i++){

				model.node[i].setCoord(coord1[i]);
				}
		
			model.scaleFactor=1;
						
			for(int i=1;i<=nEls;i++){
			model.element[i].setRegion(regMap[regNumb[i]]);
			}
			
			int[] elMapReorderd=reRegionGroupEls(model);
			
			int[] elMapReorderdRev=new int[elMapReorderd.length];
			for(int i=1;i<elMapReorderd.length;i++){
				elMapReorderdRev[elMapReorderd[i]]=i;

			}
	
	
			
			String fout=System.getProperty("user.dir")+"\\EMSol\\Hexa.txt";

			model.writeMesh(fout);
			
			String map=System.getProperty("user.dir")+"\\EMSol\\map.txt";
			
			try {
				PrintWriter pwBun = new PrintWriter(new BufferedWriter(new FileWriter(map)));
				
				for(int i=1;i<elMap.length;i++)
					if(elMap[i]>0)
						pwBun.println(i+"\t"+elMapReorderdRev[elMap[i]]);
				
				pwBun.close();
			} catch (IOException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}	
			



	
	}
	
	public void getPreHexAtlasOrig(){

		String s=util.getFile();
		
		String line;
		int max=1000000;
		Vect[] coord1=new Vect[max];
		int[][] vertNumb=new int[max][8];
		int[] nodeNumb=new int[max];
		int[] elNumb=new int[max];
		int[] regNumb=new int[max];


		int nNodes=0,nEls=0;

	
		for(int i=0;i<regNumb.length;i++)
			regNumb[i]=-1;


		try{
			File f=new File(s);
			FileReader fr=new FileReader(f);
			BufferedReader br = new BufferedReader(fr);

	
			while(true){
				line=br.readLine();
				
				if(line==null) break;
				if(line.startsWith("GRID")){
					nNodes=readAtlasNodesOrig(br,nodeNumb,coord1,nNodes);


				}
				

				
				if(line.startsWith("CONC")){	
		
					nEls=readAtlasElementsOrig(br,elNumb,regNumb,vertNumb,nEls);
				}

			}
		
		}
			catch(Exception e){System.err.println("error");	e.printStackTrace(); }

			int nNodeMax=0;
			for(int i=0;i<nNodes;i++)
				if(nodeNumb[i]>nNodeMax) nNodeMax=nodeNumb[i];
			
			int[] vertNumb0=new int[8];
			int jx=0;
			for(int i=0;i<=nNodes;i++){
				if(nodeNumb[i]>0) vertNumb0[jx++]=nodeNumb[i];
				if(jx==8) break;
			}


			int nElMax=0;
			for(int i=0;i<=nEls;i++)
				if(elNumb[i]>=nElMax) nElMax=elNumb[i];
			
	
			
			int nRegMax=0;
			for(int i=0;i<regNumb.length;i++)
				if(regNumb[i]>nRegMax) nRegMax=regNumb[i];
			
			
			int[] nRegEls=new int[nRegMax+1];
			for(int i=0;i<regNumb.length;i++)
				if(regNumb[i]>0) nRegEls[regNumb[i]]++;
						
			int[][]regEnd=new int[nRegMax+1][2];
			regEnd[0][0]=1;

					
			for(int i=0;i<nRegEls.length;i++){
		
				if(i>0)
				regEnd[i][0]=regEnd[i-1][1]+1;
				
					regEnd[i][1]=regEnd[i][0]+nRegEls[i]-1;
				}		
			

			
			Model model=new Model(nRegMax,nElMax,nNodeMax,"hexahedron");
			
			for(int i=1;i<=nElMax;i++){
				if(vertNumb[i][0]>0){
					model.element[i].setVertNumb(vertNumb[i]);
				}
				else{
					model.element[i].setVertNumb(vertNumb0);
				}
			}
				
			
		util.pr(nElMax);
		util.pr(model.numberOfElements);		
			for(int ir=1;ir<=nRegMax;ir++)
			{
				model.region[ir].setFirstEl(regEnd[ir][0]);;
				model.region[ir].setLastEl(regEnd[ir][1]);;
				model.region[ir].setName("region"+ir);
			for(int i=model.region[ir].getFirstEl();i<=model.region[ir].getLastEl();i++){
				

				model.element[i].setRegion(ir);

			}
			
		
			}
		
			util.pr(nRegMax);
			
			for(int i=1;i<=nNodeMax;i++){
				if(coord1[i]==null)
					model.node[i].setCoord(new Vect(3));
				else
				model.node[i].setCoord(coord1[i]);
				}
		
			model.scaleFactor=1;
			
			reRegionGroupEls(model);
			

			
			String fout=System.getProperty("user.dir")+"\\EMSol\\HexaOrig.txt";
			
			model.writeMesh(fout);



	}
	
	public void getPostHexAtlasOrig(){

		String s=util.getFile();
		
		String line;
		int max=100000;
		Vect[] coord1=new Vect[max];
		int[][] vertNumb=new int[max][8];
		int[] nodeNumb=new int[max];
		int[] elNumb=new int[max];
		int[] regNumb=new int[max];


		int nNodes=0,nEls=0;

	
		for(int i=0;i<regNumb.length;i++)
			regNumb[i]=-1;


		try{
			File f=new File(s);
			FileReader fr=new FileReader(f);
			BufferedReader br = new BufferedReader(fr);

	
			while(true){
				line=br.readLine();
				
				if(line==null) break;
				if(line.startsWith("GRID")){
					nNodes=readAtlasNodesOrig(br,nodeNumb,coord1,nNodes);


				}
				

				
				if(line.startsWith("CONC")){	
		
					nEls=readAtlasElementsOrig(br,elNumb,regNumb,vertNumb,nEls);
				}

			}
		
		}
			catch(Exception e){System.err.println("error");	e.printStackTrace(); }

			int nNodeMax=0;
			for(int i=0;i<nNodes;i++)
				if(nodeNumb[i]>nNodeMax) nNodeMax=nodeNumb[i];
			
			int[] vertNumb0=new int[8];
			int jx=0;
			for(int i=0;i<nNodes;i++){
				if(nodeNumb[i]>0) vertNumb0[jx++]=nodeNumb[i];
				if(jx==8) break;
			}

			int nElMax=0;
			for(int i=0;i<elNumb.length;i++)
				if(elNumb[i]>nElMax) nElMax=elNumb[i];
			
			int nRegMax=0;
			for(int i=0;i<regNumb.length;i++)
				if(regNumb[i]>nRegMax) nRegMax=regNumb[i];
			
			
			int[] nRegEls=new int[nRegMax+1];
			for(int i=0;i<regNumb.length;i++)
				if(regNumb[i]>0) nRegEls[regNumb[i]]++;
						
			int[][]regEnd=new int[nRegMax+1][2];
			regEnd[0][0]=1;

					
			for(int i=0;i<nRegEls.length;i++){
		
				if(i>0)
				regEnd[i][0]=regEnd[i-1][1]+1;
				
					regEnd[i][1]=regEnd[i][0]+nRegEls[i]-1;
				}		
			

			
			Model model=new Model(nRegMax,nElMax,nNodeMax,"hexahedron");
			
			for(int i=1;i<=nElMax;i++){
				if(vertNumb[i][0]>0){
					model.element[i].setVertNumb(vertNumb[i]);
				}
				else{
					model.element[i].setVertNumb(vertNumb0);
				}
			}
				
				
			for(int ir=1;ir<=nRegMax;ir++)
			{
				model.region[ir].setFirstEl(regEnd[ir][0]);;
				model.region[ir].setLastEl(regEnd[ir][1]);;
				model.region[ir].setName("region"+ir);
			for(int i=model.region[ir].getFirstEl();i<=model.region[ir].getLastEl();i++){
				/*for(int j=0;j<8;j++){
					int mpn=vertNumb[i][j];

					model.element[i].setVertNumb(j,mpn);
					
				}*/
				model.element[i].setRegion(ir);

			}
			
		
			}
		
	
			
			for(int i=1;i<=nNodeMax;i++){
				if(coord1[i]==null)
					model.node[i].setCoord(new Vect(3));
				else
				model.node[i].setCoord(coord1[i]);
				}
		
			model.scaleFactor=1;
			
			reRegionGroupEls(model);
			

			
			String fout=System.getProperty("user.dir")+"\\EMSol\\HexaOrig.txt";
			
			model.writeMesh(fout);



	}
	
	
	
	public int readAtlasNodes(BufferedReader br,int[] nodeMap,Vect[] coord1,int nNode){
		String line;
		String[] sp;
		try{
		line=br.readLine();
		int nn;
		int nIdent=1;
	
		while(true){
			if(line==null) break;

			sp=line.split(regex);
		
			if(sp.length<3) return nNode;

			nn=Integer.parseInt(sp[nIdent]);

			if(nn==-1) break;
			
			nodeMap[nn]=nNode;
		

			coord1[nNode]=new Vect(Double.parseDouble(sp[nIdent+1]),Double.parseDouble(sp[nIdent+2]),Double.parseDouble(sp[nIdent+3]));
	
			line=br.readLine();
			
			nNode++;

		}
		}
		catch(Exception e){System.err.println("error");	e.printStackTrace(); }
		
		return nNode;
		
	}
	
	public int readNeuNodes(BufferedReader br,int[] nodeMap,Vect[] coord1,int nNode){
		String line;
		String[] sp;
		try{
	
			int nmax1=4*coord1.length/5;
			int nmax2=coord1.length-nmax1;

		line=br.readLine();
		
		if(line==null) return nNode;

		sp=line.split(regex);
		
		int nn;
		int nIdent=0;
		
		while(!first(line).equals("-1")){
		
			nn=Integer.parseInt(sp[nIdent]);
			
			if(nn>nmax1) nn=nmax1+nn%nmax2;
			
			nodeMap[nn]=nNode;
		
			coord1[nNode]=new Vect(Double.parseDouble(sp[nIdent+11]),Double.parseDouble(sp[nIdent+12]),Double.parseDouble(sp[nIdent+13]));
	
			line=br.readLine();
			sp=line.split(regex);
			
			nNode++;


		}
		}
		catch(Exception e){System.err.println("error");	e.printStackTrace(); }
		
		return nNode;
		
	}
	
	
	
public int readNeuElements(BufferedReader br,int[] elMap,int[]regNumb,int[] nRegEls,int[][] vertNumb,int nEl){
		
	
	int nmax1=4*elMap.length/5;
	int nmax2=elMap.length-nmax1;
	
	
		String line1="",line2;
		String[] sp1,sp2;
		try{

		int nIdent=0;

		
		int ne;
		line1=br.readLine();
		line2=br.readLine();
		
		while(!first(line1).equals("-1")){
		
			sp1=line1.split(regex);
			sp2=line2.split(regex);

			int[] vertNumb1=new int[8];
			 boolean def=false;

			for(int j=0;j<8;j++){

				vertNumb1[j]=Integer.parseInt(sp2[nIdent+j]);
				
			
				
				if(vertNumb1[j]==0){
					def=true;
					break;
				}
				
				if(vertNumb1[j]>nmax1) vertNumb1[j]=nmax1+vertNumb1[j]%nmax2;
			
			}
			

		if(!def){
			
			for(int j=0;j<8;j++)
				vertNumb[nEl][j]=vertNumb1[j];
						
			
				ne=Integer.parseInt(sp1[nIdent]);
				if(ne>nmax1) ne=nmax1+ne%nmax2;
				
				
				elMap[ne]=nEl;
							
				regNumb[nEl]=Integer.parseInt(sp1[nIdent+2]);
			

				nRegEls[regNumb[nEl]]++;

			
			// change order
			
			int[] tmp=new int[8];
			for(int j=0;j<8;j++){
				tmp[j]=vertNumb[nEl][j];
			
			}
			vertNumb[nEl][0]=tmp[0];
			vertNumb[nEl][1]=tmp[3];
			vertNumb[nEl][2]=tmp[2];
			vertNumb[nEl][3]=tmp[1];
			
			vertNumb[nEl][4]=tmp[4];
			vertNumb[nEl][5]=tmp[7];
			vertNumb[nEl][6]=tmp[6];
			vertNumb[nEl][7]=tmp[5];
			
			nEl++;

		
	}

		
			br.readLine();
			br.readLine();
			br.readLine();
			br.readLine();
			br.readLine();
			
			line1=br.readLine();
			line2=br.readLine();

	
		}


		}
		
	
		catch(Exception e){System.err.println("error");	e.printStackTrace(); }
		
		return nEl;
	}
	
	public int readAtlasElements(BufferedReader br,int[] elMap,int[]regNumb,int[] nRegEls,int[][] vertNumb,int nEl){
		
		String line1,line2;
		String[] sp1,sp2;
		try{
	//	line=br.readLine();
		int nIdent=1;

		int ne;
		while(true){
			
			line1=br.readLine();
		
			sp1=line1.split(regex);
			if(sp1.length<7) {return nEl;}
			

			line2=br.readLine();
			if(line2==null) return nEl;
			
			sp2=line2.split(regex);


			ne=Integer.parseInt(sp1[nIdent]);
			elMap[ne]=nEl;
			
			regNumb[nEl]=Integer.parseInt(sp1[nIdent+2]);


			nRegEls[regNumb[nEl]]++;

			for(int j=0;j<5;j++){

				vertNumb[nEl][j]=Integer.parseInt(sp1[nIdent+j+3]);
			}
			
	
		
			for(int j=0;j<3;j++){

				vertNumb[nEl][j+5]=Integer.parseInt(sp2[nIdent+j]);
			}
			
/*
			for(int j=0;j<8;j++)
			{



					numbAddedNodes++;
			
					 vernumb[ix][j]=nnMax+numbAddedNodes;//vernumb[ix][j-1];
						if(j>=0){
					 coord1[numbAddedNodes]=coord1[vernumb[ix][j-1]].deepCopy().times(0);
						}
					//ix--;
					//break;
				}
			}*/
	
			
			// change order
			
			int[] tmp=new int[8];
			for(int j=0;j<8;j++){
				tmp[j]=vertNumb[nEl][j];
			
			}
			vertNumb[nEl][0]=tmp[0];
			vertNumb[nEl][1]=tmp[3];
			vertNumb[nEl][2]=tmp[2];
			vertNumb[nEl][3]=tmp[1];
			
			vertNumb[nEl][4]=tmp[4];
			vertNumb[nEl][5]=tmp[7];
			vertNumb[nEl][6]=tmp[6];
			vertNumb[nEl][7]=tmp[5];
			
			nEl++;

		}

		}
		
	
		catch(Exception e){System.err.println("error");	e.printStackTrace(); }
		
		return nEl;
	}
	
	public int readAtlasNodesOrig(BufferedReader br,int[] nodeNumb,Vect[] coord1, int nNode){
		String line;
		String[] sp;
		try{
		line=br.readLine();
		int nn;
		int nIdent=1;
	
		while(true){
			if(line==null) return nNode;

			sp=line.split(regex);

			

			nn=Integer.parseInt(sp[nIdent]);

			if(nn==-1) return nNode;
			
			nNode++;
			
			nodeNumb[nNode]=nn;

			coord1[nn]=new Vect(Double.parseDouble(sp[nIdent+1]),Double.parseDouble(sp[nIdent+2]),Double.parseDouble(sp[nIdent+3]));
	
			line=br.readLine();
			
		
		}
		}
		catch(Exception e){System.err.println("error");	e.printStackTrace(); }
		
		return nNode;
		
	}
	
	public int readAtlasElementsOrig(BufferedReader br,int[] elNumb,int[]regNumb,int[][] vertNumb, int nEl){
		
		String line1,line2;
		String[] sp1,sp2;
		try{
	//	line=br.readLine();
		int nIdent=1;

		int ne;
		while(true){
			
			line1=br.readLine();
		
			sp1=line1.split(regex);
			if(sp1.length<7) { return nEl;}

			line2=br.readLine();
			if(line2==null) return nEl;
			
			sp2=line2.split(regex);


			ne=Integer.parseInt(sp1[nIdent]);
			
			nEl++;
			
			elNumb[nEl]=ne;
			

			
			regNumb[ne]=Integer.parseInt(sp1[nIdent+2]);


			for(int j=0;j<5;j++){

				vertNumb[ne][j]=Integer.parseInt(sp1[nIdent+j+3]);
			}
			
	
		
			for(int j=0;j<3;j++){

				vertNumb[ne][j+5]=Integer.parseInt(sp2[nIdent+j]);
			}
			
/*
			for(int j=0;j<8;j++)
			{



					numbAddedNodes++;
			
					 vernumb[ix][j]=nnMax+numbAddedNodes;//vernumb[ix][j-1];
						if(j>=0){
					 coord1[numbAddedNodes]=coord1[vernumb[ix][j-1]].deepCopy().times(0);
						}
					//ix--;
					//break;
				}
			}*/
	
			
			// change order
			
			int[] tmp=new int[8];
			for(int j=0;j<8;j++){
				tmp[j]=vertNumb[ne][j];
			
			}
			vertNumb[ne][0]=tmp[0];
			vertNumb[ne][1]=tmp[3];
			vertNumb[ne][2]=tmp[2];
			vertNumb[ne][3]=tmp[1];
			
			vertNumb[ne][4]=tmp[4];
			vertNumb[ne][5]=tmp[7];
			vertNumb[ne][6]=tmp[6];
			vertNumb[ne][7]=tmp[5];
				
		}

		}
		
	
		catch(Exception e){System.err.println("error");	e.printStackTrace(); }
	
		return nEl;
	}
	
	
	public void getEMSolFlux(String bbf,int dim, int numb){

		String regex="[ ,\\t]+";
		try{

			File f=new File(bbf);
			FileReader fr=new FileReader(f);
			BufferedReader br = new BufferedReader(fr);
			String line="";
			String[] sp=new String[15];
		
			Mat BB=new Mat(100*1000,dim);
			
			for(int tt=0;tt<numb;tt++){

				line="";
				while(!line.startsWith("STEP")){
					line=br.readLine();
					
					}
				int[] nx=new int[dim];;

				String fout=System.getProperty("user.dir")+"\\EMSol\\flux"+tt+".txt";
				
				PrintWriter pwBun = new PrintWriter(new BufferedWriter(new FileWriter(fout)));	
				
	while(!line.startsWith("BMAG")){
	line=br.readLine();
			}
	
	for(int k=0;k<6;k++)
		line=br.readLine();
	
	int k=0;

	
			for(int i=1;i<180000;i++){
				
	
		if(line.startsWith("-1")){
			k++;
			for(int j=0;j<8;j++)
				line=br.readLine();

		}
		else

		sp=line.split(regex);
			

		double Bu=Double.parseDouble(sp[1]);
		BB.el[nx[k]][k]=Bu;

				line=br.readLine();
	
				nx[k]++;
				
			

				if(k==dim-1 && nx[k]==nx[k-1]){
					break;
				}
				
				

	}



	int Ne=nx[0];
	pwBun.println("flux");
	pwBun.println(dim);
	pwBun.println(Ne);
	
	for(int j=0;j<Ne;j++){
		for(int p=0;p<dim;p++)
		pwBun.print(BB.el[j][p]+"\t");
		
		pwBun.println();
	}
	
	System.out.println("Flux was written to "+fout);
	pwBun.close();
	
			}
		
	br.close();
	fr.close();

	
	
		}


		catch(Exception e){System.err.println("error");	e.printStackTrace(); }
		
		

	
	
	
	}
	
	
	
	
	
	public void getNeuMeshTri(){
		String regex="[ ,\\t]+";
		String s=util.getFile();
		try{
			File f=new File(s);
			FileReader fr=new FileReader(f);
			BufferedReader br = new BufferedReader(fr);
			String line="";
			String[] sp=new String[15];
			
			while(!br.readLine().startsWith("<<< End Solid Transmit <<<")){
				
			}
		
			for(int i=1;i<10000;i++)
			{

				line=br.readLine();

				sp=line.split(regex);

				if(sp.length==15 && !sp[0].equals("0")) break;
			
				
			}
			


			int nnMax=0,nn=1;
			Vect[] coord1=new Vect[1000000];
			int[] map=new int[1000000];
			int nx=0;
			for(int i=1;i<1000000;i++)
			{
				sp=line.split(regex);
		
				line=br.readLine();

				if(sp.length!=15) break;

			
				nn=Integer.parseInt(sp[0]);
				nx++;
			
				map[nn]=nx;
		

				coord1[nx]=new Vect(Double.parseDouble(sp[11]),Double.parseDouble(sp[12]),Double.parseDouble(sp[13]));
				

			}
			
			nnMax=nx;


			int[][] vernumb=new int[10*nnMax+1][3];
			int[] nReg=new int[1000000+1];
			
			for(int i=0;i<nReg.length;i++)
				nReg[i]=-1;
				
			int ix=0;
		
				line=br.readLine();
			//	line=br.readLine();
			

			sp=new String[12];
			for(int i=1;i<=vernumb.length;i++)
			{
				line=br.readLine();

				sp=line.split(regex);
				if(sp.length<5) break;
				ix++;
				nReg[ix]=Integer.parseInt(sp[2]);
				line=br.readLine();
		
				sp=line.split(regex);
				
				line=br.readLine();
				line=br.readLine();
				line=br.readLine();
				line=br.readLine();
				line=br.readLine();

				vernumb[ix][0]=map[Integer.parseInt(sp[0])];
				vernumb[ix][1]=map[Integer.parseInt(sp[1])];
				vernumb[ix][2]=map[Integer.parseInt(sp[2])];

				

				for(int j=0;j<3;j++){
					if(vernumb[ix][j]==0) {
						ix--;
						break;
					}
				}

			}
			
			int nEl=ix;


			
			List<Integer> list1=new ArrayList<Integer>();
			for(int i=1;i<nReg.length;i++){
				if(nReg[i]!=-1)
					list1.add(nReg[i]);
			}
		
				Set<Integer> set = new HashSet<Integer>(list1);
				
				ArrayList<Integer> regNumbs = new ArrayList<Integer>(set);
				
				int nRegions=regNumbs.size();
				
				
				int[] regNumber=new int[nRegions+1];
				for(int ir=1;ir<=nRegions;ir++)
					regNumber[ir]=regNumbs.get(ir-1);
				

					
			int[] elOrd=new int[nEl+1];
			int[][]regEnd=new int[nRegions+1][2];
			
			 nx=0;
			for(int ir=1;ir<=nRegions;ir++)
			{
				regEnd[ir][0]=nx+1;
				for(int i=1;i<=nEl;i++)
					if(nReg[i]==regNumber[ir])
						elOrd[++nx]=i;
				regEnd[ir][1]=nx;
				
			}


			int nNodes=nnMax;
			double scaleFactor=1000.0;
			DecimalFormat formatter;
			if(scaleFactor==1)
				formatter= new DecimalFormat("0.000000000");
			else 
				formatter= new DecimalFormat("0.000000");
			

			try{


				String fout=System.getProperty("user.dir")+"\\\\EMSol\\tri.txt";
				PrintWriter pwBun = new PrintWriter(new BufferedWriter(new FileWriter(fout)));		

				pwBun.println("triangle");
				pwBun.println("//Number_of_Node");
				pwBun.println(nNodes);

				pwBun.println("//Number_of_Element");
				pwBun.println(nEl);

				pwBun.println("//Number_of_Region");
				pwBun.println(nRegions);
				pwBun.println("//Factor");
				pwBun.println(scaleFactor);

				for(int ir=1;ir<=nRegions;ir++)
					for(int i=regEnd[ir][0];i<=regEnd[ir][1];i++){
						for(int j=0;j<3;j++){
							int mpn=vernumb[elOrd[i]][j];
							pwBun.print(mpn+",");
						}
					
						pwBun.println();
					}

			Vect v;
			for(int i=1;i<=nnMax;i++){ 
				if(coord1[i]==null)
					v=new Vect(0,0);
				else
					v=coord1[i].deepCopy();
						for(int j=0;j<2;j++){
							pwBun.print(formatter.format(v.el[j]*scaleFactor)+" ,");
						}
					pwBun.println();	

				}

				for(int ir=1;ir<=nRegions;ir++)
					pwBun.println(regEnd[ir][0]+","+regEnd[ir][1]+","+"region"+ir);
				
				pwBun.close();
				br.close();
				fr.close();
				
				System.out.println();
				System.out.println(" Bun data was written to:");
				System.out.println("    "+fout);
			}
			catch(IOException e){ System.err.println("error");}

		
		
		}
		
		

		catch(Exception e){System.err.println("error");	e.printStackTrace(); }
	}
	
	


public void getFluxAtlas(int dim, int nElMax){

	String regex="[ ,\\t]+";
	String s=util.getFile();
	try{
		
		String fout=System.getProperty("user.dir")+"\\EMSolFlux.txt";
		PrintWriter pwBun = new PrintWriter(new BufferedWriter(new FileWriter(fout)));		

		File f=new File(s);
		FileReader fr=new FileReader(f);
		BufferedReader br = new BufferedReader(fr);
		String line="";
		String[] sp=new String[15];
		
		
		Mat BB=new Mat(nElMax,dim);
		
		int[] nx=new int[dim];;
		
		
while(!line.startsWith("STEP")){
	line=br.readLine();
	pwBun.println(line);
		}
for(int k=0;k<6;k++){
	line=br.readLine();
	pwBun.println(line);
}

int k=0;


		for(int i=1;i<100000;i++){
			

/*		if(line.startsWith("-1")){
		k++;
		for(int j=0;j<8;j++){
			line=br.readLine();
			pwBun.println(line);
		}

	}
	else*/
	
	sp=line.split(regex);
		

	double Bu=Double.parseDouble(sp[1]);
	String line2=sp[0]+",\t"+sp[1]+",";
//	BB.el[nx[k]][k]=Bu;*/
	pwBun.println(line);
			line=br.readLine();
		
		//	util.pr(nx[0]+" - "+Bu);

		//	if(nx[k]<5) BB.el[nx[k]][k]=nx[k];
			
		//	nx[k]++;
			
		
/*				if(line.startsWith("-1")){
				pwBun.println(line);
				
				pwBun.println("  -1");
				break;
			}
			*/
			if(line==null || line.length()==0) break;
			

}


	
br.close();
fr.close();
pwBun.close();

System.out.println("Flux was written to "+fout);

	}


	catch(Exception e){System.err.println("error");	e.printStackTrace(); }
	

}

public int[] reRegionGroupEls(Model model){
	
	MeshFactory mf=new MeshFactory();
	
	return mf.reRegionGroupEls(model);
}

public String first(String line){
	
	String[] sp=line.split(regex);
	int b=0;
	while(b<sp.length-1 &&sp[b].equals("")){b++;}
	
	return sp[b];
}
	
	
	
}
