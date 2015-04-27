package materialData;

import java.io.FileReader;
import java.util.Scanner;

import math.Vect;
import math.util;



public class CurrentWaveForm {
	
		public double[][] TI;
		public int length;
		boolean periodic=true;

	
	public CurrentWaveForm(String f) throws Exception{
		double[][] TI1=new double[20000][2];
		
		String file = System.getProperty("user.dir") + "\\"+f;

		int i;
		// try{
		 Scanner scr=new Scanner(new FileReader(file));
		 while(scr.hasNext()){
	
		 while(!scr.next().equals("begin")){}
		 int j=0;
		 String s=scr.next();
		 while(!s.equals("end")){
		 TI1[j][0]=Double.parseDouble(s);
		 s=scr.next();
		 TI1[j++][1]=Double.parseDouble(s);
		 s=scr.next();
		 }
		this.length=j;
		 }
			
		 this.TI=new double[this.length][2];
		 for(int k=0;k<this.length;k++)
		 this.TI[k]=TI1[k];
		 
		 scr.close();

		}
	
	
	public CurrentWaveForm(double[][] ti) {
		TI=new double[ti.length][2];
		this.length=ti.length;

		 for(int k=0;k<this.length;k++){
		 this.TI[k][0]=ti[k][0];
		 this.TI[k][1]=ti[k][1];
		 }
		 

		}
			
	
	
	public static void main(String[] args) throws Exception{
		
/*		BHCurve BH=new BHCurve("H350");
		
		Curve cv=new Curve(BH,600,600);
		cv.show(true);*/

CurrentWaveForm Ix=new CurrentWaveForm("emf//Ian.txt");
/*CurrentWaveForm Iy=new CurrentWaveForm("vb0.txt");
CurrentWaveForm Iz=new CurrentWaveForm("vc0.txt");
		*/
		Curve cv2=new Curve(Ix,600,600);
		cv2.show(true);
		//util.show(Ix.TI);
		for(int k=0;k<Ix.TI.length;k++)
			util.pr(Ix.TI[k][1]);

		double[][] iy=new double[Ix.length][2];
		
		Vect  dv=new Vect(1800);

		for(int k=0;k<0*dv.length;k++){
		
			double t=.0025*0+k*5.55555e-6;
			//iy[k][1]=Ix.getI(iy[k][0])+Iy.getI(iy[k][0])+Iz.getI(iy[k][0]);
			dv.el[k]=Ix.getI(t);
			util.pr(k*5.55555e-6+"\t"+dv.el[k]);
			//util.pr(Ix.TI[k+1800][1]);
		}
		
		//dv.show();
	//util.plot(dv);
		
		//util.plot(dv);


		/*String bh = System.getProperty("user.dir") + "\\bh.xls";
		BH2.writexls(bh);
	*/
	}
	
	public  double getI(double t){


		if(t<this.TI[0][0]|| t>this.TI[length-1][0]){
			if(!this.periodic) throw new IllegalArgumentException("out of time span.");
			else{
				 t=t-Math.floor(t/this.TI[length-1][0])*this.TI[length-1][0];
			}
		}

		int j=getj(t);
		 double i=(this.TI[j][1]*(this.TI[j+1][0]-t)+this.TI[j+1][1]*(t-this.TI[j][0]))/(this.TI[j+1][0]-this.TI[j][0]);
	
		 return i;
	}
	
	public int getj(double t){
		
		int j=0;
		if(this.TI.length==1) return j;
		while(this.TI[j+1][0]<t){j++;}
		return j;
	}
	

}
