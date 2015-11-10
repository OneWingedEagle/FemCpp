package io;


import java.awt.Color;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.Arrays;

import math.util;




	public class BLAN implements ActionListener{

	Cons gui;

	private String path = System.getProperty("user.dir");
	public  boolean console=true,dated=false;

	public BLAN()
	{		
		this.gui=new Cons(this.path);
		this.gui.Run.addActionListener(this);
		this.gui.setVisible(true);
		//this.gui.Run.doClick();
	}	
	
	public static void main(String[] args){
	
		new BLAN();
	}

		
	public void run(){

		gui.inputFile=gui.tfInput.getText();
		
		int numbMesh=Integer.parseInt(gui.tfNumb_Mesh.getText());

		int ix=0;
		String input=gui.inputFile;
		String[] lines=new String[1000];
		String[] sp=null;
		int nlin=0;

		try{
			BufferedReader br = new BufferedReader(new FileReader(input));
			String line;
			line=br.readLine();
		
			
			
			lines[ix]=line;
			while(line!=null){
				line=br.readLine();
				lines[ix]=line;
				ix++;
				
			}
			
			br.close();
			System.out.println();
			System.out.println("Input file was read.");
			
			 String regex="[:; ,=\\t]+";

			 for(int i=0;i<ix-1;i++){
				if(lines[i].startsWith("* INPUT_MESH_FILE"))
				{
					
					nlin=i+1;
					sp=lines[nlin].split(regex);
			
		
					int ib=3;
					sp[ib]=Integer.toString(numbMesh);
	

				}
			}

				
		}

		catch(IOException e){
			e.printStackTrace();
		}

		
		//input=gui.inputFile+"2";

	
		try{
			
			
			PrintWriter pr=new PrintWriter(new BufferedWriter(new FileWriter(input)));		
			for(int i=0;i<ix;i++){
				if(i!=nlin)
				pr.println(lines[i]);
				else{
					pr.print("\t");
					for(int k=0;k<sp.length;k++)
						pr.print(sp[k]+"\t");
					pr.println();
				}
			}
			
			pr.close();
			
			System.out.println("Modified input file was written as "+input);
		}
		catch(IOException e){
			e.printStackTrace();
		}


		String batchString="cmd /c start "+"run-mpi"+numbMesh+".bat";

try {
	Runtime.getRuntime().exec(batchString);
} catch (IOException e) {
	// TODO Auto-generated catch block
	e.printStackTrace();
}


		
		System.exit(0);
		
		/*
		
		if(console){
		Console.redirectOutput(this.gui.progressArea);
		}
		this.model.meshFilePath=this.gui.tfMeshFile.getText();
		this.model.dataFilePath=this.gui.tfDataFile.getText();

		this.errMax=Double.parseDouble(this.gui.tfErrorMax.getText());
		this.iterMax=Integer.parseInt(this.gui.tfIterMax.getText());
		this.thread=new Thread(){
			long m1 = Runtime.getRuntime().totalMemory() - Runtime.getRuntime().freeMemory();
			@Override
			public void run(){

				
				BLAN.this.model.loadMesh(model.meshFilePath);
			
				
				prepare();
			
				model.loadData(model.dataFilePath);
	
				
				if(model.seepage){
					runGeo(); 
					}
				if(model.wavePC){
					runPC(); 
					}
				else if(model.magAnalysis && !model.mechAnalysis) {

					runMag(); 

					}
				else if(!model.magAnalysis && model.mechAnalysis) {

					runMech(); 

					}

		
				String logFilePath = model.resultFolder+ "\\log.txt";
				gui.writeLog(logFilePath);
				gui.Run.setBackground(Color.green);
				gui.Run.setEnabled(true);

			}
		};
		this.thread.start();		
	*/}

	

	@Override
	public void actionPerformed(ActionEvent e)
	{	

		 if(e.getSource()==this.gui.Run){
			this.gui.Run.setBackground(Color.gray);
			this.gui.Run.setEnabled(false);

			run();
			//loadMotor();
		}


	
	}





}

