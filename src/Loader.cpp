/*
 * Loader.cpp
 *
 *  Created on: 2015/01/28
 *      Author: Hassan
 */

#include "Loader.h"
#include "Model.h"
#include <iostream>
#include <fstream>


 Loader::Loader(){};

void Loader::loadMesh(Model *model){


	fstream fs;
	fs.open("bb.txt");

	   cout<<fs.is_open()<<endl;


	   string line;
	   string value;
	   int numline=0;
	   while(getline(fs,value,','))
	   {
	       ++numline;
	       cout << string( value, 1, value.length()-2 )<<endl;
	       cout<<line<<endl;
	   }

	/*


			FileReader fr=new FileReader(bunFilePath);
			BufferedReader br = new BufferedReader(fr);
			String line;
			String s;
			String[] sp;

			String elType=br.readLine();
			model.setElType(elType);


			br.readLine();
			line=br.readLine();
			sp=line.split(regex);
			if(!sp[0].equals(""))
				model.numberOfNodes=Integer.parseInt(sp[0]);
			else
				model.numberOfNodes=Integer.parseInt(sp[1]);


			br.readLine();
			line=br.readLine();

			sp=line.split(regex);
			if(!sp[0].equals(""))
				model.numberOfElements=Integer.parseInt(sp[0]);
			else
				model.numberOfElements=Integer.parseInt(sp[1]);

			model.element=new Element[model.numberOfElements+1];
			for(int i=1;i<=model.numberOfElements;i++){
				model.element[i]=new Element(elType);
			}


			model.node=new Node[model.numberOfNodes+1];

			for(int i=1;i<=model.numberOfNodes;i++)
				model.node[i]=new Node(model.dim);

			br.readLine();
			line=br.readLine();
			sp=line.split(regex);

			int nRegs=0;
			if(!sp[0].equals(""))
				nRegs=Integer.parseInt(sp[0]);
			else
				nRegs=Integer.parseInt(sp[1]);

			if(model.numberOfRegions<nRegs) model.numberOfRegions=nRegs;

			model.region=new Region[model.numberOfRegions+1];
			for(int i=1;i<=model.numberOfRegions;i++)
				model.region[i]=new Region(model.dim);

			br.readLine();

			line=br.readLine();


			model.scaleFactor=Double.parseDouble(line);

			double factor=1.0/model.scaleFactor;

			for(int i=1;i<=model.numberOfElements;i++){
				line=br.readLine();
				sp=line.split(regex);
				int k=0;
				for(int j=0;j<sp.length;j++){
					if(!sp[j].equals(""))
						model.element[i].setVertNumb(k++, Integer.parseInt(sp[j]));
				}
			}


			Vect z=new Vect(model.dim);



			for(int i=1;i<=model.numberOfNodes;i++){
				line=br.readLine();
				sp=line.split(regex);
				int k=0;

				for(int j=0;j<sp.length;j++)
					if(!sp[j].equals(""))
						z.el[k++]=Double.parseDouble(sp[j])*factor;


				model.node[i].setCoord(z);


				}


			model.setBounds();



				for(int i=1;i<=nRegs;i++){

					line=br.readLine();
					if(line==null) line="1,0,x"+i;
					sp=line.split(regex);
					String[] str=new String[4];
					int k=0;
					for(int j=0;j<sp.length;j++){
						if(!sp[j].equals(""))
							str[k++]=sp[j];

					}
					model.region[i].setFirstEl(Integer.parseInt(str[0]));

					model.region[i].setLastEl(Integer.parseInt(str[1]));
					model.region[i].setName(str[2]);
					model.region[i].setMaterial(str[2]);

			}


			for(int i=nRegs+1;i<=model.numberOfRegions;i++){

				model.region[i].setFirstEl(1);

					model.region[i].setLastEl(0);
					model.region[i].setName("xxx");
					model.region[i].setMaterial("xmat");

			}


				line=br.readLine();
				if(line!=null)
				model.motor=getBooleanData(line);

				if(model.motor){
					line=br.readLine();
					if(line==null) return;
					sp=line.split(regex);
					String[] str=new String[3];
					int k=0;
					for(int j=0;j<sp.length;j++){
						if(!sp[j].equals(""))
							str[k++]=sp[j];
					}

					int rotBegin,rotEnd;

					rotBegin=Integer.parseInt(str[0]);
					rotEnd=Integer.parseInt(str[1]);

					for(int ir=rotBegin;ir<=rotEnd;ir++){
						model.region[ir].rotor=true;
					model.nRotReg++;
					}
				}

				line=br.readLine();
				if(line==null|| line.equals(""))
				model.hasTwoNodeNumb=true;
				else
				model.hasTwoNodeNumb=getBooleanData(line);

				//==============

				//if(model.motor) model.hasTwoNodeNumb=true;
				//=========

			System.out.println();
			System.out.println("Loading mesh file completed.");

			br.close();
			fr.close();

			for(int ir=1;ir<=nRegs;ir++)
				for( int i=model.region[ir].getFirstEl();i<=model.region[ir].getLastEl();i++)
					model.element[i].setRegion(ir);

			model.setMaxDim();

			model.setFemCalc();



		*/}


