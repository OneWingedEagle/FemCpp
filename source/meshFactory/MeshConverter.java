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



public class MeshConverter {

	String regex="[ : ,=\\t]+";
	BoundarySet bs=new BoundarySet();
	Writer writer=new Writer();

	public static void main(String[] args){
		MeshConverter mf=new MeshConverter();
		
	}


	public void hexToPrism()
	{
		hexToPrism(2);
	}
	
	public void hexToPrism(int dir)
	
	{
		String bun=util.getFile();
	if(bun==null || bun.equals("") )return;
	Model model=new Model();
	model.loadMesh(bun);

	String prismMesh = System.getProperty("user.dir") + "//prismEl.txt";

	if(model.elCode==4)
		model.writeMeshHexaToPrism(prismMesh,dir);

	}
	
	public void hexaToPyramid()
	
	{
		String bun=util.getFile();
	if(bun==null || bun.equals("") )return;
	Model model=new Model();
	model.loadMesh(bun);

	String pyrMesh = System.getProperty("user.dir") + "//pyrEl.txt";

	if(model.elCode==4)
		writer.writeMeshHexaToPyramid(model,pyrMesh);
	}
	
public void hexaToTetra()
	
	{
		String bun=util.getFile();
	if(bun==null || bun.equals("") )return;
	Model model=new Model();
	model.loadMesh(bun);

	String pyrMesh = System.getProperty("user.dir") + "//tetEl.txt";

	if(model.elCode==4)
		writer.writeMeshHexaToTetra(model,pyrMesh);
	}
	
	public void reRegionGroupEls(Model model){
	
	MeshFactory mf=new MeshFactory();
	
	mf.reRegionGroupEls(model);
}
	
	
	
}
