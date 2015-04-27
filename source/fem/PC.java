package fem;

import static java.lang.Math.*;

import java.awt.Color;
import java.util.Arrays;

import math.SpMat;
import math.Vect;
import math.util;

public class PC {
	
	public Vect xp;
	
	public PC(){}
	
	public static void main(String[] args){
		new Main();
	}
	
	
	public void runPC(Model model, Main main){



		model.resultFolder=System.getProperty("user.dir") + "//PCResults";
		String folder=model.resultFolder;
		main.prepare();
	

			model.solver.terminate(false);


			
					String ff=folder+"\\bun.txt";
	
					//model.writeMesh(ff);
			

/*		String logFilePath = model.resultFolder+ "\\log.txt";
		main.simGui.writeLog(logFilePath);*/
		main.gui.Run.setBackground(Color.green);

	}
	
	 static double logGamma(double x) {
	      double tmp = (x - 0.5) * Math.log(x + 4.5) - (x + 4.5);
	      double ser = 1.0 + 76.18009173    / (x + 0)   - 86.50532033    / (x + 1)
	                       + 24.01409822    / (x + 2)   -  1.231739516   / (x + 3)
	                       +  0.00120858003 / (x + 4)   -  0.00000536382 / (x + 5);
	      return tmp + Math.log(ser * Math.sqrt(2 * Math.PI));
	   }
	 
	  static double gamma(double x) { 
		  return Math.exp(logGamma(x)); 
		  }

}


