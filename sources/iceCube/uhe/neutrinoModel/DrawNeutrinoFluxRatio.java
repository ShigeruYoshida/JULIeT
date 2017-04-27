package iceCube.uhe.neutrinoModel;

import iceCube.uhe.neutrinoModel.*;
import iceCube.uhe.particles.*;
import java.io.*;

public class DrawNeutrinoFluxRatio {

   public static void main(String[] args) throws IOException {

       int model1 = 1;
       int model2 = 1;
       int particleID = 0;
       double I_10PeV = 0.0;
       double epsilon = 1.0e-8;
 
       if(args.length<2){
	   System.out.println("Usage: DrawNeutrinoFlux model-parameter1  model parameter2 particleID ");
	   System.exit(0);
       }if(args.length == 3){
	   model1 = Integer.valueOf(args[0]).intValue();
	   model2 = Integer.valueOf(args[1]).intValue();
	   particleID = Integer.valueOf(args[2]).intValue();
       }

       NeutrinoFlux neutFlux1 = new NeutrinoFlux(model1);
       NeutrinoFlux neutFlux2 = new NeutrinoFlux(model2);

       System.out.println("titx Energy [GeV]");
       System.out.println("tity Ratio");
       System.out.println("scal 1.0e6 1.0e12 0.01 10.0");

       double logE = 6.0;
       while(logE<12.0){
	   double EFlux1 = neutFlux1.getEFluxwzOsci(logE,particleID);
	   double EFlux2 = neutFlux2.getEFluxwzOsci(logE,particleID);
	   double ratio = EFlux1/EFlux2;
	   double neutrinoEnergy = Math.pow(10.0,logE);

	   System.out.println("data " + neutrinoEnergy + " 0.0 " + ratio + " 0.0");
	   logE += Particle.getDeltaLogEnergy()*10.0;
       }
       System.out.println("logx");
       System.out.println("logy");
       System.out.println("join");
       System.out.println("disp");
       System.out.println("endg");

   }
}
