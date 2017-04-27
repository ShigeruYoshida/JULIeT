package iceCube.uhe.neutrinoModel;

import iceCube.uhe.neutrinoModel.*;
import iceCube.uhe.particles.*;
import java.io.*;

public class DrawNeutrinoFlux {

   public static void main(String[] args) throws IOException {

       int model = 1;
       int particleID = 0;
       double I_10PeV = 0.0;
       double epsilon = 1.0e-8;
       boolean plotDifferential = false;
       boolean normalized = false;
       boolean energyDifferential = false;
 
       if(args.length<2){
            System.out.println("Usage: DrawNeutrinoFlux model-parameter particleID (differential?(log diif yes-1 yes with normalized to be unity-2 df/de-3 no 0))");
            System.exit(0);
       }else if (args.length>=2) {
	   model = Integer.valueOf(args[0]).intValue();
	   particleID = Integer.valueOf(args[1]).intValue();
       }if(args.length == 3){
	   if(Integer.valueOf(args[2]).intValue() != 0){
	       plotDifferential = true;
	       System.err.println(" plot differential spectrum");
	       if(Integer.valueOf(args[2]).intValue() == 2){
		   normalized = true;
		   System.err.println(" normalized to be unity : like PDF");
	       }else if(Integer.valueOf(args[2]).intValue() == 3){
		   energyDifferential = true; 
		   plotDifferential = false;
	       }
       
	   }
       }

       NeutrinoFlux neutFlux = new NeutrinoFlux(model);

       if(!plotDifferential){
	   if(!energyDifferential){
	       System.out.println("titx Log E[GeV]");
	       System.out.println("tity log (Flux E^2! [GeV cm^-2! sec^-1! sr^-1!])");
	       System.out.println("scal 6.0 12.0 -14.0 -4.0");
	   }else{
	       System.out.println("titx Energy [eV]");
	       System.out.println("tity dF/dE [/(cm^2! sec^1! sr^1! GeV)])");
	       System.out.println("scal 1.0e4 1.0e20 1.0e-30 1.0e14");
	       System.out.println("gwin 0.4 1.4 0.1 1.2");
	       System.out.println("logx");
	       System.out.println("logy");
	   }

	   double logE = 6.0;
	   while(logE<12.0){
	       double EFlux = 0.0;
	       if(particleID<4) EFlux = neutFlux.getEFlux(logE,particleID);
	       else{
		   for(int id = 1; id<=3; id++) EFlux += neutFlux.getEFlux(logE,id);
	       }
	       if(!energyDifferential){
		   double logEFlux;
		   if(EFlux > 0.0){
		       logEFlux = Math.log(EFlux)/Math.log(10.0);
		   }else{
		       logEFlux = -14.0;
		   }
		   System.out.println("data " + logE + " 0.0 " + logEFlux + " 0.0");
	       }else{
		   double energyInGeV = Math.pow(10.0, logE);
		   double energyIneV = Math.pow(10.0, logE+9.0);
		   double flux = EFlux/(energyInGeV*energyInGeV);
		   System.out.format("data %e 0.0 %e 0.0\n",energyIneV,flux);

	       }
	       logE += 0.1;
	   }

	   if(!energyDifferential){
	       logE = 7.0;
	       if(particleID<4){
		   while(logE<12.0){
		       I_10PeV += neutFlux.getDFDLogE(logE,particleID)*Particle.getDeltaLogEnergy( );
		       logE += Particle.getDeltaLogEnergy( );
		   }
		   
		   System.out.println("mssg I(10PeV)=" + I_10PeV);
	       }
	   }


	   System.out.println("join");
	   System.out.println("disp");
	   System.out.println("endg");
       }else{
	   if(!normalized){
	       System.out.println("titx Energy [GeV]");
	       System.out.println("tity dF/dLogE [cm^-2! sec^-1! sr^-1!])");
	       System.out.println("scal 1.0e5 1.0e12 1.0e-20 1.0e-15");
	   }else{
	       System.out.println("titx Log(Energy [GeV])");
	       System.out.println("tity Probability");
	       System.out.println("scal 5.0 12.0 0.0 0.6");
	   }
	   if(!normalized){
	       System.out.println("gwin 0.2 0.9 0.2 0.9");
	       System.out.println("logx");
	       System.out.println("logy");
	   }

	   double logE = 5.0;
	   if(!normalized){
	       while(logE<12.0-epsilon){
		   double flux = 0.0;
		   if(particleID<4) flux = neutFlux.getDFDLogEwzOsci(logE,particleID);
		   else{
		       for(int id = 1; id<=3; id++) flux += neutFlux.getDFDLogEwzOsci(logE,id);
		   }
		   double neutrinoEnergy = Math.pow(10.0,logE);
		   System.out.format("data %e 0.0 %e 0.0\n",neutrinoEnergy,flux);
		   logE += Particle.getDeltaLogEnergy();
	       }
	   }else{
	       logE = 5.0;
	       double norm = 0.0;
	       while(logE<12.0){
		   if(particleID<4) norm += neutFlux.getDFDLogE(logE+epsilon,particleID)*Particle.getDeltaLogEnergy();
		   else{
		       for(int id = 1; id<=3; id++) norm += neutFlux.getDFDLogEwzOsci(logE+epsilon,id)*Particle.getDeltaLogEnergy();
		   }
		   logE += Particle.getDeltaLogEnergy();
	       }
	       logE = 5.0;
	       while(logE<12.0){
		   double flux = 0.0;
		   if(particleID<4) flux = neutFlux.getDFDLogE(logE+epsilon,particleID)/norm;
		   else{
		       for(int id = 1; id<=3; id++) flux += neutFlux.getDFDLogEwzOsci(logE+epsilon,id)/norm;
		   }
		   System.out.println("data " + logE + " 0.0 " + flux + " 0.0");
		   logE += Particle.getDeltaLogEnergy();
	       }
	   }

	   System.out.println("join");
	   System.out.println("disp");
	   System.out.println("endg");
       }
	    
   }
}
