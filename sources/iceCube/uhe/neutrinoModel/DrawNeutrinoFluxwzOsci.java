package iceCube.uhe.neutrinoModel;

import iceCube.uhe.neutrinoModel.*;
import iceCube.uhe.particles.*;
import java.io.*;

public class DrawNeutrinoFluxwzOsci {

   public static void main(String[] args) throws IOException {

       int model = 1;
       int particleID = 0;
       double I_10PeV = 0.0;
 
       if(args.length!=2){
            System.out.println("Usage: DrawNeutrinoFlux model-parameter particleID");
            System.exit(0);
        }else{
	    model = Integer.valueOf(args[0]).intValue();
	    particleID = Integer.valueOf(args[1]).intValue();
	}

       NeutrinoFlux neutFlux = new NeutrinoFlux(model);

        System.out.println("titx Log E[GeV]");
        System.out.println("tity log (Flux E^2 [GeV cm^-2 sec^-1 sr^-1])");
        System.out.println("scal 5.0 12.0 -14.0 -4.0");

	double logE = 5.0;
	while(logE<12.0){
	    double EFlux = 0.0;
	    if(particleID<4) EFlux = neutFlux.getEFluxwzOsci(logE,particleID);
	    else{
		for(int id = 1; id<=3; id++) EFlux += neutFlux.getEFluxwzOsci(logE,id);
	    }
	    double logEFlux;
	    if(EFlux > 0.0){
		logEFlux = Math.log(EFlux)/Math.log(10.0);
	    }else{
		logEFlux = -14.0;
	    }
	    System.out.println("data " + logE + " 0.0 " + logEFlux + " 0.0");
	    logE += 0.1;
	}

	logE = 7.0;
	while(logE<12.0){
	    if(particleID<4){
		I_10PeV += neutFlux.getDFDLogEwzOsci(logE,particleID)*Particle.getDeltaLogEnergy( );
	    }else{
		double dFlux = 0.0;
		for(int id = 1; id<=3; id++) dFlux += neutFlux.getEFluxwzOsci(logE,id);
		I_10PeV += dFlux*Particle.getDeltaLogEnergy( );
	    }
	    logE += Particle.getDeltaLogEnergy( );

	}

	System.out.println("mssg I(10PeV)=" + I_10PeV);


	System.out.println("join");
	System.out.println("disp");
	System.out.println("endg");
	    
   }
}
