package iceCube.uhe.muonModel;

import iceCube.uhe.particles.*;
import iceCube.uhe.muonModel.*;

import java.io.*;

public class DrawCosmicRayFlux {

   private static final double ln10 = Math.log(10.0);
   private static int dimension = Particle.getDimensionOfLogEnergyMatrix();

   public static void main(String[] args) throws IOException {

       String fileName = null;
       int model = 1;
       double cosTheta = 0.5;
 
       CosmicRayFlux crFlux = new CosmicRayFlux();

        System.out.println("titx Log E[GeV]");
        System.out.println("tity log (Flux E^2 [GeV cm^-2 sec^-1 sr^-1])");
        System.out.println("scal 5.0 11.0 -12.0 -3.0");


	int iLogE;
	double logE;

	// Mu
	for(iLogE=0;iLogE<dimension;iLogE++){
	    logE = Particle.getLogEnergyMinimum( ) + 
                    Particle.getDeltaLogEnergy( )*(double )iLogE;
	    double Flux = crFlux.getDFDLogE(logE);
	    double logEFlux;
	    if(Flux > 0.0){
		logEFlux = Math.log(Flux/ln10)/ln10 + logE;
	    }else{
		logEFlux = -15.0;
	    }
	    System.out.println("data " + logE + " 0.0 " + logEFlux + " 0.0");
	}

	System.out.println("join");
	System.out.println("disp");
	System.out.println("endg");
	    
   }
}
