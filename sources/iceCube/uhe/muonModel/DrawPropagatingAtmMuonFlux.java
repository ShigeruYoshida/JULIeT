package iceCube.uhe.muonModel;

import iceCube.uhe.particles.*;
import iceCube.uhe.propagation.*;
import iceCube.uhe.muonModel.*;

import java.io.*;

public class DrawPropagatingAtmMuonFlux {

   private static final double ln10 = Math.log(10.0);
   private static int dimension = Particle.getDimensionOfLogEnergyMatrix();

   public static void main(String[] args) throws IOException {

       String fileName = null;
       int model = 1;
       double cosTheta = 0.5;
 
       if(args.length!=2){
            System.out.println("Usage: DrawAtmMuonFlux MatrixFileName, ZenithAngle[deg]");
            System.exit(0);
        }else{
            fileName = args[0];
	    double zenith = Double.valueOf(args[1]).doubleValue();
	    cosTheta = Math.cos(Math.toRadians(zenith));
	    System.err.println("Zenith " + zenith + "[deg] cosine=" + cosTheta);
	}

       PropagatingAtmMuonFlux muonFlux = new PropagatingAtmMuonFlux();

        // Read the serialized object of the Neutrino Charged Interaction Matrix
       DataInputStream in = new DataInputStream(new FileInputStream(fileName));
       muonFlux.readMatrix(in);
       in.close( );

        System.out.println("titx Log E[GeV]");
        System.out.println("tity log (Flux E^2 [GeV cm^-2 sec^-1 sr^-1])");
        System.out.println("scal 6.0 12.0 -15.0 -7.0");


	int iLogE;
	double logE;

	// Mu
	for(iLogE=0;iLogE<dimension;iLogE++){
	    logE = Particle.getLogEnergyMinimum( ) + 
                    Particle.getDeltaLogEnergy( )*(double )iLogE;
	    double Flux = muonFlux.getDFMuDLogE(logE,cosTheta);
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
