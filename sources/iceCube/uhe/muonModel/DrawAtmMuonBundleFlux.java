package iceCube.uhe.muonModel;

import iceCube.uhe.particles.*;
import iceCube.uhe.muonModel.*;

import java.io.*;

public class DrawAtmMuonBundleFlux {

   private static final double ln10 = Math.log(10.0);
   private static int dimension = Particle.getDimensionOfLogEnergyMatrix();

   public static void main(String[] args) throws IOException {

       String fileName = null;
       int model = 1;
       double cosTheta = 0.5;
 
       if(args.length!=1){
            System.out.println("Usage: DrawAtmMuonBundleFlux ZenithAngle[deg]");
            System.exit(0);
        }else{
	    double zenith = Double.valueOf(args[0]).doubleValue();
	    cosTheta = Math.cos(Math.toRadians(zenith));
	    System.err.println("Zenith " + zenith + "[deg] cosine=" + cosTheta);
	}

       AtmMuonBundleFlux muonFlux = new AtmMuonBundleFlux();
       //muonFlux.setCutOffFeature(true);

       System.err.println("Threshold Energy of muons =" + muonFlux.getMuonThresholdEnergy() +
			  " [GeV]");
       System.err.println("Assumed mass number of cosmic rays : " +
			  muonFlux.getMassNumber());
       System.err.println("alpha = " +
			  muonFlux.getAlpha());
       System.err.println("CR Energy =" + muonFlux.getEffectiveEnergyOfCRs(1.0e6, cosTheta));

        System.out.println("titx Log E[GeV]");
        System.out.println("tity log (Flux E^2 [GeV cm^-2 sec^-1 sr^-1])");
        System.out.println("scal 5.0 10.0 -12.0 -4.5");


	int iLogE;
	double logE;
	muonFlux.setFluxCalculationMode(false,false);// Use the bare Elbert formula

	// with default parameters
	for(iLogE=0;iLogE<dimension;iLogE++){
	    logE = Particle.getLogEnergyMinimum( ) + 
                    Particle.getDeltaLogEnergy( )*(double )iLogE;
	    double Flux = muonFlux.getDFDLogE(logE,cosTheta);
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
	System.out.println("cont");

	// with Lower Threshold energy of muons in a bundle
	//muonFlux.setMuonThresholdEnergy(5.0e2); // 500 GeV threshold
	//muonFlux.setAlpha(1.9);                 // alpha 1.9
	muonFlux.setFluxCalculationMode(true,false); // include the fluctuation effects 
	                                            // in the average way
	for(iLogE=0;iLogE<dimension;iLogE++){
	    logE = Particle.getLogEnergyMinimum( ) + 
                    Particle.getDeltaLogEnergy( )*(double )iLogE;
	    double Flux = muonFlux.getDFDLogE(logE,cosTheta);
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
	System.out.println("cont");

	muonFlux.setFluxCalculationMode(true,true); // include the fluctuation effects 
	                                            // in the event-by-event basis
	for(iLogE=0;iLogE<dimension;iLogE++){
	    logE = Particle.getLogEnergyMinimum( ) + 
                    Particle.getDeltaLogEnergy( )*(double )iLogE;
	    double Flux = muonFlux.getDFDLogE(logE,cosTheta);
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
	System.out.println("cont");

	System.out.println("endg");
	    
   }
}
