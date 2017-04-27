package iceCube.uhe.muonModel;

import iceCube.uhe.particles.*;
import iceCube.uhe.points.*;
import iceCube.uhe.interactions.*;
import iceCube.uhe.propagation.*;
import iceCube.uhe.muonModel.*;

import java.io.*;

public class DrawPropagatingAtmMuonBundleFlux {

   private static final double ln10 = Math.log(10.0);
   private static int dimension = Particle.getDimensionOfLogEnergyMatrix();

   public static void main(String[] args) throws IOException {

       String fileName = null;
       int model = 1;
       double cosTheta = 0.5;
 
       if(args.length!=2){
            System.out.println("Usage: DrawAtmMuonBundleFlux MatrixFileName, ZenithAngle[deg]");
            System.exit(0);
        }else{
            fileName = args[0];
	    double zenith = Double.valueOf(args[1]).doubleValue();
	    cosTheta = Math.cos(Math.toRadians(zenith));
	    System.err.println("Zenith " + zenith + "[deg] cosine=" + cosTheta);
	}

       PropagatingAtmMuonFlux muonFlux = new PropagatingAtmMuonFlux();

        // Read the serialized object of the Neutrino Charged Interaction Matrix
       DataInputStream in = 
	   new DataInputStream(new FileInputStream(fileName));
       muonFlux.readMatrix(in);
       in.close( );

        System.out.println("titx Log E[GeV]");
        System.out.println("tity log (Flux E^2 [GeV cm^-2 sec^-1 sr^-1])");
        System.out.println("scal 5.0 10.0 -12.0 -6.0");


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
	System.out.println("cont");



	// Now using the analytical approximation in calculation of 
	// the muon propagation. The propagation matricis are not used.

	AtmMuonBundleFlux muonFluxModel =  new AtmMuonBundleFlux();

	muonFluxModel.setCutOffFeature(true); // with GZK cutoff

	double detectorDepth = 1.4e5;  // Detector Depth.. 1400m = 1.4e5 cm below sea level
	double sq_term = 
	Math.sqrt((ParticlePoint.REarth-detectorDepth)*(ParticlePoint.REarth-detectorDepth)
		  *cosTheta*cosTheta + 
		  2.0*ParticlePoint.REarth*detectorDepth-detectorDepth*detectorDepth);
	double trajectoryLength = sq_term - (ParticlePoint.REarth-detectorDepth)*cosTheta;
	double cos_nadir = sq_term/ParticlePoint.REarth;
	double nadirAngle = Math.acos(cos_nadir)*180.0/Math.PI;

	ParticlePoint s = new ParticlePoint(0.0,nadirAngle*Math.PI/180.0,0); // ice
	double slantDepth = trajectoryLength*s.getMediumDensity();
	System.err.println("Nadir angle at the earth surface " + nadirAngle + " [deg]");
	System.err.println("Propagation Length " + trajectoryLength + " [cm]");
	System.err.println("Medium Density " + s.getMediumDensity() + " [g/cm^3]");

	// Now drawing
	double beta7 = CELbeta.getBeta(7.0);
	for(iLogE=0;iLogE<dimension;iLogE++){
	    logE = Particle.getLogEnergyMinimum( ) + 
                    Particle.getDeltaLogEnergy( )*(double )iLogE;
	    double beta = CELbeta.getBeta(logE);
	    if(logE<7.0) beta = beta7;
	    double Flux = muonFluxModel.getDFDLogE(logE,cosTheta,beta,slantDepth);
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
