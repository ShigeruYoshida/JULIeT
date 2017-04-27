package iceCube.uhe.muonModel;

import iceCube.uhe.particles.*;
import iceCube.uhe.points.*;
import iceCube.uhe.interactions.*;
import iceCube.uhe.propagation.*;
import iceCube.uhe.muonModel.*;

import java.io.*;

public class DrawCREnergyForBundle {

   private static final double ln10 = Math.log(10.0);
   private static int dimension = Particle.getDimensionOfLogEnergyMatrix();

   public static void main(String[] args) throws IOException {

       String fileName = null;
       double alpha = 1.9; // Elbert parameter
       double muonEThreshold = 300.0; // mu Ethreshold for Bundle
       double cosTheta = 0.5; 
 
       if(args.length!=3){
            System.out.println("Usage: DrawCREnergyForBundle alpha muEthres ZenithAngle[deg]");
            System.exit(0);
        }else{
	    alpha = Double.valueOf(args[0]).doubleValue();
	    muonEThreshold = Double.valueOf(args[1]).doubleValue();
	    double zenith = Double.valueOf(args[2]).doubleValue();
	    cosTheta = Math.cos(Math.toRadians(zenith));
	    System.err.println("Zenith " + zenith + "[deg] cosine=" + cosTheta);
	    System.err.println("Muon Eth " + muonEThreshold + "[GeV] alpha=" + alpha);
	}


	// Now using the analytical approximation in calculation of 
	// the muon propagation. The propagation matricis are not used.

	AtmMuonBundleFlux muonFluxModel =  new AtmMuonBundleFlux();

	muonFluxModel.setAlpha(alpha);
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

	// Now drawing
	double beta7 = CELbeta.getBeta(7.0);
	double criticalEnergy = AtmMuonBundleFlux.criticalEnergy;

        System.out.println("titx Surface Muon Bundle Energy [GeV]");
        System.out.println("tity Cosmic Ray Energy [GeV]");
        System.out.println("scal 1.0e7 1.0e10 1.0e8 1.0e12");

	for(int iLogE=0;iLogE<dimension;iLogE++){
	    double logMuonBundleEnergy = Particle.getLogEnergyMinimum( ) + 
                    Particle.getDeltaLogEnergy( )*(double )iLogE;
	    double muonBundleEnergy = Math.pow(10.0,logMuonBundleEnergy); 

	    double beta = CELbeta.getBeta(logMuonBundleEnergy);
	    if(logMuonBundleEnergy<7.0) beta = beta7;

	    double energyEnhancedFactor = Math.exp(beta*slantDepth); 

	    double muonEThresholdAtSurface = 
		energyEnhancedFactor*(muonEThreshold+criticalEnergy)-criticalEnergy;
	    
	    muonFluxModel.setMuonThresholdEnergy(muonEThresholdAtSurface);

	    double crEnergy = 
		muonFluxModel.getEffectiveEnergyOfCRs(muonBundleEnergy,cosTheta);

	    System.out.println("data " + muonBundleEnergy + " 0.0 " + crEnergy + " 0.0");
	}

        System.out.println("logx"); 
        System.out.println("logy"); 
        System.out.println("join"); 
        System.out.println("disp"); 
        System.out.println("endg"); 

   }
}
