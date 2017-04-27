package iceCube.uhe.muonModel;

import iceCube.uhe.particles.*;
import iceCube.uhe.points.*;
import iceCube.uhe.geometry.*;
import iceCube.uhe.muonModel.*;

import java.io.*;

/**

This class implements the main method to run 
PropagatingAtmMuonFlux.getDFMuDLogCREDLogE(logCRE, logMuonE, cosZenith),
dF^2/dLogEcrDlogEmu, to make the numerical table for IceCube
downstream analysis.

 */

public class MakeCosmicRayFluxTable {

   private static final double ln10 = Math.log(10.0);
   private static int dimension = Particle.getDimensionOfLogEnergyMatrix();

   public static void main(String[] args) throws IOException {

       String fileName = null;
       // This must be changed following your environment.
       double logDistance = 0.0;
       double epsilon = 1.0e-4;
       double muonThresholdEnergy = 1.0e3;// 1TeV
       double alpha = 1.9;
       double halfWidthOfLogE = Particle.getDeltaLogEnergy()*0.5; 
 
       if(args.length<4){
            System.out.println("Usage: MakeCosmicRayFluxTable MatrixFileName logDistance alpha, Eth (halfLogEwidth)");
            System.exit(0);
       }else{
	   fileName = args[0];
	   logDistance = Double.valueOf(args[1]).doubleValue();
	   alpha = Double.valueOf(args[2]).doubleValue();
	   muonThresholdEnergy = Double.valueOf(args[3]).doubleValue();
       }

       // width of logCosmicRayEnergy bin - it determines the integral width
       if(args.length == 5) halfWidthOfLogE = Double.valueOf(args[4]).doubleValue();

       double distance = Math.pow(10.0,logDistance);

       boolean cutoffExists = true; // GZK cutoff
       PropagatingAtmMuonFlux muonFlux = 
	   new PropagatingAtmMuonFlux(alpha,muonThresholdEnergy,cutoffExists);


        // Read the serialized object of the propagation Matrix
       DataInputStream in = new DataInputStream(new FileInputStream(fileName));
       muonFlux.readMatrix(in);
       in.close( );

       //
       // Zenith angle representing the average trajectory of given distance
       //
       IceCubeCoordinate iceCube = new IceCubeCoordinate();
       double rEarth_ice = ParticlePoint.REarth + iceCube.getGlacierDepth();
       double detectorDepth = iceCube.getGlacierDepth()-iceCube.elevation;
       //880m from the icecube center

       double cosZenith = 	
	   (rEarth_ice*rEarth_ice - distance*distance - (rEarth_ice-detectorDepth)*(rEarth_ice-detectorDepth))/
	   (2.0*distance*(rEarth_ice-detectorDepth));
       if(cosZenith>1.0) cosZenith = 1.0;

       double rMC = 8.8e4;
       detectorDepth = iceCube.getGlacierDepth()-iceCube.elevation-rMC*Math.PI/4.0*cosZenith;
       cosZenith = 	
	   (rEarth_ice*rEarth_ice - distance*distance - (rEarth_ice-detectorDepth)*(rEarth_ice-detectorDepth))/
	   (2.0*distance*(rEarth_ice-detectorDepth));
       if(cosZenith>1.0) cosZenith = 1.0;

       System.err.println("Distance=" + distance + " cosZenith =" + cosZenith +
			  "MuonE threshold=" + muonThresholdEnergy);
       System.err.println("a half bin width of logEcr = " + halfWidthOfLogE);

       // Scan 
       //
       for(int iLogE=0;iLogE<dimension;iLogE+=5){ // 0.05 decade step - muE loop
	   double logEmu = Particle.getLogEnergyMinimum( ) + 
	       Particle.getDeltaLogEnergy( )*(double )iLogE;

	   double logCosmicRayEnergy = Particle.getLogEnergyMinimum( );
	   double logEmax = Particle.getLogEnergyMinimum( ) + 
	       Particle.getDeltaLogEnergy()*(double)dimension;

	   //if(iLogE%50==0) // debugging
	   //    muonFlux.listFluxes(logEmu,cosZenith);

	   while(logCosmicRayEnergy<logEmax){ // Cosmic ray energy loop
	       double flux = muonFlux.getDFMuDLogCREDLogE(logCosmicRayEnergy, halfWidthOfLogE,
							  logEmu+epsilon,cosZenith)*2.0*halfWidthOfLogE;
	       System.out.println(logDistance + " " + cosZenith + " " + 
				  logEmu + " " + logCosmicRayEnergy + " " + flux);

	       logCosmicRayEnergy += 2.0*halfWidthOfLogE;
	   }

       }
	    
   }

}
