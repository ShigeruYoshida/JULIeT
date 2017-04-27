package iceCube.uhe.muonModel;

import iceCube.uhe.particles.*;
import iceCube.uhe.points.*;
import iceCube.uhe.geometry.*;
import iceCube.uhe.muonModel.*;

import java.io.*;

public class MakeElbertFluxTable {

   private static final double ln10 = Math.log(10.0);
   private static int dimension = Particle.getDimensionOfLogEnergyMatrix();

   public static void main(String[] args) throws IOException {

       String fileName = null;
       // This must be changed following your environment.
       String referenceFileName = 
	   "/home/syoshida/java_lib/data/neutrino_earth/ice/1e5_10cm.data";
       double distanceRef = Math.pow(10.0,5.10);
       double logDistance = 0.0;
       double epsilon = 1.0e-4;
 
       if(args.length!=2){
            System.out.println("Usage: MakeElbertFluxTable MatrixFileName logDistance");
            System.exit(0);
       }else{
	   fileName = args[0];
	   logDistance = Double.valueOf(args[1]).doubleValue();
	   System.err.println("logDistance " + logDistance);
       }

       double distance = Math.pow(10.0,logDistance);

       PropagatingAtmMuonFlux muonFlux = new PropagatingAtmMuonFlux();
       PropagatingAtmMuonFlux muonFluxRef = new PropagatingAtmMuonFlux();
       muonFlux.setCutOffFeature(false); // No GZK cutoff 
       //muonFlux.setFluxCalculationMode(false,false);

        // Read the serialized object of the propagation Matrix
       DataInputStream in = new DataInputStream(new FileInputStream(fileName));
       muonFlux.readMatrix(in);
       in.close( );

       in = new DataInputStream(new FileInputStream(referenceFileName));
       muonFluxRef.readMatrix(in);
       in.close( );

       int iLogE;
       double logE;

       IceCubeCoordinate iceCube = new IceCubeCoordinate();
       double rEarth_ice = ParticlePoint.REarth + iceCube.getGlacierDepth();
       double detectorDepth = iceCube.getGlacierDepth()-iceCube.elevation;
       //double detectorDepth = iceCube.getGlacierDepth()-iceCube.elevation-8.8e4;
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

       detectorDepth = iceCube.getGlacierDepth()-iceCube.elevation;
       double cosZenithRef = 	
	   (rEarth_ice*rEarth_ice - distanceRef*distanceRef - (rEarth_ice-detectorDepth)*(rEarth_ice-detectorDepth))/
	   (2.0*distanceRef*(rEarth_ice-detectorDepth));
       if(cosZenithRef>1.0) cosZenithRef = 1.0;
       detectorDepth = iceCube.getGlacierDepth()-iceCube.elevation-rMC*Math.PI/4.0*cosZenithRef;

       cosZenithRef = 	
	   (rEarth_ice*rEarth_ice - distanceRef*distanceRef - (rEarth_ice-detectorDepth)*(rEarth_ice-detectorDepth))/
	   (2.0*distanceRef*(rEarth_ice-detectorDepth));

       AtmMuonBundleFlux bundleFlux = new AtmMuonBundleFlux();
       double iniceMuonEnergyReference = 1.0e6; // 10 PeV
       double iniceMuonEnergyThresholdReference = bundleFlux.getMuonThresholdEnergy();
       double surfaceMuonEnergyThresholdReference = 
	   (iniceMuonEnergyThresholdReference+bundleFlux.criticalEnergy)/
	   muonFluxRef.getAverageMuonEnergyLossAfterPropagation()-bundleFlux.criticalEnergy;
       bundleFlux.setMuonThresholdEnergy(surfaceMuonEnergyThresholdReference);
       double surfaceMuonEnergyReference = 
	   (iniceMuonEnergyReference+bundleFlux.criticalEnergy)/
	   muonFluxRef.getAverageMuonEnergyLossAfterPropagation()-bundleFlux.criticalEnergy;
       double energyRatioReference = surfaceMuonEnergyReference/
	   bundleFlux.getEffectiveEnergyOfCRs(surfaceMuonEnergyReference,cosZenithRef);
       System.err.println("surfaceMuonEnergyReference = " + surfaceMuonEnergyReference);
       System.err.println("surfaceMuonEnergyThreshldReference = " + 
			  surfaceMuonEnergyThresholdReference);
       System.err.println("Energy Ratio Reference = " + energyRatioReference);


       System.err.println("DistanceRef=" + distanceRef + " cosZenithRef =" + cosZenithRef);
       System.err.println("Distance=" + distance + " cosZenith =" + cosZenith);

       // Scan 
       //
       double deltaAlpha = 0.01; double deltaMuonE = 10.0;
       double alphaMin = 1.90;
       double alphaMax = 2.35;
       double alpha = alphaMin;
       int i = 0 ;
       while(alpha<=alphaMax){
	   alpha += deltaAlpha;
	   i++;
       }
       System.out.println(alphaMin + " " + i + " " + deltaAlpha);
       alpha = alphaMin;
       while(alpha<=alphaMax){
	   muonFlux.atmMuFlux.setAlpha(alpha);
	   muonFluxRef.atmMuFlux.setAlpha(alpha);
	   bundleFlux.setAlpha(alpha);
	   double energyRatio = surfaceMuonEnergyReference/
	       bundleFlux.getEffectiveEnergyOfCRs(surfaceMuonEnergyReference,cosZenithRef);
	   double surfaceMuonEnergy = 
	       energyRatio/energyRatioReference*surfaceMuonEnergyThresholdReference;
	   double iniceMuonEnergy = 
	       (surfaceMuonEnergy + bundleFlux.criticalEnergy)*
	       muonFluxRef.getAverageMuonEnergyLossAfterPropagation()-bundleFlux.criticalEnergy;
	   double muonEnergyMin = 0.7*iniceMuonEnergy;
	   double muonEnergyMax = 1.6*iniceMuonEnergy;
	   double muonThresholdEnergy = muonEnergyMin;
	   int times = 45;
	   deltaMuonE = (muonEnergyMax-muonEnergyMin)/(double)times;
	   times = 0;
	   while(muonThresholdEnergy <= muonEnergyMax){
	       times++; muonThresholdEnergy+= deltaMuonE;
	   }
	   System.out.println(muonEnergyMin + " " + times + " " + deltaMuonE);
	   muonThresholdEnergy = muonEnergyMin;
	   while(muonThresholdEnergy <= muonEnergyMax){

	       double thresholdE = (muonThresholdEnergy+muonFlux.atmMuFlux.criticalEnergy)/
		   muonFlux.getAverageMuonEnergyLossAfterPropagation()-muonFlux.atmMuFlux.criticalEnergy;
	       muonFlux.atmMuFlux.setMuonThresholdEnergy(thresholdE);
	       for(iLogE=0;iLogE<dimension;iLogE+=5){ // 0.05 decade step
		   logE = Particle.getLogEnergyMinimum( ) + 
		       Particle.getDeltaLogEnergy( )*(double )iLogE;
		   double flux = muonFlux.getDFMuDLogE(logE+epsilon,cosZenith);
		   System.out.println(logDistance + " " + alpha + " " + 
				      muonThresholdEnergy + " " + " " +
				      cosZenith + " " + logE + " " + flux);
	       }

	       muonThresholdEnergy += deltaMuonE;
	    }

	   alpha += deltaAlpha;
       }
	    
   }

}
