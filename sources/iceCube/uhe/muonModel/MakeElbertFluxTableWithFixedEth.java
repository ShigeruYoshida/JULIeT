package iceCube.uhe.muonModel;

import iceCube.uhe.particles.*;
import iceCube.uhe.points.*;
import iceCube.uhe.geometry.*;
import iceCube.uhe.muonModel.*;

import java.io.*;

public class MakeElbertFluxTableWithFixedEth {

   private static final double ln10 = Math.log(10.0);
   private static int dimension = Particle.getDimensionOfLogEnergyMatrix();

   public static void main(String[] args) throws IOException {

       String fileName = null;
       // This must be changed following your environment.
       double logDistance = 0.0;
       double epsilon = 1.0e-4;
       double muonThresholdEnergy = 1.0e3;// 1TeV
 
       if(args.length!=3){
            System.out.println("Usage: MakeElbertFluxTable MatrixFileName logDistance Eth");
            System.exit(0);
       }else{
	   fileName = args[0];
	   logDistance = Double.valueOf(args[1]).doubleValue();
	   muonThresholdEnergy = Double.valueOf(args[2]).doubleValue();
       }

       double distance = Math.pow(10.0,logDistance);

       PropagatingAtmMuonFlux muonFlux = new PropagatingAtmMuonFlux();
       //muonFlux.setCutOffFeature(false); // No GZK cutoff


        // Read the serialized object of the propagation Matrix
       DataInputStream in = new DataInputStream(new FileInputStream(fileName));
       muonFlux.readMatrix(in);
       in.close( );

       int iLogE;
       double logE;

       IceCubeCoordinate iceCube = new IceCubeCoordinate();
       double rEarth_ice = ParticlePoint.REarth + iceCube.getGlacierDepth();
       //double detectorDepth = iceCube.getGlacierDepth()-iceCube.elevation-8.8e4;
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

       AtmMuonBundleFlux bundleFlux = new AtmMuonBundleFlux();
       System.err.println("Distance=" + distance + " cosZenith =" + cosZenith +
			  "MuonE threshold=" + muonThresholdEnergy);

       // Scan 
       //
       double deltaAlpha = 0.01; 
       double alphaMin = 1.80;
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
	   bundleFlux.setAlpha(alpha);

	   System.out.println(muonThresholdEnergy + " " + "1 10.0");

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

	   alpha += deltaAlpha;
       }
	    
   }

}
