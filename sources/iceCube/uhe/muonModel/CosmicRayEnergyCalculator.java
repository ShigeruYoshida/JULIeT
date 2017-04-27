package iceCube.uhe.muonModel;

import iceCube.uhe.particles.*;
import iceCube.uhe.points.*;
import iceCube.uhe.interactions.*;
import iceCube.uhe.propagation.*;
import iceCube.uhe.muonModel.*;

import java.io.*;

public class CosmicRayEnergyCalculator {

   private static final double ln10 = Math.log(10.0);
   private static int dimension = Particle.getDimensionOfLogEnergyMatrix();
   protected static  double detectorDepth = 1.4e5;  // Detector Depth.. 1400m 

   public static void main(String[] args) throws IOException {


	// Now using the analytical approximation in calculation of 
	// the muon propagation. The propagation matricis are not used.

	AtmMuonBundleFlux muonFluxModel =  new AtmMuonBundleFlux();

	muonFluxModel.setAlpha(2.04);
	muonFluxModel.setMuonThresholdEnergy(3.73e3);
	muonFluxModel.setCutOffFeature(true); // with GZK cutoff

	ParticlePoint s = new ParticlePoint(0.0,0.0,0); // ice

	// Reading data
	DataInputStream input = new DataInputStream(System.in);

	BufferedReader  d     = new BufferedReader(new InputStreamReader(input));
	String buffer; int sep = 0; int sepstart = 0;
        char separator = ' ';

	while((buffer = d.readLine())!=null){
	    //try{

		// 1st line -- logNpe, cosZenith

		sepstart = 0;
		sep = buffer.indexOf(separator,sepstart+1);
		sepstart = sep;

		sep = buffer.indexOf(separator,sepstart+1);
		double logNpe =
		    Double.valueOf(buffer.substring(sepstart+1,sep)).doubleValue( );
		sepstart = sep;

		sep = buffer.indexOf(separator,sepstart+1);
		double cosTheta =
		    Double.valueOf(buffer.substring(sepstart+1,sep)).doubleValue( );
		sepstart = sep;

		double logEinIce = getLogEnergyFromLogNpe(logNpe);
		double slantDepth = getSlantDepth(cosTheta,s);

		double beta_loss = 4.4619776009127244E-6;
		if(logEinIce >= 7.0) beta_loss = CELbeta.getBeta(logEinIce);

		double energyEnhancementFactor = Math.exp(beta_loss*slantDepth);
		double muonBundleEnergy = Math.pow(10.0,logEinIce);
		double logEsurface = Math.log(muonBundleEnergy*energyEnhancementFactor)/ln10;
		double cosmicRayEnergy = 
  	        muonFluxModel.getEffectiveEnergyOfCRs(logEinIce,cosTheta,beta_loss,slantDepth);
		double logCosmicRayEnergy = Math.log(cosmicRayEnergy)/ln10;

		System.out.println(logNpe + " " + cosTheta + " " + logEinIce + " " +
				   logEsurface + " " + logCosmicRayEnergy);


		//}catch (EOFException e){
		//buffer = null;
		//break;
		//}

	}

   }

    protected static double getSlantDepth(double cosTheta, ParticlePoint s){
	double sq_term = 
	Math.sqrt((ParticlePoint.REarth-detectorDepth)*(ParticlePoint.REarth-detectorDepth)
		  *cosTheta*cosTheta + 
		  2.0*ParticlePoint.REarth*detectorDepth-detectorDepth*detectorDepth);
	double trajectoryLength = sq_term - (ParticlePoint.REarth-detectorDepth)*cosTheta;

	double slantDepth = trajectoryLength*s.getMediumDensity();

	return slantDepth;
    }


    protected static double getLogEnergyFromLogNpe(double logNpe){
	double logEnergy = 6.359 - 1.318*logNpe + 0.2706*logNpe*logNpe;
	return logEnergy;
    }

}
