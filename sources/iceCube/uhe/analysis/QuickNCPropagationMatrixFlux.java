package iceCube.uhe.analysis;

import iceCube.uhe.particles.*;
import iceCube.uhe.interactions.*;
import iceCube.uhe.propagation.*;
import iceCube.uhe.analysis.*;

import java.io.*;
import java.util.*;

/**
   It calculates the detectable neutrino event intensity 
   at the Earth Surface as I3ParticleFlux does
   but the calculation is made by using directly 
   the zenith angle binned propagation matrix
   and the numerically calculated effective area (I3EffectiveArea.java)
   without relying on I3Particle MC events.

   Written by S. Yoshida May 11 2008

*/

public class QuickNCPropagationMatrixFlux extends PropagationMatrixFlux {

    protected NeutrinoQuickPropagator propagator = null;

    protected InteractionsMatrix nuCCMtx = null; 
    private String nuCCMtxObjectFile = "ENeutrinoChargeMtx";
    protected InteractionsMatrix nuNCMtx = null; 
    private String nuNCMtxObjectFile = "ENeutrinoNeutralMtx";
    protected double nuNCEnhancementFactor = 1.0;
    protected static String[] intMtxPathname =
    {"iceCube/uhe/interactions/ice/","iceCube/uhe/interactions/rock/"};

    /** constroctor */
    public QuickNCPropagationMatrixFlux(double nuNCEnhancementFactor) throws IOException{
	super();
	setNeutrinoNCEnhancement(nuNCEnhancementFactor);
	// The Charged Current neutrino-nucleon
	String objectFile = intMtxPathname[0].concat(nuCCMtxObjectFile);
	InputStream in = ClassLoader.getSystemResourceAsStream(objectFile);
	nuCCMtx = InteractionsMatrixInput.inputInteractionsMatrix(in);
	// The Neutral Current neutrino-nucleon
	objectFile = intMtxPathname[0].concat(nuNCMtxObjectFile);
	in = ClassLoader.getSystemResourceAsStream(objectFile);
	nuNCMtx = InteractionsMatrixInput.inputInteractionsMatrix(in);
    }

    protected void setNeutrinoNCEnhancement(double nuNCEnhancementFactor){
	this.nuNCEnhancementFactor = nuNCEnhancementFactor;
    }

    /** Calculate the neutrino yield [cm^2 sec sr] in form of the table
        by reasing out the pre-stored propagation matrix data
	via the PropagationMatrixFactory.
	A primary neutrino bulk is assumed to be consist of
	nu_e, nu_mu, nu_tau = 1:1:1.
    */
    public void calculateYield() throws IOException {

	if(inIceFlavorList.size()<1){ // nothing stored to specify the inIceParticles!
	    System.err.println(" You have to call setInIceParticle(flavor doublet) first!");
	    System.exit(0);
	}

        for(int i=0;i<Particle.getDimensionOfLogEnergyMatrix(); i++)  
            yieldTable[i]=0.0; 
 
	// Zenith Angle Loop
       for(int upDown=1;upDown<3;upDown++){
	   for(int itheta=0;itheta<I3ParticleWeightFiller.matrixFileName[upDown].length;
	       itheta++){

	       // Read the serialized object of the Neutrino Charged Interaction Matrix
	       String fileName = 
		   I3ParticleWeightFiller.pathname[upDown%2] + 
		   I3ParticleWeightFiller.matrixFileName[upDown][itheta];
	       DataInputStream in = new DataInputStream(new FileInputStream(fileName));
	       matrix.readMatrix(in);
	       in.close( );
	       System.err.println("Reading the matrix from " + fileName + " done.");

	       // Solid angle calculation
	       double radiansUp = 0.0; double radiansDown = 0.0;
	       double cosZenith = 0.0;
	       if(upDown==1){ // rock
	          radiansUp = Math.toRadians(I3ParticleWeightFiller.rockrange[itheta]);
		  radiansDown = 
		      Math.toRadians(I3ParticleWeightFiller.rockrange[itheta+1]);
		  cosZenith = Math.cos(Math.PI-radiansUp); // nadir to Zenith
	       }else { // ice
	          radiansUp = Math.toRadians(I3ParticleWeightFiller.icerange[itheta]);
		  radiansDown = 
		      Math.toRadians(I3ParticleWeightFiller.icerange[itheta+1]);
		  cosZenith = Math.cos(radiansUp);
	       }
	       double solidAngle = 
		   2.0*Math.PI*Math.abs(Math.cos(radiansDown)-Math.cos(radiansUp));
	       System.err.println(" Solid angle " + solidAngle);


	       // Earth-surface Neutrino Energy loop
	       double epsilon = 1.0e-4; // round-off error
	       for(int iLogE = 0;iLogE<Particle.getDimensionOfLogEnergyMatrix();
		   iLogE++){

		   double logNuSurfaceEnergy = Particle.getLogEnergyMinimum()
		       + Particle.getDeltaLogEnergy()*((double )iLogE) + epsilon;

		   nuE.putLogEnergy(logNuSurfaceEnergy);
		   nuMu.putLogEnergy(logNuSurfaceEnergy);
		   nuTau.putLogEnergy(logNuSurfaceEnergy);

		   // InIceParticle loop
		   ListIterator inIceFlavorIterator = inIceFlavorList.listIterator();
		   ListIterator inIceDoubletIterator = inIceDoubletList.listIterator();
		   ListIterator inIceAreaIterator = inIceAreaList.listIterator();
		   double flux = 0.0;
		   while(inIceFlavorIterator.hasNext()){
		       int inIceFlavor = 
			   ((Integer )inIceFlavorIterator.next()).intValue();
		       int inIceDoublet = 
			   ((Integer )inIceDoubletIterator.next()).intValue();

		       Particle iceParticle = new Particle(inIceFlavor,inIceDoublet);

		       EffAreaTable areaTable = (EffAreaTable )inIceAreaIterator.next();

		       // InIce Lepton Energy loop
		       for(int jLogE= 0;jLogE<Particle.getDimensionOfLogEnergyMatrix();
			   jLogE ++) {

			   double logInIceEnergy = Particle.getLogEnergyMinimum()
			       + Particle.getDeltaLogEnergy()*((double )jLogE)+epsilon;
			   iceParticle.putLogEnergy(logInIceEnergy);

			   // effective area [cm^2]
			   double area = 
			       areaTable.getArea(logInIceEnergy,cosZenith)*1.0e10;
			   double enhanceFactor = 1.0;
			   if(inIceDoublet==0){
			       double sigmaCC = nuCCMtx.getSigmaMatrix(jLogE);
			       double sigmaNC = nuNCMtx.getSigmaMatrix(jLogE);
			       enhanceFactor = (sigmaNC*nuNCEnhancementFactor+sigmaCC)/(sigmaNC+sigmaCC);
			   }

			   flux += enhanceFactor*area*matrix.getDF(nuE,iceParticle);
			   flux += enhanceFactor*area*matrix.getDF(nuMu,iceParticle);
			   flux += enhanceFactor*area*matrix.getDF(nuTau,iceParticle);

		       }// loop over inIce Energy ends

		   }// loop over inIce Particle ends

		   yieldTable[iLogE] += flux*solidAngle*observationTime;

	       } // loop over Neutrino Energy ends

	   } // loop over Zenith Angle ends
       } // up and down

       yieldTableExists = true;

    }

}


