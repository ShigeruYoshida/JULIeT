package iceCube.uhe.analysis;

import iceCube.uhe.particles.*;
import iceCube.uhe.points.*;
import iceCube.uhe.propagation.*;
import iceCube.uhe.analysis.*;

import java.io.*;
import java.util.*;

/**
   It calculates the detectable neutrino event intensity 
   at the Earth Surface as I3ParticleFlux does
   but the calculation is made by using 
   the calculated energy distribution by NeutrinoQuickPropagator class
   and the numerically calculated effective area (I3EffectiveArea.java)
   without relying on I3Particle MC events.

   This class uses the NeutrinoQuickPropagator (in the propagation package)
   for calculating the particle propagation in the earth, while
   the superclass, PropagationMatrixFlux assigns this job to
   ProgationMatrixFactory.

   Written by S. Yoshida Feb 15 2008

*/

public class QuickPropagationMatrixFlux extends PropagationMatrixFlux {

    /** Zenith angle [deg] */
    protected static double[][] zenithAngle =
    {
	{
	0.0, 18.19, 25.84, 31.79, 36.87, 41.41, 45.57, 49.46, 53.13, 56.63,
	60.0, 63.26, 66.42, 69.51, 72.54, 75.52, 77.0, 78.46, 79.92, 81.37,
	82.82, 84.26, 85.70, 87.13, 88.56
	},

	{
       89.8, 89.5, 89.0, 88.0, 87.0, 86.0, 85.0, 84.0, 83.0, 82.0, 81.0, 80.0,
       77.5, 75.0, 72.5, 70.0, 65.0, 60.0, 50.0
	}
    };

    protected NeutrinoQuickPropagator propagator = null;

    protected double nuCCEnhancementFactor = 1.0;

    /** constroctor */
    public QuickPropagationMatrixFlux(double nuCCEnhancementFactor){
	super(false);
	setNeutrinoCCEnhancement(nuCCEnhancementFactor);
    }

    protected void setNeutrinoCCEnhancement(double nuCCEnhancementFactor){
	this.nuCCEnhancementFactor = nuCCEnhancementFactor;
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
       for(int upDown=0;upDown<2;upDown++){
	   for(int itheta=0;itheta<zenithAngle[upDown].length;itheta++){
	       double zenith = zenithAngle[upDown][itheta];

	       // propagate Neutrino
               ParticlePoint s= new ParticlePoint(0.0,Math.toRadians(zenith),upDown);
	       propagator.setParticlePoint(s);
	       propagator.propagateNeutrinoToIceCubeDepth(zenith,
							  nuCCEnhancementFactor);
	       
	       // Solid angle calculation
	       double radiansUp = 0.0; double radiansDown = 0.0;
	       double cosZenith = 0.0;
	       if(upDown==1){ // rock
		   radiansUp = Math.toRadians(zenithAngle[upDown][itheta]);
		   if(itheta==(zenithAngle[upDown].length-1)) radiansUp = 0.0;

		   if(itheta==0) {
		       radiansDown = Math.toRadians(90.0);
		   }else{
		       radiansDown = Math.toRadians(zenithAngle[upDown][itheta-1]);
		   }
		   cosZenith = Math.cos(Math.PI-Math.toRadians(zenith)); // nadir to Zenith
	       }else { // ice
		   radiansUp = Math.toRadians(zenithAngle[upDown][itheta]);
		   if(itheta<(zenithAngle[upDown].length-1)){
		       radiansDown = 
			   Math.toRadians(zenithAngle[upDown][itheta+1]);
		   }else{
		       radiansDown = Math.toRadians(90.0);
		   }
		   cosZenith = Math.cos(Math.toRadians(zenith));
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

		       EffAreaTable areaTable = (EffAreaTable )inIceAreaIterator.next();

		       // InIce Lepton Energy loop
		       for(int jLogE= 0;jLogE<Particle.getDimensionOfLogEnergyMatrix();
			   jLogE ++) {

			   double logInIceEnergy = Particle.getLogEnergyMinimum()
			       + Particle.getDeltaLogEnergy()*((double )jLogE)+epsilon;

			   // effective area [cm^2]
			   double area = 
			       areaTable.getArea(logInIceEnergy,cosZenith)*1.0e10;
			   double enhanceFactor = 1.0;
			   if(inIceDoublet==0) enhanceFactor = nuCCEnhancementFactor;

			   
			   flux += enhanceFactor*area*
			       propagator.getDF(0,logNuSurfaceEnergy, // nuE
						inIceFlavor,inIceDoublet,logInIceEnergy);
			   flux += enhanceFactor*area*
			       propagator.getDF(1,logNuSurfaceEnergy, // nuMu
						inIceFlavor,inIceDoublet,logInIceEnergy);
			   flux += enhanceFactor*area*
			       propagator.getDF(2,logNuSurfaceEnergy, // nuTau
						inIceFlavor,inIceDoublet,logInIceEnergy);

		       }// loop over inIce Energy ends

		   }// loop over inIce Particle ends

		   yieldTable[iLogE] += flux*solidAngle*observationTime;

	       } // loop over Neutrino Energy ends

	   } // loop over Zenith Angle ends
       } // up and down

       yieldTableExists = true;

    }

}

