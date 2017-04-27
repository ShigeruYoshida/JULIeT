package iceCube.uhe.analysis;

import iceCube.uhe.particles.*;
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

   Written by S. Yoshida April 11 2007

*/

public class PropagationMatrixFlux {

    protected Particle nuE;  // electron neutrino as a primary particle
    protected Particle nuMu; // muon neutrino as a primary particle
    protected Particle nuTau;// tau neutrino as a primary particle

    protected PropagationMatrixFactory matrix = null;
    double[] yieldTable;
    protected boolean yieldTableExists = false;


    /** Observation Time [sec]*/
    protected double observationTime = 1.0726387e7; // 124.148 days (IC9 data)

    /** List to store inIceParticle's flavor */
    protected List inIceFlavorList = null;
    /** List to store inIceParticle's doublet */
    protected List inIceDoubletList = null;
    /** List to store inIceParticle's effective area */
    protected List inIceAreaList = null;

    /** Constructor. Generate nu-e, nu-mu, and nu-tau objects
	as primary neutrinos at the earth surface.
	PropagationMatrixFactory is also generated.
    */
    public PropagationMatrixFlux(){
	nuE = new Particle(0,0);// electron neutrino as a primary particle
	nuMu = new Particle(1,0);// muon neutrino as a primary particle
	nuTau =new Particle(2,0);// tau neutrino as a primary particle

	inIceFlavorList = new LinkedList();
	inIceDoubletList = new LinkedList();
	inIceAreaList = new LinkedList();

	matrix = new PropagationMatrixFactory();

	yieldTable = new double[Particle.getDimensionOfLogEnergyMatrix()];
	for(int i=0;i<Particle.getDimensionOfLogEnergyMatrix(); i++) 
	    yieldTable[i]=0.0;
    }

    /** Constroctor for the subclass. Do not instance the PropagationMatrixFactory */
    protected PropagationMatrixFlux(boolean nomatrix){
	inIceFlavorList = new LinkedList();
	inIceDoubletList = new LinkedList();
	inIceAreaList = new LinkedList();

	yieldTable = new double[Particle.getDimensionOfLogEnergyMatrix()];
	for(int i=0;i<Particle.getDimensionOfLogEnergyMatrix(); i++) 
	    yieldTable[i]=0.0;
    }

    /** add flavor and doublet of inIce particle to consider.
	It also generate EffAreaTable(flavor doublet) object
	for caulculating the effecrtive area in the yield calculation.
     */
    public void addInIceParticle(int flavor, int doublet) throws IOException{
	inIceFlavorList.add(new Integer(flavor));
	inIceDoubletList.add(new Integer(doublet));
	inIceAreaList.add(new EffAreaTable(flavor,doublet));
    }

    /** Set the MC solid angle [sec] */
    public void setObservationTime(double time){
	observationTime = time;
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

			   flux += area*matrix.getDF(nuE,iceParticle);
			   flux += area*matrix.getDF(nuMu,iceParticle);
			   flux += area*matrix.getDF(nuTau,iceParticle);

		       }// loop over inIce Energy ends

		   }// loop over inIce Particle ends

		   yieldTable[iLogE] += flux*solidAngle*observationTime;

	       } // loop over Neutrino Energy ends

	   } // loop over Zenith Angle ends
       } // up and down

       yieldTableExists = true;

    }

    /** Calculate the Neutrino yeild [cm^2 sec sr] at the surface
	to give numberOfEvents you set in the argument.
	It uses the propagation matrix filled in I3Particle's.
	<pre>
	double logNeutrinoEnergy  : log10(Nu Energy at the earth surface [GeV])
	Return yield [cm^2 sec sr]
	</pre>
	Note : Yield given here is all-nu_flavor-summed value.
    */
    public double getYield(double logNeutrinoEnergy){

	int iLogE = (int )((logNeutrinoEnergy-Particle.getLogEnergyMinimum())/
			   Particle.getDeltaLogEnergy());
	// check the energy range
	if(iLogE<0 || Particle.getDimensionOfLogEnergyMatrix() <= iLogE){
	    System.err.println(" Neutrino Energy out of the range! (" +
			       iLogE + ")");
	    return 0.0;
	}
	if(!yieldTableExists){ // yield table has not been generated!
	    System.err.println(" You have to call calculateYield() first!");
	    System.exit(0);
	}

	double yieldOfIce3events = 
	    yieldTable[iLogE]/3.0; //devided by a three nu flavor

	return yieldOfIce3events;
    }

    /** Calculate the Neutrino flux at the surface
	to give numberOfEvents you set in the argument.
	It uses the propagation matrix filled in I3Particle's
	and the quasi-differential method based on dF/dLogE.
	<pre>
	double logNeutrinoEnergy  : log10(Nu Energy at the earth surface [GeV])
	double numberOfEvents     : number of events in the IceCube
	boolean averageOverDecade : true - average neutrino yield over decade of E
                                  : false - use the yield of logNeutrinoEnergy only
                                  : as the Auger/Rice people introduced.
	Return dN/dLogE [/cm^2 sec sr]
	</pre>
    */
    public double getDFDLogE(double logNeutrinoEnergy, double numberOfEvents,
			     boolean averageOverDecade){

	if(numberOfEvents<=0.0) return 0.0; // null flux

	double yieldOfIce3events;
	if(!averageOverDecade){ // The Auger/Rice method
	    yieldOfIce3events = getYield(logNeutrinoEnergy);
	}else{ // Averaging the yield over a decade of energy
	    yieldOfIce3events = 0.0;
	    double logE = logNeutrinoEnergy-0.5; // -0.5 decade
	    double logEmax = logNeutrinoEnergy+0.5; // +0.5 decade
	    double epsilon = 1.0e-4;int integralBinWidth = 0;
	    while(logE<=logEmax){
		double yield = getYield(logE+epsilon);
		if(yield>0.0){
		    yieldOfIce3events += yield;
		    integralBinWidth++;
		}
		logE += Particle.getDeltaLogEnergy();
	    }
	    yieldOfIce3events = yieldOfIce3events/(double )integralBinWidth;
	}
	double flux = numberOfEvents/yieldOfIce3events;

	return flux;
    }

    /** Calculate the Neutrino flux at the surface
	to give numberOfEvents you set in the argument.
	It uses the propagation matrix filled in I3Particle's
	and the quasi-differential method based on dF/dLogE.
	note:The energy decade average over a decade \int yield(logE) dLogE
	is involved in the calculation (averageOverDecade=true)
	<pre>
	double logNeutrinoEnergy  : log10(Nu Energy at the earth surface [GeV])
	double numberOfEvents     : number of events in the IceCube
	Return dN/dLogE [/cm^2 sec sr]
	</pre>
    */
    public double getDFDLogE(double logNeutrinoEnergy, double numberOfEvents){

	return getDFDLogE(logNeutrinoEnergy,numberOfEvents,true);
    }

    /** Calculate the Neutrino flux at the surface
	to give numberOfEvents, but the yield [cm^2 sec sr]
        given in the argument is added up to calculate the flux.
        This method may be used when you combine the estimations
        from the numerically calculated effective area.
    */
    public double getDFDLogE(double logNeutrinoEnergy, double yield, 
			     double numberOfEvents){

	if(numberOfEvents<=0.0) return 0.0; // null flux

	double yieldOfIce3events = yield + getYield(logNeutrinoEnergy);
	double flux = numberOfEvents/yieldOfIce3events;

	return flux;
    }	


}

