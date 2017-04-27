package iceCube.uhe.analysis;

import iceCube.uhe.particles.*;
import iceCube.uhe.analysis.*;

import java.io.*;
import java.util.*;
/**
   I3ParticleFlux calculates the detectable neutrino event intensity 
   at the Earth Surface using propagation matrix filled
   in each of I3Particle. You use Criteria class in the analysis
   package to difine the event criteria for "detection" by the IceCube.

   Propagation matrix must be filled by I3ParticlePropMatrixFiller
   a priori.

   A primary neutrino bulk is assumed to be consist of
   nu_e, nu_mu, nu_tau = 1:1:1. I3ParticlePropMatrixFiller object
   fills the matrix under this assumption. See I3ParticlePropMatrixFiller.java
   in the analysis package for details.

   Written by S. Yoshida April 08 2007

*/

public class I3ParticleFlux {


    /** Criteria of events */
    private Criteria criteria = null;
 
    // Container of the I3Particles
    private List i3particleList = null;
    private List  numberOfAllI3ParticlesList;
    private List  numberOfFilledI3ParticlesList;
    private int numberOfFilledI3Particles = 0;

    /** MCTruth Flag for the I3Particles.*/
    private boolean isMCTruth = false;

    /** Area where juliet particles are thrown in the MC simulation .
	Unit should be [cm^2].
     */
    protected double mcArea = 880.0*880.0*Math.PI*1.0e4; // 880*800*pi [m^2]
    /** Solid angle where juliet particles are thrown in 
	the MC simulation */
    protected double mcOmega = 4.0*Math.PI;
    /** Observation Time [sec]*/
    protected double observationTime = 1.0726387e7; // 124.148 days (IC9 data)

    /** bin widths for getInIceEffectiveArea() */
    protected double logEbinWidth = 20.0*Particle.getDeltaLogEnergy();
    protected double cosZenithBinWidth = 0.05;

    /** The lowest level Requirement of NDOMs. For saving memory.
	Only I3Particles to pass this condition will be stored
        in the container and subject to further analysis with
        higher level cuts defined by the Criteral class
     */
    protected static int minNDOMsToFill = 80;
    /** The lowest level Requirement of NPEs. For saving memory.
	Only I3Particles to pass this condition will be stored
        in the container and subject to further analysis with
        higher level cuts defined by the Criteral class
     */
    protected static double minLogNPEToFill = 4.0;


    /** Constructor */
    public I3ParticleFlux(InputStream in) throws IOException{
	i3particleList          = new LinkedList();
	numberOfAllI3ParticlesList = new LinkedList();
	numberOfFilledI3ParticlesList = new LinkedList();
	numberOfFilledI3Particles = 0;
	readI3Particles(in);
    }

    /** Constructor. Register Criteia also before reading out I3Paericles.
     This proceduer allows to save memories by removing the propagation matrix
    sotored in I3Particles that did not pass the criteria.
    */
    public I3ParticleFlux(InputStream in, Criteria criteria, boolean isMCTruth) throws IOException{
	i3particleList          = new LinkedList();
	numberOfAllI3ParticlesList = new LinkedList();
	numberOfFilledI3ParticlesList = new LinkedList();
	numberOfFilledI3Particles = 0;
	this.isMCTruth = isMCTruth;
	setCriteria(criteria);
	readI3Particles(in);
    }

    /** reading the I3Particle objects */
    public void readI3Particles(InputStream in) throws IOException{

	I3Particle iceParticle = null; 
	int numberOfI3MCdata = 0;
	while((iceParticle = I3ParticleInputStream.inputI3Particle(in))
	      !=null){

	    numberOfI3MCdata++;

	    if(isMCTruth) iceParticle.switchToMCTruth();
   	    else iceParticle.switchToReco();

	    // this event is filled with the propagation matrix
	    if(minNDOMsToFill<=iceParticle.getIceCubeData().getNDOMsLaunch() &&
	       minLogNPEToFill<=iceParticle.getIceCubeData().getLogBestNpe()){
		if(criteria!=null){
		    if(criteria.doesThisEventPass(iceParticle)){ // this event is needed
			i3particleList.add(iceParticle);
			numberOfFilledI3Particles++;
		    }
		}else{ // add this event no matter what
		    i3particleList.add(iceParticle);
		    numberOfFilledI3Particles++;
		}
	    }
	}
	numberOfAllI3ParticlesList.add(new Integer(numberOfI3MCdata));
	numberOfFilledI3ParticlesList.add(new Integer(numberOfFilledI3Particles));

	System.err.println("Number of I3Particles before applying any criteria" +
			   numberOfI3MCdata);
	System.err.println("Number of Filled I3Particles" + 
			   numberOfFilledI3Particles);
    }


    /** Set the MC area [cm^2] */
    public void setMCArea(double area){
	mcArea = area;
    }

    /** Set the MC solid angle [sr] */
    public void setMCSolidAngle(double omega){
	mcOmega = omega;
    }


    /** Set the MC solid angle [sec] */
    public void setObservationTime(double time){
	observationTime = time;
    }

    /** Switch to parameters concerned with MCtrue.*/
    public void switchToMCTruth(){
	isMCTruth = true;
    }

    /** Switch to parameters concerned with Reco results.*/
    public void switchToReco(){
	isMCTruth = false;
    }

    /** Sets the criteria on making Histogram */
    public void setCriteria(Criteria criteria){
	this.criteria = criteria;
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
	is not involved in the calculation (averageOverDecade=false)
	<pre>
	double logNeutrinoEnergy  : log10(Nu Energy at the earth surface [GeV])
	double numberOfEvents     : number of events in the IceCube
	Return dN/dLogE [/cm^2 sec sr]
	</pre>
    */
    public double getDFDLogE(double logNeutrinoEnergy, double numberOfEvents){

	return getDFDLogE(logNeutrinoEnergy,numberOfEvents,false);
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

	// Loop over I3Particles
	double  yieldOfIce3events= 0.0;
	ListIterator i3particleIterator = i3particleList.listIterator();
	int numberOfI3MCdata = 0;
	while(i3particleIterator.hasNext()){

	    numberOfI3MCdata++;

	    I3Particle iceParticle = (I3Particle)i3particleIterator.next();

	    if(isMCTruth) iceParticle.switchToMCTruth();
   	    else iceParticle.switchToReco();


	    // Criteria cut
	    if(criteria.doesThisEventPass(iceParticle)){

		double primaryFluxWeight = 
		    iceParticle.getMCPrimarySpectrumWeight();

		ListIterator numberIterator = numberOfAllI3ParticlesList.listIterator();
		ListIterator numFilled = numberOfFilledI3ParticlesList.listIterator();
		int numberOfAllI3Particles = 0;
		while(numberIterator.hasNext()){
		    numberOfAllI3Particles = ((Integer )numberIterator.next()).intValue();
		    if(numberOfI3MCdata <= 
		       ((Integer )numFilled.next()).intValue()) break;
		}

		yieldOfIce3events += iceParticle.getLogEnergyMatrix(iLogE)*
		    (1.0/(primaryFluxWeight*iceParticle.getDeltaLogEnergy()))*
		     mcArea*mcOmega*observationTime/(double )numberOfAllI3Particles;

		//if(logNeutrinoEnergy==10.0) System.out.println(numberOfI3MCdata + " " +
		//				       numberOfAllI3Particles);
		//if(iceParticle.getLogEnergyMatrix(iLogE)>0.0){
		//  System.out.println(Particle.particleName(iceParticle.getFlavor(),
		//				    iceParticle.getDoublet()) +
		//	       " logE=" + iceParticle.getLogEnergy() +
		//	       " logNuE " + logNeutrinoEnergy);
		//  double cosZ = -iceParticle.getDirectionInIceCubeCoordinate().getZ();
		//  double mtx =iceParticle.getLogEnergyMatrix(iLogE);
		//  System.out.println("   logNpe=" + 
		//       iceParticle.getIceCubeData().getLogBestNpe() +
		//       " cos(zenith)=" + cosZ +
		//       " mtx= " + mtx);
			       
		//}

	    }else{
		System.err.println("This " + Particle.particleName(iceParticle.getFlavor(),iceParticle.getDoublet()) +
				   " doesn't pass the criteria");
	    }
	}// loop of I3Particles ends

	//System.err.println(logNeutrinoEnergy + " yield [km^2 yr sr]=" + 
	//		   yieldOfIce3events/1.0e10/3.1536e7);
	//System.err.println(" number of events passed(" + 
	//	   numberOfPassedEvents + ") out of " +
	//	   numberOfAllI3Particles);

	return yieldOfIce3events/3.0; // devided by number of nu flavor.

    }

    /**
       Outputs the in-ice effective area(E-inice, cosZenith) determined by
       events passing the Criteria class. Counting number of events (I3Particle)
       satisfying the criteria and calculate the passing rate. Mutiplication with
       the MC area gives the effective area [cm^2]
       <pre>

       double logEnergy   : log10(inice Energy [GeV])
       double cosZenith   : cos(Zenith angle at the IceCube Depth)
       int flavor         : flavor of the in-ice particle. Defined by Particle class
       int doublet        : doublet of the in-ice particle. Defined by Particle class

       </pre>
       The logE binWidth is 10 times Particle.getDeltaLogEnergy() = 0.1 decade;
       The cosZenith bin is 0.1
    */

    public double getInIceEffectiveArea(double logEnergy, double cosZenith, 
					int flavor, int doublet){

	double dOmega = 2.0*Math.PI*cosZenithBinWidth/mcOmega;

	// Loop over I3Particles
	double  area = 0.0;
	ListIterator i3particleIterator = i3particleList.listIterator();
	int numberOfI3MCdata = 0;
	while(i3particleIterator.hasNext()){

	    numberOfI3MCdata++;

	    I3Particle iceParticle = (I3Particle)i3particleIterator.next();

	    iceParticle.switchToMCTruth();

	    // particle spiece
	    int inIceFlavor = iceParticle.getFlavor();
	    int inIceDoublet = iceParticle.getDoublet();
	    // energy
	    double logInIceEnergy = iceParticle.getLogEnergy();
	    // direction
	    double cosInIceZenith = -iceParticle.getDirectionInIceCubeCoordinate().getZ();
                                                // Reversed vector

	    // select the events within the given energy/zenith and the particle range
	    if((inIceFlavor == flavor && inIceDoublet == doublet) &&
	       ((cosZenith-0.5*cosZenithBinWidth)<= cosInIceZenith ) &&
	       (cosInIceZenith<(cosZenith+0.5*cosZenithBinWidth)) &&
	       ((logEnergy-0.5*logEbinWidth)<=logInIceEnergy) &&
	       (logInIceEnergy<(logEnergy+0.5*logEbinWidth))){

		// Criteria cut
		if(isMCTruth) iceParticle.switchToMCTruth();
		else iceParticle.switchToReco();

		if(criteria.doesThisEventPass(iceParticle)){

		    double primaryFluxWeight = 
			iceParticle.getMCPrimarySpectrumWeight();

		    ListIterator numberIterator = numberOfAllI3ParticlesList.listIterator();
		    ListIterator numFilled = numberOfFilledI3ParticlesList.listIterator();
		    int numberOfAllI3Particles = 0;
		    while(numberIterator.hasNext()){
			numberOfAllI3Particles = ((Integer )numberIterator.next()).intValue();
			if(numberOfI3MCdata <= 
			   ((Integer )numFilled.next()).intValue()) break;
		    }

		    area += (1.0/(primaryFluxWeight*logEbinWidth*dOmega))*
			mcArea/(double )numberOfAllI3Particles;
		}
	    }
	}

	return area;

    }

}

	
