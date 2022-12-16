package iceCube.uhe.event;

import iceCube.uhe.particles.*;
import iceCube.uhe.interactions.*;
import iceCube.uhe.decay.*;
import iceCube.uhe.event.*;
import iceCube.uhe.points.*;
import numRecipes.*;

import java.util.*;
import java.io.*;

/**
   <pre>

   This class defines the behavior of an event running
   in the madium (rock/ice). This is the fundamenthal class
   for propagating particles with the Monte Carlo Method.

   Running particles by numerically solving the transport
   equations, use PropagationMatrix class in the iceCube.uhe.propagation
   package.

   The relevant interactions are registered in the List
   inside the constructor. 

   Modified for Glashow Resonance by M.Ono December 2 2007.
   </pre>
*/

public class Event {

    /** Common logarithm */
    static final double ln10 = Math.log(10.0);

    /** Round off error */
    static final double epsilon = 10e-5; 

    /** List to contain the MonteCalroBase objects 
	invovled in the propagation */
    List MonteCarloList;
    ListIterator MonteCarloIterator = null;

    /** Step size[g/cm^2] */
    private double stepDx;

    /** Particle Point object to define the trajectory and the medium
	in the propagation */
    ParticlePoint point;

    /** MassNumber in the medium */
    double massNumber = 0.0;

    /** Propagating Particle object */
    Particle propParticle;

    /** MonteCarloBase object currently in play */
    MonteCarloBase mcBaseInPlay = null;

    /** Cascade energy */
    private double cascadeEmgEnergy;
    private double cascadeHadronEnergy;

    /** Dimension of produced log energy */
    private int expandedDim;

    /** Constructor. Allocating Lists to the MonteCarloBase,
	and register the Particle and ParticlePoint classes. */

    public Event(MonteCarloBase[] mcBases, Particle p, ParticlePoint s) {
	this.propParticle = p;
	this.point        = s;

	cascadeEmgEnergy    = 0.0;
	cascadeHadronEnergy = 0.0;

	setMassNumber();

	MonteCarloList = new LinkedList();
	registerMonteCarloBase(mcBases);

    }

    /** Calculate the mass number in the current medium */
    void setMassNumber(){
	massNumber = 0.0;
        for(int i=0;i<point.NumberOfSpecies[point.getMaterialNumber()];i++){
            massNumber += point.getNumberOfAtoms(i)*point.getAtomicNumber(i);
        }
    }

    /** Get the mass number in the current madium */
    double getMassNumber(){
	return massNumber;
    }

    /** Register the MonteCarloBase objects involved in the propagation */
    public void registerMonteCarloBase(MonteCarloBase[] mcBases){
	
	// Add MonteCarloBase objects
	for(int i = 0; i<mcBases.length; i++){
	    MonteCarloList.add(mcBases[i]);
	}
	MonteCarloIterator = MonteCarloList.listIterator();
    }
	 
    /** Calculate the stepsize determined by sampling the interaction
	points for all the interaction channels registred.
	The shortest path length is picked up and returned.
	The determined interaction is put in the class valuable
	"mcBaseInPlay".
    */
    public double getPhysicalPathLength(RandomGenerator rand){

	double logEnergy = propParticle.getLogEnergy();
	if(logEnergy<propParticle.getLogEnergyMinimum()) logEnergy = propParticle.getLogEnergyMinimum();

	int iLogE = (int )((logEnergy - propParticle.getLogEnergyMinimum())/
	                           propParticle.getDeltaLogEnergy());
	if(iLogE<0) iLogE = 0;

	double minPathLength = 1.0e25;
        mcBaseInPlay = null; // Initialization

	// Set the List Iterator to the first element
	MonteCarloIterator = MonteCarloList.listIterator();
	
	while(MonteCarloIterator.hasNext()) {
	    
	    MonteCarloBase mcBase = ((MonteCarloBase )(MonteCarloIterator.next()));
	    // Check if this interaction involves the propagating particle
	    if((mcBase.getPropFlavor()==propParticle.getFlavor()) &&
	       ((mcBase.getPropDoublet()==propParticle.getDoublet()))){
		double pathLength = 0.0;

		// Get pathLength
		if(mcBase.getTypeOfInteraction()==0) { // Interactions (not Decay)
		    if(mcBase.getPropDoublet()!=0) {  // Charged Lepton
			pathLength = massNumber/point.NA*
			    mcBase.getPathLength(iLogE, rand);
		    }else {                            // Neutrino
			/** For Glasow Resonance */
			//pathLength = mcBase.getNeutrinoPathLength(iLogE, rand)/point.NA;
			pathLength = mcBase.getNeutrinoPathLength(logEnergy, rand)/point.NA;
		    }
		}
		else if(mcBase.getTypeOfInteraction()==1){  // Decay
		    pathLength = point.getMediumDensity()*mcBase.getPathLength(iLogE,rand);
		}
		
		// compare the interaction Pathlengths
		if(pathLength<minPathLength) {
 		    minPathLength = pathLength;
		    mcBaseInPlay  = mcBase;
		}	
	    }
	}
	return minPathLength;
    }

    /** Get the Interaction's name which has just interacted with your particle */
    public String interactionsNameInPlay(){
	if(mcBaseInPlay!=null) return mcBaseInPlay.getInteractionName();
	else {
	    String message = "No interaction has occured.";
	    return message;
	}
    }

    /** Get the produced particle's Flavor */
    public int getFlavorByInteractionsInPlay(){
	if(mcBaseInPlay!=null) return mcBaseInPlay.getProducedFlavor();
	else {
	    return -1;
	}
    }

    /** The propagation particle now collides with Nuclear/Nucleon and change its energy
	via intertactions determined in GetPhysicalPathLength().
	The energy of produced particle (which may be equal to cascade energy 
	depending on the interaction) is returned. 

	if returned 0, it implies that the primary particle's kinetic energy
	becomes 0. The particle must stop at this moment.
    */
    public double collideNow(RandomGenerator rand){

	double logEnergy  = propParticle.getLogEnergy();
	double propEnergy = propParticle.getEnergy();
	int    iLogE      = (int )((logEnergy - propParticle.getLogEnergyMinimum())/
					 propParticle.getDeltaLogEnergy());
	if(iLogE<0) iLogE = 0;
	if(logEnergy<propParticle.getLogEnergyMinimum()) {
	    logEnergy = propParticle.getLogEnergyMinimum();
	}

	double producedLogEnergy = mcBaseInPlay.getProducedEnergy(logEnergy, rand);
	double producedEnergy   = Math.pow(10.0,producedLogEnergy);

    // When the particle's energy is below the minimum energy
    // we're pretending it has the minimum energy instead.
    // When sampling an energy loss were using the min energy's cross section
    // Scale down the energy of the loss to the particles actual energy
    // to improve the quality of this approximation
    if(propParticle.getLogEnergy() < propParticle.getLogEnergyMinimum()){
        producedEnergy *= propParticle.getEnergy() / Math.pow(10, propParticle.getLogEnergyMinimum());
        producedLogEnergy += propParticle.getLogEnergy() - propParticle.getLogEnergyMinimum();
    }

	// Flag for checking Neutrino interactions
	boolean compareNC = mcBaseInPlay.getInteractionName().startsWith("Neutrino-Nuclen Neutraled ");
	boolean compareCC = mcBaseInPlay.getInteractionName().startsWith("Neutrino-Nuclen Charged ");
	/** For Glashow Resonance */
	boolean compareGR = mcBaseInPlay.getInteractionName().startsWith("Glashow Resonance ");

	// change propParticle into a produced particle but its energy is not changed.
	if(mcBaseInPlay.getTypeOfInteraction() == 1){  // decay
	    if(mcBaseInPlay.getProducedFlavor() == 1){ //Tau to Mu Decay<-change particle
		changeParticle(mcBaseInPlay.getProducedFlavor(), 1); // ->muon
		//changeParticle(1, 1); // ->muon
	    }
	}else if(compareNC){   // Neutral Current
	    changeParticle(mcBaseInPlay.getPropFlavor(), 0);  // -> Neutrino
	}else if(compareCC){              // Charged Current
	    changeParticle(mcBaseInPlay.getPropFlavor(), 1);  // -> Charged Lepton

	    /** For Glashow Resonance -begin */
	}else if(compareGR){              // Glashow Resonance Leptonic
            changeParticle(mcBaseInPlay.getProducedFlavor(), 1);  // -> Charged Lepton or Pion
        }
	    /** For Glashow Resonance -end */

	
	if(producedEnergy>=propEnergy) {  // propParticle loose all the energy.
	    producedEnergy   = propEnergy - propParticle.getMass();
	    producedLogEnergy = Math.log(producedEnergy)/ln10;
	}
	double recoilEnergy    = propEnergy-producedEnergy;
	double recoilLogEnergy = Math.log(recoilEnergy)/ln10;

	if(recoilEnergy<propParticle.getMass() + epsilon){  // Propagating particle must stop 
	                                           // but it has to have its own mass as its energy.
	    recoilEnergy = propParticle.getMass();
	    recoilLogEnergy = Math.log(recoilEnergy)/ln10;
	}
	    
	propParticle.putLogEnergyMinimum(-10.0); // Allows energy down to almost 0 GeV
	propParticle.putLogEnergy(recoilLogEnergy);
	propParticle.putLogEnergyMinimum(5.0);
	propParticle.putEnergy(recoilEnergy);

	// cascade energy
	int    dim              = Particle.getDimensionOfLogEnergyMatrix();
	double logEnergyMinimum = Particle.getLogEnergyMinimum();
	expandedDim = dim + 
	    (int )((logEnergyMinimum-InteractionsBase.getLogEnergyProducedMinimum())/
			      Particle.getDeltaLogEnergy());

	/**
	   For Glashow Resonance
	   If GR occered, produced particle is neutrino which cannot produce cascades. 
	*/
	if(getFlavorByInteractionsInPlay()==0 && !compareGR) {
	    cascadeEmgEnergy += producedEnergy; // <- energy/1km
	}
	else if(getFlavorByInteractionsInPlay()==3) {
	    cascadeHadronEnergy += producedEnergy; // <- energy/1km   
	}

	return producedEnergy;
    
    } 

    public double getCascadeEmgEnergy() {
	return cascadeEmgEnergy;
    }

    public double getCascadeHadronEnergy() {
	return cascadeHadronEnergy;
    }

    public double getCascadeTotalEnergy() {
	return (cascadeEmgEnergy + cascadeHadronEnergy);
    } 

    /** Change the particle to the different kind. This method
	must be called when decay such as tau to mu takes place.
	<pre>
	newFlavor   : flavor of a new particle.
	newDoublet  : doublet of a new particle.
	</pre>
    */
    public void changeParticle(int newFlavor, int newDoublet){
	double energy    = propParticle.getEnergy();
	double logEnergy = propParticle.getLogEnergy();
	propParticle = null;
	propParticle = new Particle(newFlavor, newDoublet);
	propParticle.putEnergy(energy);
	propParticle.putLogEnergy(logEnergy);
    }

    /** Set the step size for traceParticle() */
    public void setStepDx(double dx){
	stepDx = dx;
    }

    /** Get the step size for traceParticle() */
    public double getStepDx(){
	return(stepDx);
    }

    /** Trace particle running by a given pathlength [g/cm^2]*/
    public void traceParticle(double pathLength){
	double l = point.getParticleLocation(); 
	           // The current particle location along the trajectory.
	double deltaL;
	double xSum = 0.0;
	while(xSum<pathLength){
	    deltaL = stepDx/point.getMediumDensity(); 
	    l     += deltaL; 
	    xSum  += stepDx;
	    point.setParticleLocation(l);
	}
    }

}





