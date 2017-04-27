package iceCube.uhe.interactions;

import iceCube.uhe.particles.*;
import iceCube.uhe.points.*;
import iceCube.uhe.event.*;
import iceCube.uhe.decay.*;
import iceCube.uhe.interactions.*;
import numRecipes.*;

import java.util.*;

/** <pre>
    The GlashowResonanceBase class to treat Glashow Resonance for JulietEventGenerator class and 
    RunManager class. This class inherits the MonteCarloBase class.
    </pre>
*/

public class GlashowResonanceBase extends MonteCarloBase {

    /** interactions */
    private Interactions interactions;
    
    /** ParticlePoint object to define the trajectory and the medium
        in the propagation */
    private ParticlePoint point;
    private int           materialNumber = 0;       // Ice->0, Rock->1;

    /** Ratio of anti-electron neutrino 
	Only anti-electron neutrinos can interact through Glashow Resonance, 
	although juliet can't distinguish anti-neutrino. 
    */
    private final double antiNuERatio = 0.5; 

    /** In order to save CPU time, we increase neutrino cross section by
        this factor. As long as the meanfree path is by far shorter than
        the propagation length (presumably 1km), this is equivallent
        to the case when one neutrino particle represents multiple neutrinos
        whose number is equal to this factor. You (or RunManager class) have to\ devide
        this factor later to compensate this "artificial" enhancement
        of the neutrino cross section.
    */
    public static int neutrinoFactor = 1;

    /** Propagateing particle is electron neutrino */
    private final int propFlavor  = 0;
    private final int propDoublet = 0;
	
    /** flavor of produced particle is electron neutrino */
    private int producedFlavor = 0;

    /** Default Constructor of GlashowResonanceBase. */
    public GlashowResonanceBase(){
	this(0,0); // default produced electron in ice!
    }
    
    /** Constructor of GlashowResonanceBase. */
    public GlashowResonanceBase(int flavor, int mediumID){

        this.producedFlavor = flavor;
	this.materialNumber = mediumID;
	System.err.println("GlashowResonance produced flavor is " + producedFlavor +
			   ", medium ID is " + materialNumber);

	point = new ParticlePoint(0.0, 5.0*Math.PI/180.0, materialNumber);

	if(producedFlavor != 3) interactions = new GlashowResonanceLeptonic(point, producedFlavor);
	else interactions = new GlashowResonanceHadronic(point);

        for(int iLogE=0; iLogE<Particle.getDimensionOfLogEnergyMatrix(); iLogE+=100){
	    interactions.setIncidentParticleEnergy(iLogE);

	    System.err.println("total sigma [" + iLogE + "] = " + 
			       interactions.getSigma());
	}
    }
    
    /** Get pathlength by random number. **/
    public double getPathLength(int iLogE, RandomGenerator rand) {

	double logEnergy = Particle.getLogEnergyMinimum() + iLogE*Particle.getDeltaLogEnergy();
        return getPathLength(logEnergy, rand);

    }

    /** Get pathlength by random number. **/
    public double getPathLength(double logEnergy, RandomGenerator rand) {

        double r = rand.GetRandomDouble();
	if(r> antiNuERatio){ // this is nu-e not a nu-e-bar!
	    return Double.POSITIVE_INFINITY;
	}else{ // case of anti nu-e 
	    r = rand.GetRandomDouble();
	    // Set incident particle energy
	    double energy = Math.pow(10.0,logEnergy);
	    interactions.setIncidentParticleEnergy(energy);

	    double totalsigma = interactions.getSigma(); // get total cross section
	    double path       = -Math.log(1.0-r)/totalsigma;       // path length
	    return path;
	}
	
    }


    /** Get pathlength for neutrino by given logEnergy bin. **/
    public double getNeutrinoPathLength(int iLogE, RandomGenerator rand) {

	double logEnergy = Particle.getLogEnergyMinimum() + 
	    Particle.getDeltaLogEnergy()*(double )iLogE;
        return getNeutrinoPathLength(logEnergy, rand);

    }

    /** Get pathlength for neutrino. 
	You can get exact path length from a random number
	because of simple energy dependence of Sigma.  
    **/
    public double getNeutrinoPathLength(double logEnergy, RandomGenerator rand) {

        double r = rand.GetRandomDouble();

	if(r> antiNuERatio){ // this is nu-e not a nu-e-bar!
	    return Double.POSITIVE_INFINITY;
	}else{ // case of anti nu-e 
	    r = rand.GetRandomDouble();
	    // Set incident particle energy
	    double energy = Math.pow(10.0,logEnergy);
	    interactions.setIncidentParticleEnergy(energy);

	    double totalsigma = 
		interactions.getSigma()*neutrinoFactor; //crosssection*neutrinoFactor;
	    double path       = -Math.log(1.0-r)/totalsigma;       // path length
	    return path;
	}
    }

    /** Get produced log energy. 
	You can get exact energy from a random number 
	because of simple energy dependence of dSigma/dy, 
	where y = 1 - E_{l^{-1}}/E_{\bar{\nu_e}}. 
    */
    public double getProducedEnergy(int iLogE, RandomGenerator rand){
	
	double logEnergy = Particle.getLogEnergyMinimum() + iLogE*Particle.getDeltaLogEnergy();
        return getProducedEnergy(logEnergy, rand);

    }
    
    /** Get produced log energy. 
	You can get exact energy from a random number 
	because of simple energy dependence of dSigma/dy, 
	where y = 1 - E_{l^{-1}}/E_{\bar{\nu_e}}. 
    */
    public double getProducedEnergy(double logEnergy, RandomGenerator rand){
	
	double producedLogEnergy;
	
	if(producedFlavor != 3){
	    /** Lectonic glashow resonance. 
		   dSigma/dy = 3.0 * Sigma * y^2
		therefore, 
		   random number = y^3
	    */
	    double r = rand.GetRandomDouble();
	    producedLogEnergy = logEnergy + Math.log(r)/3.0;
	    
	}else{
	    /** Hadronic glashow resonance.
		All the incident energy is immediately deposited as a hadronic cadcade.
	     */
	    producedLogEnergy = logEnergy;

	}
	
        return producedLogEnergy;

    }
    
    /** Get the flavor of the particle propagateing */
    public int getPropFlavor() { return propFlavor; }

    /** Get the doublet of the particle propagateing */
    public int getPropDoublet() { return propDoublet; }

    /** Get the flavor of the produced particle */
    public int getProducedFlavor() { return producedFlavor; }

    /** Set the flavor of the produced particle */
    public void setProducedFlavor(int flavor) {	this.producedFlavor = flavor; }

    /** Get the name of the interaction */
    public String getInteractionName() {

        String nameOfInteraction = interactions.interactionName();
        return nameOfInteraction;

    }

    /** Get type of the interaction (Interaction->0; Decay->1) */
    public int getTypeOfInteraction() {

        return 0;

    }

}
