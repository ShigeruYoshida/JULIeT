package iceCube.uhe.interactions;

import iceCube.uhe.particles.*;
import iceCube.uhe.event.*;
import iceCube.uhe.decay.*;
import iceCube.uhe.interactions.*;
import numRecipes.*;

import java.util.*;

/** <pre>
    An electron once generated by nu-e charged current interaction is subject to
    immediate electromagnetic cascades. The Event class under the current version
    is not able to simulate cascades features themselves. Instead, it simply records
    primary energy of the electron (i.e. primary energy of the emg cascades) and
    its generated location, and put an end to the particle tracing. In order to
    do so, this class provides a hypthetical "electron-to-electron interaction"  
    where all the primary energy is channeled into "produced" electron
    with pathlength of 0. By calling this class right after an electron is generated
    by interactions such as nu-e charged current interactions, all the energy 
    is deposited at the same location and the event sees its end.
    </pre>
*/

public class ElectronBase extends MonteCarloBase {

    /** Particle obejct */
    Particle p;

    /** Minimum log energy of produced particles */
    private static final double initialLogEnergyProducedMinimum = 2.0;
    private double logEnergyProducedMinimum = initialLogEnergyProducedMinimum; 
    
    /** Constructor */
    public ElectronBase(Particle p){
	if(p.getDoublet()==1 && p.getFlavor()==0){ // Requires electrons
	    this.p = p;
	}else{
            System.err.println("This particle " + 
			       p.particleName(p.getFlavor(), p.getDoublet()) + 
			       " is not an Electron!");
            System.exit(0);
        }

    }
    
    public static double getLogEnergyProducedMinimum(){
	return initialLogEnergyProducedMinimum;
    }

    /** This is the "dummy" method */
    private void setCumulativeTable(InteractionsMatrix interactMtx){
	System.err.println("This method does not have to be called (^^;");
    }

    /** Get pathlength, but this always returns 0.0  **/
    public double getPathLength(int iLogE, RandomGenerator rand) {
	double path = 0.0;  // electron deposits its energy to "produced" electron without propagating
	//System.out.println("ElectronBase path=" + path);
	return path;
    }

    /** Get pathlength, but this always returns 0.0  **/
    public double getPathLength(double logEnergy, RandomGenerator rand) {
	int iLogE = (int )((logEnergy-Particle.getLogEnergyMinimum())/Particle.getDeltaLogEnergy( ));
        return getPathLength(iLogE, rand);
    }
    
    /** This is the "dummy" method */
    public double getNeutrinoPathLength(int iLogE, RandomGenerator rand) {
	return 1.0e25;
    }

    /** This is the "dummy" method */
    public double getNeutrinoPathLength(double logEnergy, RandomGenerator rand) {
	int iLogE = (int )((logEnergy-Particle.getLogEnergyMinimum())/Particle.getDeltaLogEnergy( ));
        return getNeutrinoPathLength(iLogE, rand);
    }

    /** Get produced log energy. This method always returns same as 
        incident energy **/
    public double getProducedEnergy(int iLogE, RandomGenerator rand){

	double producedLogEnergy = (double)iLogE * Particle.getDeltaLogEnergy() + 
	                           logEnergyProducedMinimum;        // same as incident energy
	return producedLogEnergy;

    }

    /** Get produced log energy. This method always returns same as 
        incident energy in this hypothetical interaction. **/
    public double getProducedEnergy(double logEnergy, RandomGenerator rand){
        return logEnergy; 

    }

    /** Get the flavor of the particle propagating */
    public int getPropFlavor() {

	int propFlavor = p.getFlavor();  
	return propFlavor;

    }

    /** Get the doublet of the particle propagating */  
    public int getPropDoublet() {

	int propDoublet = p.getDoublet();
	return propDoublet;

    }

    /** Get the flavor of the produced particle */
    public int getProducedFlavor() {

	int producedFlavor = p.getFlavor(); // returns same flavor as getPropFlavor()
	return producedFlavor;
    
    }

    /** Get the name of the interaction */
    public String getInteractionName() {

	String nameOfInteraction = "Electron's energy deposit to electromagnetic cascade";
	return nameOfInteraction;

    }

    /** get Type Of Interaction (Interaction->0; Decay->1) */
    public int getTypeOfInteraction() {
	
	return 0;

    }


}
