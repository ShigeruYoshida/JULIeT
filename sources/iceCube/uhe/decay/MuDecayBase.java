package iceCube.uhe.decay;

import iceCube.uhe.decay.*;
import iceCube.uhe.event.*;
import iceCube.uhe.particles.*;
import iceCube.uhe.interactions.*;
import iceCube.uhe.propagation.*;
import numRecipes.*;

import java.util.*;

/** <pre>
    The MuDecayBase class to treat mu decay same as interactions for Event class.
    This class and other "Base" classes inherit the MonteCarloBase class.

    This class has MuDecayYMatrix object. It treats energy transfer as dN/dlogY,
    Y = Edecay/Emu.
    </pre>
*/

public class MuDecayBase extends MonteCarloBase {

    private MuDecayYMatrix  muDecayMtx;
    private Particle        p;
    private int             incomingLogE;
    private int             dim; 
    private double          random;
    private double[]        cumulativeTable;


    /** Constructor for making the tables of lifetime and cumulative table. */ 
    public MuDecayBase(MuDecayYMatrix muDecayMtx) {
	
	    this.p               = muDecayMtx.p;
	if(p.getDoublet()==1 && p.getFlavor()==1){ // Requires mu's
	    this.muDecayMtx      = muDecayMtx;
	    dim                  = p.getDimensionOfLogEnergyMatrix();
	    for(int kLogE=0; kLogE<dim; kLogE++){
		// Make LifeTime [sec] Matrix
		muDecayMtx.setLifeTimeMatrix(kLogE);
	    }
	    this.cumulativeTable = new double[dim];
	    this.setCumulativeTable(muDecayMtx);
	}else{
            System.err.println("This particle " + 
			       p.particleName(p.getFlavor(), p.getDoublet()) + 
			       " is not MUONs!!");
            System.exit(0);
        }
    }

   /** Make a cumulative table of differential cross section.
        The elements are normalized to 1. */
    public void setCumulativeTable(MuDecayYMatrix muDecayMtx){
	double sum =0.0;
	// MakeDecay Matrix for y=Ee/Emu 
	for(int iLogY=0; iLogY<dim; iLogY++){
		muDecayMtx.setMuDecayMatrix(iLogY);
		cumulativeTable[iLogY] = (iLogY>0?cumulativeTable[iLogY-1]:0.0) +
                                           muDecayMtx.getMuToEDecayMatrix(iLogY);
	    }

	System.out.println("Sum of cumTable for Mu Decay = " + 
			   cumulativeTable[dim-1] );
	    
    }


    /** Get pathlength by random number. **/
    public double getPathLength(int iLogE, RandomGenerator rand) {
	
	double r = rand.GetRandomDouble();
	double path = -muDecayMtx.getLifeTimeMatrix(iLogE)*PropagationMatrix.c*Math.log(1.0-r);

	return path;

    }

    public double getPathLength(double logEnergy, RandomGenerator rand) {
	
	int iLogE = (int )((logEnergy-p.getLogEnergyMinimum())/p.getDeltaLogEnergy());
        return getPathLength(iLogE, rand);

    }

    /** This is a dummy method because this class extends MonteCarloBase.class */
    public double getNeutrinoPathLength(int iLogE, RandomGenerator rand) {

	return 1.0e25;

    }	

    /** This is a dummy method because this class extends MonteCarloBase.class */
    public double getNeutrinoPathLength(double logEnergy, RandomGenerator rand) {

	return 1.0e25;

    }	
    
    /** Get produced log energy. In order to decide the value of log energy in a bin,
        use a random number **/
    public double getProducedEnergy(int iLogE, RandomGenerator rand){

        double y = 0.0;            
	double r = rand.GetRandomDouble();

	for(int kLogY=0; kLogY<dim; kLogY++){
            if(r<cumulativeTable[kLogY]){
                y = (double )kLogY;
                break;
            }
	}

	// value in a bin
	r = rand.GetRandomDouble();
	double logEnergy = (double)iLogE*Particle.getDeltaLogEnergy() + Particle.getLogEnergyMinimum();
	double producedLogEnergy = -(y+r)*Particle.getDeltaLogEnergy()+logEnergy;

	return producedLogEnergy;

    }

    /** Get produced log energy. In order to decide the value of log energy in a bin,
        use a random number **/
    public double getProducedEnergy(double logEnergy, RandomGenerator rand){

	double y = 0.0;            
	double r = rand.GetRandomDouble();
	
	for(int kLogY=0; kLogY<dim; kLogY++){
            if(r<cumulativeTable[kLogY]){
                y = (double )kLogY;
                break;
            }
	}

	// value in a bin
	r = rand.GetRandomDouble();
	double producedLogEnergy = -(y+r)*Particle.getDeltaLogEnergy()+logEnergy;

        return producedLogEnergy;

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

	return 0;  // 0 is electron flavor defined in Particle.class.
    
    }

    /** Get the name of interaction */
    public String getInteractionName() {

	String decayName = "MuDecay ";
	return decayName;

    }

    /** get Type Of Interaction (Interaction->0; Decay->1) */
    public int getTypeOfInteraction() {
	
	return 1;

    }

}

