package iceCube.uhe.decay;

import iceCube.uhe.decay.*;
import iceCube.uhe.event.*;
import iceCube.uhe.particles.*;
import iceCube.uhe.interactions.*;
import iceCube.uhe.propagation.*;
import numRecipes.*;

import java.util.*;

/** <pre>
    The TauDecayBase class to treat tau decay same as interactions for Event class.
    This class and other "Base" classes inherit the MonteCarloBase class.

    This class has TauDecayYMatrix object. It treats energy transfer as dN/dlogY,
    Y = Edecay/Etau.
    </pre>
*/

public class TauDecayBase extends MonteCarloBase {

    private TauDecayYMatrix tauDecayMtx;
    private Particle        p;
    private int             decayMode;  // 0 means tau to electron decay
                                        // 1 means tau to mu decay, and 2 means tau to hadrondecay
    private int             incomingLogE;
    private int             dim; 
    private double          random;
    private double[]        cumulativeTable;
    private double[]        cumTau2LeptonTable, cumTau2HadronTable; // one of them is choosed and copied
                                                                    //  to cumulativeTable[] in setDecayMode().

    /** Constructor for making the table of lifetime and cumulative table. */
    public TauDecayBase(TauDecayYMatrix tauDecayMtx) {

	    this.p           = tauDecayMtx.p;
	if(p.getFlavor()==2 && p.getDoublet()==1){ // Requires tau's
            this.tauDecayMtx = tauDecayMtx;
	    dim              = p.getDimensionOfLogEnergyMatrix();
	    for(int kLogE=0; kLogE<dim; kLogE++){
		tauDecayMtx.setLifeTimeMatrix(kLogE);
	    }
	    this.setCumulativeTable(tauDecayMtx);
	    
        }else{
            System.err.println("This particle " + 
			       p.particleName(p.getFlavor(), p.getDoublet()) + 
			       " is not TAUONs!!");
            System.exit(0);
        }

	
    }

    /** Make a cumulative table of differential cross section.
        The elements are normalized to branching ratio of tau decay. */
    public void setCumulativeTable(TauDecayYMatrix tauDecayMtx){
	
	cumTau2LeptonTable = new double[dim]; 
	cumTau2HadronTable = new double[dim];
	
	// MakeDecay Matrix for y=Ee/Etau or for y=Emu/Etau  BRatio=0.18
	for(int iLogY=0; iLogY<dim; iLogY++){
		tauDecayMtx.setTauDecayMatrix(iLogY);
		cumTau2LeptonTable[iLogY] += (iLogY>0?cumTau2LeptonTable[iLogY-1]:0.0)+
		    tauDecayMtx.getTauToChargedLeptonDecayMatrix(iLogY)
		/(0.5*Decay.BRatioTau2Leptons);
	}
	
	// MakeDecay Matrix for y=Eh/Etau  BRatio=0.64
	for(int iLogY=0; iLogY<dim; iLogY++){
	    cumTau2HadronTable[iLogY] += (iLogY>0?cumTau2HadronTable[iLogY-1]:0.0)+
		tauDecayMtx.getTauToHadronDecayMatrix(iLogY)
		/(Decay.BRatioTau2Pi+Decay.BRatioTau2Rho+Decay.BRatioTau2A1+Decay.BRatioTau2X);
	}

	System.out.println("Sum of cumTable for Tau Decay = " + 
			   (cumTau2LeptonTable[dim-1] + cumTau2HadronTable[dim-1]));
	System.out.println("Sum of Tau To LeptonDecay = " + 
			   cumTau2LeptonTable[dim-1] + " Sum of Tau To HadronDecay = " + cumTau2HadronTable[dim-1]);
    }

    /** Get pathlength by random number. **/
    public double getPathLength(int iLogE, RandomGenerator rand) {
	
	double r = rand.GetRandomDouble();
	double path = -tauDecayMtx.getLifeTimeMatrix(iLogE)*PropagationMatrix.c*Math.log(1.0-r);

	return path;

    }

    /** Get pathlength by random number. **/
    public double getPathLength(double logEnergy, RandomGenerator rand) {
	
	int iLogE = (int )((logEnergy-p.getLogEnergyMinimum())/p.getDeltaLogEnergy( ));
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

	setDecayMode(rand);
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
	double logEnergy        = (double)iLogE*Particle.getDeltaLogEnergy() + Particle.getLogEnergyMinimum();
	double produedLogEnergy = -(y+r)*Particle.getDeltaLogEnergy()+logEnergy;
	                            
	return produedLogEnergy;

    }

    /** Get produced log energy. In order to decide the value of log energy in a bin,
        use a random number **/
    public double getProducedEnergy(double logEnergy, RandomGenerator rand){

	setDecayMode(rand);
        double y = 0.0;            
	double r = rand.GetRandomDouble();
	
	for(int kLogY=0; kLogY<dim; kLogY++){
            if(r<cumulativeTable[kLogY]){
                y = (double)kLogY;
                break;
            }
	}

	// value in a bin
	r = rand.GetRandomDouble();
	double producedLogEnergy = -(y+r)*Particle.getDeltaLogEnergy()+logEnergy;

        return producedLogEnergy;

    }


    /** Choose the decay mode. This method is called when tau decay is chozen. 
        Tau to hadon decay is treated which produced one "hadron" particle,
        it is superposition of each hadrons. */
    public void setDecayMode(RandomGenerator rand){ 

	double   r         = rand.GetRandomDouble();
	double[] cumBRatio = new double[3];
	decayMode = 1000000;

	cumBRatio[0] = 0.5*Decay.BRatioTau2Leptons;                 // tau->e
	cumBRatio[1] = 0.5*Decay.BRatioTau2Leptons + cumBRatio[0];  // tau->mu
	cumBRatio[2] = Decay.BRatioTau2Pi+Decay.BRatioTau2Rho +
	    Decay.BRatioTau2A1  +Decay.BRatioTau2X + cumBRatio[1];  // tau->hadron  

	for(int n=0; n<cumBRatio.length; n++) {

	    if(r<cumBRatio[n]) {
		decayMode = n;
		break;
	    }
	}

	if(decayMode==0 || decayMode==1) cumulativeTable = cumTau2LeptonTable;
	else if(decayMode==2) cumulativeTable = cumTau2HadronTable;
    }


    /** Get the decay mode (0:tau to electron, 1:tau to mu, 2:tau to hadron) */
    public int getDecayMode(){
	return decayMode;
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
	
	int mode =100;
	if     (decayMode==0) mode = 0;
	else if(decayMode==1) mode = 1;
    	else if(decayMode==2) mode = 3;
	return mode; 

    }

    /** Get the name of interacion */
    public String getInteractionName() {

	String decayName = "TauDecay ";
	return decayName;

    }

    /** get Type Of Interaction (Interaction->0; Decay->1) */
    public int getTypeOfInteraction() {
	
	return 1;

    }

}

