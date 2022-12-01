package iceCube.uhe.interactions;

import iceCube.uhe.particles.*;
import iceCube.uhe.event.*;
import iceCube.uhe.decay.*;
import iceCube.uhe.interactions.*;
import numRecipes.*;

import java.util.*;

/** <pre>
    The InteractionsBase class to treat all of interactions same as decay for
    RunManager class. This class and other "Base" classes inherit the MonteCarloBase class.
    </pre>
*/

public class InteractionsBase extends MonteCarloBase {
    
    /** InteractionsMatrix obejct */
    private InteractionsMatrix interactMtx;

    /** Minimum log energy of propagating particles */
    private double logEnergyMinimum = Particle.getLogEnergyMinimum();

    /** Minimum log energy of produced particles */
    static final double initialLogEnergyProducedMinimum = 1.0;
    private double logEnergyProducedMinimum = initialLogEnergyProducedMinimum; 

    /** In order to save CPU time, we increase neutrino cross section by
        this factor. As long as the meanfree path is by far shorter than
        the propagation length (presumably 1km), this is equivallent
        to the case when one neutrino particle represents multiple neutrinos
        whose number is equal to this factor. You (or RunManager class) have to devide
        this factor later to compensate this "artificial" enhancement
        of the neutrino cross section.
    */
    public static int neutrinoFactor = 1;

    /** dimension of InteractionsMatrix */
    int dim         = Particle.getDimensionOfLogEnergyMatrix();
    int expandedDim = dim + 
                    (int )((logEnergyMinimum-logEnergyProducedMinimum)/Particle.getDeltaLogEnergy());

    /** Cumulative cross section table */
    double[][] cumulativeTable = new double[dim][expandedDim];

    /** Constructor for making the cumulative table. */
    public InteractionsBase(InteractionsMatrix intMtx){
	this.interactMtx = intMtx;
	cumulativeTable  = new double[dim][expandedDim];
	setCumulativeTable(interactMtx);
    }
    
    public static double getLogEnergyProducedMinimum(){
	return initialLogEnergyProducedMinimum;
    }

    /** Make a cumulative table of differential cross section. */
    private void setCumulativeTable(InteractionsMatrix interactMtx){

	//make corrected table
	double     correctedElement = 0.0; 
	double[][] transferSigmaTable;
	transferSigmaTable = new double[dim][expandedDim];

	for(int iLogE=0; iLogE<dim; iLogE++){
	    
	    for(int fLogE=0; fLogE<=iLogE; fLogE++){

	    	int kLogE = fLogE + expandedDim-dim;
	    	transferSigmaTable[iLogE][kLogE] = interactMtx.getTransferMatrix(iLogE,fLogE);
		    
	    }

	    // now adding the lower energy region
	    double delta = Particle.getDeltaLogEnergy();
	    Interactions interactions = interactMtx.interactions;
	    interactions.setIncidentParticleEnergy(iLogE);
	    for(int fLogE=-(expandedDim-dim); fLogE<0; fLogE++){
  	    //for(int fLogE=-(expandedDim-dim); fLogE<iLogE; fLogE++){
		int kLogE = fLogE + expandedDim-dim;
		double logY = delta*(double )(fLogE-iLogE);
		double logYUp = logY+0.5*delta;
		double yUpRange = Math.pow(10.0,logYUp);
		double logYLow = logY-0.5*delta;
		double yLowRange = Math.pow(10.0,logYLow);
		if(yLowRange>=interactions.getYmin()+2.0*interactions.roundOffError){
		    transferSigmaTable[iLogE][kLogE]=
			interactions.integralDSigmaDy(yLowRange,yUpRange);
		}else{
		    transferSigmaTable[iLogE][kLogE]=0.0;
		}
	    }

	}

	//make cumulative Table
	for(int iLogE=0; iLogE<dim; iLogE++){
	    double sum = 0.0;
	    for(int kLogE=0; kLogE<expandedDim; kLogE++){
		sum += transferSigmaTable[iLogE][kLogE];
	    }
	    if(iLogE%100==0){
		double ratio = sum/interactMtx.getSigmaMatrix(iLogE);
		System.err.println("total sigma [" + iLogE + "] = " + 
				   interactMtx.getSigmaMatrix(iLogE) + 
				   " sum = " + sum + " ratio = " + ratio);
	    }
	    // rescale the total cross section for log(produced energy) > logEnergyProducedMinimum
	    interactMtx.sigmaMtx[iLogE] = sum;

	    //Normalization
	    for(int kLogE=0; kLogE<expandedDim; kLogE++){
		transferSigmaTable[iLogE][kLogE]  = transferSigmaTable[iLogE][kLogE]/sum;
	    }

	    for(int kLogE=0; kLogE<expandedDim; kLogE++){
		for(int fLogE=0; fLogE<=kLogE; fLogE++){
		    cumulativeTable[iLogE][kLogE] += transferSigmaTable[iLogE][fLogE];
		}
	    }

	}

    }

    /** Get pathlength by random number. **/
    public double getPathLength(int iLogE, RandomGenerator rand) {
	
	double r = rand.GetRandomDouble();

	double totalsigma = interactMtx.getSigmaMatrix(iLogE); // get total cross section
	double path       = -Math.log(1.0-r)/totalsigma;       // path length
	return path ;

    }

    /** Get pathlength by random number. **/
    public double getPathLength(double logEnergy, RandomGenerator rand) {
	
	int iLogE = (int )((logEnergy-Particle.getLogEnergyMinimum())/Particle.getDeltaLogEnergy( ));
        return getPathLength(iLogE, rand);

    }
    
    /** Get pathlength for neutrino. **/
    public double getNeutrinoPathLength(int iLogE, RandomGenerator rand) {
	
	double r = rand.GetRandomDouble();

	double totalsigma = interactMtx.getSigmaMatrix(iLogE)*neutrinoFactor; //crosssection*neutrinoFactor
	double path       = -Math.log(1.0-r)/totalsigma;       // path length
	return path;

    }

    /** Get pathlength for neutrino. **/
    public double getNeutrinoPathLength(double logEnergy, RandomGenerator rand) {

	int iLogE = (int )((logEnergy-Particle.getLogEnergyMinimum())/Particle.getDeltaLogEnergy( ));
        return getNeutrinoPathLength(iLogE, rand);
    }

    /** Get produced log energy. In order to decide the value of log energy in a bin,
        use a random number **/
    public double getProducedEnergy(int iLogE, RandomGenerator rand){
	    
        double e = 0.0;
	double r = rand.GetRandomDouble();

	for(int j=0; j<expandedDim; j++){

            if(r<cumulativeTable[iLogE][j]){
                e = (double )j;
                break;
            }
	}
	// value in a bin
	r = rand.GetRandomDouble();
	double producedLogEnergy = (e + r)*Particle.getDeltaLogEnergy() + 
                        	    logEnergyProducedMinimum;
	return producedLogEnergy;

    }

    /** Get produced log energy. In order to decide the value of log energy in a bin,
        use a random number **/
    public double getProducedEnergy(double logEnergy, RandomGenerator rand){

	int iLogE = (int )((logEnergy-Particle.getLogEnergyMinimum())/Particle.getDeltaLogEnergy( ));
	if(iLogE<0) iLogE = 0;
        return getProducedEnergy(iLogE, rand);

    }

    /** Get Cumulative Probability on a given log(Primary Energy [GeV]),
        log(Produced Energy [GeV]) 
   <pre>
     double logPrimaryEnergy    : Log (Input Primary Energy [GeV])
     double logPproducedEnergy  : Log (Produced Energy [GeV] via this interaction)
   </pre>
    */
    public double getCumulativeProbability(double logPrimaryEnergy, 
					   double logProducedEnergy){

	int iLogE = 
	    (int )((logPrimaryEnergy-Particle.getLogEnergyMinimum())/
		   Particle.getDeltaLogEnergy( ));
	int jLogE = (int )((logProducedEnergy-logEnergyProducedMinimum)/
			   Particle.getDeltaLogEnergy( ));
	if(iLogE<0) iLogE = 0;
	if(iLogE>=dim) iLogE = dim-1;
	if(jLogE<0) jLogE = 0;
	if(jLogE>=expandedDim) jLogE = expandedDim-1;

	return cumulativeTable[iLogE][jLogE];
    }


    /** Get the flavor of the particle propagateing */
    public int getPropFlavor() {

	int propFlavor = interactMtx.getFlavor();
	return propFlavor;

    }

    /** Get the doublet of the particle propagateing */
    public int getPropDoublet() {

	int propDoublet = interactMtx.getDoublet();
	return propDoublet;

    }

    /** Get the flavor of the produced particle */
    public int getProducedFlavor() {

	int producedFlavor = interactMtx.getProducedFlavor();
	return producedFlavor;
    
    }

    /** Get the name of the interaction */
    public String getInteractionName() {

	String nameOfInteraction = interactMtx.interactions.interactionName();
	return nameOfInteraction;

    }

    /** Get type of the interaction (Interaction->0; Decay->1) */
    public int getTypeOfInteraction() {
	
	return 0;

    }


}

