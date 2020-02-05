package iceCube.uhe.interactions;

import iceCube.uhe.particles.*;
import iceCube.uhe.points.*;
import iceCube.uhe.interactions.*;
import numRecipes.*;

import java.util.*;
import java.io.*;

/** 
    Construct the interaction Likelihood with The InteractionsBase class.
    The interaction is separated by a stocastic interaction part and a continuous interaction part.
    The boundary of the two parts are defined by the threshold energy set with
    the method setThresholdEnergyAsStochastic(logEthreshold).

*/

public class InteractionsBaseLikelihood extends InteractionsBase implements Function{

    static double logEnergyProducedMinimum = initialLogEnergyProducedMinimum;
    static double expandedDim = Particle.getDimensionOfLogEnergyMatrix() +
	(int )((Particle.getLogEnergyMinimum()-logEnergyProducedMinimum)/Particle.getDeltaLogEnergy());
    private final static double ln10 = Math.log(10.0);
    private static double epsilon = 1.0e-7;
    private InteractionsMatrix intMtx = null;
    private ParticlePoint point = null;
    private double massNumber = 0;
    double perNucleonFactor =1.0;
    double logThresholdEnergy = logEnergyProducedMinimum;
    double yThreshold = 1.0e-5;
    double logYThreshold = -5.0;
    boolean crossSectionConstApproximation = false;
    boolean inelasticityBase = false;
    private double reductionFactor[];
    private double beta[];
    /** for interface getFunction( ). */
    double[] parameters;

    /** Consrtuctor */
    public InteractionsBaseLikelihood(InteractionsMatrix intMtx){
	super(intMtx);
	this.intMtx = intMtx;
	this.point = intMtx.interactions.s;
	parameters = new double[1];
        for(int i=0;i<point.NumberOfSpecies[point.getMaterialNumber()];i++){
            massNumber += point.getNumberOfAtoms(i)*point.getAtomicNumber(i);
        }
	if(getPropDoublet()!=0){ // charged lepton
	    perNucleonFactor = point.NA/massNumber;
	}else{ // neutrino
	    perNucleonFactor = point.NA;
	}

	reductionFactor = new double[dim];
	for(int iLogE=0;iLogE< dim;iLogE++) reductionFactor[iLogE] = 1.0;
	beta = new double[dim];
	unsetThresholdEnergyAsStochastic();
    }
    
    /** Consrtuctor */
    public InteractionsBaseLikelihood(InteractionsMatrix intMtx, double threshold, boolean inelasticityBase){
	super(intMtx);
	this.intMtx = intMtx;
	this.point = intMtx.interactions.s;
	parameters = new double[1];
        for(int i=0;i<point.NumberOfSpecies[point.getMaterialNumber()];i++){
            massNumber += point.getNumberOfAtoms(i)*point.getAtomicNumber(i);
        }
	if(getPropDoublet()!=0){ // charged lepton
	    perNucleonFactor = point.NA/massNumber;
	}else{ // neutrino
	    perNucleonFactor = point.NA;
	}

	reductionFactor = new double[dim];
	for(int iLogE=0;iLogE< dim;iLogE++) reductionFactor[iLogE] = 1.0;
	beta = new double[dim];
	this.inelasticityBase = inelasticityBase;
	setThresholdEnergyAsStochastic(threshold);
    }

    /** Consrtuctor */
    public InteractionsBaseLikelihood(InteractionsMatrix intMtx, double threshold){
	this(intMtx,threshold,false);
    }

    /** setting the threshold energy to define a stochastic interaction*/
    public void setThresholdEnergyAsStochastic(double threshold){
	
	System.err.println("...now calculating the inelasticity coefficient to handle the continuous energy loss");
	if(!inelasticityBase) logThresholdEnergy = threshold;
	else{
	    logYThreshold = threshold;yThreshold = Math.pow(10.0,threshold);
	}

	for(int iLogE=0;iLogE< dim;iLogE++){
	    if(iLogE%100==0) System.err.print(".");
	    double logPrimaryEnergy = Particle.getLogEnergyMinimum()+Particle.getDeltaLogEnergy()*(double )iLogE;
	    double logEthreshold = logThresholdEnergy;
	    if(inelasticityBase) logEthreshold = logPrimaryEnergy+logYThreshold;
	    if(logPrimaryEnergy>=logEthreshold){
		reductionFactor[iLogE] = 1.0 - getCumulativeProbability(logPrimaryEnergy+epsilon, logEthreshold+epsilon);
	    }else{
		reductionFactor[iLogE] = 0.0;
	    }
	    intMtx.interactions.setIncidentParticleEnergy(iLogE);
	    double yMax = Math.pow(10.0,logEthreshold-logPrimaryEnergy);
	    if(inelasticityBase) yMax = yThreshold;
	    if(yMax>= intMtx.interactions.getYmax()) yMax = intMtx.interactions.getYmax();
	    if(yMax > intMtx.interactions.getYmin()+intMtx.interactions.roundOffError){
		beta[iLogE] =
		    intMtx.interactions.getYDSigmaDy(intMtx.interactions.getYmin()+intMtx.interactions.roundOffError,
						     yMax-intMtx.interactions.roundOffError)*perNucleonFactor;
		//System.err.format("logPrimaryEnergy(%6.3f) reductionFac(%e) beta(%e)\n",
		//		      logPrimaryEnergy,reductionFactor[iLogE],beta[iLogE]);
	    }else{
		beta[iLogE] = 0.0;
	    }
	}
	System.err.println(".");
	System.err.format("...interactions to create cascades below %e [GeV] are calculated as continuous processes\n",
			  Math.pow(10.0, logThresholdEnergy));
    }
    
   
    /** unsetting the threshold energy to define a stochastic interaction*/
    public void unsetThresholdEnergyAsStochastic(){
	inelasticityBase = false;
	setThresholdEnergyAsStochastic(logEnergyProducedMinimum);
    }

    public void enableInelasticityBase() {inelasticityBase = true;}
    public void disableInelasticityBase() {inelasticityBase = false;}

   

    /** probability of producing produced_energy from a stocastic interaction by 
        a primary particle with energy incoming_energy **/
    public double getStochasticInteractionProbability(double logPrimaryEnergy, double logProducedEnergy){

	if(logProducedEnergy> logPrimaryEnergy) return 0.0;
	
	if(!inelasticityBase){
	    if(logProducedEnergy<= logThresholdEnergy) return 0.0;
	}else{
	    if(logProducedEnergy <= logPrimaryEnergy+logYThreshold) return 0.0;
	}
	
	int jLogE = (int )((logProducedEnergy-logEnergyProducedMinimum)/
			   Particle.getDeltaLogEnergy( ));
	double logProducedEnergyRightEdge = logProducedEnergy;
	double logProducedEnergyLeftEdge = logProducedEnergy - Particle.getDeltaLogEnergy( );

	double prob = getCumulativeProbability(logPrimaryEnergy+epsilon,logProducedEnergyRightEdge+epsilon)-
	    getCumulativeProbability(logPrimaryEnergy+epsilon,logProducedEnergyLeftEdge+epsilon);
	//double prob = getCumulativeProbability(logPrimaryEnergy,logProducedEnergyRightEdge);
	
	return prob*getTotalSigma(logPrimaryEnergy)*perNucleonFactor;
    }

    /** calculate total cross section of the stochastic interaction [cm^2]*/
    double getTotalSigma(double logPrimaryEnergy){
	int iLogE = (int )((logPrimaryEnergy+epsilon-Particle.getLogEnergyMinimum())/Particle.getDeltaLogEnergy( ));
	if(iLogE<0) iLogE =0;
	if(iLogE>=dim) iLogE = dim-1;
	double thresholdFactor = reductionFactor[iLogE];
	double totalsigma = thresholdFactor*intMtx.getSigmaMatrix(iLogE); // get total cross section
	return totalsigma;
    }


    /** Interaction Likelihood 
	<pre>
	double logPrimaryEnergy   log10(energy of propagating lepton [GeV])
	double logProducedEnergy  log10(energy of produced particles by the interaction [GeV])
	double pathLength         length of the trajectory between one stochastic interaction to another [g/cm2]
	</pre>
     */
    public double getInteractionLikelihood(double logPrimaryEnergy, double logProducedEnergy,
    					  double pathLength){
	double probability_not_to_interact =
	    getNoStocasticLossProbability(logPrimaryEnergy,pathLength);
	double logEnergyAfterPathLength = getLogPrimaryParticleEnergyAfterCELosses(logPrimaryEnergy, pathLength);
	double probability_to_interact_here =
	    getStochasticInteractionProbability(logEnergyAfterPathLength,logProducedEnergy);
	double likelihood = probability_not_to_interact*probability_to_interact_here;
	return likelihood;
    }

    
    /** return the beta, the inelsticity term of dE/dX [g/cm2]^-1 from a continuous energy loss */
    double getInelasticityTerm(double logPrimaryEnergy){
	int iLogE = (int )((logPrimaryEnergy+epsilon-Particle.getLogEnergyMinimum())/Particle.getDeltaLogEnergy( ));
	if(iLogE<0) iLogE =0;
	if(iLogE>=dim) iLogE = dim-1;
	return beta[iLogE];
    }

    /** Energy of primary particle after running pathLength [g/cm2]. CEL approximation is used.
	For sinmplicity, beta is represented by the initial logPrimaryEnergy.
     */
    public double getLogPrimaryParticleEnergyAfterCELosses(double logPrimaryEnergy, double pathLength){
	double logEnergy = logPrimaryEnergy - getInelasticityTerm(logPrimaryEnergy)*pathLength/ln10;
	return logEnergy;
    }


    /** Return the object of Particle initiating this interaction.
	This method is useful when you need to check what kind of
	particle is involved in this interaction. For example,
	this method returns muon object if the instance of this class
	is muon pair production.
     */
    public Particle getParticleInitiatingThisInteraction(){
	Particle p = intMtx.interactions.p;
	return p;
    }


    /** calculate probaility of propagatiing over pathLngth [g/cm2] without the stochastic interactions.
	It is given by exp(-\int dX Na x Sigma). 
     */
    public double getNoStocasticLossProbability(double logPrimaryEnergy, double pathLength){
	double prob = 0.0;
	if(crossSectionConstApproximation){
	    prob = Math.exp(-pathLength*perNucleonFactor*getTotalSigma(logPrimaryEnergy));
	}else{
	    int functionIndex = 1;
	    parameters[0] =  logPrimaryEnergy;
	    
	    InteractionsBaseLikelihood interactionsLikelihood = this;
	    double xIntegral = Integration.RombergIntegral(interactionsLikelihood,
			      functionIndex, parameters, 0.0, pathLength);
	    prob = Math.exp(-perNucleonFactor*xIntegral);
	}
	return prob;
    }


    /** employ the approximation that the stochastic interaction cross section is energy independent
	during the propagation over pathLngth. This is for calculation of 
	the method getNoStocasticLossProbability(double logPrimaryEnergy, double pathLength)
    */
    public void employAppoximationOfConstantCrossSection(){
	crossSectionConstApproximation = true;
    }


    /** do not employ the approximation that the stochastic interaction cross section is energy independent
	during the propagation over pathLngth. This is for calculation of 
	the method getNoStocasticLossProbability(double logPrimaryEnergy, double pathLength).
	This is a default.
    */
    public void doNotEmployAppoximationOfConstantCrossSection(){
	crossSectionConstApproximation = false;
    }
    
    

    /** 
	<pre>
	Method for interface <Function>. 

	functionIndex     1     sigma(E*exp(-beta x))
	functionIndex     2     reserved
	</pre>
    */
    public double getFunction(int functionIndex, double[] parameters, 
	double x){ 
	double dF;
	double logPrimaryEnergy = parameters[0];

        switch (functionIndex) {
        case 1 : 
	    dF = getTotalSigma(getLogPrimaryParticleEnergyAfterCELosses(logPrimaryEnergy, x));
	    break;
        default:
            dF = 0.0;
            System.err.println("Illegal parameters! Index" + functionIndex);
            System.exit(0);
        }

        return dF;

    }


    /** Main function (for debugging) */
    public static void main(String[] args) throws IOException {

        String fileName = null;
	double logEthreshold = InteractionsBaseLikelihood.logEnergyProducedMinimum;
        if(args.length<2){
            System.out.println("Usage: InteractionsBaseLikelihood file-name log(Ethreshold for stochastic interactions)");
	    System.exit(0);
        }else{
            fileName = args[0];
	    logEthreshold = Double.valueOf(args[1]).doubleValue();
        }
	System.err.format(" threshold energy for this stochastic intereactions %e [GeV]\n",Math.pow(10.0,logEthreshold));

		// Read the serialized object of the Interaction Matrix
        FileInputStream in = new FileInputStream(fileName);
	InteractionsMatrix intMtx = 
	    InteractionsMatrixInput.inputInteractionsMatrix(in);
	in.close( );

	InteractionsBaseLikelihood interactionLikelihood = new InteractionsBaseLikelihood(intMtx,logEthreshold,false);
	//interactionLikelihood.employAppoximationOfConstantCrossSection();
	
	System.out.println("zone 2 1");
	System.out.println("titx Energy [GeV]");
	System.out.println("tity probability");
	System.out.println("gwin 0.2 0.9 0.2 0.9");
	//System.out.println("scal 100.0 1.0e12 0.0 0.005");
	System.out.println("scal 1000.0 1.0e12 1.0e-8 1.0e-5");
	int iLogE,jLogE;
	for(iLogE=0;iLogE< Particle.getDimensionOfLogEnergyMatrix();iLogE+=300){
	    //for(iLogE=0;iLogE< 100;iLogE+=10){
	    double logPrimaryEnergy = Particle.getLogEnergyMinimum() + Particle.getDeltaLogEnergy()*(double )iLogE;
	    for(jLogE=0;jLogE< expandedDim;jLogE++){
		double logProducedEnergy = logEnergyProducedMinimum + 
		    Particle.getDeltaLogEnergy()*(double )jLogE;
		double producedEnergy = Math.pow(10.0,logProducedEnergy);

		double sigma = interactionLikelihood.getTotalSigma(logPrimaryEnergy);
		//double prob =  interactionLikelihood.getStochasticInteractionProbability(logPrimaryEnergy, logProducedEnergy);
		double prob =  interactionLikelihood.getInteractionLikelihood(logPrimaryEnergy, logProducedEnergy,100.0);
		System.out.format("data %e 0.0 %e 0.0\n",producedEnergy,prob);

	    }
	    System.out.println("logx");
	    System.out.println("logy");
	    System.out.println("join");
	    System.out.println("disp");
	    System.out.println("cont");
	}
	System.out.println("endg");
	
	System.out.println("titx Energy [GeV]");
	System.out.println("tity probability");
	System.out.println("gwin 0.2 0.9 0.2 0.9");
	System.out.println("scal 1000.0 1.0e12 1.0e-9 1.0e-6");
	int ix;
	double logPrimaryEnergy = 8.0; // 10^8 GeV
	for(ix=0;ix<100000;ix += 10000){
	    double pathLength = (double )ix;
	    //double sigma = interactionLikelihood.getTotalSigma(logPrimaryEnergy);
	    //double mfp = pathLength*interactionLikelihood.perNucleonFactor*interactionLikelihood.getTotalSigma(logPrimaryEnergy);
	    //System.err.format(" %d %e %e %e\n",ix,interactionLikelihood.perNucleonFactor,sigma,mfp);
	    for(jLogE=0;jLogE< expandedDim;jLogE++){
		double logProducedEnergy = logEnergyProducedMinimum + 
		    Particle.getDeltaLogEnergy()*(double )jLogE;
		double producedEnergy = Math.pow(10.0,logProducedEnergy);
		double prob =  interactionLikelihood.getInteractionLikelihood(logPrimaryEnergy, logProducedEnergy,pathLength);
		System.out.format("data %e 0.0 %e 0.0\n",producedEnergy,prob);

	    }
	    System.out.println("logx");
	    System.out.println("logy");
	    System.out.println("join");
	    System.out.println("disp");
	    System.out.println("cont");
	}

	System.out.println("endg");
    
    }

}
