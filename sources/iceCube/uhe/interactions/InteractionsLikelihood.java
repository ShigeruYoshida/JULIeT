package iceCube.uhe.interactions;

import iceCube.uhe.particles.*;
import iceCube.uhe.points.*;
import iceCube.uhe.interactions.*;
import numRecipes.*;

import java.util.*;
import java.io.*;

/** 
    Construct the interaction Likelihood using the methods in InteractionsBaseLikelihood class.
    While InteractionsLikelihood class handles a single interaction channel provided by
    the interactionMatrix object, this class handes MULTIPLE channels by a series of
    InteractionBaseLikelihood objects provided in the constructor.

    The interaction is separated by a stocastic interaction part and a continuous interaction part.
    The boundary of the two parts are defined by the threshold energy set with
    the method setThresholdEnergyAsStochastic(logEthreshold).

*/

public class InteractionsLikelihood implements Function{

    private final static double ln10 = Math.log(10.0);
    private static double epsilon = 1.0e-7;
    private List interactionsList = null;
    public LikelihoodData likelihoodData = null;
    int flavorID = 1; // particle flavor  ID of the propagating particle
    int doubletID = 1;// particle doublet ID of the propagating particle
    private boolean whetherCheckParticleID = false;
    double logThresholdEnergy = InteractionsBaseLikelihood.logEnergyProducedMinimum;
    static double yThreshold = 1.0e-5;
    boolean inelasticityBase = false;
    double perNucleonFactor = 1.0;
    boolean baseLikelihoodHasBeenAdded = false;
    boolean crossSectionConstApproximation = false;
    String interactionsMatrixDirectory  = "iceCube/uhe/interactions/ice/";
    private static String eNuCCMtxObjectZeusFile = "ENeutrinoChargeZeusNewMtx";
    private static String muNuCCMtxObjectZeusFile = "MuNeutrinoChargeZeusNewMtx";
    private static String tauNuCCMtxObjectZeusFile = "TauNeutrinoChargeZeusNewMtx";
    private static String eNuNCMtxObjectZeusFile = "ENeutrinoNeutralZeusNewMtx";
    private static String muNuNCMtxObjectZeusFile = "MuNeutrinoNeutralZeusNewMtx";
    private static String tauNuNCMtxObjectZeusFile = "TauNeutrinoNeutralZeusNewMtx";
    /** for interface getFunction( ). */
    double[] parameters;



    /** Consrtuctor */
    public InteractionsLikelihood(double threshold, boolean inelasticityBase){
	interactionsList = new LinkedList();
	parameters = new double[1];
	if(!inelasticityBase) logThresholdEnergy = threshold;
	else yThreshold = threshold;
	this.inelasticityBase = inelasticityBase;
	likelihoodData = new LikelihoodData();
    }
    
    public InteractionsLikelihood(double logEthreshold){
	this(logEthreshold,false);
    }

    /** Add the interactionBaseLikelihood object */
    public void addInteractionsBaseLikelihoodObject(InteractionsBaseLikelihood baseLikelihood){
	interactionsList.add(baseLikelihood);
	perNucleonFactor = baseLikelihood.perNucleonFactor;
	baseLikelihoodHasBeenAdded = true;
    }

    /** Utility method to add interactionBaseLikelihood objects.
	Reading the interaction matrix, generate InteractionsBaseLikelihood, and add it to the list
	if, for example, doMu2e = 1 (i.e., pair creation from a muon).
     */
    public void configureInteractionsForLikelihoodCalculation(
                                 int doCC, int doNC, int doMuBrems, int doTauBrems, 
                                 int doMuKnock, int doTauKnock, int doMu2e, int doTau2e,
                                 int doMu2mu, int doTau2mu, int doMu2tau, int doTau2tau,
                                 int doMuPN, int doTauPN) throws IOException{


        // Register Interactions and read the InteractionMatrix objects
	/** For Glashow Resonance 16->20*/
        String[] fileName = new String[20];
        String matrixName = null;
	int n = 0;

        // Charged Current Interactions [yes(1)/no(0)/allFlavor(2)]
        if(doCC == 1){
	    if(flavorID==0) {matrixName = eNuCCMtxObjectZeusFile;}
	    else if(flavorID==1) matrixName = muNuCCMtxObjectZeusFile;
	    else if(flavorID==2) matrixName = tauNuCCMtxObjectZeusFile;
	    fileName[n++]= interactionsMatrixDirectory.concat(matrixName);
        }else if(doCC == 2){
	    matrixName = eNuCCMtxObjectZeusFile;
	    fileName[n++]= interactionsMatrixDirectory.concat(matrixName);
	    matrixName = muNuCCMtxObjectZeusFile;
	    fileName[n++]= interactionsMatrixDirectory.concat(matrixName);
	    matrixName = tauNuCCMtxObjectZeusFile;
	    fileName[n++]= interactionsMatrixDirectory.concat(matrixName);
        }

        // Neutral Current Interactions [yes(1)/no(0)/allFlavor(2)]
        if(doNC == 1){
	    if(flavorID==0) matrixName      = eNuNCMtxObjectZeusFile;
	    else if(flavorID==1) matrixName = muNuNCMtxObjectZeusFile;
	    else if(flavorID==2) matrixName = tauNuNCMtxObjectZeusFile;
            fileName[n++]= interactionsMatrixDirectory.concat(matrixName);
        }else if(doNC == 2){
	    matrixName = eNuNCMtxObjectZeusFile;
	    fileName[n++]= interactionsMatrixDirectory.concat(matrixName);
	    matrixName = muNuNCMtxObjectZeusFile;
	    fileName[n++]= interactionsMatrixDirectory.concat(matrixName);
	    matrixName = tauNuNCMtxObjectZeusFile;
	    fileName[n++]= interactionsMatrixDirectory.concat(matrixName);
        }

        // Muon Bremsstrahlung [yes(1)/no(0)]
        if(doMuBrems == 1){
            matrixName   = "muBremsstrahlungMtx";
            fileName[n++]= interactionsMatrixDirectory.concat(matrixName);
        }

        // Tau Bremsstrahlung [yes(1)/no(0)]
        if(doTauBrems == 1){
            matrixName   = "tauBremsstrahlungMtx";
            fileName[n++]= interactionsMatrixDirectory.concat(matrixName);
        }

        // Muon Knock-on Electrons [yes(1)/no(0)]
        if(doMuKnock == 1){
            matrixName   = "muKnockOnElectronsMtx";
            fileName[n++]= interactionsMatrixDirectory.concat(matrixName);
        }

        // Tau Knock-on Electrons [yes(1)/no(0)]
        if(doTauKnock == 1){
            matrixName   = "tauKnockOnElectronsMtx";
            fileName[n++]= interactionsMatrixDirectory.concat(matrixName);
        }

        // Muon e+e- Pair Creation [yes(1)/no(0)]
        if(doMu2e == 1){
            matrixName   = "muToePairCreationFitMtx";
            fileName[n++]= interactionsMatrixDirectory.concat(matrixName);
        }

        // Tau e+e- Pair Creation [yes(1)/no(0)]
        if(doTau2e == 1){
            matrixName   = "tauToePairCreationFitMtx";
            fileName[n++]= interactionsMatrixDirectory.concat(matrixName);
        }

        // Muon mu+mu- Pair Creation [yes(1)/no(0)]
        if(doMu2mu == 1){
            matrixName   = "muTomuPairCreationFitMtx";
            fileName[n++]= interactionsMatrixDirectory.concat(matrixName);
        }

        // Tau mu+mu- Pair Creation [yes(1)/no(0)]
        if(doTau2mu == 1){
            matrixName   = "tauTomuPairCreationFitMtx";
            fileName[n++]= interactionsMatrixDirectory.concat(matrixName);
        }

        // Muon tau+tau- Pair Creation [yes(1)/no(0)]
        if(doMu2tau == 1){
            matrixName   = "muTotauPairCreationFitMtx";
            fileName[n++]= interactionsMatrixDirectory.concat(matrixName);
        }

        // Tau tau+tau- Pair Creation [yes(1)/no(0)]
        if(doTau2tau == 1){
            matrixName   = "tauTotauPairCreationFitMtx";
            fileName[n++]= interactionsMatrixDirectory.concat(matrixName);
        }

        // Muon Photo-nuclear interactions [yes(1)/no(0)]
        if(doMuPN == 1){
            matrixName = "muPhotoNuclearMtx";
            fileName[n++]= interactionsMatrixDirectory.concat(matrixName);
        }

        // Tau Photo-nuclear interactions [yes(1)/no(0)]
        if(doTauPN == 1){
            matrixName   = "tauPhotoNuclearMtx";
            fileName[n++]= interactionsMatrixDirectory.concat(matrixName);
        }

        // regist interactionsMtx
        System.out.println("Registered Interaction Matrix....");

	for(int i=0; i<n ;i++){       
            System.out.println(fileName[i]);
            InputStream in = ClassLoader.getSystemResourceAsStream(fileName[i]);
            InteractionsMatrix intMtx  = InteractionsMatrixInput.inputInteractionsMatrix(in);
            in.close();
	    InteractionsBaseLikelihood intBaseLikelihood = null;
	    if(!inelasticityBase) intBaseLikelihood = new InteractionsBaseLikelihood(intMtx,logThresholdEnergy,inelasticityBase);
	    else intBaseLikelihood = new InteractionsBaseLikelihood(intMtx,yThreshold,inelasticityBase);
	    addInteractionsBaseLikelihoodObject(intBaseLikelihood);
        }
    }

    /** setting the threshold energy to define a stochastic interaction*/
    public void setThresholdEnergyAsStochastic(double threshold){
	if(!inelasticityBase) logThresholdEnergy = threshold;
	else yThreshold = threshold;
	if(!baseLikelihoodHasBeenAdded){
	    printErrorMessage();
	}
	ListIterator interactionsIterator = interactionsList.listIterator();
	InteractionsBaseLikelihood intBaseLikelihood = null;
	while(interactionsIterator.hasNext()){
	    intBaseLikelihood = (InteractionsBaseLikelihood )interactionsIterator.next();
	    if(inelasticityBase) intBaseLikelihood.enableInelasticityBase();
	    else intBaseLikelihood.disableInelasticityBase();
	    intBaseLikelihood.setThresholdEnergyAsStochastic(threshold);
	}
    }
    
    /** setting the threshold energy to define a stochastic interaction*/
    public void setThresholdEnergyAsStochastic(){
	if(!baseLikelihoodHasBeenAdded){
	    printErrorMessage();
	}
	ListIterator interactionsIterator = interactionsList.listIterator();
	InteractionsBaseLikelihood intBaseLikelihood = null;
	while(interactionsIterator.hasNext()){
	    intBaseLikelihood = (InteractionsBaseLikelihood )interactionsIterator.next();
	    if(inelasticityBase) intBaseLikelihood.enableInelasticityBase();
	    else intBaseLikelihood.disableInelasticityBase();
	    if(!inelasticityBase) intBaseLikelihood.setThresholdEnergyAsStochastic(logThresholdEnergy);
	    else intBaseLikelihood.setThresholdEnergyAsStochastic(yThreshold);
	}
    }
    
   
    /** unsetting the threshold energy to define a stochastic interaction*/
    public void unsetThresholdEnergyAsStochastic(){
	setThresholdEnergyAsStochastic(InteractionsBaseLikelihood.logEnergyProducedMinimum);
	logThresholdEnergy = InteractionsBaseLikelihood.logEnergyProducedMinimum;
    }

    public void enableInelasticityBase() {inelasticityBase = true;}
    public void disableInelasticityBase() {inelasticityBase = false;}

    /** check if the propagating particle defined by flavorID and doubletID is subject to this 
	interaction */
    private boolean whetherThisInteractionActsOnThisParticle(InteractionsBaseLikelihood intBaseLikelihood){
	if(!whetherCheckParticleID) return true; // always true as we do not check particle-IDs.

	Particle particleOfThisInteraction = intBaseLikelihood.getParticleInitiatingThisInteraction();
	int flavorID_thisInteraction  = particleOfThisInteraction.getFlavor();
	int doubletID_thisInteraction  = particleOfThisInteraction.getDoublet();
	if((flavorID_thisInteraction  == flavorID) && (doubletID_thisInteraction  == doubletID)){
	    return true;
	}else{
	    return false;
	}
    }

    
    /** set the particle IDs of propagating particle which is subject to build the likelihood.
	When you set the IDs via this method, this class always check if a given interactionBaseLikeihood
	is subject to the propagating particle. For example, if propagating particle is set as muon,
	it checkes if each of the interaction is relevant to muons or not when calculating toral cross sections
	and inelasticity. If you do not set the particle IDs, all the interactions you add are used for
	the lilihood calculatuon no matter what particle is actually running.
    */
    public void setPropagatingParticleIDs(int flavor, int doublet){
	flavorID = flavor;
	doubletID = doublet;
	whetherCheckParticleID = true;
    }



    /** probability of producing produced_energy from stocastic interactions by 
        a primary particle with energy incoming_energy **/
    public double getStochasticInteractionProbability(double logPrimaryEnergy, double logProducedEnergy){
	if(!baseLikelihoodHasBeenAdded){
	    printErrorMessage();
	}
	double prob=0.0;
	ListIterator interactionsIterator = interactionsList.listIterator();
	InteractionsBaseLikelihood intBaseLikelihood = null;
	while(interactionsIterator.hasNext()){
	    intBaseLikelihood = (InteractionsBaseLikelihood )interactionsIterator.next();
	    if(whetherThisInteractionActsOnThisParticle(intBaseLikelihood)){
		prob += intBaseLikelihood.getStochasticInteractionProbability(logPrimaryEnergy, logProducedEnergy);
	    }
	}

	return prob;
    }

    /** calculate total cross section of the stochastic interactions [cm^2]*/
    public double getTotalSigma(double logPrimaryEnergy){
	if(!baseLikelihoodHasBeenAdded){
	    printErrorMessage();
	}
	double totalSigma=0.0;
	ListIterator interactionsIterator = interactionsList.listIterator();
	InteractionsBaseLikelihood intBaseLikelihood = null;
	while(interactionsIterator.hasNext()){
	    intBaseLikelihood = (InteractionsBaseLikelihood )interactionsIterator.next();
	    if(whetherThisInteractionActsOnThisParticle(intBaseLikelihood)){
		totalSigma += intBaseLikelihood.getTotalSigma(logPrimaryEnergy);
	    }
	}

	return totalSigma;
    }


    /** Interaction Likelihood 
	<pre>
	double logPrimaryEnergy   log10(energy of propagating lepton [GeV])
	double logProducedEnergy  log10(energy of produced particles by the interaction [GeV])
	double pathLength         length of the trajectory between one stochastic interaction to another [g/cm2]
	</pre>
     */
    public LikelihoodData getInteractionLikelihood(double logPrimaryEnergy, double logProducedEnergy,
    					  double pathLength){
	double probability_not_to_interact =
	    getNoStocasticLossProbability(logPrimaryEnergy,pathLength);
	double logEnergyAfterPathLength = getLogPrimaryParticleEnergyAfterCELosses(logPrimaryEnergy, pathLength);
	double probability_to_interact_here =
	    getStochasticInteractionProbability(logEnergyAfterPathLength,logProducedEnergy);
	double likelihood = probability_not_to_interact*probability_to_interact_here;
	
	double incomingEnergy = Math.pow(10.0,logEnergyAfterPathLength);
	double energyAfterThisCollision = incomingEnergy-Math.pow(10.0,logProducedEnergy);
	if(energyAfterThisCollision<0.0) energyAfterThisCollision = Interactions.roundOffError*incomingEnergy;

	// store the results in the data class.
	likelihoodData.interactionLikelihood = likelihood;
	likelihoodData.particleEnergyAfterThisInteraction = energyAfterThisCollision;
	return likelihoodData;
    }

    
    /** return the beta, the inelsticity term of dE/dX [g/cm2]^-1 from a continuous energy loss */
    double getInelasticityTerm(double logPrimaryEnergy){
	if(!baseLikelihoodHasBeenAdded){
	    printErrorMessage();
	}
	double beta=0.0;
	ListIterator interactionsIterator = interactionsList.listIterator();
	InteractionsBaseLikelihood intBaseLikelihood = null;
	while(interactionsIterator.hasNext()){
	    intBaseLikelihood = (InteractionsBaseLikelihood )interactionsIterator.next();
	    beta += intBaseLikelihood.getInelasticityTerm(logPrimaryEnergy);
	}

	return beta;
    }

    /** Energy of primary particle after running pathLength [g/cm2]. CEL approximation is used.
	For sinmplicity, beta is represented by the initial logPrimaryEnergy.
     */
    public double getLogPrimaryParticleEnergyAfterCELosses(double logPrimaryEnergy, double pathLength){
	double logEnergy = logPrimaryEnergy - getInelasticityTerm(logPrimaryEnergy)*pathLength/ln10;
	return logEnergy;
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
	    
	    InteractionsLikelihood interactionsLikelihood = this;
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


    /** print error message and exit */
    private void printErrorMessage(){
	System.err.println("No InteractionsBaseLikelihood objects have been read");
	System.err.println("call  addInteractionsBaseLikelihoodObject(InteractionsBaseLikelihood baseLikelihood) method");
	System.exit(0);
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
        if(args.length<1){
            System.out.println("Usage: InteractionsLikelihood log(Ethreshold for stochastic interactions)");
	    System.exit(0);
        }else{
	    logEthreshold = Double.valueOf(args[0]).doubleValue();
        }
	System.err.format(" threshold energy for this stochastic intereactions %e [GeV]\n",Math.pow(10.0,logEthreshold));

	InteractionsLikelihood interactionLikelihood =  new InteractionsLikelihood(logEthreshold);

	DataInputStream input = new DataInputStream(System.in);
	BufferedReader  d     = new BufferedReader(new InputStreamReader(input));
	String buffer;
		// Read the serialized object of the Interaction Matrix
	int endReadingMatrixFiles = 0;
	do{
	    System.err.print("matrix file name-> ");
	    buffer   = d.readLine();
	    fileName = buffer;
	
	    FileInputStream in = new FileInputStream(fileName);
	    InteractionsMatrix intMtx = 
		InteractionsMatrixInput.inputInteractionsMatrix(in);
	    in.close( );

	    InteractionsBaseLikelihood intBaseLikelihood = new InteractionsBaseLikelihood(intMtx,logEthreshold);
	    interactionLikelihood.addInteractionsBaseLikelihoodObject(intBaseLikelihood);

	    System.err.print("add more interaction? yes (1) or no (0)-> ");
	    buffer   = d.readLine();
	    endReadingMatrixFiles = Integer.valueOf(buffer).intValue();
	}while(endReadingMatrixFiles == 1);
	    
	System.err.print("set flavor and doublet of propagating particle ? yes (1) or no (0)-> ");
	buffer   = d.readLine();
	int setParticle = Integer.valueOf(buffer).intValue();
	if(setParticle==1){
	    System.err.print("flavor -> ");
	    buffer   = d.readLine();
	    int flavorID = Integer.valueOf(buffer).intValue();
	    System.err.print("doublet -> ");
	    buffer   = d.readLine();
	    int doubletID = Integer.valueOf(buffer).intValue();
	    interactionLikelihood.setPropagatingParticleIDs(flavorID, doubletID);
	}
	    

	System.out.println("zone 2 1");
	System.out.println("titx Energy [GeV]");
	System.out.println("tity probability");
	System.out.println("gwin 0.2 0.9 0.2 0.9");
	System.out.println("scal 1000.0 1.0e12 1.0e-10 1.0e-4");
	int iLogE,jLogE;
	for(iLogE=0;iLogE< Particle.getDimensionOfLogEnergyMatrix();iLogE+=100){
	    double logPrimaryEnergy = Particle.getLogEnergyMinimum() + Particle.getDeltaLogEnergy()*(double )iLogE;
	    for(jLogE=0;jLogE< InteractionsBaseLikelihood.expandedDim;jLogE++){
		double logProducedEnergy = InteractionsBaseLikelihood.logEnergyProducedMinimum + 
		    Particle.getDeltaLogEnergy()*(double )jLogE;
		double producedEnergy = Math.pow(10.0,logProducedEnergy);

		LikelihoodData likelihoodData 
		    =  interactionLikelihood.getInteractionLikelihood(logPrimaryEnergy, logProducedEnergy,100.0);
		double prob = likelihoodData.getLikelihoodValue();
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
	System.out.println("scal 1000.0 1.0e12 1.0e-10 1.0e-4");
	int ix;
	double logPrimaryEnergy = 9.0; // 10^8 GeV
	for(ix=0;ix<100000;ix += 10000){
	    double pathLength = (double )ix;
	    for(jLogE=0;jLogE< InteractionsBaseLikelihood.expandedDim;jLogE++){
		double logProducedEnergy = InteractionsBaseLikelihood.logEnergyProducedMinimum + 
		    Particle.getDeltaLogEnergy()*(double )jLogE;
		double producedEnergy = Math.pow(10.0,logProducedEnergy);
		LikelihoodData likelihoodData 
		    =  interactionLikelihood.getInteractionLikelihood(logPrimaryEnergy, logProducedEnergy,pathLength);
		double prob = likelihoodData.getLikelihoodValue();
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
