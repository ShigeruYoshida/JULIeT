package iceCube.uhe.analysis;

import numRecipes.*;

import java.io.*;
import java.util.*;

/** 
    A collection of the methods for the GZK (or astro) neutrino model tests
    based upon the likelihood ratio of the poisson binned likelihood.

    This class calls the methods provided by the class PoissonBinnedLikelihoodCalculator
    in the numRecipes package.

    The main data members are:
    <pre>
    PoissonBinnedLikelihoodCalculator  calBG             (for the atmospheric background)
    PoissonBinnedLikelihoodCalculator  calSinal          (for the neutrino model such as GZK to be tested)
    PoissonBinnedLikelihoodCalculator  calNuisanceSingal (for the nuisance signal - like E^-2 against the GZK test)

    PossionBinnedLikelihoodCalculator  calNull           Null Hypothesis
    PossionBinnedLikelihoodCalculator  calAlter          Alternative Hypothesis
    boolean  includeNuisance           true for including the nuisance signal contribution
    </pre>


*/
public class ModelTestByPoissonBinnedLikelihoodFactory {

    public PoissonBinnedLikelihoodCalculator calBG = null;    // Likelihood Calculator for the atmospheric backgrounds
    public PoissonBinnedLikelihoodCalculator calSignal = null;// Likelihood Calculator for the neutrino model to be tested
    public PoissonBinnedLikelihoodCalculator calNuisanceSignal = null;// Likelihood Calculator for nuisance signal

    public PoissonBinnedLikelihoodCalculator calNull = null;  // Likelihood Calculator for the null hypothesis
    public PoissonBinnedLikelihoodCalculator calAlter = null; // Likelihood Calculator for the alternative hypothesis
    public PoissonBinnedLikelihoodCalculator calHybrid = null; // Likelihood Calculator for the hybrid hypothesis

    protected boolean includeNuisance = true;
    protected double llhNull = Double.POSITIVE_INFINITY;
    protected double llhAlter = Double.POSITIVE_INFINITY;

    protected double maximizedSignalFactor = 0.0;
    protected double maximizedNuisanceFactor = 0.0;

    private boolean calNullHasBeenGenerated = false;
    private boolean builtNullHypothesisLikelihood = false;
    private boolean builtAlterHypothesisLikelihood = false;

    /** list to store ratio of likelihood -log(product of binned Poisson) from sets of many replica experiments*/
    private List llhRatioList = null;
    private List maximizedSigFacList = null;
    private boolean listHasBeenGenerated = false;


    /** The default constructor. Sets the likelihhood calculator of the BG, the neutrino signal, the nuisance signal
	hypotheses.
     */
    public ModelTestByPoissonBinnedLikelihoodFactory(PoissonBinnedLikelihoodCalculator calBG,
						     PoissonBinnedLikelihoodCalculator calSignal,
						     PoissonBinnedLikelihoodCalculator calNuisanceSignal){
	this.calBG = calBG; 
	this.calSignal = calSignal; 
	this.calNuisanceSignal = calNuisanceSignal;

	includeNuisance = true;
    }

    /** The default constructor. Sets the likelihhood calculator of the BG, and the neutrino signal.
	No nuisance signal is considered.
     */
    public ModelTestByPoissonBinnedLikelihoodFactory(PoissonBinnedLikelihoodCalculator calBG,
						     PoissonBinnedLikelihoodCalculator calSignal){
	this.calBG = calBG; 
	this.calSignal = calSignal; 

	includeNuisance = false;
    }

    /**
       Set the PoissonBinnedLikelihoodCalculator that stores the binned poisson data 
       for the nuisance signal
     */
    public void setNuisanceSignal(PoissonBinnedLikelihoodCalculator calNuisanceSignal){

	this.calNuisanceSignal = calNuisanceSignal;

	includeNuisance = true;
    }


    /**
       Set the PoissonBinnedLikelihoodCalculator that stores the binned poisson data 
       for the signal model to be tested. This method is called when you try a different
       neutrino model for testing.
     */
    public void setAnotherNeutrinoSignalModel(PoissonBinnedLikelihoodCalculator calAnotherSignal){

	this.calSignal = calAnotherSignal;
	calNullHasBeenGenerated = false;  // you have to build the null hypothesis (bg + sig) again
	builtNullHypothesisLikelihood = false; // you have to build the null hypothesis (bg + sig) again
    }

    /** Do not let the likelihood calculation include the nuisance signal */
    public void doNotIncludeNuisanceSignal(){

	includeNuisance = false;
    }


    public boolean isTheNuisanceSignalIncluded(){
	return includeNuisance;
    }


    /** 
	Build the null hypothesis : the (GZK or astro) signal + background.

	if isNuisanceSignalIncluded() is true, the PoissonBinnedLikelihoodCalculator calNuisanceSignal
        is added to maximize the lilelihood of the null hypothesis - the "profile" likelihood.

	if backgroundOnly = true, the null hypothesis contains only the background poisson data.

	if runReplicaExperiment is true, a replica experiment under the null hypothesis (the signal to tested + background)
        is performed and then calculate the likelihood.

	Return the log-likelihood -log(Products of the binned Poisson p.d.f) of the null hypothesis

     */
    public double buildLikelihoodForNullHypothesis(boolean backgroundONLY, boolean runReplicaExperiment){


	if(calNullHasBeenGenerated == false){
	    calNull = new PoissonBinnedLikelihoodCalculator();
	    calNull.useTheResultsByTheRealExperiment();
	    calBG.useTheResultsByTheRealExperiment();
	    //System.out.format("expected sum of bg hypothesis=%e\n",calBG.getSumOfExpectedValues());
	    calSignal.useTheResultsByTheRealExperiment();
	    //System.out.format("expected sum of signal hypothesis=%e\n",calSignal.getSumOfExpectedValues());
	    calNull.copyPoissonBinnedData(calBG);   // copy the Poisson binned data (obs#, expected#) stored in the BG likelihhod calculator
	    if(!backgroundONLY){
		calNull.addExpectedNumbers(calSignal);  // add the expected event rate of the signal to be tested.
		if(includeNuisance){
		    calNuisanceSignal.useTheResultsByTheRealExperiment();
		    maximizedNuisanceFactor = 0.0;
		    maximizedNuisanceFactor = calNull.getSignalFactor(1.0, calNuisanceSignal);
		    //System.err.format(" nuisance factor in the null hypothesis =%e\n",maximizedNuisanceFactor);
		    calNull.addExpectedNumbers(maximizedNuisanceFactor,calNuisanceSignal);
		}
		//System.out.format("expected sum of null hypothesis=%e\n",calNull.getSumOfExpectedValues());
		//System.out.format("observed sum of null hypothesis=%d\n",calNull.getSumOfObservedValues());
	    }

	    calNullHasBeenGenerated = true;
	}

	if(runReplicaExperiment){
	    llhNull = calNull.runReplicaExperiment();
	    calNull.useTheResultsByTheReplicaExperiment();
	    //System.out.format("observed sum of null hypothesis=%d\n",calNull.getSumOfObservedValues());
	}else{
	    calNull.useTheResultsByTheRealExperiment();
	}

	llhNull = calNull.getLikelihood();

	calNull.useTheResultsByTheRealExperiment();

	builtNullHypothesisLikelihood = true;
	return llhNull;

    }


    /** 
	Build the null hypothesis : the (GZK or astro) signal + background.

	if isNuisanceSignalIncluded() is true, the PoissonBinnedLikelihoodCalculator calNuisanceSignal
        is added to maximize the lilelihood of the null hypothesis - the "profile" likelihood.

	if runReplicaExperiment is true, a replica experiment under the null hypothesis (the signal to tested + background)
        is performed and then calculate the likelihood.

	Return the log-likelihood -log(Products of the binned Poisson p.d.f) of the null hypothesis

     */
    public double buildLikelihoodForNullHypothesis(boolean runReplicaExperiment){

	boolean backgroundONLY = false;
	return buildLikelihoodForNullHypothesis(backgroundONLY,runReplicaExperiment);
    }

    /** 
	Build the alternaive hypothesis :

	saturatedLikelihood == true:  Build the saturated likelihood from the null hypothesis calNull.
	saturatedLikelihood == false: Build the maximized likelihood with bg + model signal floating the model normalization.
	saturatedLikelihood == false
	&& nuisanceAlternative == true: 
	Build the maximized likelihood with bg + nuisance model signal floating the model normalization.


	Return the log-likelihood -log(Products of the binned Poisson p.d.f).
    */
    public double buildLikelihoodForAlternativeHypothesis(boolean saturatedLikelihood, 
							  boolean nuisanceAlternative,
							  boolean useTheResultsOfReplicaExperiment){

	if(!builtNullHypothesisLikelihood){
	    System.err.println(" First build the likelihood for null hypothesis!");
	    System.err.println(" call buildLikelihoodForNullHypothesis(boolean runReplicaExperiment)");
	    return 0.0;
	}

	calAlter = new  PoissonBinnedLikelihoodCalculator();


	//
	// calculated a saturated Poisson binned likelihood
	//
	if(saturatedLikelihood){ 

	    calAlter.useTheResultsByTheRealExperiment();
	    calAlter.copyPoissonBinnedData(calNull);    
	    // copy the Poisson binned data (obs#, expected#) stored in the null hypothesis likelihood calculator
	    if(useTheResultsOfReplicaExperiment){
		calAlter.useTheResultsByTheReplicaExperiment();
		calAlter.copyPoissonBinnedData(calNull);    
		// copy the Poisson binned data (obs#, expected#) from the replica experiment
	    }


	    calAlter.replaceExpectedNumbersWithObservedNumbers(calBG); // build a saturated likelihood
	    llhAlter = calAlter.getLikelihood();

	    calAlter.useTheResultsByTheRealExperiment();


	//
	// calculated a likelihood of bg + model signal with the model normaliztion maximizing the likelihood
	//
	}else{

	    calAlter.useTheResultsByTheRealExperiment();
	    calAlter.copyPoissonBinnedData(calNull); 
	    calAlter.copyExpectedNumbers(calBG);
	    // the expected number in each poisson bin are replaced with the BG expectation.
	    if(useTheResultsOfReplicaExperiment){
		calAlter.useTheResultsByTheReplicaExperiment();
		calAlter.copyPoissonBinnedData(calNull);    
		calAlter.copyExpectedNumbers(calBG);
		// copy the Poisson binned data (obs#, expected#) from the replica experiment
	    }

	    calSignal.useTheResultsByTheRealExperiment();
	    if(!nuisanceAlternative){
		maximizedSignalFactor = calAlter.getSignalFactor(1.0, calSignal);
		calAlter.addExpectedNumbers(maximizedSignalFactor,calSignal);
		//System.err.format(" signal floated %e\n",maximizedSignalFactor);
	    }else{
		calNuisanceSignal.useTheResultsByTheRealExperiment();
		maximizedSignalFactor = calAlter.getSignalFactor(1.0, calNuisanceSignal);
		calAlter.addExpectedNumbers(maximizedSignalFactor,calNuisanceSignal);
		//System.err.format(" nuisance floated %e\n",maximizedSignalFactor);
	    }
	    //System.out.format("expected sum of alternative hypothesis=%e\n",calAlter.getSumOfExpectedValues());
	    llhAlter = calAlter.getLikelihood();

	    calAlter.useTheResultsByTheRealExperiment();

	}

	builtAlterHypothesisLikelihood = true;
	return llhAlter;

    }

    /** 
	Build the alternaive hypothesis :

	saturatedLikelihood == true:  Build the saturated likelihood from the null hypothesis calNull.
	saturatedLikelihood == false: Build the maximized likelihood with bg + model signal floating the model normalization.
	saturatedLikelihood == false

	Return the log-likelihood -log(Products of the binned Poisson p.d.f).
    */
    public double buildLikelihoodForAlternativeHypothesis(boolean saturatedLikelihood, 
							  boolean useTheResultsOfReplicaExperiment){
	return buildLikelihoodForAlternativeHypothesis(saturatedLikelihood, false,
						       useTheResultsOfReplicaExperiment);
    }


    /**
       Build the "hybrid" hypothesis - signal with floating normalization + nuisance with floating nomalization
       to maximize the likelihood
     */
    public double buildLikelihoodForHybridHypothesis(boolean useTheResultsOfReplicaExperiment){

	if(!builtNullHypothesisLikelihood){
	    System.err.println(" First build the likelihood for null hypothesis!");
	    System.err.println(" call buildLikelihoodForNullHypothesis(boolean runReplicaExperiment)");
	    return 0.0;
	}

	if(calNuisanceSignal == null){
	    System.err.println("nuisance model is not registered!");
	    return -1.0;
	}

	calHybrid = new  PoissonBinnedLikelihoodCalculator();

	calHybrid.useTheResultsByTheRealExperiment();
	calHybrid.copyPoissonBinnedData(calNull); 
	// the expected number in each poisson bin are replaced with the BG expectation.
	if(useTheResultsOfReplicaExperiment){
	    calHybrid.useTheResultsByTheReplicaExperiment();
	    calHybrid.copyPoissonBinnedData(calNull);    
	    // copy the Poisson binned data (obs#, expected#) from the replica experiment
	}

	calSignal.useTheResultsByTheRealExperiment();
	double totalSIG = calSignal.getSumOfExpectedValues();
	long totalObs = calNull.getSumOfObservedValues();
	if(useTheResultsOfReplicaExperiment){
	    calNull.useTheResultsByTheReplicaExperiment();
	    totalObs = calNull.getSumOfObservedValues();
	    calNull.useTheResultsByTheRealExperiment();
	}
	double signalFactor = (double )totalObs/totalSIG; // the initial guess
	double signalFactorSearchRangeMin = 0.0;
	double signalFactorSearchRangeMax = 2.0*signalFactor;
	if(totalObs<=0) signalFactorSearchRangeMax = 0.1*totalSIG;
	int nSteps = 100; // 1% steps
	double deltaFactor = (signalFactorSearchRangeMax-signalFactorSearchRangeMin)/(double)nSteps;
	//System.err.format(" total obs(%d) signal range max=%f\n",totalObs,signalFactorSearchRangeMax);

	calNuisanceSignal.useTheResultsByTheRealExperiment();
	double totalNuisanceSIG = calNuisanceSignal.getSumOfExpectedValues();
	double nuisanceSignalFactor = (double )totalObs/totalNuisanceSIG; // the initial guess
	double nuisanceSignalFactorSearchRangeMin = 0.0;
	double nuisanceSignalFactorSearchRangeMax = 2.0*nuisanceSignalFactor;
	if(totalObs<=0) nuisanceSignalFactorSearchRangeMax = 0.1*totalNuisanceSIG;
	double deltaFactorNuisance = (nuisanceSignalFactorSearchRangeMax-nuisanceSignalFactorSearchRangeMin)/(double)nSteps;
	//System.err.format(" total obs(%d) nuisance range max=%f\n",totalObs,nuisanceSignalFactorSearchRangeMax);

	double llhMin = Double.POSITIVE_INFINITY;
	double signalFactorMin = 0.0;
	double nuisanceSignalFactorMin = 0.0;
	for(int i=0; i<=nSteps; i++){ // signal loop
	    double sigFactorSearch = signalFactorSearchRangeMin +deltaFactor*(double )i;
	    for(int j=0;j<=nSteps;j++){ // nuisnace signal loop
		double nuisanceSigFactorSearch = nuisanceSignalFactorSearchRangeMin +deltaFactorNuisance*(double )j;
		calHybrid.copyExpectedNumbers(calBG);
		calHybrid.addExpectedNumbers(sigFactorSearch,calSignal);
		calHybrid.addExpectedNumbers(nuisanceSigFactorSearch,calNuisanceSignal);

		double llhHybrid = calHybrid.getLikelihood();
		if(llhHybrid < llhMin){
		    llhMin = llhHybrid;
		    signalFactorMin = sigFactorSearch;
		    nuisanceSignalFactorMin = nuisanceSigFactorSearch;
		    //System.err.format("sigFac(%e) nuisFac(%e) llh=%e\n",
		    //	      sigFactorSearch,nuisanceSigFactorSearch,llhMin);
		}
	    }
	}

	calHybrid.copyExpectedNumbers(calBG);
	calHybrid.addExpectedNumbers(signalFactorMin,calSignal);
	calHybrid.addExpectedNumbers(nuisanceSignalFactorMin,calNuisanceSignal);
	//System.err.format(" signal factor=%e nuisance factor=%e\n",signalFactorMin,nuisanceSignalFactorMin);
	//System.err.format("expected sum of null hypothesis=%e\n",calNull.getSumOfExpectedValues());
	//System.err.format("observed sum of null hypothesis=%d\n",calNull.getSumOfObservedValues());
	//System.err.format("expected sum of hybrid hypothesis=%e\n",calHybrid.getSumOfExpectedValues());
	//System.err.format("observed sum of hybrid hypothesis=%d\n",calHybrid.getSumOfObservedValues());

	double llhHybrid = calHybrid.getLikelihood();

	maximizedSignalFactor = signalFactorMin;
	maximizedNuisanceFactor = nuisanceSignalFactorMin;
	builtAlterHypothesisLikelihood = true;
	calHybrid.useTheResultsByTheRealExperiment();


	return llhHybrid;
    }



    /**
       Return the normalization of the model to maximize the binned poisson pdf.
       Calculated by buildLikelihoodForAlternativeHypothesis(boolean saturatedLikelihood=false).
     */
    public double getModelNormalizationToMaximizeLikelihood(){
	if(builtAlterHypothesisLikelihood) return maximizedSignalFactor;
	else{
	    System.err.println("call buildLikelihoodForAlternativeHypothesis");
	    return -1.0;
	}
    }

    /**
       Return the nuisance normalization of the model to maximize the binned poisson pdf.
       Calculated by buildLikelihoodForHybridHypothesis(boolean useTheResultsOfReplicaExperiment)
     */
    public double getNuisanceNormalizationToMaximizeLikelihood(){
	if(builtAlterHypothesisLikelihood) return maximizedNuisanceFactor;
	else{
	    System.err.println("call buildLikelihoodForAlternativeHypothesis");
	    return -1.0;
	}
    }


    /**
       make a list collection of the log-likelihood -log(binned Poisson pdf) obtained by running the replica experiment
       by runTimes under the null hypothesis of bg + neutrino model (+ nuisance signal)

       The likeliihood ratio type is following:
       <pre>
       1   saturated likelihood/null hypothesis(bg + model signal (+ nuisance))
       2   alternative hypothesis(bg + model signal with floated normalization)/null hypothesis(bg + model signal (+ nuisance))
       3   saturated likelihood/alternative hypothesis (bg + model signal with floated normalization)
       4   (bg + nuisance model signal with floated normalization)/(bg + model signal with floated normalization)
       5   alternative hypothesis(bg + either signal or nuisance model with floated normalization, 
           select the one giving larger likelihood)/null hypothesis(bg ONLY)
       6   (bg + nuisance model signal with floated normalization + signal with floated norm)/
           (bg + model signal with floated normalization)
       </pre>
     */
    public void makeCollectionOfLogLikelihoodRatio(int type, int runTimes){
	llhRatioList = new LinkedList();
	if(type>=2 && type!=6) maximizedSigFacList = new LinkedList();
	boolean runReplicaExperiment = true;
	boolean useTheResultsOfReplicaExperiment = true;
	// run the replica experiment with the null (background + signal) hypothesis by runTimes

	boolean saturatedLikelihood = true;
	for(int run = 0 ;run<runTimes; run++){
	    double llhH0 =  Double.POSITIVE_INFINITY;
	    if(type!=5){
		llhH0 = buildLikelihoodForNullHypothesis(runReplicaExperiment); // bg + signal
	    }else{
		boolean backgroundONLY = true;
		llhH0 = buildLikelihoodForNullHypothesis(backgroundONLY,runReplicaExperiment);// bg only
	    }
	    double llhSaturated = buildLikelihoodForAlternativeHypothesis(saturatedLikelihood,useTheResultsOfReplicaExperiment);

	    double llhSignalFloated =  Double.POSITIVE_INFINITY;
	    double llhNuisanceSignalFloated =  Double.POSITIVE_INFINITY;
	    double llhHybridFloated =  Double.POSITIVE_INFINITY;
	    if(type>=2 && type!=6){
		llhSignalFloated = buildLikelihoodForAlternativeHypothesis(false,useTheResultsOfReplicaExperiment);
		double maximizedFactor = getModelNormalizationToMaximizeLikelihood();
		Double maximizedFactorObj = new Double(maximizedFactor);
		maximizedSigFacList.add(maximizedFactorObj);
		if(type>3 && type!=6){
		    llhNuisanceSignalFloated = buildLikelihoodForAlternativeHypothesis(false,true,
										       useTheResultsOfReplicaExperiment);
		}
	    }
	    if(type==6) llhHybridFloated= buildLikelihoodForHybridHypothesis(useTheResultsOfReplicaExperiment);

	    double llhH1 =  Double.POSITIVE_INFINITY;
	    switch (type){

	    case 1:
		llhH1 = llhSaturated;
		break;

	    case 2:
		llhH1 = llhSignalFloated;
		break;

	    case 3:
		llhH0 = llhSignalFloated;
		llhH1 = llhSaturated;
		break;

	    case 4:
		llhH0 = llhSignalFloated;
		llhH1 = llhNuisanceSignalFloated;
		break;

	    case 5:
		llhH1 = llhSignalFloated;
		if(llhNuisanceSignalFloated<llhSignalFloated) llhH1 = llhNuisanceSignalFloated;
		break;

	    case 6:
		llhH1 = llhHybridFloated;
		break;

	    default:
		System.err.format(" type %d you put is wrong!\n",type);
	    }


	    double llhRatio = llhH0 - llhH1;
	    Double llhRatioObj = new Double(llhRatio);
	    llhRatioList.add(llhRatioObj);
	    //System.err.format(" --runTimes(%d) llhRatio=%e\n",run,llhRatio);
	}

	// sort
	Collections.sort(llhRatioList);
	if(type>=2 && type!=6) Collections.sort(maximizedSigFacList);
    }


    public ListIterator getllhRatioIterator(){
	return llhRatioList.listIterator();
    }

    public static void outputLikelihhoRatioList(ModelTestByPoissonBinnedLikelihoodFactory testFac,
						OutputStream out)  throws IOException {

        ObjectOutputStream objectOut = new ObjectOutputStream(out);

        objectOut.writeObject(testFac.llhRatioList);
        objectOut.flush();
    }

    public static void outputMaximizedSignalFactorList(ModelTestByPoissonBinnedLikelihoodFactory testFac,
						OutputStream out)  throws IOException {

        ObjectOutputStream objectOut = new ObjectOutputStream(out);

        objectOut.writeObject(testFac.maximizedSigFacList);
        objectOut.flush();
    }

    public static LinkedList inputList(InputStream in)  throws IOException {

        LinkedList llhList = null;
        try{

            ObjectInputStream objectIn = new ObjectInputStream(in);
            llhList = (LinkedList )objectIn.readObject();

        }catch(ClassNotFoundException e){
            System.err.println("Caught ClassNotFoundException: " +
                               e.getMessage( ));
            System.exit(0);
	}catch (EOFException e){
            llhList = null;
            return llhList;
	}

	return llhList;
    }

}

