package numRecipes;

import numRecipes.*;

import java.io.*;
import java.util.*;

/** 
    Perform test based upon the ratio of the maximum likelihood ratio
    built from the binned Poisson.
    The null hypothesis (=background) and the alternative (=signal+background)
    hypothesis are describied, respectively, by
    the class PoissonBinnedLikelihoodCalculator.
    The class PoissonBinnedLikelihoodRatioFactory takes care of the relevant calculations
    <pre>

    The null hypothesis        :  bg only
    The alternative hypothesis :  bg + lambda * signal
                                (lambda>=0)

    </pre>

    Written originally by S. Yoshida for the IceCube EHE analysis.
    2013/6/24
*/
public class PoissonBinnedLikelihoodRatioTest {

    /** list to store likelihood ratio from sets of many replica experiments*/
    private List llhRatioList = null;
    private boolean listHasBeenGenerated = false;

    /**  PoissonBinnedLikelihoodRatioFactory object */
    protected PoissonBinnedLikelihoodRatioFactory llhFactory = null;

    /** The default constructor. Sets the likelihhood calculator of the null and alternative
	hypotheses.
     */
    public PoissonBinnedLikelihoodRatioTest(PoissonBinnedLikelihoodCalculator calNull,
					       PoissonBinnedLikelihoodCalculator calAlter){
	llhFactory = new PoissonBinnedLikelihoodRatioFactory(calNull, calAlter);
	llhFactory.calculateMaxLogLikelihoodRatio(1.0);

    }


    /** 
	return the maximized log likelihood ratio log(L/L0).
	if runReplicaBGexperiment sets true, run the replica MC experiment
	with the H0 (=background only) hypothesis, and calculate the ratio
	aganst the results of the replica experiment.
	This method may be used for evaluating p-value for a given 
	likelihood ratio (which can be provided by this method also setting 
	runReplicaBGexperiment=false).
	<pre> 
 	double scaleFactor:  scaleFactor which is multiplied to mu_bg_i; if set to 1.0, no fudge factor is introduced.
	boolean runReplicaBGexperiment:    true - run Replica experiment with null hypotheis
	</pre>
	This method calls the same function of the PoissonBinnedLikelihoodRatioFactory class.
     */
   public double getMaxLogLikelihoodRatio(double scaleFactor,boolean runReplicaBGexperiment){
	return llhFactory.getMaxLogLikelihoodRatio(scaleFactor, runReplicaBGexperiment);
    }

    /**
       make a list collection of the log-likelihood ratio obtained by running the replica experiment
       by runTimes under the null (background-only) hypothesis.
     */
    public void makeCollectionOfLogLikelihoodRatio(double scaleFactor, int runTimes){
	llhRatioList = new LinkedList();

	// run the replica experiment with the null (background-only) hypothesis by runTimes
	for(int run = 0 ;run<runTimes; run++){
	    double  logLikelihoodRatioReplica = getMaxLogLikelihoodRatio(scaleFactor,true);

	    Double llhRatioObj = new Double(logLikelihoodRatioReplica);
	    llhRatioList.add(llhRatioObj);
	}

	// sort
	Collections.sort(llhRatioList);
    }

    /**
       make a list collection of the log-likelihood ratio obtained by running the replica experiment
       by runTimes under the algernative (signal) hypothesis.
     */
    public void makeCollectionOfLogLikelihoodRatioAlternative(double scaleFactor, int runTimes){
	llhRatioList = new LinkedList();

	// run the replica experiment with the null (background-only) hypothesis by runTimes
	for(int run = 0 ;run<runTimes; run++){
	    llhFactory.calH1.runReplicaExperiment();
	    llhFactory.calH1.useTheResultsByTheReplicaExperiment();
	    llhFactory.calH0.useTheResultsByTheRealExperiment();
	    llhFactory.calH0.copyObservedNumbers(llhFactory.calH1);
	    double  logLikelihoodRatioReplica = getMaxLogLikelihoodRatio(scaleFactor,false);

	    Double llhRatioObj = new Double(logLikelihoodRatioReplica);
	    llhRatioList.add(llhRatioObj);
	}

	// sort
	Collections.sort(llhRatioList);
    }

    public ListIterator getListIterator(){
	ListIterator llhRatioListIterator = llhRatioList.listIterator();
	return llhRatioListIterator;
    }

    /**
       calculate the p-value for a given log-likelihood ratio (llh) value.
       Estimated by the collection of the llh obtained by calling makeCollectionOfLogLikelihoodRatio(int runTimes),
       i.e., running many replica experiment under the null hypothesis.
    */

    public double getPvalue(double scaleFactor, double precision, double observed_llhratio){

	double pValue = 0.0;
	if(precision>0.0){
	    int trialFactor = (int )Math.ceil(1.0/precision);
	    makeCollectionOfLogLikelihoodRatio(scaleFactor, trialFactor);
	}

	ListIterator llhRatioListIterator = llhRatioList.listIterator();
	int times = 0;
	while(llhRatioListIterator.hasNext()){
	    Double llhRatioObj = (Double )(llhRatioListIterator.next());
	    double llhRatio = llhRatioObj.doubleValue();
	    //System.err.format(" llhratioReplica=%e\n",llhRatio);
	    if(llhRatio>observed_llhratio) break;
	    times++;
	}
	pValue = 1.0-((double)times)/((double)llhRatioList.size());
	return pValue;
    }


    /** A simple main method */
    public static void main(String[] args) throws IOException{
 

	boolean debugMode = false;
	String sigEventRateFileName = null;
	String bgEventRateFileName = null;

	if(args.length<2){
	    System.out.println("Usage: PoissonBinnedLikelihoodRatioTest filename-to-read-BGdata filename-to-read-SIGdata (debug mode)");
	    System.exit(0);
	}else { 
	    bgEventRateFileName = args[0];
	    sigEventRateFileName = args[1];
	    if(args.length == 3) debugMode=true;
	    else if(args.length > 3){
		System.out.println("Illeagal arguments");
		System.exit(0);
	    }
	}


	// BG : background-only hypothesis
	PoissonBinnedLikelihoodCalculator calBG = new PoissonBinnedLikelihoodCalculator();
	if(debugMode) calBG.debugFlag = true;
	DataInputStream in = new DataInputStream(ClassLoader.getSystemResourceAsStream(bgEventRateFileName));
	calBG.fillData(in);
	in.close();

	// SIG : signal hypothesis
	PoissonBinnedLikelihoodCalculator calSIG = new PoissonBinnedLikelihoodCalculator();
	if(debugMode) calSIG.debugFlag = true;
	in = new DataInputStream(ClassLoader.getSystemResourceAsStream(sigEventRateFileName));
	calSIG.fillData(in);


	// Display the bin data
	System.out.println("============ Background binned data =============");
	calBG.printLogLikelihood(true);
	System.out.println("============   Signal binned data   =============");
	calSIG.printLogLikelihood(true);

	PoissonBinnedLikelihoodRatioTest tester = new PoissonBinnedLikelihoodRatioTest(calBG,calSIG);
	double llhRatioObserved = tester.getMaxLogLikelihoodRatio(1.0,false);
	System.out.println("log(LikelihoodRatio) = " + llhRatioObserved);

	double precision = 1.0e-4;
	double pValue = tester.getPvalue(1.0,precision,llhRatioObserved);
	System.out.format(" p-value=%e\n",pValue);

    }


}