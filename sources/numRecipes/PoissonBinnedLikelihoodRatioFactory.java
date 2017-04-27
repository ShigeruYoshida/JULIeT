package numRecipes;

import numRecipes.*;

import java.io.*;
import java.util.*;

/** 
    Handle the ratio of the maximum likelihood ratio
    built from the binned Poisson.
    The null hypothesis (=background) and the alternative (=signal+background)
    hypothesis are describied, respectively, by
    the class PoissonBinnedLikelihoodCalculator.
    <pre>

    The null hypothesis        :  bg only
    The alternative hypothesis :  bg + lambda * signal
                                (lambda>=0)

    </pre>

    Written originally by S. Yoshida for the IceCube EHE analysis.
    2013/6/24
*/
public class PoissonBinnedLikelihoodRatioFactory {

    protected PoissonBinnedLikelihoodCalculator calH0; // Likelihood Calculator for the null hypothesis
    protected PoissonBinnedLikelihoodCalculator calH1; // Likelihood Calculator for the alternative hypothesis

    /** log likelohood for the null hypothesis */
    protected double poissonLogLikelihoodH0 = 0.0; 
    /** log likelohood for the alternative hypothesis */
    protected  double poissonLogLikelihoodH1 = 0.0; 
    /** maximized likelihood ratio log(L/L0) */
    protected double maxLogLikelihoodRatio = 1.0;
    /** the scale factor to maximize likelihood for the alternative hypothesis */
    protected double signalFactorForMaximalLlhH1 = 0.0;

    private boolean calculatedLikelihoodRatio = false;


    /** The default constructor. Sets the likelihhood calculator of the null and alternative
	hypotheses.
     */
    public PoissonBinnedLikelihoodRatioFactory(PoissonBinnedLikelihoodCalculator calNull,
					       PoissonBinnedLikelihoodCalculator calAlter){
	calH0 = calNull; // Likelihood Calculator for the null hypothesis
	calH1 = calAlter; // Likelihood Calculator for the alternative hypothesis
    }

    /** 
	calculate the maximized log likelihood ratio log(L/L0).
    */
    public void calculateMaxLogLikelihoodRatio(double scaleFactor){

	poissonLogLikelihoodH0 = calH0.getLikelihood(scaleFactor);

	calH0.useTheResultsByTheRealExperiment();
	signalFactorForMaximalLlhH1 = calH0.getSignalFactor(scaleFactor,calH1);
	//System.err.format(" signal factor to maximize likelihood (%e)\n",signalFactorForMaximalLlhH1);

	poissonLogLikelihoodH1 = calH0.getLikelihood(scaleFactor,signalFactorForMaximalLlhH1,calH1);

	maxLogLikelihoodRatio = poissonLogLikelihoodH0-poissonLogLikelihoodH1;
	calculatedLikelihoodRatio = true;
    }

    /** return the binned Poisson likelihood for the null (background-only) hypothesis 
       <pre>
       double scaleFactor:  scaleFactor whici is multiplied to mu_bg_i; if set to 1.0, no fudge factor is introduced.
       </pre>
     */
    public double getPoissonLogLikelihoodH0(double scaleFactor){
	calculateMaxLogLikelihoodRatio(scaleFactor);  
	return  poissonLogLikelihoodH0;
    }


    /** return the binned Poisson likelihood for the alternative (background + signal) hypothesis
       <pre>
       double scaleFactor:  scaleFactor whici is multiplied to mu_bg_i; if set to 1.0, no fudge factor is introduced.
       </pre>
    */
    public double getPoissonLogLikelihoodH1(double scaleFactor){
	calculateMaxLogLikelihoodRatio(scaleFactor);  
	return  poissonLogLikelihoodH1;
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
	boolean runReplicaBGexperiment:    true - run Replica experiment with null hypotheis
	double scaleFactor:  scaleFactor which is multiplied to mu_bg_i; if set to 1.0, no fudge factor is introduced.
	</pre>
     */
    public double getMaxLogLikelihoodRatio(double scaleFactor,boolean runReplicaBGexperiment){
	if(!runReplicaBGexperiment){
	    calculateMaxLogLikelihoodRatio(scaleFactor);  
	    return maxLogLikelihoodRatio;

	}else{ // run the replica experiment with the background only hypothesis
	    double poissonLogLikelihoodH0Replica = calH0.runReplicaExperiment();
	    //System.err.format("  -- llh H0 replica=%e\n",poissonLogLikelihoodH0Replica);

	    calH0.useTheResultsByTheReplicaExperiment();
	    double replicaSignalFactorForMaximalLlhH1 = calH0.getSignalFactor(scaleFactor,calH1);
	    //System.err.format(" replica signal factor to maximize likelihood (%e)\n",replicaSignalFactorForMaximalLlhH1);

	    double poissonLogLikelihoodH1Replica = calH0.getLikelihood(scaleFactor,replicaSignalFactorForMaximalLlhH1,calH1);
	    calH0.useTheResultsByTheRealExperiment();
	    //System.err.format("  -- llh H1 replica=%e\n",poissonLogLikelihoodH1Replica);

	    double logLikelihoodRatioReplica = poissonLogLikelihoodH0Replica-poissonLogLikelihoodH1Replica;
	    if(logLikelihoodRatioReplica<0.0){
		calH0.useTheResultsByTheReplicaExperiment();
		calH0.printLogLikelihood(true);
		calH0.useTheResultsByTheRealExperiment();
		System.err.format("  -- llh H0 replica=%e\n",poissonLogLikelihoodH0Replica);
		System.err.format("  -- llh H1 replica=%e\n",poissonLogLikelihoodH1Replica);
	    }
	    return logLikelihoodRatioReplica;
	}
    }


}
