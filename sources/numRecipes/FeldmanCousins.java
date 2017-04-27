package numRecipes;

import java.io.*;
import java.util.*;

/**

Implementation of the Feldman-Cousins method to calcylate the upper limit
at a given C.L.

<UL>
  <DT> <cite>Unified Approach to the classical statistical analysis of small signals</cite>
       <TABLE CELLPADDING=2>
       <TR>
         <TD ALIGN="LEFT"><it>Gary J. Feldman and Robert D. Cousins</it></TD>
         <TD ALIGN="LEFT">
           <a href="http://link.aps.org/doi/10.1103/PhysRevD.57.3873">
            Physical Review D <b>57</b> 3873 (1998)</a></TD>
       </TR>
       </TABLE><BR>
</UL>

*/

public class FeldmanCousins {

    /** confidence level of the flux interval */
    protected static double confidenceLevel = 0.90;  //  90% C.L. default
    /** Range of searching the signal interval range decided by the Ferlman-Cousins */
    protected static double signalSearchRangeMin = 0.0;
    /** Range of searching the signal interval range decided by the Ferlman-Cousins */
    protected static double signalSearchRangeMax = 100.0;
    /** step size of searching the signal interval range decided by the Ferlman-Cousins */
    protected static double signalSearchStep = 0.0005;
    protected static double bgSearchStep = 0.00002;

    /** off set to the BG if BG=0 */
    protected static double epsilon = 1.0e-16;

    private static long observationRangeMin = (long )signalSearchRangeMax;

    /** Relative Background uncertainties */
    protected static double nBackgroundRelMin = 0.0; // (nBGMin-nBG)/nBG
    protected static double nBackgroundRelMax = 0.0; // (nBGMax-nBG)/nBG
    protected static boolean treatBGuncertainties = false;
    public static boolean inclusive = true;



    /** Set the range of signals in calculating the signal interval - min < signal < max */
    public static void setRangeOfSignalIntervalCalculation(double min, double max, double stepSize){
	if(min<=max && min>=0.0 && stepSize>0.0){
	    signalSearchRangeMin = min;
	    signalSearchRangeMax = max;
	    signalSearchStep = stepSize;
	}else{
	    System.err.println(" wrong range setting for the interval culculation");
	    System.err.println(" unchanged the range variables");
	}
    }

    public static double getMaxRangeOfSignalIntervalCalculation(){ 
	return signalSearchRangeMax;
    }

    public static double getMinRangeOfSignalIntervalCalculation(){ 
	return signalSearchRangeMin;
    }

    public static double getStepSizeOfSignalIntervalCalculation(){
	return signalSearchStep;
    }

    public static void setConfidenceLevel(double cl){
	if(0.0<=cl && cl<=1.0) confidenceLevel = cl;
	else{
	    System.err.println(" confidence level should be [0 1]");
	}
    }

    public static double getConfidenceLevel(){
	return confidenceLevel;
    }

    /** Set the bacgkround uncertainty.
	If you do not call this method, then the confidence belt calculation
        is performed assuming that background is KNOWN i.e, the exact Felman-Cousins recipe.
	<pre>
	double minusBGerr :  (BGminimum - BGmean)/BGmean
	double plusBGerr :   (BGmaximum - BGmean)/BGmean
	</pre>
    */
    public static void setRelativeBackgroundUncertainty(double minusBGerr, double plusBGerr){
	nBackgroundRelMin = minusBGerr;
	nBackgroundRelMax = plusBGerr;
	treatBGuncertainties = true;
    }


    /** Tells whether the given number of observed events are within the confidence
	belt for the given signal and expected background. Determined by 
	the FelmanCousins method based upon the likelihood ratio.
    */
    protected static boolean isNobservedWithinConfidenceBelt(double nSignal, double nBackground,
						      long nObserved){
	int timesOfUnity = 0;
	TreeMap likelihoodTable = new TreeMap();
	LinkedHashMap probabilityTable = new LinkedHashMap();
	double nBackgroundInput = nBackground;
	if(nBackgroundInput==0.0) nBackgroundInput = epsilon;

	// upper range of observed events +8 sigma
	long nMax = (long)(8.0*Math.sqrt(nBackgroundInput+nSignal)+nBackgroundInput+nSignal);
	if(nMax <= 10) nMax = 10;
	// lower range of observed events -5 sigma
	long nMin = (long)(-5.0*Math.sqrt(nBackgroundInput+nSignal)+nBackgroundInput+nSignal);
	if(nMin<0) nMin = 0;
	//System.out.format("BG(%f) sig(%f) %d < observed < %d for %d observation\n",
	//		   nBackground, nSignal, nMin, nMax, nObserved);
	if(nObserved<nMin|| nMax < nObserved){// a way out!
	    observationRangeMin = nMin;
	    return false; 
	}

	for(long i=nMin;i<=nMax;i++){
	    double bestSignal = getSignalToMaximizeProbability(nBackground,i);
	    double prob;
	    double maxProb;
	    if(!treatBGuncertainties){ // PDF is a poisson
		prob = SpecialFunctions.poisson(nSignal+nBackgroundInput,i);
		maxProb = SpecialFunctions.poisson(bestSignal+nBackgroundInput,i);
	    }else{
		prob = getPDFwithUncertainBG(nSignal,nBackgroundInput,i);
		maxProb = getPDFwithUncertainBG(bestSignal,nBackgroundInput,i);
	    }
	    double ratio = 0.0;
	    if(maxProb>0.0) ratio = maxProb/prob;
	    if(ratio == 1.0){
		timesOfUnity++;
		ratio += 1.0e5*epsilon*(double )timesOfUnity;
	    }

	    Long number = new Long(i);
	    Double r = new Double(ratio);
	    likelihoodTable.put(r,number);

	    Double probability = new Double(prob);
	    probabilityTable.put(number,probability);
	}

	long nRangeMin = ((Long )likelihoodTable.get(likelihoodTable.firstKey())).longValue();
	long nRangeMax = nRangeMin;
	Iterator iter = likelihoodTable.entrySet().iterator();
	double sumProb = 0.0;
	while(iter.hasNext()){
	    Map.Entry entry = (Map.Entry)iter.next();
	    double ratio = ((Double )(entry.getKey())).doubleValue();
	    Long numberObj = (Long)(entry.getValue());
	    Double probObj = (Double )probabilityTable.get(numberObj);

	    double prob = probObj.doubleValue();
	    sumProb += prob;

	    long number = numberObj.longValue();
	    if(nRangeMax < number) nRangeMax = number;
	    if(nRangeMin > number) nRangeMin = number;

	    //System.err.format(" -- signal(%f) BG(%f) nObs(%d) prob(%f) sumProb(%f) ratio(%f)\n",
	    //	      nSignal,nBackground,number,prob,sumProb,ratio);

	    if(sumProb> confidenceLevel) break;
	}

	observationRangeMin = nRangeMin;
	//System.out.format("nSignal=%e RangeMin=%d nObserved=%d nRangeMax=%d\n",nSignal,nRangeMin,nObserved,nRangeMax);
	if((nObserved <= nRangeMax) && (nRangeMin <= nObserved)) return true;
	else return false;

    }

    /** Return the confidence level that the given number of observed events 
	for the given signal and expected background. Determined by 
	the FelmanCousins method based upon the likelihood ratio.
    */
    public static double probabilityOfNobserved(double nSignal, double nBackground,
						   long nObserved){
	int timesOfUnity = 0;
	TreeMap likelihoodTable = new TreeMap();
	LinkedHashMap probabilityTable = new LinkedHashMap();
	double nBackgroundInput = nBackground;
	if(nBackgroundInput==0.0) nBackgroundInput = epsilon;

	// upper range of observed events +8 sigma
	long nMax = (long)(8.0*Math.sqrt(nBackgroundInput+nSignal)+nBackgroundInput+nSignal);
	if(nMax <= 10) nMax = 10;
	// lower range of observed events -5 sigma
	long nMin = (long)(-5.0*Math.sqrt(nBackgroundInput+nSignal)+nBackgroundInput+nSignal);
	if(nMin<0) nMin = 0;
	//System.out.format("BG(%f) sig(%f) %d < observed < %d for %d observation\n",
	//		   nBackground, nSignal, nMin, nMax, nObserved);
	if(nObserved<nMin|| nMax < nObserved){// a way out!
	    return 0.0; 
	}

	for(long i=nMin;i<=nMax;i++){
	    double bestSignal = getSignalToMaximizeProbability(nBackground,i);
	    double prob;
	    double maxProb;
	    if(!treatBGuncertainties){ // PDF is a poisson
		prob = SpecialFunctions.poisson(nSignal+nBackgroundInput,i);
		maxProb = SpecialFunctions.poisson(bestSignal+nBackgroundInput,i);
	    }else{
		prob = getPDFwithUncertainBG(nSignal,nBackgroundInput,i);
		maxProb = getPDFwithUncertainBG(bestSignal,nBackgroundInput,i);
	    }
	    double ratio = 0.0;
	    if(maxProb>0.0) ratio = maxProb/prob;
	    if(ratio == 1.0){
		timesOfUnity++;
		ratio += 1.0e5*epsilon*(double )timesOfUnity;
	    }

	    Long number = new Long(i);
	    Double r = new Double(ratio);
	    likelihoodTable.put(r,number);

	    Double probability = new Double(prob);
	    probabilityTable.put(number,probability);
	}

	long nRangeMin = ((Long )likelihoodTable.get(likelihoodTable.firstKey())).longValue();
	long nRangeMax = nRangeMin;
	Iterator iter = likelihoodTable.entrySet().iterator();
	double sumProb = 0.0;
	while(iter.hasNext()){
	    Map.Entry entry = (Map.Entry)iter.next();
	    double ratio = ((Double )(entry.getKey())).doubleValue();
	    Long numberObj = (Long)(entry.getValue());
	    Double probObj = (Double )probabilityTable.get(numberObj);

	    double prob = probObj.doubleValue();

	    long number = numberObj.longValue();
	    //System.err.format(" -- signal(%f) BG(%f) nObs(%d) prob(%f) ratio(%f)\n",
	    //	      nSignal,nBackground,number,prob,ratio);
	    if(inclusive) sumProb += prob;
	    if(number == nObserved) break;
	    if(!inclusive) sumProb += prob;

	}

	return sumProb;
    }

    protected static double getPDFwithUncertainBG(double nSignal, double nBackground,
						  long nObserved){
	double minBG = (1.0 + nBackgroundRelMin)*nBackground;
	double maxBG = (1.0 + nBackgroundRelMax)*nBackground;
	//System.err.format(" -- signal(%f) minBG(%f) maxBG(%f) nObs(%d)\n",
	//		  nSignal,minBG,maxBG,nObserved);
	double bgWidth = maxBG - minBG;
	double bgProb = 1.0/bgWidth;  // a flat probability [minBG maxBG]
	double bgStepWidth = bgSearchStep*bgWidth;
	double pdf = 0.0;
	double bg = minBG;
	while(bg<maxBG){
	    pdf += SpecialFunctions.poisson(nSignal+bg,nObserved)*bgProb;
	    //System.err.format(" -- signal(%f) BG(%f) nObs(%d) prob(%f)\n",
	    //	      nSignal,bg,nObserved,pdf);
	    bg += bgStepWidth;
	}
	return pdf*bgStepWidth;
    }

    protected static double getSignalToMaximizeProbability(double nBackground,long nObserved){
	if(!treatBGuncertainties){
	    double bestSignal = (double )(nObserved) - nBackground;
	    if(bestSignal<0.0) bestSignal = 0.0;
	    return bestSignal;
	}else{
	    double signal = (double )(nObserved) - nBackground - 0.2;
	    if(signal<0.0) signal = 0.0;
	    double bestSignal = signal;
	    double maxProb = 0.0;
	    boolean foundMaximum = true;
	    while(foundMaximum){
		double prob = getPDFwithUncertainBG(signal,nBackground,nObserved);
		if(maxProb < prob){
		    maxProb = prob;
		    bestSignal = signal;
		    foundMaximum = true;
		}else{
		    foundMaximum = false;
		}
		signal += signalSearchStep;
	    }

	    //System.err.format(" ---- BG(%f) nObs(%d) bestSig(%f) bestSigEstimated(%f)\n",
	    //		      nBackground,nObserved,bestSignal,
	    //		      ((double )(nObserved) - nBackground));
	    return bestSignal;
	}
	    
    }

    /**
       Calculate upper limit of number of signals for a given observed number,
       nObserved, in expectation of nBackground and return the result.
       The calculation is performed by the Feilman-Cousins method.
     */
    public static double getUpperLimit(double nBackground,long nObserved){

	double mu = (double )nObserved-nBackground;
	if(mu<0.0) mu = 0.0;
	if(mu<signalSearchRangeMin) mu =  signalSearchRangeMin;

	while(isNobservedWithinConfidenceBelt(mu,nBackground,nObserved)){
	    //if(treatBGuncertainties){ 
	    //	System.err.format(" ---- BG(%f) nObs(%d) sig(%f)\n",
	    //			  nBackground,nObserved,mu);
	    //}
	    mu += signalSearchStep;
	    if(mu>signalSearchRangeMax) break;
	}

	return mu;
    }

    /**
       Calculate lower limit of number of signals for a given observed number,
       nObserved, in expectation of nBackground and return the result.
       The calculation is performed by the Feilman-Cousins method.
    */
    public static double getLowerLimit(double nBackground,long nObserved){

	double mu =  signalSearchRangeMin;

	while(!isNobservedWithinConfidenceBelt(mu,nBackground,nObserved)){
	    if(observationRangeMin>nObserved) break; // A way out!
	    mu += signalSearchStep;
	    if(mu>signalSearchRangeMax) break;
	}

	return mu;
    }

    /**
       Calculate average upper limit of signals for expected background
       of nBackground. It integrates getUpperLimit(nBackground,nobserved)
       convoluted with Poisson(nObserved), following the recipe in
       the Gary and Catherrin's paper.
    */
    public static double getAverageUpperLimit(double nBackground){
	double mu = 0.0;
	double nBackgroundInput = nBackground;
	if(nBackgroundInput==0.0) nBackgroundInput = epsilon;
	long nMax = (long)(8.0*Math.sqrt(nBackgroundInput)+nBackgroundInput);
	if(nMax <= 10) nMax = 10;
	long nMin = (long)(-5.0*Math.sqrt(nBackgroundInput)+nBackgroundInput);
	if(nMin<0) nMin = 0;
	for(long nObserved=nMin;nObserved<=nMax;nObserved++){
	    double poissonProb=SpecialFunctions.poisson(nBackgroundInput,nObserved);
	    mu += getUpperLimit(nBackground,nObserved)*poissonProb;
	    //System.err.format(" -- BG(%f) nObs(%d) signalSum(%f)\n",
	    //	      nBackground,nObserved,mu);
	}

	return mu;
    }

    /**
       Calculate the lower signal numbers that would be discovered
       by a given confidence level (default: 90%) with significance
       of nSignificance sigma. 
       This evaluation is based on the signal discovery
       potential technique.
    */
    public static double getLeastSignalForDiscovery(double nBackground,
						    double nSignificance){

	double maxSignificance = 20.0;  // 20 sigma for the upper bound
	double probabilityOfDiscovery = 
	    SpecialFunctions.integrateGauss(0.0,1.0,nSignificance,
					    maxSignificance)*2.0; 
	                   // a factor of two comes from the two-sided integral
	//System.err.format("prob(%f sigma)=%e\n",nSignificance,
	//		  probabilityOfDiscovery);

	long nMax = (long)(1.0e2*nBackground);
	if(nMax <= 20) nMax = 20;
	long nCritical = nMax; // must be large enough

	// Calculate number of observed events required for discovery claim.
	double pValue = 0.0;
	double nBackgroundInput = nBackground;
	if(nBackgroundInput==0.0) nBackgroundInput = epsilon;
	while(pValue < probabilityOfDiscovery){
	    pValue += SpecialFunctions.poisson(nBackgroundInput,nCritical);
	    if(nCritical<=0) break;
	    nCritical--;
	}
	if(nCritical == 0) nCritical = 1;

	//System.err.format("You need %d events for discovery\n",nCritical);

	// Calculate the least signals for discovery
	double mu =  (double )nCritical - nBackground;
	double muStepSize = signalSearchStep;
	if(mu>=10.0) muStepSize = 0.001*mu;
	if(mu>=20.0) muStepSize = 0.01*mu;
	if(mu<0.0){
	    mu =0.0; muStepSize = signalSearchStep;
	}
	do{
	    pValue = 0.0; mu += muStepSize;
	    for(long i= nMax;i>=nCritical;i--) 
		pValue += SpecialFunctions.poisson(nBackgroundInput+mu,i);
	    //System.err.format("  p(%e) mu(%f)\n",pValue,mu);
	}while(pValue<confidenceLevel);

	return mu;

    }


    /**
       The simple main method - display the 90% C.L. interval for the Poisson 
       signal mean.
     */
    public static void main(String[] args) throws IOException{

	double nBackground = 0.0;
	long nObserved = 1;

        if(args.length!=2){
            System.out.println("Usage: FeldmanCousins nBackground nObserved");
	    System.exit(0);
        }else{
            nBackground = Double.valueOf(args[0]).doubleValue();
            nObserved = Long.valueOf(args[1]).longValue();
        }

	FeldmanCousins.inclusive = false;
	System.out.format(" Upper limit(%f %d)=%f\n",nBackground,nObserved,
			  FeldmanCousins.getUpperLimit(nBackground,nObserved));
	System.out.format(" Lower limit(%f %d)=%f\n",nBackground,nObserved,
			  FeldmanCousins.getLowerLimit(nBackground,nObserved));
	System.out.format(" Average Upper limit(%f)=%f\n",nBackground,
			  FeldmanCousins.getAverageUpperLimit(nBackground));

	System.out.format(" least signal for 5sigma discovery=%f\n",
			  FeldmanCousins.getLeastSignalForDiscovery(nBackground,5.0));

        DataInputStream input = new DataInputStream(System.in);  
        BufferedReader  d     = new BufferedReader(new InputStreamReader(input));  
        String buffer; 

	// Upperlimit with uncertain BGs
	System.err.print("calculate upperlimit with uncertain BG? [yes(1)/no(0)] ->");  
	buffer   = d.readLine();  
	if(Integer.valueOf(buffer).intValue()==1){
	    System.err.print("relative plus BG error ->");  
	    buffer   = d.readLine(); 
	    double maxBGerr = Double.valueOf(buffer).doubleValue();
	    System.err.print("relative minus BG error ->");  
	    buffer   = d.readLine(); 
	    double minBGerr = Double.valueOf(buffer).doubleValue();
	    FeldmanCousins.setRelativeBackgroundUncertainty(minBGerr,maxBGerr);

	    System.out.format(" Upper limit(%f+%f perc. %f perc. %d)=%f\n",
			      nBackground,maxBGerr,minBGerr,nObserved,
			      FeldmanCousins.getUpperLimit(nBackground,nObserved));

	}

	// Confidence level 
	System.err.print("calculate C.L. of nObserved for nBG and nSignal? [yes(1)/no(0)] ->");  
	buffer   = d.readLine();  
	if(Integer.valueOf(buffer).intValue()==1){
	    treatBGuncertainties = false;
	    System.err.print("nSignal ->");  
	    buffer   = d.readLine(); 
	    double nSignal = Double.valueOf(buffer).doubleValue();
	    double prob = FeldmanCousins.probabilityOfNobserved(nSignal, nBackground, nObserved);
	    System.out.format("  Conf. Level BG(%f) SIG(%f) Obs(%d) =%e\n",
			      nBackground,nSignal,nObserved,prob);

	    System.err.print("calculate C.L.  with uncertain BG? [yes(1)/no(0)] ->");  
	    buffer   = d.readLine();  
	    if(Integer.valueOf(buffer).intValue()==1){
		System.err.print("relative plus BG error ->");  
		buffer   = d.readLine(); 
		double maxBGerr = Double.valueOf(buffer).doubleValue();
		System.err.print("relative minus BG error ->");  
		buffer   = d.readLine(); 
		double minBGerr = Double.valueOf(buffer).doubleValue();
		FeldmanCousins.setRelativeBackgroundUncertainty(minBGerr,maxBGerr);

		prob = FeldmanCousins.probabilityOfNobserved(nSignal, nBackground, nObserved);
		System.out.format("  Conf. Level BG(%f) SIG(%f) Obs(%d) =%e\n",
				  nBackground,nSignal,nObserved,prob);
	    }
	}

    }

}
