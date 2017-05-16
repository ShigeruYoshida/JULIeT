package numRecipes;

import numRecipes.*;

import java.io.*;
import java.util.*;

/** 
    Calculate the binned Poisson Likelihood and its p-value
    for a series of (N_i, \mu_i) where N_i is the observed number
    and \mu_i is the predicted number at the i-th bin.

    You should call the method fillData(in) to make map of (N_i, \mu_i)

    Written originally by S. Yoshida for the IceCube EHE analysis.
    2012/8/15
*/
public class PoissonBinnedLikelihoodCalculator implements Serializable{

    private static final long serialVersionUID = -7112334109753025742L;
    private boolean isDataFilled = false;
    private boolean mapsHaveBeenGenerated = false;
    protected boolean useTheReplicaResults = false;
    private static double epsilon = 1.0e-12;
    private double expectedSum = 0.0;
    private long observedSum = 0;
    private long observedSumInTheReplicaExpetiment = 0;
    public static final double default_precision = 1.0e-6; // precision

    /** Map container for (N_i, \mu_i). N_i is the observed number in the real experiment*/
    protected Map realDataMap = null;
    protected Map replicaDataMap = null;
    protected double poissonLogLikelihood = 0.0;
    public boolean debugFlag = false;

    private RandomGenerator rand = null;


    /** Constructor. 
    */
    public PoissonBinnedLikelihoodCalculator(){
	realDataMap = new LinkedHashMap();
	mapsHaveBeenGenerated = true;
	rand = new RandomGenerator();
    }

    /** Calculate -log(Poisson) for a given set of (nObserved,expectedValue).
    */
    public static double getLogPoissonProbability(long nObserved, double expectedValue){
	double prob = SpecialFunctions.poisson(expectedValue, nObserved);
	double llh = Double.POSITIVE_INFINITY;
	if(prob >0.0) llh =  -Math.log(prob);

	return llh;
    }

    /**
       Reading a series of (N_i, \mu_i)  from the input stream.
       The format is
       <pre>
       nObserved_0 Expected-value_0 
       nObserved_1 Expected-value_1 
       .....
       </pre>
       You need a "space" before each end of line.
     */
    public void fillData(DataInputStream in) throws IOException{

	if(!mapsHaveBeenGenerated){
	    System.err.println(" maps to store the data have not been generated - must call the proper constructor");
	    System.exit(0);
	}

	poissonLogLikelihood = 0.0;
	expectedSum = 0.0;
	observedSum = 0;
	int binNumber = 0;
	// Reading data
	BufferedReader  d = new BufferedReader(new InputStreamReader(in));
	String buffer; int sep = 0; int sepstart = 0;
	char separator = ' ';
	while((buffer = d.readLine())!=null){
	    // line -- condition number, signal rate, iron rate, iron-prton-rate
	    sepstart = 0;
	    sep = buffer.indexOf(separator,sepstart+1);
	    long observedNumber =
		Long.valueOf(buffer.substring(sepstart+1,sep)).longValue( );
	    Long observedNumberObj = new Long(observedNumber);
	    sepstart = sep;

	    sep = buffer.indexOf(separator,sepstart+1);
	    double expected =
		Double.valueOf(buffer.substring(sepstart+1,sep)).doubleValue( );
	    Double expectedObject = new Double(expected);
	    sepstart = sep;

	    Map binDataMap = new LinkedHashMap();
	    binDataMap.put(observedNumberObj,expectedObject);

	    Integer binNumberObj = new Integer(binNumber);
	    realDataMap.put(binNumberObj,binDataMap);

	    poissonLogLikelihood += getLogPoissonProbability(observedNumber, expected);
	    expectedSum += expected;
	    observedSum += observedNumber;
	    binNumber++;

	}

	isDataFilled = true;

    }

    /**
       Reading a series of (N_i, \mu_i)  from the standard input by the interactive methof
     */
    public void fillData() throws IOException{
	DataInputStream input = new DataInputStream(System.in); 
	BufferedReader  d     = new BufferedReader(new InputStreamReader(input)); 
	String buffer; 

	poissonLogLikelihood = 0.0;
	expectedSum = 0.0;
	observedSum = 0;
	int binNumber = 0;
	boolean interactionEnd = false;
	do{
	    System.err.print("observed number ->"); 
	    buffer   = d.readLine(); 
	    long observedNumber = Long.valueOf(buffer).longValue();
	    Long observedNumberObj = new Long(observedNumber);

	    System.err.print("expected value ->"); 
	    buffer   = d.readLine(); 
	    double expected = Double.valueOf(buffer).doubleValue();
	    Double expectedObject = new Double(expected);

	    LinkedHashMap binDataMap = new LinkedHashMap();
	    binDataMap.put(observedNumberObj,expectedObject);

	    Integer binNumberObj = new Integer(binNumber);
	    realDataMap.put(binNumberObj,binDataMap);

	    poissonLogLikelihood += getLogPoissonProbability(observedNumber, expected);
	    expectedSum += expected;
	    observedSum += observedNumber;
	    binNumber++;

	    System.err.print("No more data? [yes(1)/no(0)] ->"); 
	    buffer   = d.readLine(); 
	    if(Integer.valueOf(buffer).intValue()==1) interactionEnd = true; 
	}while(!interactionEnd);

	isDataFilled = true;
    }

    /**
       print out the calculated Poisson binned log likelihood to the standard output
     */
    public void printLogLikelihood(boolean printData){
	if(isDataFilled){

	    if(printData){
		Iterator dataIterator = realDataMap.entrySet().iterator();
		if(useTheReplicaResults) dataIterator = replicaDataMap.entrySet().iterator();
		while(dataIterator.hasNext()){
		    Map.Entry entryData = (Map.Entry )(dataIterator.next());
		    LinkedHashMap binDataMap = (LinkedHashMap )(entryData.getValue());
		    Iterator binIterator = binDataMap.entrySet().iterator();
		    while(binIterator.hasNext()){
			Map.Entry entry = (Map.Entry )(binIterator.next());
			Long numberObj = (Long )(entry.getKey());
			Double expectedObj = (Double )(entry.getValue());
			
			System.out.format(" %5d %8.5f\n",numberObj.longValue(),expectedObj.doubleValue());
		    }
		}
	    }
	    System.out.format(" -log(PoissonLikelihood)=%e\n",poissonLogLikelihood);
	}
    }


    public void printLogLikelihood(){
	printLogLikelihood(false);
    }

    /**
       Return the sum of the expected values \mu_i.
     */
    public double getSumOfExpectedValues(){
	return    expectedSum;
    }


    /**
       Return the sum of the observed numbers \N_i.
     */
    public long getSumOfObservedValues(){
	if(!useTheReplicaResults) return    observedSum;
	else{
	    return observedSumInTheReplicaExpetiment;
	}
    }


    /**
       copy the Poisson bin data (observed number, expected number) stored in PoissonBinnedLikelihoodCalculator cal.
     */
    public void copyPoissonBinnedData(PoissonBinnedLikelihoodCalculator cal){
	if(useTheReplicaResults){
	    if(replicaDataMap == null) replicaDataMap = new LinkedHashMap();
	    replicaDataMap.clear();
	    replicaDataMap.putAll(cal.replicaDataMap);
	    observedSumInTheReplicaExpetiment = cal.observedSumInTheReplicaExpetiment;
	}else{
	    realDataMap.clear();
	    realDataMap.putAll(cal.realDataMap);
	    observedSum = cal.observedSum;
	}
    }

    /**
      copy the Map containing the Poisson binned data.
      You can call this method and put the data map generated by any external object
      so that this PoissonBinnedLikelihoodCalculator object gets functioning
      without calling the method fillData().

      The map structure should be:
      <pre>
      Map((Integer)BinNumber, LinkedHashMap((Long)Observed Number, (Double)Expected Number))
      </pre>
     */
    public void copyPoissonBinnedDataMap(Map poissonDataMap){
	realDataMap.clear();
	realDataMap.putAll(poissonDataMap);

	observedSum = 0;
	expectedSum = 0.0;
	Iterator dataIterator = realDataMap.entrySet().iterator();
	while(dataIterator.hasNext()){
	    Map.Entry entryData = (Map.Entry )(dataIterator.next());
	    LinkedHashMap binDataMap = (LinkedHashMap )(entryData.getValue());
	    Iterator binIterator = binDataMap.entrySet().iterator();
	    while(binIterator.hasNext()){
		Map.Entry entry = (Map.Entry )(binIterator.next());
		Long numberObj = (Long )(entry.getKey());
		long observedNumber = numberObj.longValue();
		observedSum += observedNumber;
		Double expectedObj = (Double )(entry.getValue());
		double expectedNumber = expectedObj.doubleValue();
		expectedSum += expectedNumber;
	    }
	}

    }


    /**
       Run replica expetiment one time following the Poissonian to the series of (N_i, mu_i).
       Return log likelihood -log(PoissonLikelihood).
       You must fill the data of (N_i, mu_i) first by, for example, calling the method fillData(DataInputStream in).
       <pre>
       double scaleFactor:  scaleFactor which is multiplied to mu_i; if set to 1.0, no fudge factor is introduced.
       </pre>
    */
    public double runReplicaExperiment(double scaleFactor){
	replicaDataMap = new LinkedHashMap();
	double llh = 0.0;
	observedSumInTheReplicaExpetiment = 0;
	Iterator dataIterator = realDataMap.entrySet().iterator();
	while(dataIterator.hasNext()){
	    Map.Entry entryData = (Map.Entry )(dataIterator.next());
	    Integer binNumberObj = (Integer )(entryData.getKey());

	    LinkedHashMap binDataMap = (LinkedHashMap )(entryData.getValue());
	    Iterator binIterator = binDataMap.entrySet().iterator();
	    while(binIterator.hasNext()){
		Map.Entry entry = (Map.Entry )(binIterator.next());
		Double expectedObj = (Double )(entry.getValue());
		double expectedValue = expectedObj.doubleValue();
		long observedNumber = rand.GetPoissonian(expectedValue*scaleFactor);
		observedSumInTheReplicaExpetiment += observedNumber;
		//System.err.format(" expected %f observed %d\n",expectedValue,observedNumber);
		llh += getLogPoissonProbability(observedNumber, expectedValue*scaleFactor);

		// Store the results to the replicalDataMap(binNumbr, map(observed,expected))
		LinkedHashMap binReplicaDataMap = new LinkedHashMap();
		Long observedNumberObj = new Long(observedNumber);
		binReplicaDataMap.put(observedNumberObj,expectedObj);
		replicaDataMap.put(binNumberObj,binReplicaDataMap);
	    }
	}
	return llh;
    }

    /**
       Run replica expetiment one time following the Poissonian to the series of (N_i, mu_i).
       Return log likelihood -log(PoissonLikelihood). No fudge factor, i.e., scaleFactor=1.0.
    */
    public double runReplicaExperiment(){
	return runReplicaExperiment(1.0);
    }

    /** call this method before calling the getLikelihood method if you want to calculate the log likelihood by
     the results by the replica experiment you had run withe the method runReplicaExperiment() */
    public void useTheResultsByTheReplicaExperiment(){
	useTheReplicaResults = true;
    }

    /** call this method before calling the getLikelihood method if you want to calculate the log likelihood by
     the real results read by the method fillData(). This is the default setting. */
    public void useTheResultsByTheRealExperiment(){
	useTheReplicaResults = false;
    }

    /** This method copies the map to store the results frmo the replica experiment to the map for storing the real results.
     Mostly for debugging purposes, but you may need to call this method to study the likelihood ratio distribution with
    the alternative (signal+bg) hypothesis.*/
    public void copyTheReplicaResults(){
	replicaDataMap = new LinkedHashMap();
	replicaDataMap.putAll(realDataMap);
    }


    /**
       Calculate the likelihood to the series of (N_i, mu_i*scale_factor).
       Return log likelihood -log(PoissonLikelihood).
       You must fill the data of (N_i, mu_i) first by, for example, calling the method fillData(DataInputStream in).
       <pre>
       double scaleFactor:  scaleFactor which is multiplied to mu_i; if set to 1.0, no fudge factor is introduced.
       </pre>
    */
    public double getLikelihood(double scaleFactor){
	double llh = 0.0;
	Iterator dataIterator = realDataMap.entrySet().iterator();
	if(useTheReplicaResults) dataIterator = replicaDataMap.entrySet().iterator();
	while(dataIterator.hasNext()){
	    Map.Entry entryData = (Map.Entry )(dataIterator.next());
	    LinkedHashMap binDataMap = (LinkedHashMap )(entryData.getValue());
	    Iterator binIterator = binDataMap.entrySet().iterator();

	    while(binIterator.hasNext()){
		Map.Entry entry = (Map.Entry )(binIterator.next());
		Double expectedObj = (Double )(entry.getValue());
		double expectedValue = expectedObj.doubleValue();
		Long numberObj = (Long )(entry.getKey());
		long observedNumber = numberObj.longValue();
		if(debugFlag) System.err.format(" expected %e observed %d\n",expectedValue,observedNumber);
		llh += getLogPoissonProbability(observedNumber, expectedValue*scaleFactor);
	    }
	}
	return llh;
    }

    /**
       Calculate the likelihood to the series of (N_i, bg_mu_i+sig_mu_i*signalFactor).
       Return log likelihood -log(PoissonLikelihood).
       You must fill the data of (N_i, mu_i) first by, for example, calling the method fillData(DataInputStream in).
       PoissonBinnedLikelihoodCalculator calSignal stores the binned poisson distribution for the signal hypothesis.
       <pre>
       double scaleFactor:  scaleFactor which is multiplied to bg_mu_i; if set to 1.0, no fudge factor is introduced.
       </pre>
    */
    public double getLikelihood(double scaleFactor, double signalFactor,PoissonBinnedLikelihoodCalculator calSignal){

	double llh = 0.0;
	Iterator dataIterator = realDataMap.entrySet().iterator();
	if(useTheReplicaResults) dataIterator = replicaDataMap.entrySet().iterator();

	Map realDataMapSignal = calSignal.realDataMap;
	Iterator signalDataIterator = realDataMapSignal.entrySet().iterator();
	while(dataIterator.hasNext()){
	    Map.Entry entryData = (Map.Entry )(dataIterator.next());
	    LinkedHashMap binDataMap = (LinkedHashMap )(entryData.getValue());
	    Integer binNumberObj = (Integer )(entryData.getKey());
	    int binNumber = binNumberObj.intValue();
	    Iterator binIterator = binDataMap.entrySet().iterator();

	    Map.Entry entrySignalData = (Map.Entry )(signalDataIterator.next());
	    LinkedHashMap binSignalDataMap = (LinkedHashMap )(entrySignalData.getValue());
	    Integer signalBinNumberObj = (Integer )(entryData.getKey());
	    int signalBinNumber = signalBinNumberObj.intValue();
	    Iterator binSignalIterator = binSignalDataMap.entrySet().iterator();
	    if(binNumber != signalBinNumber){
		System.err.format(" binNumber(%d) and signal binNumber(%d) is different!",
				  binNumber,signalBinNumber);
		return -1.0;
	    }
	    while(binIterator.hasNext()){
		Map.Entry entry = (Map.Entry )(binIterator.next());
		Double expectedObj = (Double )(entry.getValue());
		double expectedValue = expectedObj.doubleValue();
		Long numberObj = (Long )(entry.getKey());
		long observedNumber = numberObj.longValue();

		Map.Entry entrySignal = (Map.Entry )(binSignalIterator.next());
		Double expectedSignalObj = (Double )(entrySignal.getValue());
		double expectedSignalValue = expectedSignalObj.doubleValue();

		if(debugFlag) System.err.format(" expected BG %e expeted SIG %e observed %d\n",
						expectedValue,expectedSignalValue,observedNumber);
		llh += getLogPoissonProbability(observedNumber, scaleFactor*expectedValue+expectedSignalValue*signalFactor);
	    }
	}
	return llh;
    }

    /**
       Calculate the signalFactor to maximize likelihood to the series of (N_i, bg_mu_i+sig_mu_i*signalFactor).
       You must fill the data of (N_i, mu_i) first by, for example, calling the method fillData(DataInputStream in).
       PoissonBinnedLikelihoodCalculator calSignal stores the binned poisson distribution for the signal hypothesis.
       <pre>
       double scaleFactor:  scaleFactor which is multiplied to bg_mu_i; if set to 1.0, no fudge factor is introduced.
       </pre>
    */
    public double getSignalFactor(double scaleFactor, PoissonBinnedLikelihoodCalculator calSignal){

	double totalSIG = calSignal.getSumOfExpectedValues();

	long totalObs = getSumOfObservedValues();
	double signalFactor = (double )totalObs/totalSIG; // the initial guess
	double signalFactorSearchRangeMin = 0.0;
	double signalFactorSearchRangeMax = 2.0*signalFactor;
	if(totalObs<=0) signalFactorSearchRangeMax = 0.1*totalSIG;
	int nSteps = 100; // 1% steps
	double deltaFactor = (signalFactorSearchRangeMax-signalFactorSearchRangeMin)/(double)nSteps;

	double llhMin = Double.POSITIVE_INFINITY;
	double signalFactorMin = 0.0;
	for(int i=0; i<=nSteps; i++){
	    double sigFactorSearch = signalFactorSearchRangeMin +deltaFactor*(double )i;
	    double llh = getLikelihood(scaleFactor, sigFactorSearch, calSignal);
	    //System.err.format("  -- sigFactorSearch %e llh %e\n",sigFactorSearch,llh);
	    if(llh<llhMin){
		llhMin = llh;
		signalFactorMin = sigFactorSearch;
	    }
	}

	return signalFactorMin;
    }

    /**
       Calculate the likelihood to the series of (N_i, mu_i). No fudge factor.
    */
    public double getLikelihood(){
	return getLikelihood(1.0);
    }

    /** 
	Multiply the scale factor to the expected number in each of the bins.
     */
    public void multiplyFactor(double scaleFactor){

	expectedSum = 0.0;

	Iterator dataIterator = realDataMap.entrySet().iterator();
	if(useTheReplicaResults) dataIterator = replicaDataMap.entrySet().iterator();
	Map newDataMap = new LinkedHashMap();

	while(dataIterator.hasNext()){
	    Map.Entry entryData = (Map.Entry )(dataIterator.next());
	    LinkedHashMap binDataMap = (LinkedHashMap )(entryData.getValue());
	    Integer binNumberObj = (Integer )(entryData.getKey());
	    int binNumber = binNumberObj.intValue();
	    Iterator binIterator = binDataMap.entrySet().iterator();

	    while(binIterator.hasNext()){
		Map.Entry entry = (Map.Entry )(binIterator.next());
		Double expectedObj = (Double )(entry.getValue());
		double expectedValue = expectedObj.doubleValue()*scaleFactor;
		Long numberObj = (Long )(entry.getKey());
		long observedNumber = numberObj.longValue();
		Double expectedNewObj = new Double(expectedValue);

		expectedSum += expectedValue;

		LinkedHashMap newBinDataMap = new LinkedHashMap();
		newBinDataMap.put(numberObj,expectedNewObj);
		newDataMap.put(binNumberObj,newBinDataMap);
	    }
	}

	if(useTheReplicaResults){
	    replicaDataMap.clear();
	    replicaDataMap.putAll(newDataMap);
	    newDataMap.clear();
	}else{
	    realDataMap.clear();
	    realDataMap.putAll(newDataMap);
	    newDataMap.clear();
	}

    }

    /** 
	Copy the observed number in calSignal to the internal map to this object.
     */
    public void copyObservedNumbers(PoissonBinnedLikelihoodCalculator calSignal){

	long observedSumInCalSignal = 0;

	Iterator dataIterator = realDataMap.entrySet().iterator();
	if(useTheReplicaResults) dataIterator = replicaDataMap.entrySet().iterator();
	Map newDataMap = new LinkedHashMap();

	Map realDataMapSignal = calSignal.realDataMap;
	if(calSignal.useTheReplicaResults) realDataMapSignal = calSignal.replicaDataMap;
	Iterator signalDataIterator = realDataMapSignal.entrySet().iterator();
	while(signalDataIterator.hasNext()){
	    Map.Entry entryData = (Map.Entry )(dataIterator.next());
	    LinkedHashMap binDataMap = (LinkedHashMap )(entryData.getValue());
	    Integer binNumberObj = (Integer )(entryData.getKey());
	    int binNumber = binNumberObj.intValue();
	    Iterator binIterator = binDataMap.entrySet().iterator();

	    Map.Entry entrySignalData = (Map.Entry )(signalDataIterator.next());
	    LinkedHashMap binSignalDataMap = (LinkedHashMap )(entrySignalData.getValue());
	    Integer signalBinNumberObj = (Integer )(entryData.getKey());
	    int signalBinNumber = signalBinNumberObj.intValue();
	    Iterator binSignalIterator = binSignalDataMap.entrySet().iterator();
	    if(binNumber != signalBinNumber){
		System.err.format(" binNumber(%d) and signal binNumber(%d) is different!",
				  binNumber,signalBinNumber);
	    }
	    while(binSignalIterator.hasNext()){
		Map.Entry entry = (Map.Entry )(binIterator.next());
		Double expectedObj = (Double )(entry.getValue());
		double expectedValue = expectedObj.doubleValue();
		Long numberObj = (Long )(entry.getKey());
		long observedNumber = numberObj.longValue();

		Map.Entry entrySignal = (Map.Entry )(binSignalIterator.next());
		Double expectedSignalObj = (Double )(entrySignal.getValue());
		double expectedSignalValue = expectedSignalObj.doubleValue();
		Long numberSignalObj = (Long )(entrySignal.getKey());
		long observedSignalNumber = numberSignalObj.longValue();
		observedSumInCalSignal += observedSignalNumber;

		LinkedHashMap newBinDataMap = new LinkedHashMap();
		newBinDataMap.put(numberSignalObj,expectedObj);
		newDataMap.put(binNumberObj,newBinDataMap);
	    }
	}

	if(useTheReplicaResults){
	    replicaDataMap.clear();
	    replicaDataMap.putAll(newDataMap);
	    newDataMap.clear();
	    observedSumInTheReplicaExpetiment = observedSumInCalSignal;
	}else{
	    realDataMap.clear();
	    realDataMap.putAll(newDataMap);
	    newDataMap.clear();
	    observedSum = observedSumInCalSignal;
	}

    }

    /** 
	Copy the expected number in calSignal to the internal map to this object.
	This method is used when the null pypothesis is NOT background only
	but signalFactor*signal +  background. 
     */
    public void copyExpectedNumbers(double signalFactor, PoissonBinnedLikelihoodCalculator calSignal){

	expectedSum = 0.0;
	Iterator dataIterator = realDataMap.entrySet().iterator();
	if(useTheReplicaResults) dataIterator = replicaDataMap.entrySet().iterator();
	Map newDataMap = new LinkedHashMap();

	Map realDataMapSignal = calSignal.realDataMap;
	if(calSignal.useTheReplicaResults) realDataMapSignal = calSignal.replicaDataMap;
	Iterator signalDataIterator = realDataMapSignal.entrySet().iterator();
	while(signalDataIterator.hasNext()){
	    Map.Entry entryData = (Map.Entry )(dataIterator.next());
	    LinkedHashMap binDataMap = (LinkedHashMap )(entryData.getValue());
	    Integer binNumberObj = (Integer )(entryData.getKey());
	    int binNumber = binNumberObj.intValue();
	    Iterator binIterator = binDataMap.entrySet().iterator();

	    Map.Entry entrySignalData = (Map.Entry )(signalDataIterator.next());
	    LinkedHashMap binSignalDataMap = (LinkedHashMap )(entrySignalData.getValue());
	    Integer signalBinNumberObj = (Integer )(entryData.getKey());
	    int signalBinNumber = signalBinNumberObj.intValue();
	    Iterator binSignalIterator = binSignalDataMap.entrySet().iterator();
	    if(binNumber != signalBinNumber){
		System.err.format(" binNumber(%d) and signal binNumber(%d) is different!",
				  binNumber,signalBinNumber);
	    }
	    while(binSignalIterator.hasNext()){
		Map.Entry entry = (Map.Entry )(binIterator.next());
		Double expectedObj = (Double )(entry.getValue());
		double expectedValue = expectedObj.doubleValue();
		Long numberObj = (Long )(entry.getKey());
		long observedNumber = numberObj.longValue();

		Map.Entry entrySignal = (Map.Entry )(binSignalIterator.next());
		Double expectedSignalObj = (Double )(entrySignal.getValue());
		double expectedSignalValue = expectedSignalObj.doubleValue()*signalFactor;
		expectedSum += expectedSignalValue;
		expectedSignalObj = new Double(expectedSignalValue);
		Long numberSignalObj = (Long )(entrySignal.getKey());
		long observedSignalNumber = numberSignalObj.longValue();

		LinkedHashMap newBinDataMap = new LinkedHashMap();
		newBinDataMap.put(numberObj,expectedSignalObj);
		newDataMap.put(binNumberObj,newBinDataMap);
	    }
	}

	if(useTheReplicaResults){
	    replicaDataMap.clear();
	    replicaDataMap.putAll(newDataMap);
	    newDataMap.clear();
	}else{
	    realDataMap.clear();
	    realDataMap.putAll(newDataMap);
	    newDataMap.clear();
	}

    }

    /** 
	Copy the expected number in calSignal to the internal map to this object.
	This method is used when the null pypothesis is NOT background only
	but signal +  background. 
     */
    public void copyExpectedNumbers(PoissonBinnedLikelihoodCalculator calSignal){
	copyExpectedNumbers(1.0,calSignal);
    }

    /** 
	Add the expected number in calSignal to the internal map to this object.
	This method is used when the null pypothesis is NOT background only
	but signalFactor*signal +  background. 
     */
    public void addExpectedNumbers(double signalFactor, PoissonBinnedLikelihoodCalculator calSignal){

	expectedSum = 0.0;
	Iterator dataIterator = realDataMap.entrySet().iterator();
	if(useTheReplicaResults) dataIterator = replicaDataMap.entrySet().iterator();
	Map newDataMap = new LinkedHashMap();

	Map realDataMapSignal = calSignal.realDataMap;
	if(calSignal.useTheReplicaResults) realDataMapSignal = calSignal.replicaDataMap;
	Iterator signalDataIterator = realDataMapSignal.entrySet().iterator();
	while(signalDataIterator.hasNext()){
	    Map.Entry entryData = (Map.Entry )(dataIterator.next());
	    LinkedHashMap binDataMap = (LinkedHashMap )(entryData.getValue());
	    Integer binNumberObj = (Integer )(entryData.getKey());
	    int binNumber = binNumberObj.intValue();
	    Iterator binIterator = binDataMap.entrySet().iterator();

	    Map.Entry entrySignalData = (Map.Entry )(signalDataIterator.next());
	    LinkedHashMap binSignalDataMap = (LinkedHashMap )(entrySignalData.getValue());
	    Integer signalBinNumberObj = (Integer )(entryData.getKey());
	    int signalBinNumber = signalBinNumberObj.intValue();
	    Iterator binSignalIterator = binSignalDataMap.entrySet().iterator();
	    if(binNumber != signalBinNumber){
		System.err.format(" binNumber(%d) and signal binNumber(%d) is different!",
				  binNumber,signalBinNumber);
	    }
	    while(binSignalIterator.hasNext()){
		Map.Entry entry = (Map.Entry )(binIterator.next());
		Double expectedObj = (Double )(entry.getValue());
		double expectedValue = expectedObj.doubleValue();
		Long numberObj = (Long )(entry.getKey());
		long observedNumber = numberObj.longValue();

		Map.Entry entrySignal = (Map.Entry )(binSignalIterator.next());
		Double expectedSignalObj = (Double )(entrySignal.getValue());
		double expectedSignalValue = expectedSignalObj.doubleValue()*signalFactor;
		expectedSignalObj = new Double(expectedSignalValue+expectedValue);
		expectedSum += expectedSignalValue+expectedValue;
		Long numberSignalObj = (Long )(entrySignal.getKey());
		long observedSignalNumber = numberSignalObj.longValue();

		LinkedHashMap newBinDataMap = new LinkedHashMap();
		newBinDataMap.put(numberObj,expectedSignalObj);
		newDataMap.put(binNumberObj,newBinDataMap);
	    }
	}

	if(useTheReplicaResults){
	    replicaDataMap.clear();
	    replicaDataMap.putAll(newDataMap);
	    newDataMap.clear();
	}else{
	    realDataMap.clear();
	    realDataMap.putAll(newDataMap);
	    newDataMap.clear();
	}

    }

    /** 
	Add the expected number in calSignal to the internal map to this object.
	This method is used when the null pypothesis is NOT background only
	but signal +  background. 
     */
    public void addExpectedNumbers(PoissonBinnedLikelihoodCalculator calSignal){
	addExpectedNumbers(1.0,calSignal);
    }


    /** 
	Replace the expected number with the observed number so that
        the expected number vector exactly equals to the observed number vector.
        This operation is used for building a saturated likelihood.
     */
    public void replaceExpectedNumbersWithObservedNumbers(){

	Map newDataMap = new LinkedHashMap();

	Iterator dataIterator = realDataMap.entrySet().iterator();
	if(useTheReplicaResults) dataIterator = replicaDataMap.entrySet().iterator();
	while(dataIterator.hasNext()){
	    Map.Entry entryData = (Map.Entry )(dataIterator.next());

	    Integer binNumberObj = (Integer )(entryData.getKey());
	    LinkedHashMap binDataMap = (LinkedHashMap )(entryData.getValue());
	    Iterator binIterator = binDataMap.entrySet().iterator();

	    while(binIterator.hasNext()){
		Map.Entry entry = (Map.Entry )(binIterator.next());
		Double expectedObj = (Double )(entry.getValue());
		double expectedValue = expectedObj.doubleValue();
		Long numberObj = (Long )(entry.getKey());
		long observedNumber = numberObj.longValue();
		double newExpectedNumber = (double )observedNumber;
		if(newExpectedNumber <= 0.0) newExpectedNumber = 1.0e-20;
		Double newExpectedObj = new Double(newExpectedNumber);
		
		LinkedHashMap newBinDataMap = new LinkedHashMap();
		newBinDataMap.put(numberObj,newExpectedObj);
		newDataMap.put(binNumberObj,newBinDataMap);
	    }
	}

	if(useTheReplicaResults){
	    replicaDataMap.clear();
	    replicaDataMap.putAll(newDataMap);
	    newDataMap.clear();
	}else{
	    realDataMap.clear();
	    realDataMap.putAll(newDataMap);
	    newDataMap.clear();
	}
    }

    /** 
	Replace the expected number with the observed number - expected number in calBase so that
        the expected number vector maximizes the likelihood under the conditaion that
	the minimum expected numbers are bounded by the expected numbers stored in calBase.
	calBase usually stored the numbers expected in the pure background hypothesis.

        This operation is used for building a saturated likelihood.
     */
    public void replaceExpectedNumbersWithObservedNumbers(PoissonBinnedLikelihoodCalculator calBase){

	expectedSum = 0.0;
	Map newDataMap = new LinkedHashMap();

	Iterator dataIterator = realDataMap.entrySet().iterator();
	if(useTheReplicaResults) dataIterator = replicaDataMap.entrySet().iterator();

	Map realDataMapBase = calBase.realDataMap;
	if(calBase.useTheReplicaResults) realDataMapBase = calBase.replicaDataMap;
	Iterator baseDataIterator = realDataMapBase.entrySet().iterator();
	while(baseDataIterator.hasNext()){
	    Map.Entry entryData = (Map.Entry )(dataIterator.next());
	    LinkedHashMap binDataMap = (LinkedHashMap )(entryData.getValue());
	    Integer binNumberObj = (Integer )(entryData.getKey());
	    int binNumber = binNumberObj.intValue();
	    Iterator binIterator = binDataMap.entrySet().iterator();

	    Map.Entry entryBaseData = (Map.Entry )(baseDataIterator.next());
	    LinkedHashMap binBaseDataMap = (LinkedHashMap )(entryBaseData.getValue());
	    Integer baseBinNumberObj = (Integer )(entryData.getKey());
	    int baseBinNumber = baseBinNumberObj.intValue();
	    Iterator binBaseIterator = binBaseDataMap.entrySet().iterator();
	    if(binNumber != baseBinNumber){
		System.err.format(" binNumber(%d) and base binNumber(%d) is different!",
				  binNumber,baseBinNumber);
	    }
	    while(binBaseIterator.hasNext()){
		Map.Entry entry = (Map.Entry )(binIterator.next());
		Double expectedObj = (Double )(entry.getValue());
		double expectedValue = expectedObj.doubleValue();
		Long numberObj = (Long )(entry.getKey());
		long observedNumber = numberObj.longValue();

		Map.Entry entryBase = (Map.Entry )(binBaseIterator.next());
		Double expectedBaseObj = (Double )(entryBase.getValue());
		double expectedBaseValue = expectedBaseObj.doubleValue();
		double newExpectedSignalValue = (double)observedNumber - expectedBaseValue;
		if(newExpectedSignalValue <=0.0) newExpectedSignalValue = 0.0;
		expectedBaseObj = new Double(newExpectedSignalValue+expectedBaseValue);
		expectedSum += newExpectedSignalValue+expectedBaseValue;
		Long numberBaseObj = (Long )(entryBase.getKey());
		long observedBaseNumber = numberBaseObj.longValue();

		LinkedHashMap newBinDataMap = new LinkedHashMap();
		newBinDataMap.put(numberObj,expectedBaseObj);
		newDataMap.put(binNumberObj,newBinDataMap);
	    }
	}

	if(useTheReplicaResults){
	    replicaDataMap.clear();
	    replicaDataMap.putAll(newDataMap);
	    newDataMap.clear();
	}else{
	    realDataMap.clear();
	    realDataMap.putAll(newDataMap);
	    newDataMap.clear();
	}
    }

    /**
       Calculate the p-value of the experiment to yield the series of bins (N_obs, \mu_i).
       You must fill the data of (N_i, mu_i) first by, for example, calling the method fillData(DataInputStream in).
       The MC method with the binned Poisson likelihood computes the p-value.
       <pre>
       double precision :  precision you want to caluclate p-value. Running replica expetiments 1/precision times
                           to compute p-value.
       </pre>
     */
    public double getPvalue(double precision, double scaleFactor){
	if(!isDataFilled){
	    System.err.println(" You must fill the map to store the bins of (N_obs_i, mu_i) first!");
	    return 0.0;
	}
	double pValue = 0.0;
	poissonLogLikelihood = getLikelihood(scaleFactor);
	if(precision>0.0){
	    long trialFactor = (long )Math.ceil(1.0/precision);
	    long runNumber = 0;
	    while(runNumber<trialFactor){
		double llhOfThisReplicaExperiment = runReplicaExperiment(scaleFactor);
		//useTheResultsByTheReplicaExperiment();
		//double llhOfThisReplicaExperiment2 = getLikelihood(scaleFactor);
		//useTheResultsByTheRealExperiment();
		//if(llhOfThisReplicaExperiment!=llhOfThisReplicaExperiment2) System.err.println("something is wrong");
		if(llhOfThisReplicaExperiment> poissonLogLikelihood) pValue += 1.0;
		runNumber ++;
	    }
	    pValue = pValue/(double )trialFactor;
	}
	return pValue;
    }

    /**
       Calculate the p-value of the experiment to yield the series of bins (N_obs, \mu_i) with
       the default_precision.
    */
    public double getPvalue(){
	return getPvalue(default_precision,1.0);
    }

    /**
       This method is purely for debugging purposes.
       Run an expetiment to yield N_i from the expected value \mu_i following Poissonian.
       Then you can check if the calculated p-value makes sense by calling getPvalue().

       You must fill the data of (N_i, mu_i) first by, for example, calling the method fillData(DataInputStream in).
    */
    protected void runPoissonExperiment(){
	if(!isDataFilled){
	    System.err.println(" You must fill the map to store the bins of (N_obs_i, mu_i) first!");
	}
	poissonLogLikelihood = runReplicaExperiment();
    }


    public static void outputPoissonBinnedLikelihoodCalculator(PoissonBinnedLikelihoodCalculator cal, 
						OutputStream out)  throws IOException {

	ObjectOutputStream objectOut = new ObjectOutputStream(out);

	objectOut.writeObject(cal);
	objectOut.flush();
    }

    public static PoissonBinnedLikelihoodCalculator inputPoissonBinnedLikelihoodCalculator(
						InputStream in)  throws IOException {

        PoissonBinnedLikelihoodCalculator cal = null;
        try{

            ObjectInputStream objectIn = new ObjectInputStream(in);
            cal = (PoissonBinnedLikelihoodCalculator )objectIn.readObject();

        }catch(ClassNotFoundException e){
            System.err.println("Caught ClassNotFoundException: " + 
                               e.getMessage( ));
            System.exit(0);
        }catch (EOFException e){
            //System.err.println("Caught EOFException: " + e.getMessage());
            cal = null;
            return cal;
        }

        return cal;

    }


    /** A simple main method */
    public static void main(String[] args) throws IOException{
 

	boolean readDataFromFileName = false;
	boolean debugMode = false;
	String eventRateFileName = null;

	if(args.length==0){
	    System.out.println("Usage: PoissonBinnedLikelihoodCalculator interactive?(yes 0 no 1) filename-to-read-data (debug mode)");
	    System.exit(0);
	}else { 
	    if(Integer.valueOf(args[0]).intValue() == 1) readDataFromFileName = true;
	    if(args.length >= 2) eventRateFileName = args[1];
	    if(args.length == 3) debugMode=true;
	    else if(args.length > 3){
		System.out.println("Illeagal arguments");
		System.exit(0);
	    }
	}

	PoissonBinnedLikelihoodCalculator cal = new PoissonBinnedLikelihoodCalculator();
	if(debugMode) cal.debugFlag = true;

	if(!readDataFromFileName){// calculate the llh by the interactive way
	    cal.fillData();

	}else{ // calculate the llh for the data from the file
	    DataInputStream in = new DataInputStream(ClassLoader.getSystemResourceAsStream(eventRateFileName));
	    cal.fillData(in);
	}

	cal.printLogLikelihood(true);
	double pValue = cal.getPvalue();
	System.out.format(" p-value=%e\n",pValue);

	//System.out.println(" Run the poisson experiment");
	//cal.runPoissonExperiment();
	//cal.printLogLikelihood(false);
	//pValue = cal.getPvalue();
	//System.out.format(" p-value=%e\n",pValue);

    }
}
