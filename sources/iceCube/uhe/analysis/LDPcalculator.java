package iceCube.uhe.analysis;

import iceCube.uhe.analysis.*;
import numRecipes.*;

import java.io.*;
import java.util.*;

/** 
    Calculate either the model rejection factor (MRF)
    or the least signal discovery potential (LDP). Rely on numRecipes.FeldmanCousins.java

    If you read the method FillRateData(in), you can calculate
    the MRF/LDP from all the series of (signalRate, bgRate) and sort them out
    in order of MRF/LDP

    Written originally by S. Yoshida for the IceCube EHE analysis.
    2012/3/7
*/
public class LDPcalculator {

    private boolean isSignalFilled = false;
    private boolean isBGFilled = false;
    private boolean mapsHaveBeenGenerated = false;
    private static double epsilon = 1.0e-12;

    protected Map ldpMap = null;
    protected Map sigMap = null;
    protected Map bgMap = null;

    public boolean isLDP = true;   // if false, calculate MRF instead.
    public static double significance = 4.0;  // 4 sigma for claiming "discovery"


    /** Constructor. 
    */
    public LDPcalculator(boolean isLDP, boolean generateMap){
	this.isLDP = isLDP;
	if(generateMap){
	    ldpMap = new TreeMap();
	    sigMap = new LinkedHashMap();
	    bgMap = new LinkedHashMap();
	    mapsHaveBeenGenerated = true;
	}
    }

    public LDPcalculator(boolean isLDP){
	this(isLDP, false);
    }

    /** Calculate the model rejection factor for a given set of (nSignal,nBG).
    */
    public static double getModelRejectionFactor(double nSignal, double nBG){
	// average Feldman-Cousins upper limit
	double mu = FeldmanCousins.getAverageUpperLimit(nBG);
	if(nSignal<=epsilon) nSignal = epsilon;
	double modelRejectionFactor = mu/nSignal;

	return modelRejectionFactor;
    }

    /** Calculate the signal discovery potential for a given set of (nSignal,nBG) at significanceOfDiscovery sigma.
    */
    public static double getLeastDiscoveryPotential(double nSignal, double nBG, double significanceOfDiscovery){


	// least signal for discovery - 5sigma significance
	FeldmanCousins.setRangeOfSignalIntervalCalculation(0.0,100.0,0.00005);
	double mu = FeldmanCousins.getLeastSignalForDiscovery(nBG,significanceOfDiscovery);
	if(nSignal<=epsilon) nSignal = epsilon;
	double signalDiscoveryPotential = mu/nSignal;

	return  signalDiscoveryPotential;
    }

    /**
       Reading the signal and bg rates from the input stream.
       The format is
       <pre>
       index% signal-rate bg-rate bg-proton-rate 
       </pre>
       You need a "space" before each end of line.
     */
    public void readAndFillLDP(DataInputStream in) throws IOException{

	if(!mapsHaveBeenGenerated){
	    System.err.println(" maps to store the data have not been generated - must call the proper constructor");
	    System.exit(0);
	}

	// Reading data
	BufferedReader  d = new BufferedReader(new InputStreamReader(in));
	String buffer; int sep = 0; int sepstart = 0;
	char separator = ' ';
	while((buffer = d.readLine())!=null){
	    // line -- condition number, signal rate, iron rate, iron-prton-rate
	    sepstart = 0;
	    sep = buffer.indexOf(separator,sepstart+1);
	    int conditionNumber =
		Integer.valueOf(buffer.substring(sepstart+1,sep)).intValue( );
	    Integer conditionNumberObj = new Integer(conditionNumber);
	    sepstart = sep;

	    sep = buffer.indexOf(separator,sepstart+1);
	    double sigRate =
		Double.valueOf(buffer.substring(sepstart+1,sep)).doubleValue( );
	    Double sigObject = new Double(sigRate);
	    sigMap.put(conditionNumberObj,sigObject);
	    sepstart = sep;

	    sep = buffer.indexOf(separator,sepstart+1);
	    double bgRate =
		Double.valueOf(buffer.substring(sepstart+1,sep)).doubleValue( );
	    Double bgObject = new Double(bgRate);
	    bgMap.put(conditionNumberObj,bgObject);
	    sepstart = sep;

	    sep = buffer.indexOf(separator,sepstart+1);
	    double bgProtonRate =
		Double.valueOf(buffer.substring(sepstart+1,sep)).doubleValue( );
	    sepstart = sep;

	    //System.err.format(" condition=%5d sig(%f) BG(%f)\n",conditionNumber,sigRate,bgRate);

	    if(!isLDP){ // calculate the model rejection factor
		double mrf = getModelRejectionFactor(sigRate,bgRate);
		Double mrfObject = new Double(mrf);

		ldpMap.put(mrfObject, conditionNumberObj);
	    }else{ // calcuate the Least Discovery Potential
		double ldp = getLeastDiscoveryPotential(sigRate,bgRate,significance);
		Double ldpObject = new Double(ldp);

		ldpMap.put(ldpObject,conditionNumberObj);
	    }


	}

	isSignalFilled = true; isBGFilled = true;

    }

    /**
       print out the calculated LDPs to the standard output
     */
    public void printResults(){
	if(isSignalFilled && isBGFilled){

	    Iterator ldpIterator = ldpMap.entrySet().iterator();
	    while(ldpIterator.hasNext()){
		Map.Entry entry = (Map.Entry )(ldpIterator.next());
		Double ldpObj = (Double )(entry.getKey());
		Integer numberObj = (Integer )(entry.getValue());

		Double sigRateObj = (Double )sigMap.get(numberObj);
		Double bgRateObj = (Double )bgMap.get(numberObj);

		System.out.format(" %5d %8.5f %11.9f %11.9f\n",numberObj.intValue(),ldpObj.doubleValue(),
				  sigRateObj.doubleValue(),bgRateObj.doubleValue());
	    }
	}
    }


    /** A simple main method */
    public static void main(String[] args) throws IOException{
 

	boolean readDataFromFileName = false;
	boolean isLDP = false;
	double bgRate = 1.0;
	double sigRate = 1.0;
	String eventRateFileName = null;

	if(args.length==0){
	    System.out.println("Usage: LDPcalculator filename  LDP?(0 or 1)");
	    System.out.println("Usage: LDPcalculator sigRate bgRate LDP?(0 or 1)");
	    System.exit(0);
	}else if(args.length==2) { // calculate the LDPs for the data from the file
	    eventRateFileName = args[0];
	    if(Double.valueOf(args[1]).doubleValue() == 0) isLDP = true;
	    readDataFromFileName = true;
	}else if(args.length==3){ // calculate the LDPs for the rates given in the arguments
	    sigRate = Double.valueOf(args[0]).doubleValue();
	    bgRate = Double.valueOf(args[1]).doubleValue();
	    if(Double.valueOf(args[2]).doubleValue() == 0) isLDP = true;
	}else{
	    System.out.println("Illeagal arguments");
	    System.exit(0);
	}


	if(!readDataFromFileName){// calculate the LDPs for the rates given in the arguments

	    if(isLDP){ // calculate the Least Discovery Potential
		double ldp = LDPcalculator.getLeastDiscoveryPotential(sigRate,bgRate,significance); // 5sigma
		System.out.format(" LDP %8.4f sig(%10.5f ) bg(%10.5f )\n",ldp,sigRate,bgRate);
	    }else { //  calculate the Model Rejection Factor
		double mrf = LDPcalculator.getModelRejectionFactor(sigRate,bgRate);
		System.out.format(" MRF %8.4f sig(%10.5f ) bg(%10.5f )\n",mrf,sigRate,bgRate);
	    }

	}else{ // calculate the LDPs for the data from the file
	    if(isLDP) System.err.println("sorted by the Least Discovery Potential");
	    else System.err.println("sorted by the Model Rejection Factor");
	    LDPcalculator cal = new LDPcalculator(isLDP,true);

	    DataInputStream in = new DataInputStream(ClassLoader.getSystemResourceAsStream(eventRateFileName));
	    cal.readAndFillLDP(in);

	    cal.printResults();
	
	}
    }
}
