package iceCube.uhe.analysis;

import iceCube.uhe.analysis.*;
import numRecipes.*;

import java.io.*;
import java.util.*;

/**

   Build the Map to contain the PoissonBinnedLikelihoodCalculator
   reading from the f2k data via standard input.
   The f2k data can be created by DumpNuRateHistogramFromTree.C.
   By this way, you can leave the Root framework and move any downstream
   analysis in the pure java world.

   The map format is
   <pre>
   Map(Map("parameter name", (Double) parameter value), PoissonBinnedLikelihoodCalculator object)
   </pre>

   The input f2k file format should be:
   <pre>
    event event# maximal-redshift  evoltion parameter
    event event# cosZ  energyProxy expected-event-number
    event event# cosZ  energyProxy expected-event-number
    .....

   </pre>
*/


public class PoissonBinnedLikelihoodMapBuilder {

    protected int numberOfBinsNPE = 1240;
    protected int numberOfBinsRecoE = 1344;
    protected int numberOfBins = numberOfBinsNPE;

    /** map container for the poisson binned likelihood calculator objects */
    protected Map pblMap = null;  

    public final static  String parameterRedshift = "zmax";
    public final static  String parameterEvolution = "evolution";
    public final static  String parameterPowerLawIndex = "powerLawIndex";
    public final static  String parameterEnergyFluxAtBase = "eFluxAtBase";
    public final static  String parameterEnergyBase = "energyBase";
    public final static  String parameterCutoffEnergy = "cutoffEnergy";

    /** the flag to tell the input file is GZK or power law*/
    private boolean isGZK = true;

    private boolean rebinned = false;

    private boolean noObservedEvent = false;

    /** number of observed events  */
    public final static int numberOfObservedEvents = 2;

    /** the observed event cosZ and log(energy(or NPE))*/
    public final static double[][] observedEventsNPE = { 
	{4.59945035846730244e-02,4.87099}, // Event Calliandra
	{-0.2085,5.11394}, // Event Kloppo
    };
    public final static double[][] observedEventsRecoE = { 
	{-0.17,6.2367891}, // Event Calliandra
	{-0.2085,6.719745}, // Event Kloppo
    };

    protected static double[][] observedEvents = null;

    /** cos z bin width */
    public static double cosZbinWidth = 0.025; // half width

    /** energy bin width */
    public final static double npeBinWidth = 0.05; // half width
    public final static double recoEBinWidth = 0.2; // half width
    protected static double energyBinWidth = npeBinWidth;

    private boolean hasBeenProcessed = false;


    /** constructor */
    public PoissonBinnedLikelihoodMapBuilder(boolean useNPEbin, boolean rebinned){
	pblMap = new LinkedHashMap();
	observedEvents = new double[2][2];
	if(useNPEbin){
	    numberOfBins = numberOfBinsNPE;
	    energyBinWidth = npeBinWidth;
	    for(int i=0;i<2;i++){
		for(int j=0;j<2;j++) observedEvents[i][j]=observedEventsNPE[i][j];
	    }
	}else{
	    numberOfBins = numberOfBinsRecoE;
	    energyBinWidth = recoEBinWidth;
	    for(int i=0;i<2;i++){
		for(int j=0;j<2;j++) observedEvents[i][j]=observedEventsRecoE[i][j];
	    }
	}

	if(rebinned){
	    //numberOfBins = numberOfBins/4;
	    //energyBinWidth = energyBinWidth*2.0;
	    //cosZbinWidth = cosZbinWidth*2.0;
	    numberOfBins = numberOfBins/8;
	    energyBinWidth = energyBinWidth*4.0;
	    cosZbinWidth = cosZbinWidth*2.0;
	}
    }


    /** return the number of observed events in this bin */
    protected static long getObservedNumber(double binCenterCosZ, double binCenterEnergyProxy){
	long numberOfEventsInThisBin = 0;
	for(int i=0;i<numberOfObservedEvents;i++){
	    double cosZOfThisEvent = observedEvents[i][0];
	    double energyOfThisEvent = observedEvents[i][1];
	    double deltaCosZ = Math.abs(binCenterCosZ-cosZOfThisEvent);
	    if(deltaCosZ< cosZbinWidth){
		double deltaEnergy = Math.abs(binCenterEnergyProxy-energyOfThisEvent);
		if(deltaEnergy<energyBinWidth) numberOfEventsInThisBin++;
	    }
	}
	return numberOfEventsInThisBin;
    }

    /** tell this object that the binned data is concerned with the GZK neutrinos with m and zmax varied. 
	This is default.
     */
    public void isGZK(){
	isGZK = true;
    }

    /** tell this object that the binned data is concerned with the power law fluxes  with cutoff energy etc  varied. 
     */
    public void isPowerLaw(){
	isGZK = false;
    }


    public void process(DataInputStream in) 
       throws IOException {

	int eventNumber = 0;

        // Reading data
        BufferedReader  d     = new BufferedReader(new InputStreamReader(in));
        String buffer; int sep = 0; int sepstart = 0;
        char separator = ' ';
        while((buffer = d.readLine())!=null){
            try{

		Map parametersMap = new LinkedHashMap();

		if(isGZK){
		    // 1st line -- eventNumber, redshift, evolution
		    sepstart = 0;
		    sep = buffer.indexOf(separator,sepstart+1);
		    sepstart = sep;

		    sep = buffer.indexOf(separator,sepstart+1);
		    eventNumber =
			Integer.valueOf(buffer.substring(sepstart+1,sep)).intValue( );
		    sepstart = sep;

		    sep = buffer.indexOf(separator,sepstart+1);
		    double zmax =
			Double.valueOf(buffer.substring(sepstart+1,sep)).doubleValue( );
		    Double zmaxObj = new Double(zmax);
		    parametersMap.put(parameterRedshift,zmaxObj);
		    sepstart = sep;

		    sep = buffer.indexOf(separator,sepstart+1);
		    double evolution  =
			Double.valueOf(buffer.substring(sepstart+1,sep)).doubleValue( );
		    Double evolutionObj = new Double(evolution);
		    parametersMap.put(parameterEvolution,evolutionObj);
		    sepstart = sep;

		    System.err.format("zmax=%f evolution=%f\n",zmax,evolution);
		}else{
		    // 1st line -- eventNumber, powerLawIndex, energyFlux
		    sepstart = 0;
		    sep = buffer.indexOf(separator,sepstart+1);
		    sepstart = sep;

		    sep = buffer.indexOf(separator,sepstart+1);
		    eventNumber =
			Integer.valueOf(buffer.substring(sepstart+1,sep)).intValue( );
		    sepstart = sep;

		    sep = buffer.indexOf(separator,sepstart+1);
		    double powerLaw =
			Double.valueOf(buffer.substring(sepstart+1,sep)).doubleValue( );
		    Double powerLawObj = new Double(powerLaw);
		    parametersMap.put(parameterPowerLawIndex,powerLawObj);
		    sepstart = sep;

		    sep = buffer.indexOf(separator,sepstart+1);
		    double eFlux  =
			Double.valueOf(buffer.substring(sepstart+1,sep)).doubleValue( );
		    Double eFluxObj = new Double(eFlux);
		    parametersMap.put(parameterEnergyFluxAtBase,eFluxObj);
		    sepstart = sep;

		    System.err.format("power law=%f Energy Flux=%e\n",powerLaw,eFlux);

		    // 2nd line -- eventNumber, energyBase, cut-off energy
		    buffer = d.readLine();  		
		    sepstart = 0;
		    sep = buffer.indexOf(separator,sepstart+1);
		    sepstart = sep;

		    sep = buffer.indexOf(separator,sepstart+1);
		    eventNumber =
			Integer.valueOf(buffer.substring(sepstart+1,sep)).intValue( );
		    sepstart = sep;

		    sep = buffer.indexOf(separator,sepstart+1);
		    double energyBase =
			Double.valueOf(buffer.substring(sepstart+1,sep)).doubleValue( );
		    Double energyBaseObj = new Double(energyBase);
		    parametersMap.put(parameterEnergyBase,energyBaseObj);
		    sepstart = sep;

		    sep = buffer.indexOf(separator,sepstart+1);
		    double cutoffEnergy  =
			Double.valueOf(buffer.substring(sepstart+1,sep)).doubleValue( );
		    Double cutoffEnergyObj = new Double(cutoffEnergy);
		    parametersMap.put(parameterCutoffEnergy,cutoffEnergyObj);
		    sepstart = sep;

		    System.err.format("base Energy=%e [GeV] cut-off Energy =%e [GeV]\n",energyBase,cutoffEnergy);
		}

		//2nd line eventNumber, cosZ, log(NPE)(or energy proxy), expected rate
		Map realDataMap = new LinkedHashMap();
		for(int i= 0; i<numberOfBins;i++){
		    buffer = d.readLine();  		
		    sepstart = 0;
		    sep = buffer.indexOf(separator,sepstart+1);
		    sepstart = sep;

		    sep = buffer.indexOf(separator,sepstart+1);
		    eventNumber =
			Integer.valueOf(buffer.substring(sepstart+1,sep)).intValue( );
		    sepstart = sep;

		    sep = buffer.indexOf(separator,sepstart+1);
		    double cosZ =
                    Double.valueOf(buffer.substring(sepstart+1,sep)).doubleValue( );
		    sepstart = sep;

		    sep = buffer.indexOf(separator,sepstart+1);
		    double energyProxy  =
			Double.valueOf(buffer.substring(sepstart+1,sep)).doubleValue( );
		    sepstart = sep;

		    sep = buffer.indexOf(separator,sepstart+1);
		    double expectedNumber  =
                    Double.valueOf(buffer.substring(sepstart+1,sep)).doubleValue( );
		    if(expectedNumber<=0.0) expectedNumber=1.0e-20;
		    Double expectedObject = new Double(expectedNumber);
		    sepstart = sep;

		    long observedNumber=0;
		    if(!noObservedEvent) observedNumber=getObservedNumber(cosZ, energyProxy);
		    Long observedNumberObj = new Long(observedNumber);

		    //if(i%100==0 || observedNumber>0) 
		    if(observedNumber>0) 
			System.err.format("cosZ=%f energyProxy=%f rate=%e obs rate=%d\n",
					  cosZ,energyProxy,expectedNumber,observedNumber);

		    Map binDataMap = new LinkedHashMap();
		    binDataMap.put(observedNumberObj,expectedObject);

		    Integer binNumberObj = new Integer(i);
		    realDataMap.put(binNumberObj,binDataMap);

		}

		PoissonBinnedLikelihoodCalculator cal = new PoissonBinnedLikelihoodCalculator();
		cal.copyPoissonBinnedDataMap(realDataMap);

		pblMap.put(cal,parametersMap);

            }catch (EOFException e){
                buffer = null;
                break;
            }
        }

	hasBeenProcessed = true;
    }


    public void outputCreatedMap(OutputStream out) throws IOException{
	if(!hasBeenProcessed){
	    System.err.println("you have to call process() first!");
	}else{
	    ObjectOutputStream objectOut = new ObjectOutputStream(out);

	    objectOut.writeObject(pblMap);
	    objectOut.flush();
	}
    }


    public void printSummaryOfCreatedMap(){
	if(hasBeenProcessed){
	    Iterator dataIterator = pblMap.entrySet().iterator();
	    while(dataIterator.hasNext()){
		Map.Entry entryData = (Map.Entry )(dataIterator.next());
		PoissonBinnedLikelihoodCalculator calInMap = 
		    (PoissonBinnedLikelihoodCalculator)(entryData.getKey());
		Map paraMap = (Map )(entryData.getValue());
		long observedNumberOfEvents = calInMap.getSumOfObservedValues();
		double expectedNumberOfEvents = calInMap.getSumOfExpectedValues();

		if(isGZK){
		    double zmax = ((Double)(paraMap.get(parameterRedshift))).doubleValue();
		    double m = ((Double )(paraMap.get(parameterEvolution))).doubleValue();

		    System.out.format("(m,zmax)= (%f %f) expectedEventRate=%e observedNumber=%d\n",
				      m,zmax,expectedNumberOfEvents,observedNumberOfEvents);
		}else{
		    double powerLaw = ((Double)(paraMap.get(parameterPowerLawIndex))).doubleValue();
		    double eFlux = ((Double )(paraMap.get(parameterEnergyFluxAtBase))).doubleValue();
		    double energyBase = ((Double )(paraMap.get(parameterEnergyBase))).doubleValue();
		    double cutoffEnergy = ((Double )(paraMap.get(parameterCutoffEnergy))).doubleValue();

		    System.out.format("(powerLawIndex,energyBase, energyFlux, cutoffEnergy)= (%f %e %e %e) expectedEventRate=%e observedNumber=%d\n",
				      powerLaw,energyBase,eFlux,cutoffEnergy,expectedNumberOfEvents,observedNumberOfEvents);
		}
	    }
	}
    }
		


    public static void main(String[] args) throws IOException{

	String outFileName = null;
	boolean isGZK = true;
	boolean useNPEbin = true;
	boolean rebinned = false;
	boolean noObservedEvent = false;
        if(args.length<2){
            System.out.println("Usage: PoissonBinnedLikelihoodMapBuilder NPEbin(yes 1 no 0 rebin 2) filename-to-output-BinnedLikelohoodMap (0 if power law) (0 if no observed events - debugging)");
            System.exit(0);
        }else {
	    int index = Integer.valueOf(args[0]).intValue();
	    if(index==1){
		System.err.println(" reading NPE-binned TH2D");
	    }else{
		System.err.println(" reading recoE-binned TH2D");
		useNPEbin=false;
		if(index==2){
		    rebinned = true;
		    System.err.println(" rebinned TH2D");
		}
	    }
            outFileName = args[1];
	    if(args.length>=3) isGZK = false;
	    if(args.length==4) noObservedEvent = true;
	}
 	PoissonBinnedLikelihoodMapBuilder builder = new PoissonBinnedLikelihoodMapBuilder(useNPEbin,rebinned);
	if(!isGZK){
	    System.err.println("input file is power law");
	    builder.isPowerLaw();
	    if(noObservedEvent){
		System.err.println("null events observed case");
		builder.noObservedEvent = noObservedEvent;
	    }
	}

	DataInputStream input = new DataInputStream(System.in); 
	builder.process(input);
	builder.printSummaryOfCreatedMap();

	FileOutputStream out = new FileOutputStream(outFileName); 
	builder.outputCreatedMap(out);
	out.close();

   }
}
