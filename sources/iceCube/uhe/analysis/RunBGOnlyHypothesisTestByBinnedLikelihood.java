package iceCube.uhe.analysis;

import numRecipes.*;
import iceCube.uhe.analysis.*;
import java.io.*;
import java.util.*;

/**

   Run the bg only hypothesis test based on the binned Poisson likelihood. 
   The PoissonBinnedLikelihoodCalculator :
   <pre>
   PoissonBinnedLikelihoodCalculator  calBG     (for the atmospheric background)
   PoissonBinnedLikelihoodCalculator  calSinal   (for the neutrino model such as GZK to be tested)
   PoissonBinnedLikelihoodCalculator  calNuisanceSingal (for the nuisance signal - like E^-2 against the GZK test) 
   </pre>
   are generated by the correponding binned data from the files and then executes the methods
   provied by ModelTestByPoissonBinnedLikelihoodFactory.

The likeliihood ratio type 5:
  <pre> 
       5   alternative hypothesis(bg + either signal or nuisance model with floated normalization, 
           select the one giving larger likelihood)/null hypothesis(bg ONLY)
  </pre>
  is calculated by pseudo-experiments with the method 
  void makeCollectionOfLogLikelihoodRatio(int type=5, int runTimes)
  provied by ModelTestByPoissonBinnedLikelihoodFactory.

*/
public class RunBGOnlyHypothesisTestByBinnedLikelihood {

    public static void main(String[] args) throws IOException{

        PoissonBinnedLikelihoodCalculator calBG = null;
        PoissonBinnedLikelihoodCalculator calSignal = null;
        PoissonBinnedLikelihoodCalculator calNuisanceSignal = null;
	ModelTestByPoissonBinnedLikelihoodFactory testFactory = null;
        boolean debugMode = false;
	boolean includeNuisance = false;
	boolean runReplicaExperiment = false;
	boolean DataIsMap = false;
	boolean nuisanceAlternative = false;

        String sigEventRateFileName = null;
        String nuisanceSigEventRateFileName = null;
        String bgEventRateFileName = null;

        if(args.length<3){
            System.out.println("Usage: RunBGOnlyHypothesisTestByBinnedLikelihood dataIsMap(yes 1 no 0) filename-to-read-BGdata filename-to-read-SIGdata filename-to-read-NuisanceData");
            System.exit(0);
        }else {
	    int index = Integer.valueOf(args[0]).intValue();
	    if(index == 1)  DataIsMap = true;
	    if(DataIsMap) System.err.println(" Binned Data in Map format is read out");
            bgEventRateFileName = args[1];
            sigEventRateFileName = args[2];
	    nuisanceSigEventRateFileName = args[3];
	}

        // BG : background-only hypothesis
	if(!DataIsMap){
	    calBG = new PoissonBinnedLikelihoodCalculator();
	    if(debugMode) calBG.debugFlag = true;
	    DataInputStream in = 
		new DataInputStream(ClassLoader.getSystemResourceAsStream(bgEventRateFileName));
	    calBG.fillData(in);
	    in.close();
	}else{
	    FileInputStream in = new FileInputStream(bgEventRateFileName);
	    ObjectInputStream objectIn = new ObjectInputStream(in);
	    Map bgBinnedDataMap = null;
	    try{
		bgBinnedDataMap = (Map )objectIn.readObject();
	    }catch(ClassNotFoundException e){
		System.err.println("Caught ClassNotFoundException: " + 
				   e.getMessage( ));
		System.exit(0);
	    }
	    in.close();
	    calBG = PoissonBinnedLHExtracter.getPoisonBinnedLikelihoodCalculator(bgBinnedDataMap);
	}


	// SIG : GZK model
	if(!DataIsMap){
	    calSignal = new PoissonBinnedLikelihoodCalculator();
	    if(debugMode) calSignal.debugFlag = true;
	    DataInputStream in = 
		new DataInputStream(ClassLoader.getSystemResourceAsStream(sigEventRateFileName));
	    calSignal.fillData(in);
	    in.close();
	}else{
	    FileInputStream in = new FileInputStream(sigEventRateFileName);
	    ObjectInputStream objectIn = new ObjectInputStream(in);
	    Map modelBinnedDataMap = null;
	    try{
		modelBinnedDataMap = (Map )objectIn.readObject();
	    }catch(ClassNotFoundException e){
		System.err.println("Caught ClassNotFoundException: " + 
				   e.getMessage( ));
		System.exit(0);
	    }
	    in.close();
	    calSignal = 
		PoissonBinnedLHExtracter.getPoisonBinnedLikelihoodCalculator(modelBinnedDataMap);
	}

        // SIG: nuisance model, probably E^-2
	if(!DataIsMap){
	    calNuisanceSignal = new PoissonBinnedLikelihoodCalculator();
	    if(debugMode) calNuisanceSignal.debugFlag = true;
	    DataInputStream in = new DataInputStream(ClassLoader.getSystemResourceAsStream(nuisanceSigEventRateFileName));
	    calNuisanceSignal.fillData(in);
	    in.close();
	}else{
	    FileInputStream in = new FileInputStream(nuisanceSigEventRateFileName);
	    ObjectInputStream objectIn = new ObjectInputStream(in);
	    Map nuisanceModelBinnedDataMap = null;
	    try{
		nuisanceModelBinnedDataMap = (Map )objectIn.readObject();
	    }catch(ClassNotFoundException e){
		System.err.println("Caught ClassNotFoundException: " + 
				   e.getMessage( ));
		System.exit(0);
	    }
	    in.close();
	    calNuisanceSignal = 
		PoissonBinnedLHExtracter.getPoisonBinnedLikelihoodCalculator(nuisanceModelBinnedDataMap);
	}


	// generate the ModelTestByPoissonBinnedLikelihoodFactory.
	testFactory = new ModelTestByPoissonBinnedLikelihoodFactory(calBG, calSignal, calNuisanceSignal);
	testFactory.doNotIncludeNuisanceSignal();

	// the lilelihood of null hypothesis (model + bg)
	boolean bgOnly = true;
	double llhNull = testFactory.buildLikelihoodForNullHypothesis(bgOnly,runReplicaExperiment);
	double llhSignalFloated = testFactory.buildLikelihoodForAlternativeHypothesis(false,runReplicaExperiment);
	double llhNuisanceSignalFloated = testFactory.buildLikelihoodForAlternativeHypothesis(false,true,runReplicaExperiment);
	double  maximizedFactor = testFactory.getModelNormalizationToMaximizeLikelihood();
	double llhH1 = llhSignalFloated;
	if(llhNuisanceSignalFloated<llhSignalFloated) llhH1 = llhNuisanceSignalFloated;
	double llhRatioObserved =  llhNull - llhH1;

	System.out.format("llh Null = %e\n",llhNull);
	System.out.format("llh signal Floated = %e\n",llhSignalFloated);
	System.out.format("llh nuisance signal Floated = %e\n",llhNuisanceSignalFloated);
	System.out.format("normalization to maximize llh = %e\n",maximizedFactor);


        DataInputStream input = new DataInputStream(System.in); 
        BufferedReader  d     = new BufferedReader(new InputStreamReader(input)); 
        String buffer; 
	System.err.print("File name to store lists of the llh ratio->"); 
	buffer   = d.readLine(); 
	String outputLLHRatioFileName = buffer;
	System.err.format(" OK will store the results in the file %s\n",outputLLHRatioFileName); 

	int runTimes = 500000;
	int likelihoodRatioType = 5;
	testFactory.makeCollectionOfLogLikelihoodRatio(likelihoodRatioType,runTimes);
	FileOutputStream out = new FileOutputStream(outputLLHRatioFileName); 
	ModelTestByPoissonBinnedLikelihoodFactory.outputLikelihhoRatioList(testFactory,out);
	out.close();

	ListIterator llhRatioListIterator = testFactory.getllhRatioIterator();
	double significance = 1.0;
	int times = 0;
	while(llhRatioListIterator.hasNext()){
	    times++;
	    Double llhRatioObj = (Double )(llhRatioListIterator.next());
	    double llhRatio = llhRatioObj.doubleValue();
	    double pValue = 1.0-((double )times)/((double )runTimes);
	    if(llhRatio< llhRatioObserved) significance = pValue;
	}
	System.out.format("p-value %e\n",significance);



    }

   
}