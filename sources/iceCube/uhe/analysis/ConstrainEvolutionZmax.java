package iceCube.uhe.analysis;

import numRecipes.*;
import iceCube.uhe.analysis.*;
import java.io.*;
import java.util.*;

/**
   Calculate the p-value of GZK Yoshida-Ishihara in the plane of m zmax.

   The Map to contain the PoissonBinnedLikelihoodCalculator of the model predicion
   with various m and zmax, which is generated by PoissonBinnedLikelihoodMapBuilder,
   is fed into the p-value caluclator factory,
   ModelParaPvalueCalculatorByPoissonBinnedLH class.

 */

public class ConstrainEvolutionZmax {

   public static void main(String[] args) throws IOException{

	String inSigFileName = null;
	String inBGFileName = null;
	String inNuisanceFileName = null;
	boolean quickScan = false;
	boolean includeNuisance = false;

        if(args.length<3){
            System.out.println("Usage: ConstrainEvolutionZmax quickScan(yes=1) filename-of-GZK-YI-BinnedLikelohoodMap filename-BG-Map (filename-E2-Map as nuisance)");
            System.exit(0);
        }else {
	    if(Integer.valueOf(args[0]).intValue()==1) quickScan = true;
            inSigFileName = args[1];
            inBGFileName = args[2];
	    if(args.length == 4){
		inNuisanceFileName = args[3];
		includeNuisance = true;
	    }
	}
	if(quickScan) System.err.println("Do quick scan to search for the model to maximize llh");


	//
	// Reading the Map of the GZK Yoshida Ishihara
	//
        FileInputStream in = new FileInputStream(inSigFileName);
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

	//
	// Reading the Map of the atmospheric backgrounds
	//
        in = new FileInputStream(inBGFileName);
	objectIn = new ObjectInputStream(in);

 	Map bgBinnedDataMap = null;
        try{
            bgBinnedDataMap = (Map )objectIn.readObject();
        }catch(ClassNotFoundException e){
            System.err.println("Caught ClassNotFoundException: " + 
                               e.getMessage( ));
            System.exit(0);
        }
	in.close();

	//
	// extract the PoissonBinnedLikelihoodCalculator of the atrmospheric
	// background from the map
	//

	PoissonBinnedLikelihoodCalculator calBG = 
	    PoissonBinnedLHExtracter.getPoisonBinnedLikelihoodCalculator(bgBinnedDataMap);

	//
	// Reading the Map of the nuisance E2 signal model
	//
	PoissonBinnedLikelihoodCalculator calNuisance = null;
	if(includeNuisance){
	    System.err.println("include Nuisance E2 model");
	    in = new FileInputStream(inNuisanceFileName);
	    objectIn = new ObjectInputStream(in);

	    Map nuisanceBinnedDataMap = null;
	    try{
		nuisanceBinnedDataMap = (Map )objectIn.readObject();
	    }catch(ClassNotFoundException e){
		System.err.println("Caught ClassNotFoundException: " + 
				   e.getMessage( ));
		System.exit(0);
	    }
	    in.close();

	    calNuisance = 
		PoissonBinnedLHExtracter.getPoisonBinnedLikelihoodCalculator(nuisanceBinnedDataMap);
	}

	//=========================================================
	// Main Loop
	//
	// loop over the model predictions of the poisson binned data 
	//
	//=========================================================
	Iterator dataIterator = modelBinnedDataMap.entrySet().iterator();
	while(dataIterator.hasNext()){
	    Map.Entry entryData = (Map.Entry )(dataIterator.next());

	    //
	    // The signal to be tested
	    //
	    PoissonBinnedLikelihoodCalculator calSignal = 
		(PoissonBinnedLikelihoodCalculator)(entryData.getKey());
	    Map paraMap = (Map )(entryData.getValue());
	    Double mInThisModelObj = 
		(Double)paraMap.get(PoissonBinnedLikelihoodMapBuilder.parameterEvolution);
	    double mInThisModel = mInThisModelObj.doubleValue();
	    Double zmaxInThisModelObj = 
		(Double)paraMap.get(PoissonBinnedLikelihoodMapBuilder.parameterRedshift);
	    double zmaxInThisModel = zmaxInThisModelObj.doubleValue();


	    //
	    // put the model and bg predictions to the calculator factory
	    //
	    ModelParaPvalueCalculatorByPoissonBinnedLH pvalueCalculator = null;
	    if(!includeNuisance){
		pvalueCalculator = 
		    new ModelParaPvalueCalculatorByPoissonBinnedLH(calBG,calSignal,modelBinnedDataMap);
	    }else{
		pvalueCalculator = 
		    new ModelParaPvalueCalculatorByPoissonBinnedLH(calBG,calSignal,calNuisance,modelBinnedDataMap);
	    }
	    if(quickScan) pvalueCalculator.useTheQuickSearchToMaximizeLLH();

	    //
	    // The likelihood ratio calculation
	    //
	    boolean runReplicaExperiment = false;
	    double llhRatioObserved = pvalueCalculator.getLLHRatio(runReplicaExperiment);

	    //
	    // Extract the model to mazimize the llh
	    //
	    PoissonBinnedLikelihoodCalculator modelToMaxLH =  pvalueCalculator.getModelToMaximizeLH();
	    Map paraToMaxLH =  (Map )(modelBinnedDataMap.get(modelToMaxLH));
	    Double mInTheMostProbableModelObj = 
		(Double)paraToMaxLH.get(PoissonBinnedLikelihoodMapBuilder.parameterEvolution);
	    double mInTheMostProbableModel = mInTheMostProbableModelObj.doubleValue();
	    Double zmaxInTheMostProbableModelObj = 
		(Double)paraToMaxLH.get(PoissonBinnedLikelihoodMapBuilder.parameterRedshift);
	    double zmaxInTheMostProbableModel = zmaxInTheMostProbableModelObj.doubleValue();

	    //
	    // output the results in the stderr (for debugging)
	    //
	    System.err.format("(m zmax)=(%f %f) llhRatioObserved = %e\n",mInThisModel,zmaxInThisModel,
			      llhRatioObserved);
	    System.err.format(" the model to max llh (m zmax)=(%f %f)\n",
			      mInTheMostProbableModel,zmaxInTheMostProbableModel);

	    double pValue = 0.0;
	    if(!includeNuisance) pValue = pvalueCalculator.getPvalue(llhRatioObserved);
	    else  pValue = pvalueCalculator.getPvalue(llhRatioObserved,1000);
	    System.out.format("%f %f %e\n",mInThisModel,zmaxInThisModel,pValue);
	    //double pValueByChi2 = SpecialFunctions.gammaQ(1.0,llhRatioObserved);
	    //System.out.format("%f %f %e chi\n",mInThisModel,zmaxInThisModel,pValueByChi2);
	}

   }

}