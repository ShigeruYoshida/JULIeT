package iceCube.uhe.analysis;

import iceCube.uhe.particles.*;
import iceCube.uhe.interactions.*;
import iceCube.uhe.points.*;
import iceCube.uhe.propagation.*;
import iceCube.uhe.neutrinoModel.*;
import iceCube.uhe.analysis.*;

import java.io.*;
import java.util.*;

/** Calculate the neutrino detection yield by running QuickProrpagationMatrixFlux
and Make a Map(log(Neutrino Energy) Map(cross section, yield)).
The map is written to ObjectOutputStream.
An enhancement factor
of the neutrino-nucleon CC interaction ("nuCCenhance") is
given by an interactive way in execution of this program. 


This class is used for constraining
the neutrino-nucleon cross section possibly increaced by
the extra-dimension mechanism.

Written by S.Yoshida 2008 July 13th
*/

public class MakeQuickPropagationYieldMap {

    static final double ln10 = Math.log(10.0);
    // Path to Intertaction Matrix file
    private static String[] intMtxPathname =
    {"iceCube/uhe/interactions/ice/","iceCube/uhe/interactions/rock/"};
    protected static String nuCCMtxObjectFile = "ENeutrinoChargeMtx";


    public static void main(String[] args) throws IOException{

	String fileName = null;
        if(args.length<1){
            System.out.println("MakeQuickPropagationYieldMap outputFileName");
            System.exit(0);
        }else{
            fileName = args[0];
	}

	// Maps to store the Yield Table
	//LinkedHashMap<Double,LinkedHashMap<Double,Double>> yieldMap = null;
	Map yieldMap = null;

	// Neutrino InteractionMatrix object 
	//    to calculate the STANDARD cross section
	InteractionsMatrix nuCCMtx = null; 

	boolean calTheory = false;
	double numberOfEvents = 0.0;
	double nuCCEnhancementFactor = 1.0;

	String objectFile = 
	    intMtxPathname[0].concat(nuCCMtxObjectFile); // ice
	InputStream inMtx = ClassLoader.getSystemResourceAsStream(objectFile);
	nuCCMtx = InteractionsMatrixInput.inputInteractionsMatrix(inMtx);
 

	QuickPropagationMatrixFlux iceFlux = new QuickPropagationMatrixFlux(nuCCEnhancementFactor);
	ParticlePoint s = new ParticlePoint(0.0,0.0,0);
	iceFlux.propagator = new NeutrinoQuickPropagator(s);

	double time =  365.0*24.0*3600.0*1.0; // 1year
	iceFlux.setObservationTime(time);

	// Interactive session to set the inIceParticle flavor and doublets
        DataInputStream input = new DataInputStream(System.in); 
        BufferedReader  d     = new BufferedReader(new InputStreamReader(input)); 

	// Set the inice particle spieces.
	iceFlux.addInIceParticle(0,0); // nuE
	iceFlux.addInIceParticle(1,0); // nuMu
	iceFlux.addInIceParticle(1,1); // Muon
	iceFlux.addInIceParticle(2,0); // nuTau
	iceFlux.addInIceParticle(2,1); // Tau

	// Set the nuCCEnhancementFactor range by an interactive way
        String buffer; 

	System.err.print("Log10(nuCCenhance) lower bound ->");
	buffer   = d.readLine(); 
	double logCCenhanceLow = Double.valueOf(buffer).doubleValue();

	System.err.print("Log10(nuCCenhance) higher bound ->");
	buffer   = d.readLine(); 
	double logCCenhanceHigh = Double.valueOf(buffer).doubleValue();

	System.err.print("Log10(nuCCenhance) bin size ->");
	buffer   = d.readLine(); 
	double binSize = Double.valueOf(buffer).doubleValue();


	// display nuCCenhancefactor range you will survey
	double epsilon = 1.0e-4; // round-off margin for binning
	double logCCenhanceFactor = logCCenhanceLow;
	while((logCCenhanceFactor+epsilon)<logCCenhanceHigh){
	    nuCCEnhancementFactor = Math.pow(10.0,logCCenhanceFactor);
	    System.err.println("nuCCEnhancementFactor=" + nuCCEnhancementFactor);
	    logCCenhanceFactor += binSize;
	}

        //yieldMap = new LinkedHashMap<Double,LinkedHashMap<Double,Double>>();
        yieldMap = new LinkedHashMap();

	// Now calculate the yield for a given enhancement factor
	// and write it to the standard output
	logCCenhanceFactor = logCCenhanceLow;
	while((logCCenhanceFactor+epsilon)<logCCenhanceHigh){
	    nuCCEnhancementFactor = Math.pow(10.0,logCCenhanceFactor);
	    System.err.println("Now Caluculating.. " + nuCCEnhancementFactor);
	    iceFlux.setNeutrinoCCEnhancement(nuCCEnhancementFactor);
	    iceFlux.calculateYield();

	    double logNeutrinoEnergy = Particle.getLogEnergyMinimum();
	    double logNeutrinoEnergyMax = Particle.getLogEnergyMinimum()
		+ Particle.getDeltaLogEnergy()*(double )(Particle.getDimensionOfLogEnergyMatrix());
	    while(logNeutrinoEnergy < logNeutrinoEnergyMax){ // E < 10^12 GeV

		// Calculate the relevant cross section
		int iLogE = (int )((logNeutrinoEnergy+epsilon-Particle.getLogEnergyMinimum())/
				   Particle.getDeltaLogEnergy());
		double sigmaCC = nuCCMtx.getSigmaMatrix(iLogE)*nuCCEnhancementFactor;

		// Add the values to the map
		logNeutrinoEnergy = Particle.getLogEnergyMinimum()+
		    Particle.getDeltaLogEnergy()*(double )iLogE;
		Double logEobject = new Double(logNeutrinoEnergy);
		Double sigmaObject = new Double(sigmaCC);
		Double yieldObject = new Double(iceFlux.getYield(logNeutrinoEnergy + epsilon));

		//LinkedHashMap<Double,Double> sigmaYieldMap = null;
		LinkedHashMap sigmaYieldMap = null;
		if(yieldMap.containsKey(logEobject)){
		    //sigmaYieldMap = yieldMap.get(logEobject);
		    sigmaYieldMap = (LinkedHashMap )yieldMap.get(logEobject);
		}else{
		    //sigmaYieldMap = new LinkedHashMap<Double,Double>();
		    sigmaYieldMap = new LinkedHashMap();
		}
		sigmaYieldMap.put(sigmaObject,yieldObject);
		yieldMap.put(logEobject,sigmaYieldMap);


		logNeutrinoEnergy += Particle.getDeltaLogEnergy();
	    }

	    logCCenhanceFactor += binSize;

	}

	FileOutputStream fos = new FileOutputStream(fileName);
	ObjectOutputStream objectOut = new ObjectOutputStream(fos);

	objectOut.writeObject(yieldMap);
	objectOut.flush();
	fos.close();

    }


}
