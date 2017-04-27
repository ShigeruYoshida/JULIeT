package iceCube.uhe.analysis;

import iceCube.uhe.particles.*;
import iceCube.uhe.points.*;
import iceCube.uhe.propagation.*;
import iceCube.uhe.neutrinoModel.*;
import iceCube.uhe.analysis.*;

import java.io.*;
import java.util.*;

/** Calculate the neutrino detection yield by running ProrpagationMatrixFlux
and write it to the standard output. An enhancement factor
of the neutrino-nucleon CC interaction ("nuCCenhance") is
given by an interactive way in execution of this program. 


This class is used for constraining
the neutrino-nucleon cross section possibly increaced by
the extra-dimension mechanism.

Written by S.Yoshida 2008 April 4th
*/

public class DumpQuickNCPropagationMatrixYield {

    static final double ln10 = Math.log(10.0);

    public static void main(String[] args) throws IOException{

	boolean calTheory = false;
	double numberOfEvents = 0.0;
	double nuNCEnhancementFactor = 1.0;


	QuickNCPropagationMatrixFlux iceFlux = new QuickNCPropagationMatrixFlux(nuNCEnhancementFactor);

	//double time =  365.0*24.0*3600.0*1.0; // 1year
	double time =  2.091744e7; //  242.1 days in [sec]
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

	// Set the nuNCEnhancementFactor range by an interactive way
        String buffer; 

	System.err.print("Log10(nuNCenhance) lower bound ->");
	buffer   = d.readLine(); 
	double logNCenhanceLow = Double.valueOf(buffer).doubleValue();

	System.err.print("Log10(nuNCenhance) higher bound ->");
	buffer   = d.readLine(); 
	double logNCenhanceHigh = Double.valueOf(buffer).doubleValue();

	System.err.print("Log10(nuNCenhance) bin size ->");
	buffer   = d.readLine(); 
	double binSize = Double.valueOf(buffer).doubleValue();


	// display nuCCenhancefactor range you will survey
	double epsilon = 1.0e-4; // round-off margin for binning
	double logNCenhanceFactor = logNCenhanceLow;
	while((logNCenhanceFactor+epsilon)<logNCenhanceHigh){
	    nuNCEnhancementFactor = Math.pow(10.0,logNCenhanceFactor);
	    System.err.println("nuNCEnhancementFactor=" + nuNCEnhancementFactor);
	    logNCenhanceFactor += binSize;
	}


	// Now calculate the yield for a given enhancement factor
	// and write it to the standard output
	logNCenhanceFactor = logNCenhanceLow;
	while((logNCenhanceFactor+epsilon)<logNCenhanceHigh){
	    nuNCEnhancementFactor = Math.pow(10.0,logNCenhanceFactor);
	    System.err.println("Now Caluculating.. " + nuNCEnhancementFactor);
	    iceFlux.setNeutrinoNCEnhancement(nuNCEnhancementFactor);
	    iceFlux.calculateYield();

	    double logNeutrinoEnergy = Particle.getLogEnergyMinimum();
	    double logNeutrinoEnergyMax = Particle.getLogEnergyMinimum()
		+ Particle.getDeltaLogEnergy()*(double )(Particle.getDimensionOfLogEnergyMatrix());
	    while(logNeutrinoEnergy < logNeutrinoEnergyMax){ // E < 10^12 GeV
		//System.out.println(" " + nuNCEnhancementFactor + " " + logNeutrinoEnergy +
		//		   " " + iceFlux.getYield(logNeutrinoEnergy + epsilon) +
		//		   " ");
		System.out.format(" %17.13f %5.2f %15.10e \n",nuNCEnhancementFactor,
				  logNeutrinoEnergy,
				  iceFlux.getYield(logNeutrinoEnergy + epsilon));
		logNeutrinoEnergy += Particle.getDeltaLogEnergy();
	    }

	    logNCenhanceFactor += binSize;

	}

    }


}
