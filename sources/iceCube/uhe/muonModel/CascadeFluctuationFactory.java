package  iceCube.uhe.muonModel;

import iceCube.uhe.muonModel.*;
import numRecipes.*;
import java.io.*;

/**

This class describes the intrinsic flucuation of atmospheric muon energies
due to extensive airshower process. The Elbert model that is used
in AtmMuonBunfleFlux.java associates muon energies to their primary
cosmic ray energy by one on one relation. The CascadeFlucuationFactory
simulates their flucuation in hanndy manner using inputs from the Corsika
full simulation. AtmMuonBundleFlux can use this factory to calculate
the flux taking int account the fluctuation if the option is set.

Written originally for the IceCube EHE background analysis
by Shigeru Yoshida 2007-06-08

*/

public class CascadeFluctuationFactory {

    //double[][] muBundleSurfaceCRtable;
    //double[][] muBundleInIceCRtable;

    // Gaussian sigma of log((E_muon/E0)/ Bar(E_muon/E0)) in ice


    final double sigmaLogRInIce_6_LogE_7 = 0.450658;
    final double sigmaLogRInIce_7_LogE_8 = 0.372203;
    final double sigmaLogRInIce_8_LogE_9 = 0.300118;
    final double sigmaLogRInIce_9_LogE_10 = 0.240683;
    final double sigmaLogRInIce_10_LogE_11 = 0.195287;

    /** default constructor */
    public CascadeFluctuationFactory(){
    }

    /** Return the probability of Log10((E_muon/E0)/ Bar(E_muon/E0)), a relative depature
    from the mean of the energy ratio.

    <pre>
    double logR           : log10((E_muon/E0)/ Bar(E_muon/E0))
    double logPrimaryEnergy: log10(Primary Cosmic Ray Energy [GeV])
    boolean asInIce       : true for inIce muon energies, false for those at the earth surface
    </pre>
    */
    public double getProbability(double logR, double logPrimaryEnergy, boolean asInIce){
	double prob = 0.0;
	double sigmaLogR = sigmaLogRInIce_6_LogE_7;
	if(logPrimaryEnergy < 7.0){
	    sigmaLogR = sigmaLogRInIce_6_LogE_7;
	} else if (logPrimaryEnergy < 8.0){
	    sigmaLogR = sigmaLogRInIce_7_LogE_8;
	} else if (logPrimaryEnergy < 9.0){
	    sigmaLogR = sigmaLogRInIce_8_LogE_9;
	} else if (logPrimaryEnergy < 10.0){
	    sigmaLogR = sigmaLogRInIce_9_LogE_10;
	}else{
	    sigmaLogR = sigmaLogRInIce_10_LogE_11;
	}
	if(asInIce){ // against in-ice muon energy
	    prob = SpecialFunctions.gauss(0.0,sigmaLogR,logR);
	}else{ // at surface
	    prob = SpecialFunctions.gauss(0.0,sigmaLogR,logR);
	}
	return prob;
    }

    /** Returns minimum of Log((E_muon/E0)/ Bar(E_muon/E0)) */
    public double getLogMuOverCREnergyMin(double logPrimaryEnergy, boolean asInIce) {
	double sigmaLogR = sigmaLogRInIce_6_LogE_7;
	if(logPrimaryEnergy < 7.0){
	    sigmaLogR = sigmaLogRInIce_6_LogE_7;
	} else if (logPrimaryEnergy < 8.0){
	    sigmaLogR = sigmaLogRInIce_7_LogE_8;
	} else if (logPrimaryEnergy < 9.0){
	    sigmaLogR = sigmaLogRInIce_8_LogE_9;
	} else if (logPrimaryEnergy < 10.0){
	    sigmaLogR = sigmaLogRInIce_9_LogE_10;
	}else{
	    sigmaLogR = sigmaLogRInIce_10_LogE_11;
	}
	return -5.0*sigmaLogR;
    }

    /** Returns maximum of Log((E_muon/E0)/ Bar(E_muon/E0)) */
    public double getLogMuOverCREnergyMax(double logPrimaryEnergy, boolean asInIce) {
	double sigmaLogR = sigmaLogRInIce_6_LogE_7;
	if(logPrimaryEnergy < 7.0){
	    sigmaLogR = sigmaLogRInIce_6_LogE_7;
	} else if (logPrimaryEnergy < 8.0){
	    sigmaLogR = sigmaLogRInIce_7_LogE_8;
	} else if (logPrimaryEnergy < 9.0){
	    sigmaLogR = sigmaLogRInIce_8_LogE_9;
	} else if (logPrimaryEnergy < 10.0){
	    sigmaLogR = sigmaLogRInIce_9_LogE_10;
	}else{
	    sigmaLogR = sigmaLogRInIce_10_LogE_11;
	}
	return 5.0*sigmaLogR;
    }

    /** Returns Log((E_muon/E0)/ Bar(E_muon/E0)) for a given confidenceLevel.
	<pre>
	double confidenceLevel >0:   
               integral_(logR)^(infinity) probablity =0.5*(1.0-confidenceLevel)
	double confidenceLevel <0:   
               integral_(-infinity)^(logR) probablity =0.5*(1.0+confidenceLevel)

	double logPrimaryEnergy: log10(Primary Energy [GeV])
	boolean asInIce  : true for inIce muon energies, false for those at the earth surface
       return logR
       </pre>
     */
    public double getLogMuOverCREnergy(double confidenceLevel, 
				       double logPrimaryEnergy, boolean asInIce){
	double sigmaLogR = sigmaLogRInIce_6_LogE_7;
	if(logPrimaryEnergy < 7.0){
	    sigmaLogR = sigmaLogRInIce_6_LogE_7;
	} else if (logPrimaryEnergy < 8.0){
	    sigmaLogR = sigmaLogRInIce_7_LogE_8;
	} else if (logPrimaryEnergy < 9.0){
	    sigmaLogR = sigmaLogRInIce_8_LogE_9;
	} else if (logPrimaryEnergy < 10.0){
	    sigmaLogR = sigmaLogRInIce_9_LogE_10;
	}else{
	    sigmaLogR = sigmaLogRInIce_10_LogE_11;
	}

	double alpha = 0.5*(1.0-Math.abs(confidenceLevel));
	if(alpha>=0.4999999999) return 0.0;

	double z_alpha = 0.0;
	if(alpha<=1.0e-5) z_alpha = 4.2649;
	else if(alpha<=1.0e-4) z_alpha = 3.7190;
	else if(alpha<=2.0e-4) z_alpha = 3.49;
	else if(alpha<=3.0e-4) z_alpha = 3.39;
	else if(alpha<=4.0e-4) z_alpha = 3.33;
	else if(alpha<=5.0e-4) z_alpha = 3.2905;
	else if(alpha<=6.0e-4) z_alpha = 3.22;
	else if(alpha<=7.0e-4) z_alpha = 3.18;
	else if(alpha<=8.0e-4) z_alpha = 3.14;
	else if(alpha<=9.0e-4) z_alpha = 3.11;
	else if(alpha<=1.0e-3) z_alpha = 3.0902;
	else if(alpha<=1.1e-3) z_alpha = 3.05;
	else if(alpha<=1.5e-3) z_alpha = 2.96;
	else if(alpha<=5.0e-3) z_alpha = 2.5758;
	else if(alpha<=1.0e-2) z_alpha = 2.3263;
	else if(alpha<=2.5e-2) z_alpha = 1.96;
	else if(alpha<=5.0e-2) z_alpha = 1.6449;
	else if(alpha<=7.5e-2) z_alpha = 1.44;
	else if(alpha<=1.0e-1) z_alpha = 1.2816;
	else if(alpha<=1.25e-1) z_alpha = 1.15;
	else if(alpha<=1.5e-1) z_alpha = 1.035;
	else if(alpha<=1.6e-1) z_alpha = 1.00;
	else if(alpha<=1.7e-1) z_alpha = 0.95;
	else z_alpha = 0.67;

	if(confidenceLevel>=0.0) return z_alpha*sigmaLogR;
	else return -z_alpha*sigmaLogR;
    }


    /** sample and returns Log((E_muon/E0)/ Bar(E_muon/E0)) with MC method */
    public double sampleLogEnergyRatioFactor(RandomGenerator rand, 
					     double logPrimaryEnergy, boolean asInIce){
	double sigmaLogR = sigmaLogRInIce_6_LogE_7;
	if(logPrimaryEnergy < 7.0){
	    sigmaLogR = sigmaLogRInIce_6_LogE_7;
	} else if (logPrimaryEnergy < 8.0){
	    sigmaLogR = sigmaLogRInIce_7_LogE_8;
	} else if (logPrimaryEnergy < 9.0){
	    sigmaLogR = sigmaLogRInIce_8_LogE_9;
	} else if (logPrimaryEnergy < 10.0){
	    sigmaLogR = sigmaLogRInIce_9_LogE_10;
	}else{
	    sigmaLogR = sigmaLogRInIce_10_LogE_11;
	}
	if(asInIce){ // against in-ice muon energy
	    return rand.GetGaussianDouble(0.0,sigmaLogR);
	}else{
	    return rand.GetGaussianDouble(0.0,sigmaLogR);
	}
    }


}
