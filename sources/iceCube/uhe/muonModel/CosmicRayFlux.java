package iceCube.uhe.muonModel;

import java.io.*;
import java.util.*;

//import numRecipes.*;

/**
<pre>
UHE Cosmic Ray Flux

</pre>
*/


public class CosmicRayFlux extends ParticleFlux{

    private static final double ln10 = Math.log(10.0);
    private double energyBase = 4.0e8;   // 4x10^8 GeV where the 2nd knee emerges
    private double energyAnkle = 6.3e9;  // 6.3x10^9 GeV where the angle starts
    private double energyKnee = 4.0e6;   // 4x10^6 GeV where the 1st knee point
    private double powerIndexBeforeKnee = 2.66894; // power law index before knee
    private double powerIndexLow = 3.0;  // power law index below energyBase
    private double powerIndexHigh = 3.2; // power law index above energyBase
    private double powerIndexAnkle = 2.75; // power law index above the ankle energy
    private double fluxAtBase = 6.26e-24;     // Flux [/cm^2 sec sr GeV] at energyBase
    private double fluxAtAnkle = 9.23e-28;    // Flux [/cm^2 sec sr GeV] at ankle
    private double logELowestBound = 5.0;    // This spectrum should have E> 10^5 GeV
    private double logEHighestBound = 10.8;  // The upper bound of E (10^11 GeV)
    private double logEmaximumBound = 15.0;  // The solod upper bound of E 
    private double logFluxAtBase;
    private double logEnergyBase;
    private double logEnergyAnkle;
    private double logFluxAnkle;
    private double logEnergyKnee;
    private double logFluxKnee = -17.1715;    // common logarithm of Flux [/cm^2 sec sr GeV] at knee
    private boolean cutoffExists = true; // No assumption of the (GZK) cutoff

    /** Constructor: Calculate log(dJ/dE) @ E = energyBase
    */
    public CosmicRayFlux() {
	logFluxAtBase = Math.log(fluxAtBase)/ln10;
	logFluxAnkle = Math.log(fluxAtAnkle)/ln10;
	logEnergyKnee = Math.log(energyKnee)/ln10;
	logEnergyBase = Math.log(energyBase)/ln10;
	logEnergyAnkle = Math.log(energyAnkle)/ln10;
    }


    /** 
	Tell if the CR spectrum calculated in this object
        has a cutoff feature at the highest energy end
    */
    public boolean doesCutOffExists(){
	return cutoffExists;
    }


    /** 
	Sets on/off the cutoff feature at the highest energy end
    */
    public void setCutOffFeature(boolean cutoffExists){
	this.cutoffExists = cutoffExists;
    }

    /** 
	calculate the differential Energy Flux [GeV /cm^2 sec sr] 

	<pre>
	logEnergy [GeV]
	</pre>

    */
    public double getEFlux(double logEnergy){

	double energy = Math.pow(10.0,logEnergy);
	double EFlux = getDFDLogE(logEnergy)*energy/ln10;
	
	return(EFlux);
    }


    /** 
	calculate the differential Energy Flux [GeV /cm^2 sec sr] 

	<pre>
	logEnergy [GeV] cos(zenith angle)
	</pre>

	This method is indentical to getEFlux(double logEnergy)
        because the cosmic ray flux is isotropic.
    */
    public double getEFlux(double logEnergy, double cosTheta){
	return getEFlux(logEnergy);
    }


    /** 
	<pre>
	calculate the log differential Flux dF/dLogE [/cm^2 sec sr] 

	logEnergy [GeV]
	</pre>

    */
    public double getDFDLogE(double logEnergy){

	if(isValidEnergy(logEnergy)){ // Energy is withinth the valid range
	    double logFlux = 0.0;
	    if(logEnergy <= logEnergyKnee){
		logFlux = logFluxKnee-powerIndexBeforeKnee
		    *(logEnergy-logEnergyKnee)+logEnergy;
	    }else if(logEnergy <= logEnergyBase){
		logFlux = logFluxAtBase-powerIndexLow*(logEnergy-logEnergyBase)+
		    logEnergy;
	    }else if (logEnergy <= logEnergyAnkle){
		logFlux = logFluxAtBase-powerIndexHigh*(logEnergy-logEnergyBase)+
		    logEnergy;
	    }else{
		logFlux = logFluxAnkle-powerIndexAnkle*(logEnergy-logEnergyAnkle)+
		    logEnergy;
	    }
	    double flux = Math.pow(10.0,logFlux)*ln10;
	    return(flux);
	}else{
	    return 0.0;
	}

    }


    /** 
	<pre>
	calculate the log differential Flux dF/dLogE [/cm^2 sec sr] 

	logEnergy [GeV] cos(zenith angle)
	</pre>


	This method is indentical to getDFDLogE(double logEnergy)
        because the cosmic ray flux is isotropic.
    */
    public double getDFDLogE(double logEnergy, double cosTheta){
	return getDFDLogE(logEnergy);
    }



    /** 
	Check if the energy is in the range where this power law model is valid.
        It depends on whether the spectrum has a cutoff feature which you
        can set by calling the method void setCutOffFeature(boolean ).
    */
    public boolean isValidEnergy(double logEnergy){
	if(logELowestBound <= logEnergy && logEnergy <= logEmaximumBound){
	    if(cutoffExists &&  logEnergy > logEHighestBound) return false;
	    else return true;
	}else return false;
    }

}

