package iceCube.uhe.neutrinoModel;

import numRecipes.*;
import iceCube.uhe.particles.*;

import java.io.*;

/** Calculate the source neutrino flux by the analytical functions */

public class NeutrinoFluxFromSource extends NeutrinoFluxFunction {

    public static double optDepth = 0.1; // optical depth of source
    /** The reference energy for the optical depth [GeV]*/
    public static double energyRef = 1.0e6;  // 1.0e6 GeV = 1 PeV
    /** The power law index of the energy spectrum of the photon field */
    public static double gamma = 2.0; 
    /** The Lorentz factor of the jet */
    public static double lorentzBulkFactor = 1.0;
    /** The maximum neutrino energy [GeV] */
    public static double neutrinoEnergyMax = 1.0e8; // 10^7 GeV = 10 PeV
    /** The minumum enery of the target photon [GeV] */
    public static double targetPhotonEnergyMin = getTargetPhotonEnergyMinimum();
    /** The maximal enery of the target photon [GeV] */
    public static double targetPhotonEnergyMax = 1.0e5*targetPhotonEnergyMin;
    /** The default CR flux relative to that extrpolated from the GZK regime*/
    public static double relativeCRFlux = 1.0;
    /** The reference energy of cosmic ray spectrum at which the relative intensity to that extrapolated
    from the GZK regime is calculated */
    public static double crEnergyRef = 1.0e7; // 1.0e7 GeV = 10 PeV
    /** Whether the expenential dumping factor exp(tau) is applied in calculation of CR flux*/
    public static boolean applyExpDumingToCRFlux = false;
    /** maximal CR proton energy determined from the maximum neutrino energy [GeV] */
    public static double crEnergyMax = neutrinoEnergyMax/(piEnergyMinus*(1.0-rPi));

    private static final double epsilon = 1.0e-7;


    /** The constructor: set the opetical depth */
    public NeutrinoFluxFromSource(double optDepth){
	super();
	this.optDepth = optDepth;
    }

    /** The constructor: set the refeerence energy and opetical depth */
    public NeutrinoFluxFromSource(double energyRef,double optDepth){
	super();
	this.optDepth = optDepth;
	this.energyRef = energyRef;
    }

    public NeutrinoFluxFromSource(){
	super();
    }

    public void setOpticalDepth(double depth){
	optDepth = depth;
    }

    public void setEnergyReference(double energy){
	energyRef = energy;
    }

    public void setCREnergyReference(double energy){
	crEnergyRef = energy;
    }

    public double getCREnergyReference(){
	return crEnergyRef;
    }

    /** set the relative normalization relativeCRFlux so that the normalization of the cosmic ray spectrum
    gives the eFlux at the cosmic ray reference energy */
    public void setCRFluxAtReferenceEnergy(double zMax, double m, double alpha, double eFlux){
	relativeCRFlux = 1.0; // back to the original default value
	double crFluxFromGZK =  getCRDFDE(zMax, m, alpha, crEnergyRef)*crEnergyRef*crEnergyRef;
	relativeCRFlux = eFlux/crFluxFromGZK;
    }

    /** set the relative normalization relativeCRFlux so that the normalization of the cosmic ray spectrum
    gives the eFlux at the cosmic ray reference energy. The case of the constant evolution beyond zConst. */
    public void setCRFluxAtReferenceEnergy(double zMax, double zConst, double m, double alpha, double eFlux){
	relativeCRFlux = 1.0; // back to the original default value
	double crFluxFromGZK =  getCRDFDE(zMax, zConst,m, alpha, crEnergyRef)*crEnergyRef*crEnergyRef;
	relativeCRFlux = eFlux/crFluxFromGZK;
    }


    public void returnToTheCRFluxNormalizationEstimatedFromGZK(){
	relativeCRFlux = 1.0; // back to the original default value
    }

    public void applyDumpingFactorToCRFlux(){
	applyExpDumingToCRFlux = true;
    }
    public void doNOTapplyDumpingFactorToCRFlux(){
	applyExpDumingToCRFlux = false;
    }

    /** set the power law index gamma. The radiation spectrum follows E<sup>-gamma</sup> */
    public void setPowerLawIndexOfRadiation(double gamma){
	this.gamma = gamma;
    }

    /** calculate the target photon energy corresponsing to a given neutrinoEnergy [GeV] */
    public static double getTargetPhotonEnergy(double neutrinoEnergy, double piEnergyRatio){
	double energy = (sRes-Mp*Mp)*lorentzBulkFactor*lorentzBulkFactor/(4.0*neutrinoEnergy)*piEnergyRatio*(1.0-rPi);
	return energy;
    }

    /** calculate the minimum target photon energy for the setted neutrinoEnergyMax */
    public static double getTargetPhotonEnergyMinimum(){
	double energy = getTargetPhotonEnergy(neutrinoEnergyMax,piEnergyMinus);
	return energy;
    }

    /** set the maximum energy of generated neutrinos [GeV] */
    public static void setMaximumNeutrinoEnergy(double energy){
	if(energy>energyRef){
	    neutrinoEnergyMax = energy;
	    targetPhotonEnergyMin = getTargetPhotonEnergyMinimum();
	    targetPhotonEnergyMax = 1.0e5*targetPhotonEnergyMin;
	    crEnergyMax = energy/(piEnergyMinus*(1.0-rPi));
	}else{
	    System.err.println(" You must set the maximal energy larger than the reference energy " +
			      energyRef + " [GeV]");
	    System.err.println(" Otherwise call the method setEnergyReference to reset the reference energy");
	}
    }

    /** return the manimum energy of the target photons [GeV].  
	This energy is reponsible for maximal energy of neutrinos emitted by p-gamma collision.*/
    public static double getMinumumTargetPhotonEnergy(){
	return targetPhotonEnergyMin;
    }

    /** return the maximum energy of the target photons [GeV].  */
    public static double getMaximumTargetPhotonEnergy(){
	return targetPhotonEnergyMax;
    }

    /** set the minimum energy of target photons [GeV]. You can specify this value
	to determine the maximum neutrino energy. The other option is
	to set the maximum neutrino energy by calling the method of setMaximumNeutrinoEnergy(double energy),
	which then automatically determines the minimum energy of target photons.
     */
    public static void setMinimumTargetPhotonEnergy(double energy){
	targetPhotonEnergyMin = energy;
	neutrinoEnergyMax = 
	    (sRes-Mp*Mp)*lorentzBulkFactor*lorentzBulkFactor/(4.0*targetPhotonEnergyMin)*piEnergyMinus*(1.0-rPi);
    }

    /** set the maximal energy of target photons [GeV]*/
    public static void setMaximumTargetPhotonEnergy(double energy){
	targetPhotonEnergyMax = energy;
    }

    public double getFunction(int functionIndex, double[] parameters, 
			      double x){
	if(functionIndex==0){
	    double alpha = parameters[0];
	    double zMax = parameters[1];
	    double m = parameters[2];
	    return getDFDE(zMax,m,alpha,x);
	}else{
	    double alpha = parameters[0];
	    double zMax = parameters[1];
	    double zConst = parameters[2];
	    double m = parameters[3];
	    return getDFDE(zMax,zConst,m,alpha,x);
	}

    }

    /** return the differential flux dF/dE [/cm^2 sec sr GeV] for the emission density per co-moving volume (1+z)^m upto zMax.*/
    public double getDFDE(double zMax, double m, double alpha, double neutrinoEnergy){
	double gzkEnergy = 1.0e11;   // 10^20 eV = 10^11 GeV
	//double gzkEnergy = gzkEnergyThreshold;
	double crossSectionTerm = 1.0/(piEnergyPlus-piEnergyMinus);
	double neutrinoIntensityTerm = 3.0/(1.0-rPi);
	double zTerm = getRedShiftTermSource(zMax,m,alpha,neutrinoEnergy);

	double flux = (alpha-1.0)*(alpha-1.0)/((alpha+1.0-gamma)*(alpha+1.0-gamma))*Math.pow(gzkEnergy,alpha-1.0)/Math.pow(energyRef, gamma-1.0)*getGZKCRFlux(gzkEnergy)*relativeCRFlux/(gzkRatio)*crossSectionTerm*neutrinoIntensityTerm*(1.0-Math.exp(-optDepth))*zTerm;

	return flux;

    }

    /** return the differential flux dF/dE [/cm^2 sec sr GeV] for the emission density per co-moving volume 
	(1+z)^m upto zConst and no evolution above upto zMax.*/
    public double getDFDE(double zMax, double zConst, double m, double alpha, double neutrinoEnergy){
	if(zMax<= zConst) return getDFDE(zMax,m,alpha,neutrinoEnergy);

	double fluxUptoZConst = getDFDE(zConst,m,alpha,neutrinoEnergy);
	if(fluxUptoZConst<=0.0) return 0.0;
	double zTerm = getRedShiftTermSource(zConst,m,alpha,neutrinoEnergy);
	double additionalzTerm = Math.pow((1+zConst),m)*getRedShiftTermSource(zMax,zConst,0.0,alpha,neutrinoEnergy);
	//System.out.format("zTerm=%e additionalzTerm=%e\n",zTerm,additionalzTerm);
	double flux = fluxUptoZConst + (additionalzTerm/zTerm)*fluxUptoZConst;

	return flux;
    }


    private static double getRedShiftTermSource(double zMax, double zMin, double m, double alpha, double neutrinoEnergy){
	double aa = m+1.0-alpha+gamma;
	double massDensityTerm = Math.pow(omegaM,-(aa-1.0)/3.0);
	if(aa==2.5) aa = epsilon+2.5;
	double powerIndexTerm = 2.0/(2.0*aa-5.0);
	double distancePowerIndex = 2.0*(aa/3.0-5.0/6.0);

	double zUp = getZSourceLimit(zMax,zMin,piEnergyMinus,targetPhotonEnergyMin,neutrinoEnergy);
	double zDown = getZSourceLimit(zMax,zMin,piEnergyPlus,targetPhotonEnergyMax,neutrinoEnergy);
	double zDownPrime = getZSourceLimit(zMax,zMin,piEnergyMinus,targetPhotonEnergyMax,neutrinoEnergy);
	//System.err.format(" zUp=%e zDown=%e zDownPrim=%e neutrinoEnergy=%e\n",zUp,zDown,zDownPrime,neutrinoEnergy);

	double distanceUp = getDistance(zUp);
	double distanceDown = getDistance(zDown);
	double distanceDownPrime = getDistance(zDownPrime);
	double distanceMin = getDistance(zMin);
	//System.err.format(" distanceUp=%e distanceDown=%e distanceDownPrim=%e\n",distanceUp,distanceDown,distanceDownPrime);

	double termA = powerIndexTerm*massDensityTerm*(Math.pow(distanceUp,distancePowerIndex)-Math.pow(distanceDown,distancePowerIndex));
	double termB = powerIndexTerm*massDensityTerm*(Math.pow(distanceDown,distancePowerIndex)-Math.pow(distanceDownPrime,distancePowerIndex));

	double energyTermA = Math.pow(neutrinoEnergy/(piEnergyPlus*(1.0-rPi)),-(alpha+1.0-gamma));
	double energyTermB = Math.pow(neutrinoEnergy/(piEnergyMinus*(1.0-rPi)),-(alpha+1.0-gamma));

	double energyPlusInLogTerm = getTargetPhotonEnergy(neutrinoEnergy,piEnergyPlus)/(targetPhotonEnergyMax*(1.0+zDown));
	double energyMinusInLogTerm = getTargetPhotonEnergy(neutrinoEnergy,piEnergyMinus)/(targetPhotonEnergyMax*(1.0+zDownPrime));

	double distancePowerIndexC = 2.0*((m+1.0)/3.0-0.5);
	double term = termA*energyTermA - termB*energyTermB;
	if(Math.abs(distancePowerIndexC)>1.0e-3 && zUp>0.0){
	    double dd = 2.0/(3.0*distancePowerIndexC);
	    double protonEnergyMin = (sRes-Mp*Mp)*lorentzBulkFactor*lorentzBulkFactor/(4.0*targetPhotonEnergyMax);
	    double distanceDownTermC = Math.pow(distanceDown,distancePowerIndexC);
	    double distanceDownPrimeTermC = Math.pow(distanceDownPrime,distancePowerIndexC);

	    double termC = Math.pow(distanceDown,distancePowerIndexC)+(alpha+1.0-gamma)*Math.pow(distanceDown,distancePowerIndexC)*(Math.log(energyPlusInLogTerm)+dd)-Math.pow(distanceDownPrime,distancePowerIndexC)-(alpha+1.0-gamma)*Math.pow(distanceDownPrime,distancePowerIndexC)*(Math.log(energyMinusInLogTerm)+dd)-(alpha+1.0-gamma)*Math.pow(distanceMin,distancePowerIndexC)*logFactor;
	    term +=  dd*Math.pow(omegaM,-(m+1.0)/3.0)*Math.pow(protonEnergyMin,-(alpha+1.0-gamma))*termC;
	    double dlog = Math.log(energyPlusInLogTerm)-Math.log(energyMinusInLogTerm);
	}

	return term;
    }

    private static double getRedShiftTermSource(double zMax, double m, double alpha, double neutrinoEnergy){
	return getRedShiftTermSource(zMax, 0.0, m, alpha, neutrinoEnergy);
    }


    private static double getZSourceLimit(double zMax, double zMin, double piEnergyRatio, double targetPhotonEnergy, double neutrinoEnergy){
	double xLimit = (sRes-Mp*Mp)*lorentzBulkFactor*lorentzBulkFactor/(4.0*targetPhotonEnergy*neutrinoEnergy)*piEnergyRatio*(1.0-rPi);
	double zLimit = xLimit-1;
	if(zLimit>zMax) zLimit = zMax;
	if(zLimit<zMin) zLimit = zMin;
	return zLimit;
    }	

    public double getCRDFDE(double zMax, double m, double alpha, double cosmicRayEnergy){
	double gzkEnergy = 1.0e11;   // 10^20 eV = 10^11 GeV
	double distanceMax = getDistance(zMax);
	double distancePowerIndexC = 2.0*((m-alpha)/3.0-1.0/6.0);
	double redShiftTerm = 2.0/(2.0*(m-alpha)-1.0)*Math.pow(omegaM,-(m-alpha+1.0)/3.0)*(Math.pow(distanceMax,distancePowerIndexC)-1.0);
	double flux = (alpha-1.0)*(alpha-1.0)*getGZKCRFlux(gzkEnergy)*relativeCRFlux/(gzkRatio)*Math.pow(cosmicRayEnergy,-alpha)*Math.pow(gzkEnergy,alpha-1.0)*redShiftTerm;

	double dumpFactor = 1.0;
	if(applyExpDumingToCRFlux){
	    if(cosmicRayEnergy< crEnergyMax) dumpFactor = Math.exp(-optDepth);
	}
	return flux*dumpFactor;
    }

    public double getCRDFDE(double zMax, double zConst, double m, double alpha, double cosmicRayEnergy){
	double gzkEnergy = 1.0e11;   // 10^20 eV = 10^11 GeV
	if(zMax<= zConst) return getCRDFDE(zMax,m,alpha,cosmicRayEnergy);
	double fluxUptoZConst = getCRDFDE(zConst,m,alpha,cosmicRayEnergy);
	if(fluxUptoZConst<=0.0) return 0.0;

	double distanceMax = getDistance(zMax);
	double distanceConst = getDistance(zConst);
	double distancePowerIndexC = 2.0*((-alpha)/3.0-1.0/6.0);
	double redShiftTerm = 2.0/(3.0*distancePowerIndexC)*Math.pow(omegaM,-(alpha+1.0)/3.0)*(Math.pow(distanceMax,distancePowerIndexC)-Math.pow(distanceConst,distancePowerIndexC))*Math.pow((1.0+zConst),m);
	double flux = (alpha-1.0)*(alpha-1.0)*getGZKCRFlux(gzkEnergy)*relativeCRFlux/(gzkRatio)*Math.pow(cosmicRayEnergy,-alpha)*Math.pow(gzkEnergy,alpha-1.0)*redShiftTerm + fluxUptoZConst;

	double dumpFactor = 1.0;
	if(applyExpDumingToCRFlux){
	    if(cosmicRayEnergy< crEnergyMax) dumpFactor = Math.exp(-optDepth);
	}
	return flux*dumpFactor;
    }

    public static void main(String[] args){
	if(args.length < 5) {
	    System.out.println("Usage: NeutrinoFluxFromSource zMax powerIndex m opticalDepth neutrinoEnergy (zConst)");
	    System.exit(0);
	}
	double zMax = 1.0;
	zMax = Double.valueOf(args[0]).doubleValue();
	double alpha = 3.0;
	alpha = Double.valueOf(args[1]).doubleValue();
	double m = 2.0;
	m = Double.valueOf(args[2]).doubleValue();
	double opticalDepth = 0.1;
	opticalDepth = Double.valueOf(args[3]).doubleValue();
	double neutrinoEnergy = 1.0e8;
	neutrinoEnergy = Double.valueOf(args[4]).doubleValue();
	boolean assumeConstEvolution = false;
	double zConst = 1.0;
	if(args.length==6){
	    zConst = Double.valueOf(args[5]).doubleValue();
	    assumeConstEvolution = true;
	    System.err.format(" constant evolution above %f\n",zConst);
	}

	NeutrinoFluxFromSource neutrinoFlux = new NeutrinoFluxFromSource( );
	neutrinoFlux.setOpticalDepth(opticalDepth);
	double[] parameters = new double[5];
	parameters[0] = alpha;
	parameters[1] = zMax;
	if(!assumeConstEvolution) parameters[2] = m;
	else{
	    parameters[2]= zConst;
	    parameters[3] = m;
	}


	neutrinoFlux.setPowerLawIndexOfRadiation(1.0);
	double crEnergyRef = 1.0e9; // [GeV]
	neutrinoFlux.setCREnergyReference(crEnergyRef);
	double eCRFlux = 2.0e-7; // [GeV/cm2 sec sr]
	if(!assumeConstEvolution) neutrinoFlux.setCRFluxAtReferenceEnergy(zMax, m, alpha, eCRFlux);
	else  neutrinoFlux.setCRFluxAtReferenceEnergy(zMax, zConst, m, alpha, eCRFlux);

	double energyFlux = 0.0;
	if(!assumeConstEvolution) energyFlux = neutrinoFlux.getDFDE(zMax,m, alpha,neutrinoEnergy)*neutrinoEnergy*neutrinoEnergy;
	else energyFlux = neutrinoFlux.getDFDE(zMax, zConst, m, alpha,neutrinoEnergy)*neutrinoEnergy*neutrinoEnergy;
	System.out.format(" E^2dF/dE (normalized at the 1EeV CR Flux) = %10.7e [GeV /cm^2 sec sr]\n",energyFlux);

	neutrinoFlux.returnToTheCRFluxNormalizationEstimatedFromGZK();
	if(!assumeConstEvolution) energyFlux = neutrinoFlux.getDFDE(zMax,m, alpha,neutrinoEnergy)*neutrinoEnergy*neutrinoEnergy;
	else energyFlux = neutrinoFlux.getDFDE(zMax, zConst, m, alpha,neutrinoEnergy)*neutrinoEnergy*neutrinoEnergy;
	System.out.format(" E^2dF/dE (normalized by the GZK intensity) = %10.7e [GeV /cm^2 sec sr]\n",energyFlux);

	// integral flux
	//double maxNeutrinoEnergy = 1.0e11;   // 1.0e11 GeV
	//double integralFlux = 0.0;
	//if(!assumeConstEvolution) integralFlux = Integration.RombergIntegral(neutrinoFlux, 
	//		  0, parameters, neutrinoEnergy, maxNeutrinoEnergy);
	//else integralFlux = Integration.RombergIntegral(neutrinoFlux, 
	//		1, parameters, neutrinoEnergy, maxNeutrinoEnergy);
	//System.out.format(" F(>%7.3e GeV) = %10.7e [/cm^2 sec sr]\n",neutrinoEnergy,integralFlux);
    }

}
