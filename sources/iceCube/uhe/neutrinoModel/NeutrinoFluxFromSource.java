package iceCube.uhe.neutrinoModel;

import numRecipes.*;
import iceCube.uhe.particles.*;

import java.io.*;

/** Calculate the source neutrino flux by the analytical functions */

public class NeutrinoFluxFromSource extends NeutrinoFluxFunction {

    public static double optDepth = 0.1; // optical depth of source
    /** The reference energy for the optical depth [GeV]*/
    public static double energyRef = 1.0e7;  // 1.0e7 GeV = 10 PeV
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
    /** the threshold CR energy for the GZK photopion production processa */
    public static double crGZKthresholdEnergy = 6.0e10;  // 6x10^10 [GeV]
    /** the GZK CR energy in calculation of spectrum normalization */
    public static double crGZKnormalizationEnergy = 1.0e11;  // 10^11 [GeV]
    /** the threshold CR energy for the GZK BH processa */
    public static double crBHthresholdEnergy = 2.0e9;  // 2x 10^9 [GeV]

    /** parameters regaring synchrotron cooling **/
    public static double sigmaT = 6.66524e-25;  // Thompson scattering cross section [cm2]
    public static double magneticB = 1.0e3;  // Magnetic Field at source [gauss]
    public static double lorentzBoostFactor = 1.0;  // lorentz factor of jet as
                                                    // the magnetic field is defined in the shock co-moving frame
    public static boolean synchrotronCooling = true; // Take into account Sync.cooling in neutrino spectrum calculation

    /** whether calculation of UHECR spectrum is half-numerical or full-analytical */
    private static boolean calcUHECRfluxByNumerical = false;


    /** constants related to the ultra-high energy cosmic-ray propagation 
	involving BH pair creation */
    static double bhSphere = 2.0e3*Mpc;
    static double bhRatio = bhSphere*H0/c;


    private static final double epsilon = 1.0e-7;

    private static final double yr2sec = 365.0*24.0*3600.0; // [sec]
    private static final double erg = 6.242e2;  // [GeV]


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

    public static void calculateUHECRSpectrumByHalfNumericalMethod(){
	calcUHECRfluxByNumerical = true;	
    }

    public static void calculateUHECRSpectrumByFullAnalyticalMethod(){
	calcUHECRfluxByNumerical = false;	
    }


    /**
       set the relative normalization relativeCRFlux from a source photon luminosity L_gamma_source [erg/yr].
       <pre>
       double L_gamma_source  [erg/sec]     :  photon luminosity
       double cr_energy_min [GeV]           :  minimal cosmic energy emitted from a source
       double cr_energy_max [GeV]           :  maximal cosmic energy emitted from a source
       double alpha                         :  powerlaw index of CR spectrum
       double soucce_number_density [/Mpc3] :  source number density 
       double baryon_loading_factor         :  relative cr luminosity per gamma luminosity
       </pre>
     */
    public static void setCRFluxFromSourcePhotonLuminosity(double L_gamma_source, double cr_energy_min, double cr_energy_max,
							   double alpha,
							   double source_number_density, double baryon_loading_factor){
	relativeCRFlux = 1.0; // back to the original default value
	double x_min = cr_energy_min/energyRef;
	double x_max = cr_energy_max/energyRef;
	double kappa_cr = 0.0; // [/GeV sec sr]
	if(alpha!=2.0){
	    kappa_cr = (alpha-2.0)*L_gamma_source/(4.0*Math.PI*energyRef*energyRef)/
		(Math.pow(x_min,-alpha+2.0)-Math.pow(x_max,-alpha+2.0))*erg;
	}else{
	    kappa_cr = L_gamma_source/(4.0*Math.PI*energyRef*energyRef*Math.log(x_max/x_min))*erg;
	}
	double nokappa_cr = baryon_loading_factor*source_number_density/(Mpc*Mpc*Mpc)*kappa_cr; // [/GeV cm3 sec sr]
	relativeCRFlux = nokappa_cr/getCRSpectrumNormalization(crGZKnormalizationEnergy,alpha);
    }
	
    

    public static void setTheCRFluxNormalizationEstimatedFromGZK(double relativeFlux){
	relativeCRFlux = relativeFlux;
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


    
    /** assume Synchrotroton cooling of muons and pions */
    public static void assumeSynchrotronCooling(){
	synchrotronCooling = true;
    }

    /** assume NO Synchrotroton cooling of muons and pions */
    public static void assumeNoSynchrotronCooling(){
	synchrotronCooling = false;
    }

    /** Calculate the magnetic energy density [erg/cm3]
	<pre>
	B_gauss : Magnetic field [gauss]
	</pre>
    */
    public static double getMagneticEnergyDensity(double B_gauss){
	double rho_B = B_gauss*B_gauss/(8.0*Math.PI);
	return rho_B;
    }

    /** set the source magnetic field [gauss] */
    public static void setSourceMagneticField(double B_gauss){
	magneticB = B_gauss;  // [Gauss]
    }

    /** set the source magnetic field [gauss] from a given magnetic energy density [erg/cm3] */
    public static void setSourceMagneticFieldFromMagneticDensity(double rho_B){
	magneticB = Math.sqrt(8.0*Math.PI*rho_B);  // [Gauss]
    }

    /** set the lorentz boost factor for synchrotron cooling */
    public static void setLorentzBoostFactorForSynchrotronLoss(double lorentzFactor){
	lorentzBoostFactor = lorentzFactor;
    }

    
    /** get the muon critical energy for the synchrotron cooling [GeV] 
    */
    public static double getCriticalEnergyOfMuon(){
	double muonLifeTime =  2.19703e-6; // [sec]
	double Melectron = Particle.particleMasses[0][1]; // [GeV] flavor=0 doublet=1
	double energyScale = muonLifeTime*c*getMagneticEnergyDensity(magneticB)*sigmaT*erg; // [GeV]
	double inelasticityFactor = 2.0/3.0*rPi;
	double numFactor = 3.0/4.0;
	double criticalEnergySquared = numFactor*inelasticityFactor*inelasticityFactor*
	    Mmuon*Mmuon*Mmuon*Mmuon*Mmuon/(Melectron*Melectron*energyScale);
	double criticalEnergyMuon = Math.sqrt(criticalEnergySquared);
	return criticalEnergyMuon*lorentzBoostFactor;
    }

    /** get the base energy of muon synchrotron [GeV] */
    public static double getBaseEnergyOfMuon(){
	double Melectron = Particle.particleMasses[0][1]; // [GeV] flavor=0 doublet=1
	double inelasticityFactor = 2.0/3.0*rPi;
	double baseEnergy = Math.sqrt(3.0/4.0)*inelasticityFactor*Mmuon*Mmuon/Melectron;
	return baseEnergy;
    }

    /** get the effective muon synchrotron volume [cm3] */
    public static double getEffectiveMuonVolume(){
	double muonLifeTime =  2.19703e-6; // [sec]
	double v_sync = muonLifeTime*c*sigmaT;
	return v_sync;
    }


    /** get the pion critical energy for the synchrotron cooling [GeV] 
    */
    public static double getCriticalEnergyOfPion(){
	double pionLifeTime =  2.603e-8; // [sec]
	double Melectron = Particle.particleMasses[0][1]; // [GeV] flavor=0 doublet=1
	double energyScale = pionLifeTime*c*getMagneticEnergyDensity(magneticB)*sigmaT*erg; // [GeV]
	double inelasticityFactor = 1.0-rPi;
	double numFactor = 3.0/4.0;
	double criticalEnergySquared = numFactor*inelasticityFactor*inelasticityFactor*
	    Mpi*Mpi*Mpi*Mpi*Mpi/(Melectron*Melectron*energyScale);
	double criticalEnergyPion = Math.sqrt(criticalEnergySquared);
	return criticalEnergyPion*lorentzBoostFactor;
    }

    /** get the base energy of pion synchrotron [GeV] */
    public static double getBaseEnergyOfPion(){
	double Melectron = Particle.particleMasses[0][1]; // [GeV] flavor=0 doublet=1
	double inelasticityFactor = 1.0-rPi;
	double baseEnergy = Math.sqrt(3.0/4.0)*inelasticityFactor*Mpi*Mpi/Melectron;
	return baseEnergy;
    }

    /** get the effective pion synchrotron volume [cm3] */
    public static double getEffectivePionVolume(){
	double pionLifeTime =  2.603e-8; // [sec]
	double v_sync = pionLifeTime*c*sigmaT;
	return v_sync;
    }


    /** 
	Calculate the CR spectrum normalization kappa_cr*n_0.
	<pre>
	dN_CR/dE_p = kappa_cr*(E_p/E_0)**(-alpha):  E_0 is set by setEnergyReference(double energy)
                                                    (public static double energyRef)
        n_0 : the source number density at the present epoch z=0
        return kappa_cr*n0 for the pre-setted relative CR flux at energy = gzkEnergy [GeV]
	</pre>
     */
    public static double getCRSpectrumNormalization(double gzkEnergy, double alpha){
	double gzkCRFlux = getGZKCRFlux(gzkEnergy);
	double n0kappa_CR = (alpha-1.0)*(alpha-1.0)/gzkSphere*gzkCRFlux*Math.pow(gzkEnergy,alpha-1.0)*Math.pow(energyRef,-alpha)*relativeCRFlux;
	return n0kappa_CR;
    }

    /**
       Calculate the CR souce budget in unit of [erg/Mpc3 yr] with energies above crEnergyThreshold
     */
    public static double getCRSourceEnergyBudget(double crEnergyThreshold, double alpha){
	double n0kappa_CR = getCRSpectrumNormalization(crGZKnormalizationEnergy,alpha);
	double relativeEnergy = crEnergyThreshold/energyRef;
	double qCR = 0.0;
	if(alpha!=2.0){
	    qCR = 4.0*Math.PI*n0kappa_CR*energyRef*energyRef*Math.pow(relativeEnergy,-(alpha-2.0))/(alpha-2.0)*Mpc*Mpc*Mpc*yr2sec/erg;
	}else{ // differential energy budget
	    qCR = 4.0*Math.PI*n0kappa_CR*energyRef*energyRef*Mpc*Mpc*Mpc*yr2sec/erg;
	}
	return qCR;

    }

    public double getFunction(int functionIndex, double[] parameters, 
			      double x){
	if(functionIndex==0){
	    double alpha = parameters[0];
	    double zMax = parameters[1];
	    double m = parameters[2];
	    return getDFDE(zMax,m,alpha,x);
	}else if(functionIndex==1){
	    double alpha = parameters[0];
	    double zMax = parameters[1];
	    double zConst = parameters[2];
	    double m = parameters[3];
	    return getDFDE(zMax,zConst,m,alpha,x);
	}else if(functionIndex==2){
	    double alpha = parameters[0];
	    double zMax = parameters[1];
	    double m = parameters[2];
	    if(!calcUHECRfluxByNumerical){
		return getCRDFDE(zMax,m,alpha,x);
	    }else{
		return getCRDFDEHalfNumerical(zMax,m,alpha,x);
	    }
	}else if(functionIndex==3){
	    double alpha = parameters[0];
	    double zMax = parameters[1];
	    double zConst = parameters[2];
	    double m = parameters[3];
	    if(!calcUHECRfluxByNumerical){
		return getCRDFDE(zMax,zConst,m,alpha,x);
	    }else{
		return getCRDFDEHalfNumerical(zMax,zConst,m,alpha,x);
	    }
	}else  if(functionIndex==4){
            double alpha = parameters[0];
            double m = parameters[1];
	    return (1.0/getDistance(x))*Math.pow(1.0+x,-alpha+m);
	}else{
            double alpha = parameters[0];
            double m = parameters[1];
            double sphereRatio = parameters[2];
	    double zBound = parameters[3];
            double expLogTerm =  (alpha-1.0)/sphereRatio*2.0/(3.0*omegaM)*(getDistance(x)-getDistance(zBound));
	    double expTerm = Math.exp(-expLogTerm);
	    return (1.0/getDistance(x))*Math.pow(1.0+x,-alpha+m)*expTerm;
	}

    }

    /** return the differential flux dF/dE [/cm^2 sec sr GeV] for the emission density per co-moving volume (1+z)^m upto zMax.
	This is the method based on the retired normalization definition. Just for  debugging pruposes
    */
    public double getDFDEold(double zMax, double m, double alpha, double neutrinoEnergy){
	double crossSectionTerm = 1.0/(piEnergyPlus-piEnergyMinus);
	double neutrinoIntensityTerm = 3.0/(1.0-rPi);
	double zTerm = getRedShiftTermSource(zMax,m,alpha,neutrinoEnergy);

	double flux = (alpha-1.0)*(alpha-1.0)/((alpha+1.0-gamma)*(alpha+1.0-gamma))*Math.pow(crGZKnormalizationEnergy,alpha-1.0)/Math.pow(energyRef, gamma-1.0)*getGZKCRFlux(crGZKnormalizationEnergy)*relativeCRFlux/(gzkRatio)*crossSectionTerm*neutrinoIntensityTerm*(1.0-Math.exp(-optDepth))*zTerm;

	return flux;

    }

    /** return the differential flux dF/dE [/cm^2 sec sr GeV] for the emission density per co-moving volume (1+z)^m upto zMax.*/
    public double getDFDE(double zMax, double m, double alpha, double neutrinoEnergy){
	double crossSectionTerm = 1.0/(piEnergyPlus-piEnergyMinus);
	double neutrinoIntensityTerm = 3.0/(1.0-rPi);
	double zTerm = getRedShiftTermSource(zMax,m,alpha,neutrinoEnergy);

	double flux = getCRSpectrumNormalization(crGZKnormalizationEnergy,alpha)/((alpha+1.0-gamma)*(alpha+1.0-gamma))*Math.pow(energyRef, alpha-gamma+1.0)*c/H0*crossSectionTerm*neutrinoIntensityTerm*(1.0-Math.exp(-optDepth))*zTerm;

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


    /** return the differential flux dF/dE [/cm^2 sec sr GeV] from UHECR nuclei 
	for the emission density per co-moving volume (1+z)^m upto zMax.*/
    public double getDFDEFromNuclei(double zMax, double m, double alpha,
				    double neutrinoEnergy, double massNumber){
	double flux =  getDFDE(zMax, m, alpha, neutrinoEnergy);
	double reductionFactor = 1.0;
	if(massNumber>1.0){
	    reductionFactor = Math.pow(massNumber,2.0-alpha);
	}

	return flux*reductionFactor;
    }

    /** return the differential flux dF/dE [/cm^2 sec sr GeV] from UHECR nuclei
	for the emission density per co-moving volume 
	(1+z)^m upto zConst and no evolution above upto zMax.*/
    public double getDFDEFromNuclei(double zMax, double zConst, double m, double alpha,
				    double neutrinoEnergy, double massNumber){
	double flux = getDFDE(zMax, zConst, m, alpha, neutrinoEnergy);
	double reductionFactor = 1.0;
	if(massNumber>1.0){
	    reductionFactor = Math.pow(massNumber,2.0-alpha);
	}

	return flux*reductionFactor;
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
	double zMuon = getZMuonSynchrotronCooling(zMax, zMin, neutrinoEnergy);
	double zPion = getZPionSynchrotronCooling(zMax, zMin, neutrinoEnergy);

	double distanceUp = getDistance(zUp);
	double distanceDown = getDistance(zDown);
	double distanceDownPrime = getDistance(zDownPrime);
	double distanceMin = getDistance(zMin);
	//System.err.format(" distanceUp=%e distanceDown=%e distanceDownPrim=%e\n",distanceUp,distanceDown,distanceDownPrime);
	double distanceMax = getDistance(zMax);
	double distanceMuon = getDistance(zMuon);
	double distancePion = getDistance(zPion);

	double termA = 0.0;
	if(!synchrotronCooling){ // no synchrotron. The original formula
	    termA = powerIndexTerm*massDensityTerm*(Math.pow(distanceUp,distancePowerIndex)-Math.pow(distanceDown,distancePowerIndex));
	}else{
	    termA = powerIndexTerm*massDensityTerm*(Math.pow(distanceMuon,distancePowerIndex)-Math.pow(distanceDown,distancePowerIndex));
	    termA += powerIndexTerm*massDensityTerm*
		(Math.pow(distancePion,distancePowerIndex)-Math.pow(distanceMuon,distancePowerIndex))/3.0;
	    double aaB = m-1.0-alpha+gamma;
	    double massDensityTermB = Math.pow(omegaM,-(aaB-1.0)/3.0);
	    if(aaB==2.5) aaB = epsilon+2.5;
	    double powerIndexTermB = 2.0/(2.0*aaB-5.0);
	    double distancePowerIndexB = 2.0*(aaB/3.0-5.0/6.0);
	    double relativeEnergyMu = neutrinoEnergy/getCriticalEnergyOfMuon();
	    termA += powerIndexTermB*massDensityTermB/(relativeEnergyMu*relativeEnergyMu)*
	    (Math.pow(distanceMax,distancePowerIndexB)-Math.pow(distanceMuon,distancePowerIndexB))*2.0/3.0;
	    double relativeEnergyPi = neutrinoEnergy/getCriticalEnergyOfPion();
	    termA += powerIndexTermB*massDensityTermB/(relativeEnergyPi*relativeEnergyPi)*
	    	(Math.pow(distanceMax,distancePowerIndexB)-Math.pow(distancePion,distancePowerIndexB))/3.0;
	}
	
	double termB = powerIndexTerm*massDensityTerm*(Math.pow(distanceDown,distancePowerIndex)-Math.pow(distanceDownPrime,distancePowerIndex));

	double energyTermA = Math.pow(neutrinoEnergy/(piEnergyPlus*(1.0-rPi)),-(alpha+1.0-gamma));
	double energyTermB = Math.pow(neutrinoEnergy/(piEnergyMinus*(1.0-rPi)),-(alpha+1.0-gamma));

	double energyPlusInLogTerm = getTargetPhotonEnergy(neutrinoEnergy,piEnergyPlus)/(targetPhotonEnergyMax*(1.0+zDown));
	double energyMinusInLogTerm = getTargetPhotonEnergy(neutrinoEnergy,piEnergyMinus)/(targetPhotonEnergyMax*(1.0+zDownPrime));

	double distancePowerIndexC = 2.0*((m+1.0)/3.0-0.5);
	double term = termA*energyTermA - termB*energyTermB;
	double zIndex = zUp;
	if(synchrotronCooling) zIndex = zMuon;
	if(Math.abs(distancePowerIndexC)>1.0e-3 && zIndex>zMin){
	    double dd = 2.0/(3.0*distancePowerIndexC);
	    double protonEnergyMin = (sRes-Mp*Mp)*lorentzBulkFactor*lorentzBulkFactor/(4.0*targetPhotonEnergyMax);
	    double distanceDownTermC = Math.pow(distanceDown,distancePowerIndexC);
	    double distanceDownPrimeTermC = Math.pow(distanceDownPrime,distancePowerIndexC);

	    double termC = Math.pow(distanceDown,distancePowerIndexC)+(alpha+1.0-gamma)*Math.pow(distanceDown,distancePowerIndexC)*(Math.log(energyPlusInLogTerm)+dd)-Math.pow(distanceDownPrime,distancePowerIndexC)-(alpha+1.0-gamma)*Math.pow(distanceDownPrime,distancePowerIndexC)*(Math.log(energyMinusInLogTerm)+dd)-(alpha+1.0-gamma)*Math.pow(distanceMin,distancePowerIndexC)*logFactor;
	    term +=  dd*Math.pow(omegaM,-(m+1.0)/3.0)*Math.pow(protonEnergyMin,-(alpha+1.0-gamma))*termC;
	    double dlog = Math.log(energyPlusInLogTerm)-Math.log(energyMinusInLogTerm);
	    //	System.err.format("nuE %e zMuon %f zPion %f termA %e termB %e termC %e term %e\n",neutrinoEnergy,zMuon,zPion,termA, termB, termC,term);
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

    private static double getZMuonSynchrotronCooling(double zMax, double zMin, double neutrinoEnergy){
	double criticalEnergy = getCriticalEnergyOfMuon();
	double xLimit = criticalEnergy/neutrinoEnergy;
	double zLimit = xLimit-1;
	if(zLimit>zMax) zLimit = zMax;
	if(zLimit<zMin) zLimit = zMin;
	return zLimit;
    }

    private static double getZPionSynchrotronCooling(double zMax, double zMin, double neutrinoEnergy){
	double criticalEnergy = getCriticalEnergyOfPion();
	double xLimit = criticalEnergy/neutrinoEnergy;
	double zLimit = xLimit-1;
	if(zLimit>zMax) zLimit = zMax;
	if(zLimit<zMin) zLimit = zMin;
	return zLimit;
    }


    protected static double getZCRSourceGZK2BHLimit(double sphereRatio, double gzkThresholdZ,
						    double zMax, double zMin, double gzkEnergy, double cosmicRayEnergy){
	double xLimit = Math.sqrt(gzkEnergy/cosmicRayEnergy);
	double xLimitLow = gzkThresholdZ+1.0;
	double energyDistance = getDistance(gzkThresholdZ);

	while(xLimit>=xLimitLow){
	    double rightTerm = energyDistance+getLogTermOnZCRSourceLimit(sphereRatio,gzkEnergy, xLimit-1.0, cosmicRayEnergy);
	    double zDistance = getDistance(xLimit-1.0);
	    //	    System.err.format("sphereRatio(%f) z_s(%f) z_s_high((%f) zDistance=%e rightTerm=%e\n",
	    //	    		  sphereRatio,xLimit-1.0,xLimitLow-1.0,zDistance,rightTerm);
	    if(zDistance<=rightTerm) break;
	    xLimit -= 1.0e-3;
	}
	double zLimit = xLimit-1.0;
	if(zLimit>zMax) zLimit = zMax;
	if(zLimit<zMin) zLimit = zMin;
	return zLimit;
    }

    
    protected static double getZCRSourceGZK2BHLimit(double sphereRatio, double gzkThresholdZ, double zMax, double gzkEnergy, double cosmicRayEnergy){
	return getZCRSourceGZK2BHLimit(sphereRatio, gzkThresholdZ, zMax, 0.0, gzkEnergy, cosmicRayEnergy);
    }

    private static double getLogTermOnZCRSourceLimit(double sphereRatio, double gzkEnergy, double redshift, double cosmicRayEnergy){
	double energyTerm = gzkEnergy/(cosmicRayEnergy*(1.0+redshift)*(1.0+redshift));
	double term = sphereRatio*3.0*omegaM/2.0*Math.log(energyTerm);
	//System.err.format("   redshift=%f energyTerm=%e logTerm=%e\n", redshift,energyTerm,term);
	return term;
    }


    protected static double getZCRSourceLimit(double sphereRatio,double zMax, double zMin,double gzkEnergy, double cosmicRayEnergy){
	double xLimit = Math.sqrt(gzkEnergy/cosmicRayEnergy);
	double zLimit = xLimit-1.0;
	if(zLimit>zMax) zLimit = zMax;
	if(zLimit<zMin) zLimit = zMin;
	return zLimit;
    }

    protected static double getZCRSourceLimit(double sphereRatio, double zMax, double gzkEnergy, double cosmicRayEnergy){
	return getZCRSourceLimit(sphereRatio, zMax, 0.0, gzkEnergy, cosmicRayEnergy);
    }
    

    public double getCRDFDE(double zMax, double m, double alpha, double cosmicRayEnergy){
	double flux = 0.0;

	// flux from sources at distances that only matters redshift loss
	double zUp = getZCRSourceLimit(bhRatio, zMax, crBHthresholdEnergy, cosmicRayEnergy);
	if(zUp>0.0){
		double distanceMax = getDistance(zUp);
		double distancePowerIndexC = 2.0*((m-alpha)/3.0-1.0/6.0);
		double redShiftTerm = 2.0/(2.0*(m-alpha)-1.0)*Math.pow(omegaM,-(m-alpha+1.0)/3.0)*(Math.pow(distanceMax,distancePowerIndexC)-1.0);
		double relativeEnergy = cosmicRayEnergy/energyRef;
		flux = getCRSpectrumNormalization(crGZKnormalizationEnergy,alpha)*c/H0*
		    Math.pow(relativeEnergy,-alpha)*redShiftTerm;
	}

	// flux from sources at distances that matters redshift-evolved GZK BH  process
	if(cosmicRayEnergy<=crBHthresholdEnergy){
	    double zDown = getZCRSourceLimit(bhRatio, zMax, crBHthresholdEnergy, cosmicRayEnergy);
	    if(zDown<zMax){
		double fluxBHprocess = getCRFluxFromBHprocess(bhRatio, crBHthresholdEnergy,
								zDown, zDown, zMax, m, alpha, false);
		flux += fluxBHprocess;
	    }
	}else{  // flux from sources within GZK BH sphere
	    double thresholdZ = 0.0;
	    double zUpBH = getZCRSourceGZK2BHLimit(bhRatio, thresholdZ, zMax, crGZKthresholdEnergy, cosmicRayEnergy);
	    double fluxWithinBHSphere = getCRFluxWithinGZKsphere(bhSphere, bhRatio,
								  cosmicRayEnergy, 
								  zUpBH, m, alpha, true);
	    flux += fluxWithinBHSphere;
	}
	
	// flux from sources at distances that matters redshift-evolved GZK photopion process
	if(cosmicRayEnergy<=crGZKthresholdEnergy){
	    double thresholdZ = 0.0;
	    double zDown = getZCRSourceGZK2BHLimit(bhRatio, thresholdZ, zMax, crGZKthresholdEnergy, cosmicRayEnergy);
	    if(zDown<zMax){
		double fluxGZKprocess = getCRFluxFromGZKprocess(zDown, zDown, zMax, m, alpha);
		flux += fluxGZKprocess;
	    }
	}else{ // flux from sources within GZK sphere
	    double fluxWithinGZKSphere = getCRFluxWithinGZKsphere(gzkSphere, gzkRatio,
								  cosmicRayEnergy, 
								  zMax, m, alpha, false);

	    flux += fluxWithinGZKSphere;
	}

	double dumpFactor = 1.0;
	if(applyExpDumingToCRFlux){
	    if(cosmicRayEnergy< crEnergyMax) dumpFactor = Math.exp(-optDepth);
	}
	return flux*dumpFactor;
    }

    

    protected static double getCRFluxFromBHprocess(double sphereRatio, double gzkEnergy,
						    double thresholdZ,
						    double zDown, double zUp, double m, double alpha,
						    boolean include2ndOrder){
	double relativeRefEnergy = gzkEnergy/energyRef;
	double energyDistance = getDistance(thresholdZ);
	double distanceDown = getDistance(zDown);
	double distanceZmax =  getDistance(zUp);
	double expPowerDown = (alpha-1.0)/sphereRatio*2.0/(3.0*omegaM)*(distanceDown-energyDistance);
	double expPowerMax = (alpha-1.0)/sphereRatio*2.0/(3.0*omegaM)*(distanceZmax-energyDistance);
	double redShiftTermGZK = 1.0/(alpha-1.0)*
		(Math.pow(1.0+zDown,m-alpha-2.0)*Math.exp(-expPowerDown)-Math.pow(1.0+zUp,m-alpha-2.0)*Math.exp(-expPowerMax));
	if(include2ndOrder){
	    redShiftTermGZK += (m-alpha-2.0)/(alpha-1.0)*sphereRatio*Math.pow(1.0+zDown,m-alpha-5.0)*Math.exp(-expPowerDown)-
		(m-alpha-2.0)/(alpha-1.0)*sphereRatio*Math.pow(1.0+zUp,m-alpha-5.0)*Math.exp(-expPowerMax);
	}
	double fluxGZKprocess = 0.5*getCRSpectrumNormalization(crGZKnormalizationEnergy,alpha)*c/H0/energyDistance*
	    Math.pow(relativeRefEnergy,-alpha)*Math.pow(1.0+thresholdZ,(2.0*alpha+3.0))*redShiftTermGZK;
	//System.err.format(" flux from BH process (%e)\n", fluxGZKprocess);

	return fluxGZKprocess;
    }


    protected static double getCRFluxFromBHprocess(double sphereRatio, double gzkEnergy,
						    double thresholdZ,
						    double zDown, double zUp, double zConst, double m, double alpha,
						    boolean include2ndOrder){
	double zDownPrime = zDown;
	if(zDown < zConst) zDownPrime = zConst;

	double relativeRefEnergy = gzkEnergy/energyRef;
	double energyDistance = getDistance(thresholdZ);
	double distanceDown = getDistance(zDownPrime);
	double distanceZmax =  getDistance(zUp);
	double expPowerDown = (alpha-1.0)/sphereRatio*2.0/(3.0*omegaM)*(distanceDown-energyDistance);
	double expPowerMax = (alpha-1.0)/sphereRatio*2.0/(3.0*omegaM)*(distanceZmax-energyDistance);
	double redShiftTermGZK = 1.0/(alpha-1.0)*
		(Math.pow(1.0+zDownPrime,-alpha-2.0)*Math.exp(-expPowerDown)-Math.pow(1.0+zUp,-alpha-2.0)*Math.exp(-expPowerMax))*Math.pow(1.0+zConst,m);
	if(include2ndOrder){
	    redShiftTermGZK += (-alpha-2.0)/(alpha-1.0)*sphereRatio*Math.pow(1.0+zDownPrime,-alpha-5.0)*Math.exp(-expPowerDown)-
		(-alpha-2.0)/(alpha-1.0)*sphereRatio*Math.pow(1.0+zUp,-alpha-5.0)*Math.exp(-expPowerMax)*Math.pow(1.0+zConst,m);
	}
	double fluxGZKprocess = 0.5*getCRSpectrumNormalization(crGZKnormalizationEnergy,alpha)*c/H0/energyDistance*
	    Math.pow(relativeRefEnergy,-alpha)*Math.pow(1.0+thresholdZ,(2.0*alpha+3.0))*redShiftTermGZK;

	return fluxGZKprocess;
    }
    



    protected static double getCRFluxFromGZKprocess(double thresholdZ,
						    double zDown, double zUp, double m, double alpha){

	double relativeRefEnergy = crGZKthresholdEnergy/energyRef;
	double energyDistance = getDistance(thresholdZ);
	double distanceDown = getDistance(zDown);
	double distanceZmax =  getDistance(zUp);
	double expPowerDown = (alpha-1.0)/gzkRatio*2.0/(3.0*omegaM)*(distanceDown-energyDistance);
	double expPowerMax = (alpha-1.0)/gzkRatio*2.0/(3.0*omegaM)*(distanceZmax-energyDistance);
	double expFactor =  1.0/bhRatio*2.0/(3.0*omegaM)*(energyDistance-1.0);
	double redShiftTermGZK = 1.0/(alpha-1.0)*
		(Math.pow(1.0+zDown,m-alpha-2.0)*Math.exp(-expPowerDown)-Math.pow(1.0+zUp,m-alpha-2.0)*Math.exp(-expPowerMax));
	double distanceFactor = 2.0*energyDistance/Math.pow(1.0+thresholdZ,5.0)+1.0/(bhRatio*Math.pow(1.0+thresholdZ,2.0));
	double fluxGZKprocess = getCRSpectrumNormalization(crGZKnormalizationEnergy,alpha)*Math.exp(expFactor)*c/H0/distanceFactor*
	    Math.pow(relativeRefEnergy,-alpha)*Math.pow(1.0+thresholdZ,(2.0*alpha-2.0))*redShiftTermGZK;
	//System.err.format("norm(%e)  zterm (%e) distance factor(%e) flux from GZK process (%e)\n",
	//		  relativeCRFlux,redShiftTermGZK, distanceFactor,fluxGZKprocess);
	//System.err.format("flux from GZK process (%e)\n",fluxGZKprocess);

	return fluxGZKprocess;
    }


    protected static double getCRFluxFromGZKprocess(double thresholdZ,
						    double zDown, double zUp, double zConst, double m, double alpha){

	double relativeRefEnergy = crGZKthresholdEnergy/energyRef;
	double energyDistance = getDistance(thresholdZ);
	double distanceDown = getDistance(zDown);
	double distanceZmax =  getDistance(zUp);
	double expPowerDown = (alpha-1.0)/gzkRatio*2.0/(3.0*omegaM)*(distanceDown-energyDistance);
	double expPowerMax = (alpha-1.0)/gzkRatio*2.0/(3.0*omegaM)*(distanceZmax-energyDistance);
	double expFactor =  1.0/bhRatio*2.0/(3.0*omegaM)*(energyDistance-1.0);
	double redShiftTermGZK = 1.0/(alpha-1.0)*
	    (Math.pow(1.0+zDown,-alpha-2.0)*Math.exp(-expPowerDown)-Math.pow(1.0+zUp,-alpha-2.0)*Math.exp(-expPowerMax))*Math.pow(1.0+zConst,m);
	double distanceFactor = 2.0*energyDistance/Math.pow(1.0+thresholdZ,5.0)+1.0/(bhRatio*Math.pow(1.0+thresholdZ,2.0));
	double fluxGZKprocess = getCRSpectrumNormalization(crGZKnormalizationEnergy,alpha)*Math.exp(expFactor)*c/H0/distanceFactor*
	    Math.pow(relativeRefEnergy,-alpha)*Math.pow(1.0+thresholdZ,(2.0*alpha-2.0))*redShiftTermGZK;
	//System.err.format("zth(%e) zdown(%e) zconst(%e) zup(%e) zterm(%e) flux from GZK process (%e)\n",
	//		  thresholdZ, zDown,zConst,zUp,redShiftTermGZK,fluxGZKprocess);

	return fluxGZKprocess;
    }

    

    protected static double getCRFluxWithinGZKsphere(double sphereRadius, double sphereRatio,
						     double cosmicRayEnergy,
						     double zUp, double m, double alpha,
						     boolean include2ndOrder){
	double distanceZmax =  getDistance(zUp);
	double expPowerMax = (alpha-1.0)/sphereRatio*2.0/(3.0*omegaM)*(distanceZmax-1.0);
	double relativeEnergy = cosmicRayEnergy/energyRef;
	double redshiftTerm = (1.0-Math.pow(1+zUp,m-alpha-2.0)*Math.exp(-expPowerMax));
	if(include2ndOrder){
	    redshiftTerm += (m-alpha-2.0)/(alpha-1.0)*sphereRatio*(1.0-Math.pow(1.0+zUp,m-alpha-5.0)*Math.exp(-expPowerMax));
	}
	double fluxWithinGZKSphere = getCRSpectrumNormalization(crGZKnormalizationEnergy,alpha)*sphereRadius*
	    Math.pow(relativeEnergy,-alpha)*redshiftTerm/(alpha-1.0);
	//System.err.format(" flux within GZK sphere (%e)\n", fluxWithinGZKSphere);
	return fluxWithinGZKSphere;
    }    
    
    protected static double getCRFluxWithinGZKsphere(double sphereRadius, double sphereRatio,
						     double cosmicRayEnergy,
						     double zUp, double zConst, double m, double alpha,
						     boolean include2ndOrder){
	    
	double distanceZmax =  getDistance(zUp);
	double expPowerMax = (alpha-1.0)/sphereRatio*2.0/(3.0*omegaM)*(distanceZmax-1.0);
	double distanceConst = getDistance(zConst);
	double expPowerConst = (alpha-1.0)/sphereRatio*2.0/(3.0*omegaM)*(distanceConst-1.0);
	double relativeEnergy = cosmicRayEnergy/energyRef;
	double redshiftTerm =  (Math.pow(1.0+zConst,-alpha-2.0)*Math.exp(-expPowerConst)-Math.pow(1+zUp,-alpha-2.0)*Math.exp(-expPowerMax))*Math.pow(1.0+zConst,m);
	if(include2ndOrder){
	    redshiftTerm +=  (m-alpha-2.0)/(alpha-1.0)*sphereRatio*(Math.pow(1.0+zConst,-alpha-5.0)*Math.exp(-expPowerConst)-Math.pow(1+zUp,-alpha-5.0)*Math.exp(-expPowerMax))*Math.pow(1.0+zConst,m);
	}
	double fluxWithinGZKSphere = getCRSpectrumNormalization(crGZKnormalizationEnergy,alpha)*sphereRadius*
		Math.pow(relativeEnergy,-alpha)*redshiftTerm/(alpha-1.0);
	return fluxWithinGZKSphere;
    }

    
    
    public double getCRDFDE(double zMax, double zConst, double m, double alpha, double cosmicRayEnergy){
	if(zMax<= zConst) return getCRDFDE(zMax,m,alpha,cosmicRayEnergy);
	double fluxUptoZConst = getCRDFDE(zConst,m,alpha,cosmicRayEnergy);
	if(fluxUptoZConst<=0.0) return 0.0;

	double flux = 0.0;

	// flux from sources at distances that only matters redshift loss
	double zUp = getZCRSourceLimit(bhRatio, zMax, zConst, crBHthresholdEnergy, cosmicRayEnergy);
	if(zUp>zConst){
	    double distanceMax = getDistance(zUp);
	    double distanceConst = getDistance(zConst);
	    double distancePowerIndexC = 2.0*((-alpha)/3.0-1.0/6.0);
	    double redShiftTerm = 2.0/(3.0*distancePowerIndexC)*Math.pow(omegaM,(alpha-1.0)/3.0)*(Math.pow(distanceMax,distancePowerIndexC)-Math.pow(distanceConst,distancePowerIndexC))*Math.pow((1.0+zConst),m);
	    double relativeEnergy = cosmicRayEnergy/energyRef;
	    flux = getCRSpectrumNormalization(crGZKnormalizationEnergy,alpha)*c/H0*Math.pow(relativeEnergy,-alpha)*redShiftTerm;
	}
	
	// flux from sources at distances that matters redshift-evolved GZK BH process
	if(cosmicRayEnergy<=crBHthresholdEnergy){
	    double zDown = getZCRSourceLimit(bhRatio, zMax, crBHthresholdEnergy, cosmicRayEnergy);
	    if(zConst<zDown){
		double fluxBHprocess = getCRFluxFromBHprocess(bhRatio, crBHthresholdEnergy, zDown,
								zDown, zMax, zConst, m, alpha, false);
		flux += fluxBHprocess;
	    }
	}else{ // flux from sources within GZK BH sphere
	    double thresholdZ = getZCRSourceGZK2BHLimit(bhRatio, 0.0, zMax, crGZKthresholdEnergy, cosmicRayEnergy);
	    if(thresholdZ>=zConst){
		double fluxWithinBHSphere = getCRFluxWithinGZKsphere(bhSphere, bhRatio,
								     cosmicRayEnergy,
								     thresholdZ, zConst, m, alpha, true);
		flux += fluxWithinBHSphere;
	    }
	}

	// flux from sources at distances that matters redshift-evolved GZK photopion process
	if(cosmicRayEnergy<=crGZKthresholdEnergy){
	    double thresholdZ = getZCRSourceGZK2BHLimit(bhRatio, 0.0, zMax, crGZKthresholdEnergy, cosmicRayEnergy);
	    if(thresholdZ>=zConst){
		double fluxGZKprocess = getCRFluxFromGZKprocess(thresholdZ, thresholdZ, zMax, zConst, m, alpha);
		flux += fluxGZKprocess;
	    }
	}else{ // flux from sources within GZK sphere
	    double fluxWithinGZKSphere = getCRFluxWithinGZKsphere(gzkSphere, gzkRatio,
								  cosmicRayEnergy,
								  zMax, zConst, m, alpha, false);
	    flux += fluxWithinGZKSphere;
	}
	
	flux += fluxUptoZConst;

	double dumpFactor = 1.0;
	if(applyExpDumingToCRFlux){
	    if(cosmicRayEnergy< crEnergyMax) dumpFactor = Math.exp(-optDepth);
	}
	return flux*dumpFactor;
    }

    /** Calculate CR spectrum with numerical redshift integration instead of the analytical approximations 
     */
    public double getCRDFDEHalfNumerical(double zMax, double m, double alpha, double cosmicRayEnergy){
	double flux = 0.0;
        double[] parameters = new double[4];

	// flux from sources at distances that only matters redshift loss
	double zUp = getZCRSourceLimit(bhRatio, zMax, crBHthresholdEnergy, cosmicRayEnergy);
	if(zUp>0.0){
                parameters[0] = alpha;
                parameters[1] = m;
		double redShiftTerm = Integration.RombergIntegral(this, 4, parameters, 0.0, zUp);
		double relativeEnergy = cosmicRayEnergy/energyRef;
		flux = getCRSpectrumNormalization(crGZKnormalizationEnergy,alpha)*c/H0*
		    Math.pow(relativeEnergy,-alpha)*redShiftTerm;
		//System.err.format(" flux from redshift only (%e)\n", flux);
	}

	// flux from sources at distances that matters redshift-evolved GZK BH  process
	if(cosmicRayEnergy<=crBHthresholdEnergy){
	    double zDown = getZCRSourceLimit(bhRatio, zMax, crBHthresholdEnergy, cosmicRayEnergy);
	    if(zDown<zMax){
		double relativeRefEnergy = crBHthresholdEnergy/energyRef;
		double energyDistance = getDistance(zDown);
		parameters[0] = alpha;
		parameters[1] = m;
		parameters[2] = bhRatio;
		parameters[3] = zDown;
		double redshiftTerm = Integration.RombergIntegral(this,5,parameters,zDown,zMax);
		double fluxBHprocess = 0.5*getCRSpectrumNormalization(crGZKnormalizationEnergy,alpha)*c/H0/(bhRatio*energyDistance)*
		    Math.pow(relativeRefEnergy,-alpha)*Math.pow(1.0+zDown, (2.0*alpha+3.0))*redshiftTerm;
		//System.err.format(" flux from BH process (%e)\n", fluxBHprocess);

		flux += fluxBHprocess;
	    }
	}else{  // flux from sources within GZK BH sphere
	    double thresholdZ = 0.0;
	    double zUpBH = getZCRSourceGZK2BHLimit(bhRatio, thresholdZ, zMax, crGZKthresholdEnergy, cosmicRayEnergy);
	    double distanceZmax =  getDistance(zUpBH);
	    parameters[0] = alpha;
	    parameters[1] = m;
	    parameters[2] = bhRatio;
	    parameters[3] = 0.0;
	    double redshiftTerm = Integration.RombergIntegral(this,5,parameters,0.0,zUpBH);
	    double relativeEnergy = cosmicRayEnergy/energyRef;
	    double fluxWithinBHSphere = getCRSpectrumNormalization(crGZKnormalizationEnergy,alpha)*c/H0*
		Math.pow(relativeEnergy,-alpha)*redshiftTerm;
	    //System.err.format(" flux within BH process (%e)\n", fluxWithinBHSphere);
	    flux += fluxWithinBHSphere;
	}
	
	// flux from sources at distances that matters redshift-evolved GZK photopion process
	if(cosmicRayEnergy<=crGZKthresholdEnergy){
	    double thresholdZ = 0.0;
	    double zDown = getZCRSourceGZK2BHLimit(bhRatio, thresholdZ, zMax, crGZKthresholdEnergy, cosmicRayEnergy);
	    if(zDown<zMax){
		double relativeRefEnergy = crGZKthresholdEnergy/energyRef;
		double distanceDown = getDistance(zDown);
		double distanceZmax =  getDistance(zMax);
		double expFactor =  1.0/bhRatio*2.0/(3.0*omegaM)*(distanceDown-1.0);
		double distanceFactor = 2.0*distanceDown/Math.pow(1.0+zDown,5.0)+1.0/(bhRatio*Math.pow(1.0+zDown,2.0));
		parameters[0] = alpha;
		parameters[1] = m;
		parameters[2] = gzkRatio;
		parameters[3] = zDown;
		double redshiftTerm = Integration.RombergIntegral(this,5,parameters,zDown,zMax);
		double fluxGZKprocess = getCRSpectrumNormalization(crGZKnormalizationEnergy,alpha)*Math.exp(expFactor)/gzkRatio/distanceFactor*
		    Math.pow(relativeRefEnergy,-alpha)*Math.pow(1.0+zDown,(2.0*alpha-2.0))*redshiftTerm;
		//System.err.format("norm(%e)  zterm (%e) distance factor(%e)\n",
		//		  Math.exp(expFactor),redshiftTerm, distanceFactor);
		//System.err.format("flux from GZK process (%e)\n",fluxGZKprocess);
		flux += fluxGZKprocess;
	    }
	}else{ // flux from sources within GZK sphere
	    double distanceZmax =  getDistance(zMax);
	    double relativeEnergy = cosmicRayEnergy/energyRef;
	    parameters[0] = alpha;
	    parameters[1] = m;
	    parameters[2] = gzkRatio;
	    parameters[3] = 0.0;
	    double redshiftTerm = Integration.RombergIntegral(this,5,parameters,0.0,zMax);
	    double fluxWithinGZKSphere = getCRSpectrumNormalization(crGZKnormalizationEnergy,alpha)*c/H0*
		Math.pow(relativeEnergy,-alpha)*redshiftTerm;
	    //System.err.format(" flux within GZK sphere (%e)\n", fluxWithinGZKSphere);

	    flux += fluxWithinGZKSphere;
	}

	double dumpFactor = 1.0;
	if(applyExpDumingToCRFlux){
	    if(cosmicRayEnergy< crEnergyMax) dumpFactor = Math.exp(-optDepth);
	}
	return flux*dumpFactor;
    }

    
    public double getCRDFDEHalfNumerical(double zMax, double zConst, double m, double alpha, double cosmicRayEnergy){
	if(zMax<= zConst) return getCRDFDEHalfNumerical(zMax,m,alpha,cosmicRayEnergy);
	double fluxUptoZConst = getCRDFDEHalfNumerical(zConst,m,alpha,cosmicRayEnergy);
	if(fluxUptoZConst<=0.0) return 0.0;

        double[] parameters = new double[4];
	double flux = 0.0;

	// flux from sources at distances that only matters redshift loss
	double zUp = getZCRSourceLimit(bhRatio, zMax, zConst, crBHthresholdEnergy, cosmicRayEnergy);
	if(zUp>zConst){
	    parameters[0] = alpha;
	    parameters[1] = 0.0;
	    double redShiftTerm = Integration.RombergIntegral(this, 4, parameters, zConst, zUp);
	    double relativeEnergy = cosmicRayEnergy/energyRef;
	    flux = getCRSpectrumNormalization(crGZKnormalizationEnergy,alpha)*c/H0*
		Math.pow(relativeEnergy,-alpha)*redShiftTerm*Math.pow((1.0+zConst),m);
	    //System.err.format(" flux from redshift only (%e)\n", flux);
	}
	
	// flux from sources at distances that matters redshift-evolved GZK BH process
	if(cosmicRayEnergy<=crBHthresholdEnergy){
	    double zDown = getZCRSourceLimit(bhRatio, zMax, crBHthresholdEnergy, cosmicRayEnergy);
	    if((zDown<=zMax) && (zConst<=zDown)){
		double relativeRefEnergy = crBHthresholdEnergy/energyRef;
		double energyDistance = getDistance(zDown);
		parameters[0] = alpha;
		parameters[1] = 0.0;
		parameters[2] = bhRatio;
		parameters[3] = zDown;
		double redshiftTerm = Integration.RombergIntegral(this,5,parameters,zDown,zMax);
		double fluxBHprocess = 0.5*getCRSpectrumNormalization(crGZKnormalizationEnergy,alpha)*c/H0/(bhRatio*energyDistance)*
		    Math.pow(relativeRefEnergy,-alpha)*Math.pow(1.0+zDown, (2.0*alpha+3.0))*redshiftTerm*Math.pow((1.0+zConst),m);
		//System.err.format(" flux from BH process (%e)\n", fluxBHprocess);

		flux += fluxBHprocess;

	    }
	}else{ // flux from sources within GZK BH sphere
	    double thresholdZ = getZCRSourceGZK2BHLimit(bhRatio, 0.0, zMax, crGZKthresholdEnergy, cosmicRayEnergy);
	    if(thresholdZ>=zConst){
		double distanceZmax =  getDistance(thresholdZ);
		parameters[0] = alpha;
		parameters[1] = 0.0;
		parameters[2] = bhRatio;
		parameters[3] = 0.0;
		double redshiftTerm = Integration.RombergIntegral(this,5,parameters,zConst,thresholdZ);
		double relativeEnergy = cosmicRayEnergy/energyRef;
		double fluxWithinBHSphere = getCRSpectrumNormalization(crGZKnormalizationEnergy,alpha)*c/H0*
		    Math.pow(relativeEnergy,-alpha)*redshiftTerm*Math.pow((1.0+zConst),m);
		//System.err.format(" flux within BH process (%e)\n", fluxWithinBHSphere);
		flux += fluxWithinBHSphere;
	    }
	}

	// flux from sources at distances that matters redshift-evolved GZK photopion process
	if(cosmicRayEnergy<=crGZKthresholdEnergy){
	    double thresholdZ = 0.0;
	    double zDown = getZCRSourceGZK2BHLimit(bhRatio, thresholdZ, zMax, crGZKthresholdEnergy, cosmicRayEnergy);
	    if((zDown<=zMax) && (zConst <= zDown)){
		double relativeRefEnergy = crGZKthresholdEnergy/energyRef;
		double distanceDown = getDistance(zDown);
		double distanceZmax =  getDistance(zMax);
		double expFactor =  1.0/bhRatio*2.0/(3.0*omegaM)*(distanceDown-1.0);
		double distanceFactor = 2.0*distanceDown/Math.pow(1.0+zDown,5.0)+1.0/(bhRatio*Math.pow(1.0+zDown,2.0));
		parameters[0] = alpha;
		parameters[1] = 0.0;
		parameters[2] = gzkRatio;
		parameters[3] = zDown;
		double redshiftTerm = Integration.RombergIntegral(this,5,parameters,zDown,zMax);
		double fluxGZKprocess = getCRSpectrumNormalization(crGZKnormalizationEnergy,alpha)*Math.exp(expFactor)/gzkRatio/distanceFactor*
		    Math.pow(relativeRefEnergy,-alpha)*Math.pow(1.0+zDown,(2.0*alpha-2.0))*redshiftTerm*Math.pow((1.0+zConst),m);
		//System.err.format("norm(%e)  zterm (%e) distance factor(%e)\n",
		//		  Math.exp(expFactor),redshiftTerm, distanceFactor);
		//System.err.format("flux from GZK process (%e)\n",fluxGZKprocess);
		flux += fluxGZKprocess;
	    }

	}else{ // flux from sources within GZK sphere
	    double distanceZmax =  getDistance(zMax);
	    double relativeEnergy = cosmicRayEnergy/energyRef;
	    parameters[0] = alpha;
	    parameters[1] = 0.0;
	    parameters[2] = gzkRatio;
	    parameters[3] = 0.0;
	    double redshiftTerm = Integration.RombergIntegral(this,5,parameters,zConst,zMax);
	    double fluxWithinGZKSphere = getCRSpectrumNormalization(crGZKnormalizationEnergy,alpha)*c/H0*
		Math.pow(relativeEnergy,-alpha)*redshiftTerm*Math.pow((1.0+zConst),m);
	    //System.err.format(" flux within GZK sphere (%e)\n", fluxWithinGZKSphere);

	    flux += fluxWithinGZKSphere;
	}
	
	flux += fluxUptoZConst;

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
	System.err.format("alpha(%f) m(%f) zMax(%f) nuE=%e opticalDepth(%f)\n",alpha,m,zMax,neutrinoEnergy,neutrinoFlux.optDepth);
	double[] parameters = new double[5];
	parameters[0] = alpha;
	parameters[1] = zMax;
	if(!assumeConstEvolution) parameters[2] = m;
	else{
	    parameters[2]= zConst;
	    parameters[3] = m;
	}

	double zBoundBH = neutrinoFlux.getZCRSourceLimit(bhRatio, zMax, crBHthresholdEnergy, neutrinoEnergy);
	double zBoundGZK2BH = neutrinoFlux.getZCRSourceGZK2BHLimit(bhRatio, 0.0, zMax, 0.0, crGZKthresholdEnergy, neutrinoEnergy);
	System.err.format("zboundBH (%e) zBoundGZK2BH(%e) for E=%e\n",
			  zBoundBH,zBoundGZK2BH, neutrinoEnergy);
	//System.exit(0);

	neutrinoFlux.setPowerLawIndexOfRadiation(1.0);
	double crEnergyRef = 1.0e10; // [GeV]
	neutrinoFlux.setCREnergyReference(crEnergyRef);
	double L_gamma = 1.0e43;
	neutrinoFlux.setCRFluxFromSourcePhotonLuminosity(L_gamma, 0.1*crEnergyRef, 10*crEnergyRef,
							 alpha, 3.0e-7, 1.0);
	System.out.format("relative CR flux at GZK energy %e for L_g(%e)\n",relativeCRFlux,L_gamma);
	neutrinoFlux.returnToTheCRFluxNormalizationEstimatedFromGZK();
	
	double eCRFlux = 1.6e-8; // [GeV/cm2 sec sr]
	if(!assumeConstEvolution) neutrinoFlux.setCRFluxAtReferenceEnergy(zMax, m, alpha, eCRFlux);
	else  neutrinoFlux.setCRFluxAtReferenceEnergy(zMax, zConst, m, alpha, eCRFlux);
	System.out.format("relative CR flux at GZK energy %e\n",relativeCRFlux);

	if(!assumeConstEvolution){
	    System.out.format("CR flux at reference energy(%e)= %e\n",crEnergyRef,neutrinoFlux.getCRDFDE(zMax,m,alpha,crEnergyRef));
	}else{
	    System.out.format("CR flux at reference energy(%e)= %e\n",crEnergyRef,neutrinoFlux.getCRDFDE(zMax,zConst,m,alpha,crEnergyRef));
	}

	double energyFlux = 0.0;
	if(!assumeConstEvolution) energyFlux = neutrinoFlux.getDFDE(zMax,m, alpha,neutrinoEnergy)*neutrinoEnergy*neutrinoEnergy;
	else energyFlux = neutrinoFlux.getDFDE(zMax, zConst, m, alpha,neutrinoEnergy)*neutrinoEnergy*neutrinoEnergy;
	System.out.format(" E^2dF/dE (normalized at the 1EeV CR Flux) = %e [GeV /cm^2 sec sr]\n",energyFlux);

	double relativeCRFluxBackup = relativeCRFlux;
	neutrinoFlux.returnToTheCRFluxNormalizationEstimatedFromGZK();
	if(!assumeConstEvolution) energyFlux = neutrinoFlux.getDFDE(zMax,m, alpha,neutrinoEnergy)*neutrinoEnergy*neutrinoEnergy;
	else energyFlux = neutrinoFlux.getDFDE(zMax, zConst, m, alpha,neutrinoEnergy)*neutrinoEnergy*neutrinoEnergy;
	System.out.format(" E^2dF/dE (normalized by the GZK intensity) = %10.7e [GeV /cm^2 sec sr]\n",energyFlux);
	neutrinoFlux.setTheCRFluxNormalizationEstimatedFromGZK(relativeCRFluxBackup);

	double gzkEnergy = crBHthresholdEnergy; //[GeV]
	for(int i=0;i<12;i++){
	    double cosmicRayEnergy = Math.pow(10.0,9.0+0.2*(double )i);
	    zBoundBH = neutrinoFlux.getZCRSourceLimit(bhRatio, zMax, crBHthresholdEnergy, cosmicRayEnergy);
	    zBoundGZK2BH = neutrinoFlux.getZCRSourceGZK2BHLimit(bhRatio, 0.0, zMax, crGZKthresholdEnergy, cosmicRayEnergy);
	    double flux = 0.0;
	    if(!assumeConstEvolution) flux = neutrinoFlux.getCRDFDE(zMax,m,alpha,cosmicRayEnergy);
	    else flux =  neutrinoFlux.getCRDFDE(zMax,zConst,m,alpha,cosmicRayEnergy);
	    double fluxNumerical = 0.0;
	    if(!assumeConstEvolution) fluxNumerical = neutrinoFlux.getCRDFDEHalfNumerical(zMax,m,alpha,cosmicRayEnergy);
	    else fluxNumerical =  neutrinoFlux.getCRDFDEHalfNumerical(zMax,zConst,m,alpha,cosmicRayEnergy);
	    
	    double intensity = 0.0; double crEndEnergy = 2.0e11;
	    //if(!assumeConstEvolution)
	    //	intensity = Integration.RombergIntegral(neutrinoFlux, 
	    //						2, parameters, cosmicRayEnergy, crEndEnergy);
	    //else
	    //	intensity = Integration.RombergIntegral(neutrinoFlux, 
	    //						3, parameters, cosmicRayEnergy, crEndEnergy);
		
	    System.out.format("CR Energy(%e) zBoundBH (%f) zBoundGZK2BH (%f) Flux(%e) FluxHalfNumerical(%e)\n",
			      cosmicRayEnergy,zBoundBH,zBoundGZK2BH,flux,fluxNumerical);
	}

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
