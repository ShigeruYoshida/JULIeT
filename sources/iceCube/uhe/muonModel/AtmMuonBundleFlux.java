package iceCube.uhe.muonModel;

import iceCube.uhe.interactions.*;
import iceCube.uhe.geometry.*;
import geometry.*;
import numRecipes.*;

import java.util.*;


/**
<pre>
UHE Atmospheric Muon fluxes based on the Elbert model
taking into account the muon bundles from EAS cascades.

</pre>
<UL>
  <DT> <cite>Cosmic Rays and Particle Physics</cite> 
       <TABLE CELLPADDING=2>
       <TR>
         <TD ALIGN="LEFT"><it> Tom Gaisser ,</it></TD>
         <TD ALIGN="LEFT">
         Cambridge University Press, Cambridge, England, 1990
       </TR>
       </TABLE><BR>
</UL>

*/


public class AtmMuonBundleFlux extends ParticleFlux implements Function{

    private static final double ln10 = Math.log(10.0);

    /** 
	a power law index of the Elbert Formula. The model-predicted value
    */
    private final static double alpha_model = 1.757;
    /**
	a power law index of the Elbert Formula. The IceCube favored value
    */
    private final static double alphaByIceCube = 2.04;

    /** 
	a power law index of the hard term of the Elbert formula.
	Its nominal value is 5.25 and NEVER changed.
    */
    private final static double beta = 5.25;

    /** effective energy of the Elbert formula */
    private final static double effEnergy = 14.5; // 14.5 GeV

    /** critical energy of ionization loss */
    public final static double criticalEnergy = 500.0; // 500 GeV

    /** mass number of Primary cosmic rays */
    private double massNumber;

    /** 
	a power law index of the Elbert Formula. 
    */
    private double alpha;

    /**
       A lowest threshold energy of muons in a bundle [GeV].
    */
    private double muonThresholdE;

    /**
       Primary Cosmic ray flux object
    */
    private CosmicRayFlux primaryFlux = null;

    /**
       Cascade flucutation calculator object 
    */
    public CascadeFluctuationFactory cascade = null;

    public boolean includeFluctuationEffects = true; 
    // if false, calculate only the bare number determined by the Elbert formula 

    /** Options to calculate flux with/without taking into
        account fluxcuation of muon energies due to
	the EAS cascade process

	<pre>
	true :  Flux is calculated such that an event by event flucuation is included.
                Useful to evaluate flucuation of number of background events
        false: the flucuation effects is included in average way. i.e.,
               The systematic factor of the flux shift due the cascade flucuation
               is calculated and multiplied in the flux calculation.
               Useful to the MC-data fitting, calculation of average number of events
       </pre>
    */
    public boolean fluctuateEventByEvent= false;

    /** Confidence level of the Energy rato parameter 
	R = (E_muon/E0)/ Bar(E_muon/E0). If this value is non zero (zero in its default)
	and (includeFluctuationEffect fluctuateEventByEvent) = (true, true)
	then all the relevant fluxes are calculated based on a fixed value of R 
        given by this confidence level.
    */
    public double confidenceLevel_of_R = 0.0;

    /** The typical height of the muon production in the atmosphere */
    private double muon_production_height = 1.0e6; // 10 km
    protected boolean considerAtmosphericCurvature = true;

    private RandomGenerator rand = null;
    private double epsilon = 0.01; // the relative accuracy in the integration



    /**
       Default Constructor. Fill alpha, mass number, and muon energy threshold
       with their moninal vlues. It also instansizes the cosmic ray flux class.
    */
    public AtmMuonBundleFlux(){
	setMassNumber(1.0); // A = 1
	//setAlpha(alpha_model);
	setAlpha(alphaByIceCube);
	setMuonThresholdEnergy(3.73e3); // 3730 GeV - IceCube favorite
	primaryFlux = new CosmicRayFlux();
	cascade = new CascadeFluctuationFactory();
	Integration.setRelativeAccuracy(epsilon);
    }


    /**
       This constuctor sets alpha, mass number, and muon energy threshold
       by the values given in its arguments. The cosmic ray flux object is also 
       generated.
    */
    public AtmMuonBundleFlux(double alpha, double mass_number, 
			     double energy_threshold){
	setMassNumber(mass_number);
	setAlpha(alpha);
	setMuonThresholdEnergy(energy_threshold);
	primaryFlux = new CosmicRayFlux();
	cascade = new CascadeFluctuationFactory();
	Integration.setRelativeAccuracy(epsilon);
    }



    /** 
	Check if the alpha value in setting is within the valid
        range of the Elbert model.
    */
    public boolean isValidAlpha(double alpha){
	if(1.0<alpha && alpha < 2.5) return true;
	else return false;
    }


    /**
       Set the alpha paramerter
    */
    public void setAlpha(double alpha){
	if(isValidAlpha(alpha)){
	       this.alpha = alpha;
	}else{
	    System.err.println("You tried setting alpha out of range");
	    System.err.println("alpha has remained same");
	}
    }


    /**  Get the alpha parameter */
    public double getAlpha(){
	return alpha;
    }


    /** set mass number of the primary cosmic rays */
    public void setMassNumber(double mass_number){
	massNumber = mass_number;
    }


    /**  Get the mass number */
    public double getMassNumber(){
	return massNumber;
    }

    /**
       Set the threshold energy [GeV] of muons contained in a bundle
    */
    public void setMuonThresholdEnergy(double energy_threshold){
	muonThresholdE = energy_threshold;
    }

    /**  Get the threshold energy of muons [GeV] */
    public double getMuonThresholdEnergy(){
	return muonThresholdE;
    }


    /** Set on/off the GZK feature of the primary cosmic ray spectrum */
    public void setCutOffFeature(boolean cutoffExists){
	primaryFlux.setCutOffFeature(cutoffExists);
    }

    /** set the mode on  the flux caulation 
	<pre>
	boolean includeFluctuationEffects : 
                    true  default. calculate the flux taking into account flucuation of 
                          muon energies due to the EAS cascading. 
                          CascadeFluctuationFactory class does this part of calculation.
                    false calculate the bare flux given by the Elbert formula.
                          an energy of muon bundles is associated with the energy of primary
                          cosmic ray by one-on-one relation.

	boolean fluctuateEventByEvent : Valid when includeFluctuationEffects = true

	   true :  Flux is calculated such that an event by event flucuation is included.
                   Useful to evaluate flucuation of number of background events
           false: default. the flucuation effects is included in average way. i.e.,
                  The systematic factor of the flux shift due the cascade flucuation
                  is calculated and multiplied in the flux calculation.
                  Useful to the MC-data fitting, calculation of average number of events
	</pre>
     **/
    public void setFluxCalculationMode(boolean includeFluctuationEffects,
				       boolean fluctuateEventByEvent){
	this.includeFluctuationEffects = includeFluctuationEffects;
	this.fluctuateEventByEvent = fluctuateEventByEvent;

	if(fluctuateEventByEvent && rand == null) rand =  new RandomGenerator();
    }


    /** Set Confidence level of the Energy rato parameter 
	R = (E_muon/E0)/ Bar(E_muon/E0). If this value is non zero (zero in its default)
	and (includeFluctuationEffect fluctuateEventByEvent) = (true, true)
	then all the relevant fluxes are calculated based on a fixed value of R 
        given by this confidence level. If zero (default), the fluxes are given
	by the integral over R which correponds to taking the central value of R-distribution.
    */
    public void setConfidenceLevelOfFluctuation(double prob){
	confidenceLevel_of_R = prob;
    }


    /** Calculate the zenith angle at the earth surface.
	The zenith angle at the IceCube depth is
	transformed to that at the surface accounting
	the earth curvature.
    */
    public double getCosineOfZenithAngleAtEarthSurface(double cos_theta_ice3){

	IceCubeCoordinate ice3Coordinate =  new IceCubeCoordinate();
	double ice3elevation =
	    EarthCenterCoordinate.REarth+IceCubeCoordinate.elevation;
	double radius = 
	    EarthCenterCoordinate.REarth+ice3Coordinate.getGlacierDepth();

	double a = ice3elevation*cos_theta_ice3;
	double trackLength = -a +
	    Math.sqrt(a*a + radius*radius - ice3elevation*ice3elevation);
	double cos_theta_center = 
	    (trackLength*trackLength + trackLength*ice3elevation*cos_theta_ice3)/
	    (trackLength*radius);

	return cos_theta_center;
    }

    /** Calculate the zenith angle at the muon production height up in the
	atmosphere. This angle is equal to the value returned by
	getCosineOfZenithAngleAtEarthSurface(double cos_theta_ice3)
	IF YOU NEGLET the EARTH CURVATURE effect in the atomosphere.
	Here this method considers the curvature effect.

	<pre>
	double cos_theta_earth: cosZenith at the earth surface
	</pre>
    */

    public double getCosineOfZenithAngleAtAtmosphericMuonHeight(double cos_theta_earth){
	IceCubeCoordinate ice3Coordinate =  new IceCubeCoordinate();
	double radius = 
	    EarthCenterCoordinate.REarth+ice3Coordinate.getGlacierDepth();
	double a = radius*cos_theta_earth;
	double trackLength = -a +
	    Math.sqrt(a*a - radius*radius + 
		      (radius + muon_production_height)*(radius + muon_production_height));
	double cos_theta_muHeight = muon_production_height/trackLength;
	return cos_theta_muHeight;
    }

    /** If true, then the Elbert formula takes into account the atmosphere density
	profile is curved following the earth curvature. This is a newly added
	function. The older version is equal to the case of false. */
    public void considerAtmosphericCurvature(boolean consider){
	considerAtmosphericCurvature = consider;
    }



    /** Calculate the effective energy [GeV] of the parent cosmic ray
	for a given energy of the EHE muon bundle and its zenith angle
        at the IceCube depth.
    */
    public double getEffectiveEnergyOfCRs(double muonBundleEnergy, 
					  double cosTheta_ice3){
	double cosTheta_earth = 1.0;
	if(considerAtmosphericCurvature){ // The curvature in the atmosphere considered
	    cosTheta_earth = 
		getCosineOfZenithAngleAtAtmosphericMuonHeight(
			  getCosineOfZenithAngleAtEarthSurface(cosTheta_ice3));
	}else{
	    cosTheta_earth = getCosineOfZenithAngleAtEarthSurface(cosTheta_ice3);
	}
	double powerIndex = 1.0/(alpha-1.0);
	double eFactor = cosTheta_earth/massNumber*((alpha-1.0)/alpha)*
	    muonBundleEnergy/effEnergy;

	double cosmicRayEnergy = Math.pow(eFactor,powerIndex)*massNumber*
	    muonThresholdE;
	return cosmicRayEnergy;
    }

    /** Calculate the effective energy [GeV] of the parent cosmic ray
	for a given IN-ICE energy of the EHE muon bundle and its zenith angle
        at the IceCube depth.
        The muon propagation effects are accounted by the simple 
        approximation exp[-beta*x] ignoring the ionization loss.
    */
    public double getEffectiveEnergyOfCRs(double logEnergy, double cosTheta, 
		      double beta_loss, double slantDepth){

	double energyEnhancedFactor = Math.exp(beta_loss*slantDepth);
	double muonBundleEnergy = Math.pow(10.0,logEnergy);

	double muonThresholdE_org = muonThresholdE;
	muonThresholdE = energyEnhancedFactor*(muonThresholdE+criticalEnergy) 
	    -criticalEnergy;
	double cosmicRayEnergy =
	    getEffectiveEnergyOfCRs(energyEnhancedFactor*
				    muonBundleEnergy,cosTheta);
	muonThresholdE = muonThresholdE_org;
	return cosmicRayEnergy;
    }


    /** Calculate dF/dLogE [/cm^2 sec sr] at the earth surface as a function
        of log E [GeV] and cosine of zenith angle */
    public double getDFDLogE(double logEnergy, double cosTheta){

	double muonBundleEnergy = Math.pow(10.0,logEnergy);
	double cosmicRayEnergy = 
	    getEffectiveEnergyOfCRs(muonBundleEnergy,cosTheta);
	double logCosmicRayEnergy = Math.log(cosmicRayEnergy)/ln10;
	double dMuonBundleDLogE = (1.0/(alpha-1.0))*
	    primaryFlux.getDFDLogE(logCosmicRayEnergy);

	if(!includeFluctuationEffects)	return dMuonBundleDLogE; // the bere flux

	else if (!fluctuateEventByEvent){ // includes the fluctuation 
	                                   // in average bases

	    if(Math.abs(confidenceLevel_of_R) <= 1.0e-4){
                                          // Integral to obtain the central value of 
                        		  //(E_muon/E0)/ Bar(E_muon/E0)
		double[] param = new double[2];
		param[0] = logCosmicRayEnergy;
		param[1]= cosTheta;
		double fluxShiftFactor = 
		    Integration.RombergIntegral(this,2,param,
					cascade.getLogMuOverCREnergyMin(logCosmicRayEnergy,false),
					cascade.getLogMuOverCREnergyMax(logCosmicRayEnergy,false));
	    // integral dLog10R F_R/F_1 x rho(log10R), 
	    // R = (E_muon/E0)/ Bar(E_muon/E0) where E_muon is the one at surface
  	    //System.err.println(" " +  logEnergy + " r=" + fluxShiftFactor);
		return fluxShiftFactor*dMuonBundleDLogE;
	    }else{
		double logR =           // (E_muon/E0)/ Bar(E_muon/E0) at surface
		    cascade.getLogMuOverCREnergy(confidenceLevel_of_R,logCosmicRayEnergy,false);
		dMuonBundleDLogE =  (1.0/(alpha-1.0))*primaryFlux.getDFDLogE(logCosmicRayEnergy - logR);
		return dMuonBundleDLogE;
	    }
	}else{
	    double logR =               // (E_muon/E0)/ Bar(E_muon/E0)
	    cascade.sampleLogEnergyRatioFactor(rand,logCosmicRayEnergy,false); 
	    dMuonBundleDLogE =  (1.0/(alpha-1.0))*primaryFlux.getDFDLogE(logCosmicRayEnergy - logR);
	    return dMuonBundleDLogE;
	}
    }



    /** Calculate in-ice dF/dLogE [/cm^2 sec sr] at the slant depth x [g/cm^2] 
        as a function  of log E [GeV] in ice and cosine of zenith angle.
        The muon propagation effects are accounted by the simple 
        approximation exp[-beta*x] ignoring the ionization loss.
    */
    public double getDFDLogE(double logEnergy, double cosTheta, 
		      double beta_loss, double slantDepth){

	double logCosmicRayEnergy = 
	    Math.log(getEffectiveEnergyOfCRs(logEnergy,cosTheta,
					     beta_loss,slantDepth))/ln10;

	double dMuonBundleDLogE = (1.0/(alpha-1.0))*
	    primaryFlux.getDFDLogE(logCosmicRayEnergy);


	if(!includeFluctuationEffects)	return dMuonBundleDLogE; // the bere flux

	else if (!fluctuateEventByEvent){ // includes the fluctuation 
	                                   // in average bases

	    if(Math.abs(confidenceLevel_of_R) <= 1.0e-4){
                                          // Integral to obtain the central value of 
                        		  //(E_muon/E0)/ Bar(E_muon/E0)
		double[] param = new double[2];
		param[0] = logCosmicRayEnergy;
		param[1]= cosTheta;
		double fluxShiftFactor = 
		    Integration.RombergIntegral(this,1,param,
					cascade.getLogMuOverCREnergyMin(logCosmicRayEnergy,true),
					cascade.getLogMuOverCREnergyMax(logCosmicRayEnergy,true));
	    // integral dLog10R F_R/F_1 x rho(log10R), 
	    // R = (E_muon/E0)/ Bar(E_muon/E0) where E_muon is the one in ince
		return fluxShiftFactor*dMuonBundleDLogE;
	    }else{
		double logR =           // (E_muon/E0)/ Bar(E_muon/E0) in ice
		    cascade.getLogMuOverCREnergy(confidenceLevel_of_R,logCosmicRayEnergy,true);
		dMuonBundleDLogE =  (1.0/(alpha-1.0))*primaryFlux.getDFDLogE(logCosmicRayEnergy - logR);
		return dMuonBundleDLogE;
	    }
	}else{                         // (E_muon/E0)/ Bar(E_mu/E0) inice
	    double logR = cascade.sampleLogEnergyRatioFactor(rand,logCosmicRayEnergy,true);
	    dMuonBundleDLogE =  (1.0/(alpha-1.0))*primaryFlux.getDFDLogE(logCosmicRayEnergy - logR);
	    return dMuonBundleDLogE;
	}
    }

    /** Calculate in-ice dF/dLogE [/cm^2 sec sr] at the slant depth x [g/cm^2] 
        as a function  of log E [GeV] in ice and cosine of zenith angle.
        The muon propagation effects are accounted by the simple 
        approximation exp[-beta*x] ignoring the ionization loss.
	The "beta" term used in this method is caluclated inside this function
        by the value given by CELbeta static class in the interaction package.
    */
    public double getDFDLogE(double logEnergy, double cosTheta, 
			     double slantDepth){

	double beta_loss = 4.4619776009127244E-6;
	if(logEnergy >= 7.0) beta_loss = CELbeta.getBeta(logEnergy);

	return getDFDLogE(logEnergy,cosTheta,beta_loss,slantDepth);

    }


    /** Calculate the net fluctuation of the flux due
        to the cascade process. Returns the relative flux, i.e.,
        multiplying this factor to the bare intrinsic flux 
	gives the flux including the fluctuation effect.
        Depending on the option set by  setFluxCalculationMode(),
       this method calculate the value in the average basis or with MC samping.

       <pre>
       double logEnergy       : Log(Muon Bundle Energy in ice [GeV])
       </pre>
        The muon propagation effects are accounted by the simple 
        CEL approximation.
    */
    public double getFluctuatedRelativeFlux(double logEnergy, double cosTheta, 
					    double slantDepth){

	double beta_loss = 4.4619776009127244E-6;
	if(logEnergy >= 7.0) beta_loss = CELbeta.getBeta(logEnergy);

	double logCosmicRayEnergy = 
	    Math.log(getEffectiveEnergyOfCRs(logEnergy,cosTheta,
					     beta_loss,slantDepth))/ln10;

	double fluxShiftFactor;
	if (!fluctuateEventByEvent){ // includes the fluctuation 
	                            // in average bases

	    if(Math.abs(confidenceLevel_of_R) <= 1.0e-4){
                                          // Integral to obtain the central value of 
                        		  //(E_muon/E0)/ Bar(E_muon/E0)
		double[] param = new double[2];
		param[0] = logCosmicRayEnergy;
		param[1]= cosTheta;
		fluxShiftFactor = 
		    Integration.RombergIntegral(this,1,param,
					cascade.getLogMuOverCREnergyMin(logCosmicRayEnergy,true),
					cascade.getLogMuOverCREnergyMax(logCosmicRayEnergy,true));
	    // integral dLog10R F_R/F_1 x rho(log10R), 
	    // R = (E_muon/E0)/ Bar(E_muon/E0) where E_muon is the one in ince
	    }else{
		double logR =           // (E_muon/E0)/ Bar(E_muon/E0) in ice
		    cascade.getLogMuOverCREnergy(confidenceLevel_of_R,logCosmicRayEnergy,true);
		fluxShiftFactor =  primaryFlux.getDFDLogE(logCosmicRayEnergy - logR)/
		 primaryFlux.getDFDLogE(logCosmicRayEnergy);
	    }
	}else{                         // (E_muon/E0)/ Bar(E_mu/E0) inice
	    double logR = cascade.sampleLogEnergyRatioFactor(rand,logCosmicRayEnergy,true); 
	    fluxShiftFactor =  primaryFlux.getDFDLogE(logCosmicRayEnergy - logR)/
		 primaryFlux.getDFDLogE(logCosmicRayEnergy);
	}
	return fluxShiftFactor;
    }

    /** Calculate energy flux E^2 dF/dE [GeV/cm^2 sec sr] */
    public double getEFlux(double logEnergy, double cosTheta){
	double energy = Math.pow(10.0,logEnergy);
	double EFlux = getDFDLogE(logEnergy,cosTheta)*energy/ln10;
	
	return(EFlux);
    }

    /** Calculate energy flux E^2 dF/dE [GeV/cm^2 sec sr] 
	at the slant depth x [g/cm^2].
        The muon propagation effects are accounted by the simple 
        approximation exp[-beta*x] ignoring the ionization loss.
    */
    public double getEFlux(double logEnergy, double cosTheta,
		    double beta, double slantDepth){
	double energy = Math.pow(10.0,logEnergy);
	double EFlux = getDFDLogE(logEnergy,cosTheta,beta,slantDepth)*energy/ln10;
	
	return(EFlux);
    }


    /**
       Calculate d^2F/dLogEcrdLogEmuon [/cm^2 sec sr] as a function of
       cosmic ray energy and muon (bundle) energy at the earth surface.
       When neglect the airshower cascade flucutation, these two energies are
       related on one-on-one, but introducing the fluctuation leads to cosmic ray
       energy distribution that is responsieble for bundles with a given energy.

       The method  getDFDLogE(double logEnergy, double cosTheta) with
       includeFluctuationEffects =true, and fluctuateEventByEvent =false
       corresponds to integral of this flux over the cosmic ray energy.

       <pre>
       double logCosmicRayEneergy  : log10(Cosmic Ray   Energy [GeV])
       double logMuonBundleEneergy : log10(Muone Bundle Energy [GeV])
       double cosTheta             : cos(Zenith) of the bundle
       </pre>
    */
    public double getDFDLogCREDLogE(double logCosmicRayEnergy, double logMuonBundleEnergy,
					 double cosTheta){

	double muonBundleEnergy = Math.pow(10.0,logMuonBundleEnergy);

	// Mean cosmic ray energy responsible for mu-bundles with energy of muonBundleEnergy 
	// determined by the Elbert formula parameterization
	double meanCosmicRayEnergy = 
	    getEffectiveEnergyOfCRs(muonBundleEnergy,cosTheta);
	double logMeanCosmicRayEnergy = Math.log(meanCosmicRayEnergy)/ln10;
	double meanCosmicFlux = primaryFlux.getDFDLogE(logMeanCosmicRayEnergy);
	if(meanCosmicFlux<=0.0) return 0.0;

	double logR = logMeanCosmicRayEnergy - logCosmicRayEnergy;
	double prob_R=0.0;
	if(cascade.getLogMuOverCREnergyMin(logMeanCosmicRayEnergy,false)<= logR &&
	   logR <= cascade.getLogMuOverCREnergyMax(logMeanCosmicRayEnergy,false)){
	    prob_R = cascade.getProbability(logR,logMeanCosmicRayEnergy,false); 
	}else{
	    return 0.0;
	}

	// 
	double dMuonBundleDLogE = (1.0/(alpha-1.0))*
	    primaryFlux.getDFDLogE(logCosmicRayEnergy)*prob_R;

	return dMuonBundleDLogE;

    }

    public double integralDFDLogCREDLogEOverCREnergy(double lowerLogCosmicRayEnergy,
						     double upperLogCosmicRayEnergy,
						     double logMuonBundleEnergy,
						     double cosTheta){
	double[] param = new double[2];
	param[0] = logMuonBundleEnergy;
	param[1]= cosTheta;
	// integral from lowerLogCosmicRayEnergy to upperLogCosmicRayEnergy
	double dFdLogEbyIntegral = 
	    Integration.RombergIntegral(this,3,param,lowerLogCosmicRayEnergy,upperLogCosmicRayEnergy); 
	return dFdLogEbyIntegral;
    }

    /** integrate getDFDLogCREDLogE(logEcr,logEmu,cosTheta) over logEcr
	and compare with getDFDLogE(logEmu,cosTheta) for the consistency check
    */
    protected void listFluxes(double logEnergy, double cosTheta){

	double dFdLogEbyIntegral = 
	    integralDFDLogCREDLogEOverCREnergy(4.0,12.0,logEnergy,cosTheta); // integral from 100TeV to 1000 EeV
	System.err.println(" muon bundle flux by integral(logE=" + logEnergy + ", cosZen=" + cosTheta +
			   ") = " + dFdLogEbyIntegral);

	double dFDLogE = getDFDLogE(logEnergy, cosTheta);
	System.err.println(" muon bundle flux(logE=" + logEnergy + ", cosZen=" + cosTheta +
			   ") = " + dFDLogE);
    }

    /** Implementation of interface numRecipes.Function
	used for numerical integration of the probability
	of (E_muon/E0)/ Bar(E_muon/E0) multiplied by CosmicRayFlux (Energy x R)

	<pre>
	int functionIndex  :  1 calculation with inice muon, 2 with the surface
	double x = R = (E_muon/E0)/ Bar(E_muon/E0)
	
	functionIndex =3 then x = logCosmicRayEnergy

	parameter[0]   : log10(primary Cosmic Ray Energy [GeV]) to be responsible
                         for the muon bundle in AVERAGE i.e. R=1
	parameter[1]   : cosine(Zenith angle)

	</pre>

    */
    public double getFunction(int functionIndex, double[] parameters, 
	double x){

	double logEnergy = parameters[0];
	double cosZenith = parameters[1];

	if(functionIndex<=2){ // concerns the bundle flux with logR
	    double referenceCosmicFlux = primaryFlux.getDFDLogE(logEnergy); // flux at R=1;
	    if(referenceCosmicFlux<=0.0) return 0.0;

	    double logR = x;
	    double r_Emin; double r_Emax;
	    if(functionIndex<=1) r_Emin = cascade.getLogMuOverCREnergyMin(logEnergy,true) ;
	    else r_Emin = cascade.getLogMuOverCREnergyMin(logEnergy,false) ;
	    if(logR< r_Emin) logR = r_Emin;
	    double cosmicFlux = primaryFlux.getDFDLogE(logEnergy - logR);

	    double prob_R;
	    if(functionIndex<=1) prob_R = cascade.getProbability(x,logEnergy,true); // inice
	    else  prob_R = cascade.getProbability(x,logEnergy,false); // surface
	

	    return cosmicFlux/referenceCosmicFlux*prob_R;
	}else{ // concerns the bundle flux with log(cosmic ray energy)
	    double logCosmicRayEnergy = x;
	    return getDFDLogCREDLogE(logCosmicRayEnergy,logEnergy,cosZenith);
	}
    }

}

