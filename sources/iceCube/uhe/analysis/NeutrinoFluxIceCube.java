package iceCube.uhe.analysis;

import numRecipes.*;

import iceCube.uhe.neutrinoModel.*;
import java.io.*;

/** Calculate neutrino fluxes at the IceCube depth as a function
of neutrno energies at EARTH SURFACE. The propagation effects and 
the detection efficiency is taken care of by NeutrinoEffAreaTable class
in the analysis package - it just grabbed the results of the IceCube
IC40 and IC86 analysis by Aya.

Three modes exist:

Mode 1 : Use analtical model for GZK, no in-situ neutrino emission, i.e., pure GZK
Mode 2 : Use analtical model for GZK + in-situ neutrino emission
Mode 3 : Use numericaly calculated GZK model fluxes using NeutrinoFlux object.

Written by S.Yoshida 2011 May 1st
*/

class NeutrinoFluxIceCube implements Function {

    static final double ln10 = Math.log(10.0);
    static double solidAngle = 4.0*Math.PI; // 4pi, the neutrino observation solid angle

    static final double day_to_sec = 8.64e4;
    double livetime = 5.0*365.0*day_to_sec; // 365 days = 1 year


    NeutrinoFlux neutFlux = null;
    NeutrinoFluxFunction neutFluxFunction = null;
    NeutrinoFluxFromSource neutFluxSource = null;
    NeutrinoEffAreaTable areaTable = null;
    NeutrinoExposureTable exposureTable = null;
    NeutrinoExposureTableHESE exposureTableHESE = null;

    boolean analyticalFunction = true; // Whether you use the NeutrinoFluxFuntion class
    boolean addInSituEmission = false; // Whether you use the NeutrinoFluxFromSource class
    boolean isExposure = true;
    boolean addHESE = false;

    boolean assumeConstEvolution = false;
    double zConst = 2.0;

    int modelID = 4; // ID to specify GZK model in NeutrinoFlux class

    /**
       Constructor
       Mode 1 : Use analtical model for GZK, no in-situ neutrino emission, i.e., pure GZK
       Mode 2 : Use analtical model for GZK + in-situ neutrino emission
       Mode 3 : Use numericaly calculated GZK model fluxes using NeutrinoFlux object.
    */
    public NeutrinoFluxIceCube(int mode, boolean isExposure) throws IOException {

	this.isExposure = isExposure;

	switch (mode) {
	case 1:
	    analyticalFunction = true;
	    addInSituEmission = false;
	    break;
	case 2:
	    analyticalFunction = true;
	    addInSituEmission = true;
	    break;
	case 3:
	    analyticalFunction = false;
	    addInSituEmission = false;
	    break;
	default:
	    System.err.println("Illegal parameters! mode" + mode);
	    System.exit(0);
	}

	DataInputStream input = new DataInputStream(System.in); 
	BufferedReader  d     = new BufferedReader(new InputStreamReader(input)); 
	String buffer; 

	if(analyticalFunction){ // Go to the analytical model
	    neutFluxFunction = new NeutrinoFluxFunction();
	    System.err.print("Assume const evolution above a given redshift? (yes 1 no 0) ->"); 
	    buffer   = d.readLine(); 
	    if(Integer.valueOf(buffer).intValue()==1) {
		assumeConstEvolution = true;
		System.err.print("threshold redshift ->"); 
		buffer   = d.readLine(); 
		zConst = Double.valueOf(buffer).doubleValue();
		System.err.format("Setted zConst=%f\n",zConst);
	    }
	    if(addInSituEmission){
		neutFluxSource = new NeutrinoFluxFromSource();
		System.err.print("maximum neutrino energy from source [GeV] ->"); 
		buffer   = d.readLine(); 
		double neutrinoEnergyMax = Double.valueOf(buffer).doubleValue();
		neutFluxSource.setMaximumNeutrinoEnergy(neutrinoEnergyMax);
		System.err.println("minimum target photon energy [GeV] " + neutFluxSource.getTargetPhotonEnergyMinimum()); 
		System.err.print("power law index of the target photon field spectrum ->"); 
		buffer   = d.readLine(); 
		double gamma = Double.valueOf(buffer).doubleValue();
		neutFluxSource.setPowerLawIndexOfRadiation(gamma);
		System.err.print("optical depth ->"); 
		buffer   = d.readLine(); 
		double opticalDepth = Double.valueOf(buffer).doubleValue();
		neutFluxSource.setOpticalDepth(opticalDepth);
		System.err.print("Apply exponential duming factor to the CR flux? (yes 1 no 0) ->"); 
		buffer   = d.readLine(); 
		if(Integer.valueOf(buffer).intValue()==1) neutFluxSource.applyDumpingFactorToCRFlux();
	    }
	}else{ // Go to the numerical model
	    System.err.print("neutrino flux model ID ->"); 
	    buffer   = d.readLine(); 
	    modelID = Integer.valueOf(buffer).intValue();
	    neutFlux = new NeutrinoFlux(modelID);
	    if(neutFlux.isPowerLaw()){ // Set the logEmax for the power law
		System.err.print("log(EnergyMax/GeV) ->"); 
		buffer   = d.readLine(); 
		double logEmax = Double.valueOf(buffer).doubleValue();
		neutFlux.setLogEmax(logEmax);
		System.err.format("log(EnergyMax/GeV) set as %f\n",neutFlux.getLogEmax());
	    }

	}
	    //System.err.println(" flux= " + neutFlux.getDFDLogE(8.0,1));

	if(!isExposure){
	    System.err.print("Full IceCube configuration? [yes(1 - IC86 offline)/no(0 - IC86 onlie)] ->"); 
	    buffer   = d.readLine(); 
	    if(Integer.valueOf(buffer).intValue()==1) {
		areaTable = new NeutrinoEffAreaTable(true); // IC86
		System.err.println("Considering the full IC86 array");
	    }else{
		//areaTable = new NeutrinoEffAreaTable(false);// IC40
		areaTable = new NeutrinoEffAreaTable(false);// IC86 online alert area
		solidAngle = 2.0*Math.PI; // 2pi, upward-going only
		//livetime = 333.5*day_to_sec;   // 333.5 days IC40 livetime (exclud. burn sample)
		System.err.println("Considering the IC40 array");
	    }
	}else{
	    System.err.print("HESE + IC40+79+86? yes(1) no(0)->"); 
	    buffer   = d.readLine(); 
	    if(Integer.valueOf(buffer).intValue()==1) {
		addHESE = true;
		exposureTableHESE = new NeutrinoExposureTableHESE(true); // HESE + EHE
	    }
	    else{
		System.err.print("Full IceCube configuration? [yes(1 - IC86)/no(excluded IC40)] ->"); 
		buffer   = d.readLine(); 
		if(Integer.valueOf(buffer).intValue()==1) {
		    exposureTable = new NeutrinoExposureTable(true);
		    System.err.println("Calculate besed upon the IC40 + IC79-86 combined exposure table");
		}else{
		    exposureTable = new NeutrinoExposureTable(false);
		    System.err.println("Calculate besed upon the IC79-86 combined exposure table");
		}
	    }
	}
    }

    /**
       Default Constructor: Use the NeutrinoExposureTable class to calculate the flux with IC79 and IC86-1 combined
       Mode 1 : Use analtical model for GZK, no in-situ neutrino emission, i.e., pure GZK
       Mode 2 : Use analtical model for GZK + in-situ neutrino emission
       Mode 3 : Use numericaly calculated GZK model fluxes using NeutrinoFlux object.
    */
    public NeutrinoFluxIceCube(int mode) throws IOException {
	this(mode,true);
    }

    /**

       Setting the optical depth for in-situ neutrino emission
     */
    public void setOpticalDepth(double depth){
	if(addInSituEmission) neutFluxSource.setOpticalDepth(depth);
	else{
	    System.err.println("You are in the wrong mode to set the depth");
	    System.err.println("Nothing has been changed");
	}

    }

    /**
       Whether the neutrino flux assumes the constant evolution beyond zConst.
       The constant evolution option can be set via the constructor.
       This method enables any other object to tell if this option is set or not.
     */
    public boolean isConstantEvolution(){
	return assumeConstEvolution;
    }

    /**

       get the CRflux respnsible for in-situ neutrino emission
     */
    public double getCRDFDE(double zMax, double m, double alpha, double cosmicRayEnergy){
	if(addInSituEmission){
	    if(!assumeConstEvolution) return neutFluxSource.getCRDFDE(zMax,m,alpha,cosmicRayEnergy);
	    else return neutFluxSource.getCRDFDE(zMax,zConst,m,alpha,cosmicRayEnergy);
	}else{
	    System.err.println("You are in the wrong mode to set the depth");
	    System.err.println("Nothing has been changed");
	    return 0.0;
	}

    }

    /**

       set the CRflux reference energy
     */
    public void setCREnergyReference(double energy){
	if(addInSituEmission){
	    neutFluxSource.setCREnergyReference(energy);
	}else{
	    System.err.println("You are in the wrong mode to set the depth");
	    System.err.println("Nothing has been changed");
	}

    }
    /**

       set the CRflux respnsible for in-situ neutrino emission
     */
    public void setCRFluxAtReferenceEnergy(double zMax, double m, double alpha, double eFlux){
	if(addInSituEmission){
	    if(!assumeConstEvolution) neutFluxSource.setCRFluxAtReferenceEnergy(zMax, m, alpha, eFlux);
	    else neutFluxSource.setCRFluxAtReferenceEnergy(zMax, zConst, m, alpha, eFlux);
	}else{
	    System.err.println("You are in the wrong mode to set the depth");
	    System.err.println("Nothing has been changed");
	}

    }

    /**

       set the neutrino cutoff energy [GeV] for the simple power law model (model # 24)
     */
    public void setCutoffEnergy(double energy){
	if(!analyticalFunction){ 
	    if(neutFlux.isPowerLaw()){
		    double logEnergyMax = Math.pow(10.0,energy);
		    neutFlux.setLogEmax(logEnergyMax);
	    }
	}else{
	    System.err.println("You are in the wrong mode to set the maximal energy");
	    System.err.println("Nothing has been changed");
	}

    }

    /**

       set the neutrino energy flux [GeV/cm2 sec sr] for the simple power law model (model # 24)
     */
    public void setEflux(double eFlux){
	if(!analyticalFunction){ 
	    if(neutFlux.isPowerLaw()){
		    neutFlux.setEflux(eFlux);
	    }
	}else{
	    System.err.println("You are in the wrong mode to set the maximal energy");
	    System.err.println("Nothing has been changed");
	}

    }

    /** 
	number of events IceCube would see as a function of log10(NeutrinoEnergy[GeV]).
	All three neutrino flavors adds up.
     */
    public double getDNDLogE(double zMax, double m, double alpha, double logNeutrinoEnergy){
	double intensity = 0.0;
	if(analyticalFunction){ // anatitical function
	    double neutrinoEnergy = Math.pow(10.0,logNeutrinoEnergy);
	    if(!assumeConstEvolution){
		intensity = neutFluxFunction.getDFDE(zMax,m,alpha,neutrinoEnergy);
		if(addInSituEmission)  intensity += neutFluxSource.getDFDE(zMax,m,alpha,neutrinoEnergy);
	    }else{
		intensity = neutFluxFunction.getDFDE(zMax,zConst,m,alpha,neutrinoEnergy);
		if(addInSituEmission)  intensity += neutFluxSource.getDFDE(zMax,zConst,m,alpha,neutrinoEnergy);
	    }
	    // Now IceCube neutrino effective area - all three flavors add up
	    double exposure;
	    if(!isExposure){
		exposure = (areaTable.getArea(0,logNeutrinoEnergy)+areaTable.getArea(1,logNeutrinoEnergy)
			    +areaTable.getArea(2,logNeutrinoEnergy))*solidAngle*livetime;
	    }else{
		if(!addHESE){
		    exposure = exposureTable.getExposure(0,logNeutrinoEnergy)+exposureTable.getExposure(1,logNeutrinoEnergy)
			+exposureTable.getExposure(2,logNeutrinoEnergy);
		}else{
		    exposure = exposureTableHESE.getExposure(0,logNeutrinoEnergy)+exposureTableHESE.getExposure(1,logNeutrinoEnergy)
			+exposureTableHESE.getExposure(2,logNeutrinoEnergy);
		}
	    }
	    intensity = (1.0/3.0)*intensity*neutrinoEnergy*ln10*exposure;

	}else{ // numerical model. Doesn't matter zMax and m. 
	    for(int pID = 1; pID <=3; pID++){ // add all three neutrino flavors
		if(!isExposure){
		    intensity += neutFlux.getDFDLogEwzOsci(logNeutrinoEnergy,pID)*areaTable.getArea(pID-1,logNeutrinoEnergy)*
		    solidAngle*livetime;
		}else{
		    if(!addHESE){
			intensity += neutFlux.getDFDLogEwzOsci(logNeutrinoEnergy,pID)*exposureTable.getExposure(pID-1,logNeutrinoEnergy);
		    }else{
			intensity += neutFlux.getDFDLogEwzOsci(logNeutrinoEnergy,pID)*exposureTableHESE.getExposure(pID-1,logNeutrinoEnergy);
		    }
		}
	    }
	    //System.err.format("   logE(%f) intensity(%e)\n",logNeutrinoEnergy,intensity);
	}

	return intensity;
    }

    /** 
	method for number of events in the space of log(neutrino Energy [GeV])
	<pre>
	parameters[0]  :: alpha cosmic ray spectrum power law index
	parameters[1]  :: zMax  maximum redshift
	parameters[2]  :: m  source evolution parameter
	</pre>

     */
    public double getFunction(int functionIndex, double[] parameters, 
			      double x){
	double alpha = parameters[0];
	double zMax = parameters[1];
	double m = parameters[2];
	return getDNDLogE(zMax,m,alpha,x);
    }



    /** A simple main method for test */
    public static void main(String[] args) throws IOException{

	NeutrinoFluxIceCube flux;
	int mode = 0;
	double m = 2.0; // evolution
	double zMax = 2.0;
	double alpha = 2.5;
	double depth = 0.0;
	double[] parameters = new double[5];

	if(args.length!=1){
	    System.out.println("Usage:NeutrinoFluxIceCube mode(1 analytical GZK 2 plus in-situ 3 numerical model)");
	    System.exit(0);
	}else{
	    mode = Integer.valueOf(args[0]).intValue();
	}

	flux = new NeutrinoFluxIceCube(mode);

	// set up the parameters
	DataInputStream input = new DataInputStream(System.in); 
	BufferedReader  d     = new BufferedReader(new InputStreamReader(input)); 
	String buffer; 
	if(mode!=3){
	    System.err.print("source evolution m ->"); 
	    buffer   = d.readLine(); 
	    m = Double.valueOf(buffer).doubleValue();
	    System.err.print("Zmax ->"); 
	    buffer   = d.readLine(); 
	    zMax = Double.valueOf(buffer).doubleValue();
	    System.err.print("alpha ->"); 
	    buffer   = d.readLine(); 
	    alpha = Double.valueOf(buffer).doubleValue();

	    System.err.format(" m(%f) Zmax(%f) alpha(%f)\n",m,zMax,alpha);

	    if(mode==2){
		System.err.print("optical depth ->"); 
		buffer   = d.readLine(); 
		depth = Double.valueOf(buffer).doubleValue();
		System.err.format(" depth(%f)\n",depth);
		flux.setOpticalDepth(depth);
	    }

	   parameters[0] = alpha;
	   parameters[1] = zMax;
	   parameters[2] = m;

	}


	double logEmin = 3.0; // 10^3 GeV
	double logEmax = 11.0;// 10^11 GeV
	double numberOfEvents = Integration.RombergIntegral(flux, 0, parameters, logEmin,logEmax);
	System.out.format("Number of Events seen in IceCube %f\n",numberOfEvents);

    }  
}
