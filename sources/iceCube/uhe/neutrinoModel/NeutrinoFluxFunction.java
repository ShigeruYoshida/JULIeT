package iceCube.uhe.neutrinoModel;

import numRecipes.*;
import iceCube.uhe.particles.*;

import java.io.*;

/** Calculate the cosmic neutrino flux by the analytical functions.
    The formula is based on
<UL>
  <DT> <cite>Constraints on the origin of the ultrahigh energy cosmic rays using cosmic diffuse neutrino flux limits: An analytical approach</cite>
       <TABLE CELLPADDING=2>
       <TR>
         <TD ALIGN="LEFT"><it> Shigeru Yoshida and Aya Ishihara</it></TD>
         <TD ALIGN="LEFT">
            <a href="http://link.aps.org/doi/10.1103/PhysRevD.85.063002">
            Physical Review D <b>85</b> 063002 (2012)</a></TD>
       </TR>
       </TABLE><BR>
</UL>


*/
public class NeutrinoFluxFunction implements Function{

    /** Speed of light */
    public final static double c = 2.99792e10; // [cm/sec]

    /** Mpc unit */
    public final static double Mpc = 3.08568025e24;   // [cm]
    /** km unit */
    public final static double Km = 1.0e5;        // [cm]
    /** The plank const */
    public final static double hBar = 6.5829507e-25;  // [GeV sec]
    /** The boltzmann const */
    public final static double kB = 8.617385e-14;   // [GeV/K]


    /** constants of the cosmology */
    public static double omegaM = 0.3;        // matter density 
    public static double omegaLambda = 0.7;   // density originated in the cosmological constant
    public static double H0 =73.5*Km/Mpc;     // The habble constatnt [km/sec/Mpc] estimated by WMAP
    public static double TCMB = 2.725;        // The CMB black-body temperature [K]
    public static double ECMB = kB*TCMB;      // The CMB thermal energy [GeV]

    /** constatnts of the particle (physics) */
    public final static double Mp      = 938.272e-3; // Proton mass [GeV]
    public final static double Mpi     = Particle.particleMasses[3][1]; // Pi+ mass [GeV]
    public final static double Mmuon   = Particle.particleMasses[1][1]; // muon mass [GeV]
    public final static double rPi     = Mmuon*Mmuon/(Mpi*Mpi);         // mu/pi mass ratio
    public final static double sigmaRes= 2.1e-28;   // p(gamma pi+)n production at the Delta resonance [cm^2]
    public final static double sRes    = 1.47;      // The CMS energy square at the Delta resonance [GeV^2]
    public final static double deltaRes= 0.6;       // The width of the Delta resonance [GeV^2]
    public final static double piCMSEnergy = (sRes+Mpi*Mpi-Mp*Mp)/(2.0*sRes); // relative pi energy at CMS
    public final static double piCMSMomentum = Math.sqrt(piCMSEnergy*piCMSEnergy-Mpi*Mpi); // momentum at CMS
    public final static double piEnergyPlus = piCMSEnergy + piCMSMomentum;
    public final static double piEnergyMinus = piCMSEnergy - piCMSMomentum;
    static double logFactor = Math.log(piEnergyPlus/piEnergyMinus);
    //static double logFactor = Math.log(0.379/0.03);

    /** constants related to the ultra-high energy cosmic-ray propagation */
    static double gzkSphere = 100.0*Mpc;
    static double gzkRatio = gzkSphere*H0/c;

    /** parameters regarding the ultrah-high energy cosmic ray flux */
    //double alpha = 2.5;   // the power law index of the cosmic ray spectrum
    static double gzkEnergyThreshold = 5.6e10;   // "threshold energy" for the GZK process [GeV] 10^19.75 eV
    static double gzkIntegralFluxBase = 2.25e-20; // Flux at the GZK threshold [/cm^2 sec sr] by HiRes

    /** numerical constatnts for the integral approximation */
    static double x0 = 0.275;
    //static double x0 = 0.255;
    static double x1 = 0.16;
    //static double x1 = 0.12;
    static double eMinus2 = Math.exp(-2.0);

    /** The defult constructor */
    public NeutrinoFluxFunction(){
	System.err.format("Mp(%6.3e) Mpi(%6.3e) piEnergyPlus(%7.5f) piEnergyMinus(%7.5f) rPi(%4.3f)\n",
			  Mp,Mpi,piEnergyPlus,piEnergyMinus,rPi);
	System.err.format("cosmic distance (%7.4e Mpc) gzkRatio(%7.5f)\n",(c/H0)/Mpc, gzkRatio);

    }

    /** Return the UHECR flux [/cm^2 sec sr] above the GZK cutoff energy*/
    public static double getGZKCRFlux(double cosmicrayEnergy){
	if(cosmicrayEnergy<=gzkEnergyThreshold) return gzkIntegralFluxBase;
	else{
	    double flux = gzkIntegralFluxBase*Math.pow((cosmicrayEnergy/gzkEnergyThreshold),-3.5);
	    return flux;
	}
    }

    protected  static double getDistance(double z){
	double x = 1.0+z;
	double distance = Math.sqrt(omegaM*x*x*x+omegaLambda);
	return distance;
    }

    private static double getExpDumpFactor(double zFrom, double zTo, double alpha){
	double distanceFrom = getDistance(zFrom);
	double distanceTo = getDistance(zTo);
	double expTerm = Math.exp(-(alpha+1.0)*(1.0/gzkRatio)*2.0/(3.0*omegaM)*(distanceFrom-distanceTo));
	return expTerm;
    }

    private static double getExpDumpFactor(double z, double alpha){
	double distance = getDistance(z);
	double expTerm = Math.exp(-(alpha+1.0)*(1.0/gzkRatio)*2.0/(3.0*omegaM)*distance);
	return expTerm;
    }

    private static double getLogEnergyTerm(double z, double neutrinoEnergy){
	double b = getEnergyBaseTerm(z,neutrinoEnergy);
	return Math.log(x1/b);
    }

    private static double getEnergyBaseTerm(double z, double neutrinoEnergy){
	double b= 4.0*ECMB*neutrinoEnergy*(1.0+z)*(1.0+z)/((sRes-Mp*Mp)*piEnergyMinus*(1.0-rPi));
	//double b = beta/energyA*neutrinoEnergy*(1.0+z)*(1.0+z);
	return b;
    }

    private static double getZLowLimit(double zMax,double neutrinoEnergy){
	double xLimit = Math.sqrt(x1/getEnergyBaseTerm(0.0,neutrinoEnergy));
	double zLimit = xLimit-1.0;
	if(zLimit>zMax) zLimit = zMax;
	if(zLimit<0.0) zLimit = 0.0;
	return zLimit;
    }

    private static double getZLowLimit(double zMax,double zMin, double neutrinoEnergy){
	double xLimit = Math.sqrt(x1/getEnergyBaseTerm(0.0,neutrinoEnergy));
	double zLimit = xLimit-1.0;
	if(zLimit>zMax) zLimit = zMax;
	if(zLimit<zMin) zLimit = zMin;
	return zLimit;
    }

    private static double getZUpLimit(double zMax,double neutrinoEnergy){
	double b= ((sRes-Mp*Mp)*piEnergyPlus*(1.0-rPi))/(4.0*ECMB*neutrinoEnergy);
	double xLimit = Math.sqrt(x1*b);
	double zLimit = xLimit-1.0;
	if(zLimit>zMax) zLimit = zMax;
	if(zLimit<0.0) zLimit = 0.0;
	return zLimit;
    }

    private static double getZUpLimit(double zMax,double zMin,double neutrinoEnergy){
	double b= ((sRes-Mp*Mp)*piEnergyPlus*(1.0-rPi))/(4.0*ECMB*neutrinoEnergy);
	double xLimit = Math.sqrt(x1*b);
	double zLimit = xLimit-1.0;
	if(zLimit>zMax) zLimit = zMax;
	if(zLimit<zMin) zLimit = zMin;
	return zLimit;
    }

    private static double getEpsilon0(double z, double neutrinoEnergy){
	double b= x1*((sRes-Mp*Mp)*piEnergyPlus*(1.0-rPi))/(4.0*ECMB*neutrinoEnergy*(1.0+z)*(1.0+z));
	if(b>=1.0) return 1.0;
	else return 0.0;
    }

    private static double getEpsilon1(double z, double neutrinoEnergy){
	double bPlus= x1*((sRes-Mp*Mp)*piEnergyPlus*(1.0-rPi))/(4.0*ECMB*neutrinoEnergy*(1.0+z)*(1.0+z));
	double bMinus= x1*((sRes-Mp*Mp)*piEnergyMinus*(1.0-rPi))/(4.0*ECMB*neutrinoEnergy*(1.0+z)*(1.0+z));
	if(bMinus>=1.0) return 0.0;
	else if(bPlus>=1.0) return 1.0;
	else return 0.0;
    }

    public double getFunction(int functionIndex, double[] parameters, 
			      double x){
	if(functionIndex == 0){
	    double alpha = parameters[0];
	    double zMax = parameters[1];
	    double m = parameters[2];
	    return getDFDE(zMax,m,alpha,x);
	}else{
	    double alpha = parameters[0];
	    double zMax = parameters[1];
	    double zConst = parameters[2];
	    double m = parameters[3];
	    //System.err.format(" alpha=%f zMax=%f zConst=%f energy=%f\n",alpha,zMax,zConst,x);
	    return getDFDE(zMax,zConst,m,alpha,x);
	}

    }

    public double getDFDE(double zMax, double m, double alpha, double neutrinoEnergy){

	double gzkEnergy = 1.0e11;   // 10^20 eV = 10^11 GeV
	//double gzkEnergy = gzkEnergyThreshold;
	//System.err.format(" gzk flux=%e\n",getGZKCRFlux(gzkEnergy));
	double crossSectionTerm = sigmaRes/(piEnergyPlus-piEnergyMinus)*deltaRes;
	double resonanceEnergy = (sRes-Mp*Mp)/(4.0*ECMB);
	double resonanceEnergyTerm = Math.pow(resonanceEnergy,-alpha-2.0)*(sRes-Mp*Mp);
	double neutrinoIntensityTerm = 3.0/(1.0-rPi);
	double cosmicRadius = c/H0;
	double cmbDensityTerm = ECMB/(8.0*Math.PI*Math.PI*hBar*hBar*hBar*c*c*c);
	double zTerm = getRedShiftTerm(zMax,m,alpha,neutrinoEnergy);
	//System.err.format(" resE=%7.4e\n",resonanceEnergy);

	double flux = (alpha-1.0)*Math.pow(gzkEnergy,alpha-1.0)*getGZKCRFlux(gzkEnergy)*cosmicRadius*cmbDensityTerm*resonanceEnergyTerm*crossSectionTerm*neutrinoIntensityTerm*zTerm;

	return flux;

    }

    public double getDFDE(double zMax, double zConst, double m, double alpha, double neutrinoEnergy){
	if(zMax<= zConst) return getDFDE(zMax,m,alpha,neutrinoEnergy);

	double fluxUptoZConst = getDFDE(zConst,m,alpha,neutrinoEnergy);
	if(fluxUptoZConst<=0.0) return 0.0;
	double zTerm = getRedShiftTerm(zConst,m,alpha,neutrinoEnergy);
	double additionalzTerm = getRedShiftTerm(zMax,zConst,m,alpha,neutrinoEnergy);
	double flux = fluxUptoZConst + (additionalzTerm/zTerm)*fluxUptoZConst;

	return flux;
    }

    public static double getRedShiftTerm(double zMax, double m, double alpha, double neutrinoEnergy){

	double numericalTerm = Math.pow(x0,-alpha-1.0)*eMinus2;
	double numericalTermEdge = Math.pow(x1,-alpha-3.0)*Math.exp(-1.0/x1)*eMinus2;

	double powerIndexTerm = 2.0/(2.0*(alpha+m)-3.0);
	//double powerIndexTerm = 2.0/(2.0*(alpha-1.0+m)-3.0);

	double distanceMax = getDistance(zMax);
	double expDumpMaxToPresent = getExpDumpFactor(zMax,0.0,alpha);

	double zUp = getZUpLimit(zMax,neutrinoEnergy);
	double kinematicsTermUp = x1/getEnergyBaseTerm(zUp,neutrinoEnergy);
	double logKinematicsTermUp = Math.log(kinematicsTermUp);
	double distanceUp = getDistance(zUp);
	double expDumpMaxToUp = getExpDumpFactor(zMax,zUp,alpha);

	double zMin = getZLowLimit(zMax,neutrinoEnergy);
	double kinematicsTermMin = x1/getEnergyBaseTerm(zMin,neutrinoEnergy);
	double logKinematicsTermMin = Math.log(kinematicsTermMin);
	double distanceMin = getDistance(zMin);
	double expDumpMaxToMin = getExpDumpFactor(zMax,zMin,alpha);

	double distance0 = getDistance(0.0);

	double distancePowerIndex = 2.0*(alpha+m)/3.0-1.0;
	//double distancePowerIndex = 2.0*(alpha-1.0+m)/3.0-1.0;
	double massDensityTerm = Math.pow(omegaM,-(alpha+m)/3.0);
	//double massDensityTerm = Math.pow(omegaM,-(alpha-1.0+m)/3.0);

	double term1 = powerIndexTerm*massDensityTerm
	    *(Math.pow(distanceUp,distancePowerIndex)*(numericalTerm*logFactor+numericalTermEdge*(logKinematicsTermUp+2.0*powerIndexTerm))
	      -Math.pow(distance0,distancePowerIndex)*numericalTerm*logFactor
	      -Math.pow(distanceMin,distancePowerIndex)*numericalTermEdge*(logKinematicsTermMin+2.0*powerIndexTerm));

	//double term2 = -1.0/(alpha+1.0)*gzkRatio*Math.pow(1.0+zMax,m-3.0)*Math.pow(1.0+zUp,alpha)*expDumpMaxToUp*(numericalTerm*logFactor+numericalTermEdge*logKinematicsTermUp);
	//double term3 = 1.0/(alpha+1.0)*gzkRatio*Math.pow(1.0+zMax,m-3.0)
	//    *(numericalTerm*logFactor*expDumpMaxToPresent+Math.pow(1.0+zMin,alpha)*numericalTermEdge*logKinematicsTermMin*expDumpMaxToMin);

	//System.err.format(" term1=%7.4e term2=%8.4e term3=%8.4e\n",term1,term2,term3);
	return term1;

    }


    /** 
	additional contribution when the evolution goes const when redshift is beyond zConst
     */
    public static double getRedShiftTerm(double zMax, double zConst, double m, double alpha, double neutrinoEnergy){

	double numericalTerm = Math.pow(x0,-alpha-1.0)*eMinus2;
	double numericalTermEdge = Math.pow(x1,-alpha-3.0)*Math.exp(-1.0/x1)*eMinus2;

	double powerIndexTerm = 2.0/(2.0*alpha-3.0);

	double zUp = getZUpLimit(zMax,zConst,neutrinoEnergy);
	double kinematicsTermUp = x1/getEnergyBaseTerm(zUp,neutrinoEnergy);
	double logKinematicsTermUp = Math.log(kinematicsTermUp);
	double distanceUp = getDistance(zUp);

	double zMin = getZLowLimit(zMax,zConst,neutrinoEnergy);
	double kinematicsTermMin = x1/getEnergyBaseTerm(zMin,neutrinoEnergy);
	double logKinematicsTermMin = Math.log(kinematicsTermMin);
	double distanceMin = getDistance(zMin);

	double distanceConst = getDistance(zConst);

	double distancePowerIndex = 2.0*alpha/3.0-1.0;
	double massDensityTerm = Math.pow(omegaM,-alpha/3.0);

	double term1 = powerIndexTerm*massDensityTerm
	    *(Math.pow(distanceUp,distancePowerIndex)*(numericalTerm*logFactor+numericalTermEdge*(logKinematicsTermUp+2.0*powerIndexTerm))
	      -Math.pow(distanceConst,distancePowerIndex)*numericalTerm*logFactor
	      -Math.pow(distanceMin,distancePowerIndex)*numericalTermEdge*(logKinematicsTermMin+2.0*powerIndexTerm));

	return term1*Math.pow((1.0+zConst),m);

    }

    public static double getNeutrinoYieldFromAsource(double zSource, double zNeutrino, double alpha, double neutrinoEnergy){
	double gzkEnergy = 1.0e11;   // 10^20 eV = 10^11 GeV
	//double gzkEnergy = gzkEnergyThreshold;
	//System.err.format(" gzk flux=%e\n",getGZKCRFlux(gzkEnergy));
	double crossSectionTerm = sigmaRes/(piEnergyPlus-piEnergyMinus)*deltaRes;
	double resonanceEnergy = (sRes-Mp*Mp)/(4.0*ECMB);
	double resonanceEnergyTerm = Math.pow(resonanceEnergy,-alpha-2.0)*(sRes-Mp*Mp);
	double neutrinoIntensityTerm = 3.0/(1.0-rPi);
	double cosmicRadius = c/H0;
	double cmbDensityTerm = ECMB/(8.0*Math.PI*Math.PI*hBar*hBar*hBar*c*c*c);

	double distance = getDistance(zNeutrino);
	double expDump = getExpDumpFactor(zSource,zNeutrino,alpha);
	double zTerm = Math.pow((1.0+zNeutrino),alpha+3.0)*expDump/(distance*(1.0+zNeutrino));

	double numericalTerm = Math.pow(x0,-alpha-1.0)*eMinus2;
	double numericalTermEdge = Math.pow(x1,-alpha-3.0)*Math.exp(-1.0/x1)*eMinus2;
	double kinematicsTerm = x1/getEnergyBaseTerm(zNeutrino,neutrinoEnergy);
	double logKinematicsTerm = Math.log(kinematicsTerm);

	double energyTerm = getEpsilon0(zNeutrino, neutrinoEnergy)*numericalTerm*logFactor+
	    getEpsilon1(zNeutrino,neutrinoEnergy)*numericalTermEdge*logKinematicsTerm;

	double yield = (alpha-1.0)*(alpha-1.0)*Math.pow(gzkEnergy,alpha-1.0)*getGZKCRFlux(gzkEnergy)*cosmicRadius*cmbDensityTerm*resonanceEnergyTerm*crossSectionTerm*neutrinoIntensityTerm*zTerm*energyTerm;

	return yield;
    }


    public static void main(String[] args){

	if(args.length < 4) {
	    System.out.println("Usage: NeutrinoFluxFunction zMax powerIndex m neutrinoEnergy");
	    System.exit(0);
	}
	double zConst = 2.0;
	double zMax = 1.0;
	zMax = Double.valueOf(args[0]).doubleValue();
	double alpha = 3.0;
	alpha = Double.valueOf(args[1]).doubleValue();
	double m = 2.0;
	m = Double.valueOf(args[2]).doubleValue();
	double neutrinoEnergy = 1.0e8;
	neutrinoEnergy = Double.valueOf(args[3]).doubleValue();
	double zLimit = getZLowLimit(zMax,neutrinoEnergy);
	double zUpLimit = getZUpLimit(zMax,neutrinoEnergy);
	System.out.format(" alpha(%5.2f) zMax(%5.2f) zLowLimit(%5.2f) zUpLimit(%5.2f) m(%5.2f) neutrinoEnergy(%7.3e)\n",alpha,zMax,zLimit,zUpLimit,m, neutrinoEnergy);

	double[] parameters = new double[5];
	parameters[0] = alpha;
	parameters[1] = zMax;
	parameters[2] = m;

	NeutrinoFluxFunction neutrinoFlux = new NeutrinoFluxFunction( );
	System.out.format(" analitical Integral= %10.7f\n",getRedShiftTerm(zMax,m, alpha,neutrinoEnergy));
	if(zMax>zConst) 
	    System.out.format(" additional analitical Integral= %10.7f\n",getRedShiftTerm(zMax,zConst,m, alpha,neutrinoEnergy));

	double energyFlux = neutrinoFlux.getDFDE(zMax,m, alpha,neutrinoEnergy)*neutrinoEnergy*neutrinoEnergy;
	System.out.format(" E^2dF/dE= %10.7e [GeV /cm^2 sec sr]\n",energyFlux);
	if(zMax>zConst){
	    energyFlux = neutrinoFlux.getDFDE(zMax,zConst,m, alpha,neutrinoEnergy)*neutrinoEnergy*neutrinoEnergy;
	    System.out.format(" E^2dF/dE= %10.7e [GeV /cm^2 sec sr] for zConst=%f\n",energyFlux,zConst);
	}

	// integral flux
	double maxNeutrinoEnergy = 1.0e11;   // 1.0e11 GeV
	double integralFlux = Integration.RombergIntegral(neutrinoFlux, 
					  0, parameters, neutrinoEnergy, maxNeutrinoEnergy);
	System.out.format(" F(>%7.3e GeV) = %10.7e [/cm^2 sec sr]\n",neutrinoEnergy,integralFlux);
	if(zMax>zConst){
	    parameters[2] = zConst;
	    parameters[3] = m;
	    integralFlux = Integration.RombergIntegral(neutrinoFlux, 
					  1, parameters, neutrinoEnergy, maxNeutrinoEnergy);
	    System.out.format(" F(>%7.3e GeV) = %10.7e [/cm^2 sec sr] for zConst=%f\n",neutrinoEnergy,integralFlux,zConst);
	}
    }

}

