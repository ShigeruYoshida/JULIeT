package iceCube.uhe.neutrinoModel;

import numRecipes.*;
import iceCube.uhe.particles.*;

import java.io.*;

/** Calculate the GZK neutrino flux by the energetics argument */

public class NeutrinoFluxFunctionEnergetics extends NeutrinoFluxFunction implements Function{

    public static final double R_nu = 0.45;  // Fraction of proton energy channeling into neutrinos
    public static final double energy_scale = 2.0e11; // GZK energy scale 2.0x10^11 GeV
    public static final double xTh = 1.0/3.0; 
    public static final double xMax = 5.0; 
    public static final double gzkMinimumEnergy = energy_scale*xTh; // 6.0x10^10 GeV
    public static final double gzkMaximumEnergy = energy_scale*xMax; // 1.0x10^11 GeV, maximal energy at source

    public NeutrinoFluxFunctionEnergetics(){
	super();
    }

    public double getFunction(int functionIndex, double[] parameters, 
			      double x){
	if(functionIndex == 0){
	    double alpha = parameters[0];
	    double zMax = parameters[1];
	    double m = parameters[2];
	    return getDFDE(zMax,m,alpha,x);
	}else if (functionIndex < 2){
	    double alpha = parameters[0];
	    double zMax = parameters[1];
	    double zConst = parameters[2];
	    double m = parameters[3];
	    return getDFDE(zMax,zConst,m,alpha,x);
	}else{
	    double m = parameters[0];
	    double neutrinoEnergy = parameters[1];
	    return getZFunction(m,x,neutrinoEnergy); // integral over redshift z
	}

    }

    public double getDFDE(double zMax, double m, double alpha, double neutrinoEnergy){
	double gzkEnergy = 1.0e11;   // 10^20 eV = 10^11 GeV
	//double gzkEnergy = gzkEnergyThreshold;
	double crossSectionTerm = 1.0/(piEnergyPlus-piEnergyMinus);
	double neutrinoIntensityTerm = 1.0/(1.0-rPi);
	double[] parameters = new double[2];
	parameters[0] = m;
	parameters[1] = neutrinoEnergy;
	double zTerm = Integration.RombergIntegral(this, 2, parameters, 0.0, zMax);

	double eflux = getGZKCRFlux(gzkEnergy)/(gzkRatio)*(gzkEnergy/energy_scale)*
	    crossSectionTerm*neutrinoIntensityTerm*R_nu*logFactor*zTerm;

	return eflux/neutrinoEnergy;

    }

    public double getDFDE(double zMax, double zConst, double m, double alpha, double neutrinoEnergy){
	if(zMax<= zConst) return getDFDE(zMax,m,alpha,neutrinoEnergy);

	double fluxUptoZConst = getDFDE(zConst,m,alpha,neutrinoEnergy);
	if(fluxUptoZConst<=0.0) return 0.0;
	double[] parameters = new double[2];
	parameters[0] = m;
	parameters[1] = neutrinoEnergy;
	double zTerm = Integration.RombergIntegral(this, 2, parameters, 0.0, zConst);

	parameters[0] = 0.0;
	double additionalzTerm = Integration.RombergIntegral(this,
				   2, parameters, zConst, zMax)*Math.pow((1.0+zConst),m);
	double flux = fluxUptoZConst + (additionalzTerm/zTerm)*fluxUptoZConst;

	return flux;
    }

    private static double getZLowLimit(double zMax, double neutrinoEnergy){
	double zLimit = getZLowLimit(zMax,0.0,neutrinoEnergy);
	return zLimit;
    }

    private static double getZLowLimit(double zMax, double zMin, double neutrinoEnergy){
	double xLimit = Math.sqrt(gzkMinimumEnergy*piEnergyMinus*(1.0-rPi)/neutrinoEnergy);
 	double zLimit = xLimit-1.0;
	if(zLimit>zMax) zLimit = zMax;
	if(zLimit<zMin) zLimit = zMin;
	return zLimit;
    }

    private static double getZUpLimit(double zMax, double neutrinoEnergy){
	double zLimit = getZUpLimit(zMax,0.0,neutrinoEnergy);
	return zLimit;
    }

    private static double getZUpLimit(double zMax, double zMin, double neutrinoEnergy){
	double xLimit = Math.sqrt(gzkMinimumEnergy*piEnergyPlus*(1.0-rPi)/neutrinoEnergy);
 	double zLimit = xLimit-1.0;
	if(zLimit>zMax) zLimit = zMax;
	if(zLimit<zMin) zLimit = zMin;
	return zLimit;
    }

    private static double getZFunction(double m, double z, double neutrinoEnergy){
	double sourceTerm = Math.pow(1.0+z,m+1)/getDistance(z);
	double xPlus = neutrinoEnergy*(1.0+z)*(1.0+z)/(energy_scale*piEnergyPlus*(1.0-rPi));
	double xMinus =neutrinoEnergy*(1.0+z)*(1.0+z)/(energy_scale*piEnergyMinus*(1.0-rPi));
	double tanMaxTerm = Math.atan(xMax*(1.0+z));
	double tanThTerm = Math.atan(xTh);
	double tanTerm = 0.0;
	if(xMinus<= xTh){
	    tanTerm = tanMaxTerm - tanThTerm;
	    //System.out.println(" category 0");
	}else if(xPlus <= xTh){
	    tanTerm = tanMaxTerm - 0.5*(Math.atan(xMinus) + tanThTerm);
	    //System.out.println(" category 1");
	}else if(xMinus <= xMax*(1.0+z)){
	    tanTerm = tanMaxTerm - 0.5*(Math.atan(xMinus) + Math.atan(xPlus));
	    //System.out.println(" category 2");
	}else if(xPlus <= xMax*(1.0+z)){
	    tanTerm = tanMaxTerm - 0.5*(tanMaxTerm + Math.atan(xPlus));
	    //System.out.println(" category 3");
	}else{
	    tanTerm = 0.0;
	    //System.out.println(" category 4");
	}

	double zFunc = sourceTerm*tanTerm;
	//System.err.format("E=%8.3e z=%6.3f m=%6.3f zfunc=%10.7e\n",neutrinoEnergy,z,m,zFunc);
	return zFunc;
    }

    public static void main(String[] args){
	if(args.length < 3) {
	    System.out.println("Usage: NeutrinoFluxFromSource zMax m neutrinoEnergy");
	    System.exit(0);
	}
	double zMax = 1.0;
	zMax = Double.valueOf(args[0]).doubleValue();
	double m = 2.0;
	m = Double.valueOf(args[1]).doubleValue();
	double neutrinoEnergy = 1.0e8;
	neutrinoEnergy = Double.valueOf(args[2]).doubleValue();

	double alpha = 2.0;
	NeutrinoFluxFunctionEnergetics neutrinoFlux = new NeutrinoFluxFunctionEnergetics( );

	// check zfunction
	//double logE = 7.0;
	//while(logE<=11.0){
	//    double energy = Math.pow(10.0,logE);
	//    double zFunc = neutrinoFlux.getZFunction(m,zMax,energy);
	//    System.out.format(" logE=%5.2f m=%5.3f zMax = %6.3f zFunc = %e\n",logE,m,zMax,zFunc);
	//    logE += 0.1;
	//}
	//System.exit(0);

	double[] parameters = new double[5];
	parameters[0] = alpha;
	parameters[1] = zMax;
	parameters[2] = m;

	double energyFlux = neutrinoFlux.getDFDE(zMax,m, alpha,neutrinoEnergy)*neutrinoEnergy*neutrinoEnergy;
	System.out.format(" E^2dF/dE= %10.7e [GeV /cm^2 sec sr]\n",energyFlux);

	// integral flux
	double maxNeutrinoEnergy = 1.0e11;   // 1.0e11 GeV
	double integralFlux = Integration.RombergIntegral(neutrinoFlux, 
			  0, parameters, neutrinoEnergy, maxNeutrinoEnergy);
	System.out.format(" F(>%7.3e GeV) = %10.7e [/cm^2 sec sr]\n",neutrinoEnergy,integralFlux);

	// Numerical integral for check
	double energy = neutrinoEnergy;
	double dE = 0.1*neutrinoEnergy;
	integralFlux = 0.0;
	while(energy<=maxNeutrinoEnergy){
	    integralFlux += neutrinoFlux.getDFDE(zMax,m, alpha, energy)*dE;
	    energy += dE;
	}
	System.out.format(" F(>%7.3e GeV) = %10.7e [/cm^2 sec sr]\n",neutrinoEnergy,integralFlux);
    }

}
