package iceCube.uhe.neutrinoModel;

import iceCube.uhe.particles.*;
import iceCube.uhe.neutrinoModel.*;
import java.io.*;


public class MakeNeutrinoFluxTable {

    public static void main(String[] args) throws IOException {

	NeutrinoFlux neutFlux = null;
	int modelNumber = 4;
	double binWidth = 0.05; // 0.05 decade
	if(args.length<1){
	    System.out.println("Usage: MakeNeutrinoFluxTable model-name");
	    System.exit(0);
	}else{
	    modelNumber= Integer.valueOf(args[0]).intValue();
	}

	neutFlux = new NeutrinoFlux(modelNumber);

	double logEnergy = Particle.getLogEnergyMinimum();
	double logMaxEnergy = Particle.getLogEnergyMinimum()
	    + Particle.getDeltaLogEnergy()*(double )(Particle.getDimensionOfLogEnergyMatrix());


	while(logEnergy<logMaxEnergy){
	    double cummulativeDeltaLogE = 0.0;

	    double[] flux = new double[3];
	    for(int iP=0;iP<3;iP++) flux[iP]=0.0;

	    while(cummulativeDeltaLogE < binWidth){

		// nu-e
		flux[0] += neutFlux.getDFDLogEwzOsci(logEnergy+cummulativeDeltaLogE,1); 
		// nu-mu
		flux[1] += neutFlux.getDFDLogEwzOsci(logEnergy+cummulativeDeltaLogE,2); 
		// nu-tau
		flux[2] += neutFlux.getDFDLogEwzOsci(logEnergy+cummulativeDeltaLogE,3); 

		cummulativeDeltaLogE +=  Particle.getDeltaLogEnergy();
	    }

	    System.out.format("%5.2f %11.5e %11.5e %11.5e\n",
			      logEnergy,flux[0],flux[1],flux[2]);

	    logEnergy += binWidth;
	}
	    
    }
}

