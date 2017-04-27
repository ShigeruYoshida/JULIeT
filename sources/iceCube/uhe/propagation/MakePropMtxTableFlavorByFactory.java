package  iceCube.uhe.propagation;

import iceCube.uhe.particles.*;
import iceCube.uhe.propagation.*;
import java.io.*;

/**

This class is to generate propagation matrix table in F2K format
for offline IceCube analysis. By using PropagationMatrixFactory,
it writes (nu-e to iniceParticle) +(nu-mu to iniceParticle) +
(nu-tau to iniceParticle). They constitute a neutrino yield
in calculating the sensitivity. Refer to I3ParticleFlux class
in the analysis package to see how these matrix elements
contribute to the yield calculation.

This class has only a main method.

Written by S. Yoshida Feb 5th 2009.

*/

public class MakePropMtxTableFlavorByFactory {

    public static void main(String[] args) throws IOException {


	int iniceFlavor = 1;
	int iniceDoublet = 1;
	String propMtxFileName = "1e5_00cm.data";
	double binWidth = 0.05; // 0.02 decade 
        if(args.length<3){
            System.out.println(
		       "Usage: MakePropMtxTableByFactory inice-flavor inice-doublet propMatrixFileName");
	    System.exit(0);
	}else{
	    iniceFlavor = Integer.valueOf(args[0]).intValue();
	    iniceDoublet = Integer.valueOf(args[1]).intValue();
	    propMtxFileName = args[2];
	}

	String iniceParticleName = Particle.particleName(iniceFlavor,iniceDoublet);

	System.err.println("now generating matrix table of " + iniceParticleName);
	System.err.println("input matrix file name " + propMtxFileName);

	// PropagatingMatrixFactory object
	PropagationMatrixFactory matrix = new PropagationMatrixFactory();
	DataInputStream in = 
	    new DataInputStream(new FileInputStream(propMtxFileName));
	matrix.readMatrix(in);


	//
	// Loop over the propagation matrix
	//

	double epsilon = 1.0e-4; // round-off margin for binning
	double logInIceEnergy = Particle.getLogEnergyMinimum() + epsilon;
	double logInIceMaxEnergy = Particle.getLogEnergyMinimum()
	    + Particle.getDeltaLogEnergy()*(double )(Particle.getDimensionOfLogEnergyMatrix());
	int iLogE = 0;
	while(logInIceEnergy < logInIceMaxEnergy){// loop over inice particle

	    double logSurfaceEnergy = Particle.getLogEnergyMinimum() + epsilon;
	    double logSurfaceMaxEnergy = Particle.getLogEnergyMinimum()
		+ Particle.getDeltaLogEnergy()*(double )(Particle.getDimensionOfLogEnergyMatrix());

	    double cummulativeInIceDeltaLogE = 0.0; 
	    int cummulativeInIceBinNumber = 0;

	    int jLogE = 0; double fluxSum = 0.0;
	    while(logSurfaceEnergy < logSurfaceMaxEnergy){ // loop over surface particle

		cummulativeInIceDeltaLogE = 0.0;cummulativeInIceBinNumber = 0;
		double[] fluxToThisInIceEnergy = new double[3];
		for(int iFlavor=0;iFlavor<3;iFlavor++) fluxToThisInIceEnergy[iFlavor]=0.0;

		double cummulativeDeltaLogE = 0.0; int cummulativeBinNumber = 0;
		while(cummulativeInIceDeltaLogE < binWidth){

		    double[] flux = new double[3];
		    for(int iFlavor=0;iFlavor<3;iFlavor++) flux[iFlavor]=0.0;

		    cummulativeDeltaLogE = 0.0; cummulativeBinNumber = 0;
		    while(cummulativeDeltaLogE < binWidth){
			// from nu-e
			flux[0] += matrix.getDF(0,0,logSurfaceEnergy+cummulativeDeltaLogE,
					     iniceFlavor,iniceDoublet,
					     logInIceEnergy+cummulativeInIceDeltaLogE);
			// from nu-mu
			flux[1] += matrix.getDF(1,0,logSurfaceEnergy+cummulativeDeltaLogE,
					     iniceFlavor,iniceDoublet,
					     logInIceEnergy+cummulativeInIceDeltaLogE);
			// from nu-tau
			flux[2] += matrix.getDF(2,0,logSurfaceEnergy+cummulativeDeltaLogE,
					     iniceFlavor,iniceDoublet,
					     logInIceEnergy+cummulativeInIceDeltaLogE);

			cummulativeBinNumber ++;
			cummulativeDeltaLogE +=  Particle.getDeltaLogEnergy();
		    }
		    // average over the bins
		    for(int iFlavor=0;iFlavor<3;iFlavor++) 
			flux[iFlavor] = flux[iFlavor]/(double)cummulativeBinNumber; 

		    for(int iFlavor=0;iFlavor<3;iFlavor++) 
			fluxToThisInIceEnergy[iFlavor] += flux[iFlavor];
		    cummulativeInIceBinNumber ++;
		    cummulativeInIceDeltaLogE +=  Particle.getDeltaLogEnergy();
		}

		System.out.format("%5.2f %5.2f %11.5e %11.5e %11.5e\n",
				  Particle.getLogEnergyMinimum()+
				  Particle.getDeltaLogEnergy()*(double )iLogE,
				  Particle.getLogEnergyMinimum()+
				  Particle.getDeltaLogEnergy()*(double )jLogE,
				  fluxToThisInIceEnergy[0],
				  fluxToThisInIceEnergy[1],
				  fluxToThisInIceEnergy[2]);

		logSurfaceEnergy +=  cummulativeDeltaLogE;
		jLogE += cummulativeBinNumber;

		for(int iFlavor=0;iFlavor<3;iFlavor++) 
		    fluxSum += fluxToThisInIceEnergy[iFlavor];
	    }

	    System.err.format("%5.2f FluxSum=%11.5e\n",
			      Particle.getLogEnergyMinimum()+
			      Particle.getDeltaLogEnergy()*(double )iLogE,
			      fluxSum);
	    iLogE += cummulativeInIceBinNumber;
	    logInIceEnergy += cummulativeInIceDeltaLogE;
	}
    }
 
}
