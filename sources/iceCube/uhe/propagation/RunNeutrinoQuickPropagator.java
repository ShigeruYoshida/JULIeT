package  iceCube.uhe.propagation;

import iceCube.uhe.particles.*;
import iceCube.uhe.points.*;
import iceCube.uhe.propagation.*;
import java.io.*;

public class RunNeutrinoQuickPropagator {

    private static final double ln10 = Math.log(10.0);


    public static void main(String[] args) throws IOException {

	String fileName = null;
	int material = 0;
	double propDistance = 0.0;
	double angle = 0.0;
	int index = 0;
	double nuCCEnhancement = 1.0;

        if(args.length<4){
            System.out.println(
   "Usage: RunNeutrinoQuickPropagator medium(0 ice 1 rock) PropagationDistance[cm]/angle CCenhanceFactor switch(0 for distance, 1 for angle)");
            System.exit(0);
        }else{
            material = Integer.valueOf(args[0]).intValue();
            propDistance = Double.valueOf(args[1]).doubleValue();
	    nuCCEnhancement = Double.valueOf(args[2]).doubleValue();
            index = Integer.valueOf(args[3]).intValue();
	    if(index == 0){
		System.err.println("medium(" + material + ") Propagation Distance = " + 
				   propDistance + " [cm]");
	    }else{
		angle = Double.valueOf(args[1]).doubleValue();
		System.err.println("medium(" + material + ") angle = " + 
				   angle + " [deg]");
	    }
	}
        // Generate the ParticlePoint class.
        ParticlePoint s = new ParticlePoint(0.0, angle*Math.PI/180.0,material);

	// Generate the NeutrinoQuickPropagator class
	NeutrinoQuickPropagator propagator = new NeutrinoQuickPropagator(s);

	//
	// Propagte the Neutrino
	//

	if(index==0){
	    propagator.propagateNeutrino(propDistance*s.getMediumDensity(),nuCCEnhancement);
	}else{
	    propagator.propagateNeutrinoToIceCubeDepth(angle,nuCCEnhancement);
	}


	//
	// outputs the energy distributions after the propagation
	//

	System.out.println("titx Log E[GeV]");
	System.out.println("tity Log Count");
	System.out.println("scal 6.0 12.0 -15.0 1.0");

	double logEnergyInput = 10.0;   // 10^10 GeV
	for(int neutrinoFlavor = 1; neutrinoFlavor <3; neutrinoFlavor++){

	    // in-ice mu 
	    int color = 18+2*(neutrinoFlavor-1);
	    System.out.println("lncl " + color);
	    double logEnergyOutput = Particle.getLogEnergyMinimum();
	    while(logEnergyOutput < 12.0){ // E < 10^12 GeV 
		double count = propagator.getDF(neutrinoFlavor,logEnergyInput,
					       1,1,logEnergyOutput);
		double logCount = -15.0;
		if(count>0.0) logCount = Math.log(count)/ln10;
		System.out.println("data " + logEnergyOutput + " 0.0 " +
				   logCount + " 0.0");
		logEnergyOutput += Particle.getDeltaLogEnergy();
	    }

	    System.out.println("join");
	    System.out.println("disp");
	    System.out.println("cont");

	    // in-ice tau 
	    color = 8+2*(neutrinoFlavor-1);
	    System.out.println("lncl " + color);
	    logEnergyOutput = Particle.getLogEnergyMinimum();
	    while(logEnergyOutput < 12.0){ // E < 10^12 GeV 
		double count = propagator.getDF(neutrinoFlavor,logEnergyInput,
					       2,1,logEnergyOutput);
		double logCount = -15.0;
		if(count>0.0) logCount = Math.log(count)/ln10;
		System.out.println("data " + logEnergyOutput + " 0.0 " +
				   logCount + " 0.0");
		logEnergyOutput += Particle.getDeltaLogEnergy();
	    }

	    System.out.println("join");
	    System.out.println("disp");
	    System.out.println("cont");

	    // in-ice nu-mu 
	    color = 27+2*(neutrinoFlavor-1);
	    System.out.println("lncl " + color);
	    logEnergyOutput = Particle.getLogEnergyMinimum();
	    while(logEnergyOutput < 12.0){ // E < 10^12 GeV 
		double count = propagator.getDF(neutrinoFlavor,logEnergyInput,
					       1,0,logEnergyOutput);
		double logCount = -15.0;
		if(count>0.0) logCount = Math.log(count)/ln10;
		System.out.println("data " + logEnergyOutput + " 0.0 " +
				   logCount + " 0.0");
		logEnergyOutput += Particle.getDeltaLogEnergy();
	    }

	    System.out.println("join");
	    System.out.println("disp");
	    System.out.println("cont");
	}

	System.out.println("endg");


    }


}
