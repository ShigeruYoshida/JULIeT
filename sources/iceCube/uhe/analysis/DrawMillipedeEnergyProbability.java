package iceCube.uhe.analysis;

import java.io.*;
import java.util.*;
import iceCube.uhe.particles.*;
import iceCube.uhe.propagation.*;
import iceCube.uhe.event.*;
import geometry.*;

public class DrawMillipedeEnergyProbability {

    /** Main method. In order to run JULIET, JulietEventGenerator object 
        and DataOutPutStream are generated. */
    public static void main(String[] args) throws IOException {

	String i3dataFileName = null;
	int runNumberToDraw = 0;
	double epsilon = 1.0e-12;
	boolean logPlot = false;

        if(args.length<2){
            System.out.println("Usage: DrawMillepedeEnergyProbability input-i3particle-file-name run-number log-plot?(yes 1)");
            System.exit(0);
        }else{
             i3dataFileName = args[0];
             runNumberToDraw = Integer.valueOf(args[1]).intValue();
             if(args.length ==3) {
               if(Integer.valueOf(args[2]).intValue()==1) logPlot = true;
             }
        }

	// input stream to read I3Particle
	InputStream in = ClassLoader.getSystemResourceAsStream(i3dataFileName);

	//
	// Now reading i3particle
	//
	I3Particle iceParticle = null; 
	boolean findEvent = false;
	while((iceParticle = I3ParticleInputStream.inputI3Particle(in)) !=null){

	    int flavor = iceParticle.getFlavor();
	    int doublet = iceParticle.getDoublet();
	    int runNumber = iceParticle.getIceCubeData().getEventNumber();

	    if(runNumber == runNumberToDraw){ // OK, we select this event
		findEvent = true;
		break;
	    }
	}

	if(findEvent){ // now draw the energy probability;

	    System.out.format(" titx Log(Energy [GeV])\n");
	    System.out.format(" tity probability [arbitrary]\n");
	    System.out.format(" scal 5.0 12.0 1.0e-5 0.1\n");

	    for(int iLogE=0; iLogE<Particle.getDimensionOfLogEnergyMatrix(); iLogE++){
		double logEnergy = 
		    Particle.getDeltaLogEnergy()*(double )iLogE + Particle.getLogEnergyMinimum();
		double prob = iceParticle.getLogEnergyMatrix(iLogE);

		if(prob==0.0) prob = epsilon;
		System.out.format(" data %f 0.0 %e 0.0\n",logEnergy,prob);
	    }

	    double logEnergy = iceParticle.getLogEnergy();
	    double distance = iceParticle.getDistanceFromEarthSurfaceToIceCube();
	    System.out.format(" line %f 1.0e-5 %f 0.1\n",logEnergy,logEnergy);
	    if(logPlot) System.out.format("logy\n");
	    System.out.format(" mssg distance from the earth surface %e cm\n",distance);
	    System.out.format(" mssg NPE %e \n",iceParticle.getIceCubeData().getBestNpe());
	    System.out.format("hist\n");
	    System.out.format("disp\n");
	    System.out.format("endg\n");
	}

	in.close();
    }


}
