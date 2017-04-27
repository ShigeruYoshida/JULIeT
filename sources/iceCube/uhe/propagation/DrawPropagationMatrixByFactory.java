package  iceCube.uhe.propagation;

import iceCube.uhe.particles.*;
import iceCube.uhe.propagation.*;
import java.io.*;

public class DrawPropagationMatrixByFactory {

    static final double ln10 = Math.log(10.0);

    public static void main(String[] args) throws IOException {

        String fileName = null;
	String matrixName = null;
	int inputFlavor = 1;
	int inputDoublet = 0;
	double logEnergy= 6.0;
	int outputFlavor = 1;
	int outputDoublet = 1;
	boolean isInputE = true;

	PropagationMatrixFactory matrix = null;


        if(args.length<7){
            System.out.println(
   "Usage: DrawPropagationMatrixByFactory file-name inFlavor inDoublet logE outFlavor outDoublet inputEnergy?(yes 1 no 0) (withoutGlashow - 0)");
	    System.exit(0);
        }else{
            fileName = args[0];
	    inputFlavor = Integer.valueOf(args[1]).intValue();
	    inputDoublet = Integer.valueOf(args[2]).intValue();
	    logEnergy = Double.valueOf(args[3]).doubleValue();
	    outputFlavor = Integer.valueOf(args[4]).intValue();
	    outputDoublet = Integer.valueOf(args[5]).intValue();
	    if(Integer.valueOf(args[6]).intValue() == 1){
		isInputE = true;
	    }else{
		isInputE = false;
	    }
		
        }

	// Input particle object
	double energyInput = Math.pow(10.0,logEnergy);
	Particle inputParticle = 
	    new Particle(inputFlavor,inputDoublet,energyInput);
	System.err.println("Primary particle is " + 
			   inputParticle.particleName(inputFlavor,inputDoublet));

	// Output particle object
	double energyOutput = Math.pow(10.0,logEnergy);
	Particle outputParticle = 
	    new Particle(outputFlavor,outputDoublet,energyOutput);
	System.err.println("Post-propagation particle is " + 
			   outputParticle.particleName(outputFlavor,outputDoublet));


	// Instance of PropagationMatrixFactory
	matrix = new PropagationMatrixFactory();
	// Read the serialized object of the Neutrino Charged Interaction Matrix
	if(args.length>7) matrix.whetherPropagationMatrixWithGlashowResonance(false);
	DataInputStream in = new DataInputStream(new FileInputStream(fileName));
	matrix.readMatrix(in);
	in.close( );


	System.out.println("titx Log E[GeV]");
	System.out.println("tity Log Count");
	System.out.println("scal 5.0 12.0 -15.0 1.0");

	if(isInputE){
	    double logEnergyOutput = Particle.getLogEnergyMinimum();
	    while(logEnergyOutput <= 12.0){ // E < 10^12 GeV
		energyOutput = Math.pow(10.0,logEnergyOutput);
		outputParticle.putLogEnergy(logEnergyOutput);
		outputParticle.putEnergy(logEnergyOutput);
		double count = matrix.getDF(inputParticle,outputParticle);
		double logCount = -15.0;
		if(count>0.0) logCount = Math.log(count)/ln10;
		System.out.println("data " + logEnergyOutput + " 0.0 " +
				   logCount + " 0.0");
		logEnergyOutput += Particle.getDeltaLogEnergy();
	    }
	}else{
	    double logEnergyInput = Particle.getLogEnergyMinimum();
	    while(logEnergyInput <= 12.0){ // E < 10^12 GeV
		energyInput = Math.pow(10.0,logEnergyInput);
		inputParticle.putLogEnergy(logEnergyInput);
		inputParticle.putEnergy(logEnergyInput);
		double count = matrix.getDF(inputParticle,outputParticle);
		double logCount = -15.0;
		if(count>0.0) logCount = Math.log(count)/ln10;
		System.out.println("data " + logEnergyInput + " 0.0 " +
				   logCount + " 0.0");
		logEnergyInput += Particle.getDeltaLogEnergy();
	    }
	}

	System.out.println("join");
	System.out.println("disp");
	System.out.println("endg");

    }
}
