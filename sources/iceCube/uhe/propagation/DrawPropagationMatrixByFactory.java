package  iceCube.uhe.propagation;

import iceCube.uhe.particles.*;
import iceCube.uhe.propagation.*;
import java.io.*;

public class DrawPropagationMatrixByFactory {

    static final double ln10 = Math.log(10.0);
    static double logEnergyBase = Particle.getLogEnergyMinimum();
    static double energyBase = Math.pow(10.0,logEnergyBase);

    public static void main(String[] args) throws IOException {

        String fileName = null;
	String matrixName = null;
	int inputFlavor = 1;
	int inputDoublet = 0;
	double logEnergy= 6.0;
	int outputFlavor = 1;
	int outputDoublet = 1;
	boolean isInputE = true;
	boolean withGlashowResonance = true;
	boolean powerLawWeighted = false;
	double gamma = 1.0;

	PropagationMatrixFactory matrix = null;


        if(args.length<7){
            System.out.println(
   "Usage: DrawPropagationMatrixByFactory file-name inFlavor inDoublet logE outFlavor outDoublet inputEnergy?(yes 1 no 0) (withoutGlashow - 0) (power law index)");
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
	    if(args.length>7){
		if(Integer.valueOf(args[7]).intValue() == 0){
		    withGlashowResonance = false;
		    System.err.println(" matrix data with no Glasshow Resonance channels");
		}
	    }
	    if(args.length>8){
		gamma = Double.valueOf(args[8]).doubleValue();
		powerLawWeighted = true;
		System.err.println(" weighted by powerlaw spectrum with index of " + gamma);
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
	if(!withGlashowResonance) matrix.whetherPropagationMatrixWithGlashowResonance(withGlashowResonance);
	DataInputStream in = new DataInputStream(new FileInputStream(fileName));
	matrix.readMatrix(in);
	in.close( );


	System.out.println("titx Log E[GeV]");
	System.out.println("tity Log Count");
	System.out.println("scal 5.0 12.0 -15.0 1.0");

	double sum_count = 0.0;
	if(isInputE){
	    double logEnergyOutput = Particle.getLogEnergyMinimum();
	    while(logEnergyOutput <= 12.0){ // E < 10^12 GeV
		energyOutput = Math.pow(10.0,logEnergyOutput);
		outputParticle.putLogEnergy(logEnergyOutput);
		outputParticle.putEnergy(logEnergyOutput);
		double count = matrix.getDF(inputParticle,outputParticle);
		sum_count += count;
		if(powerLawWeighted) count = count*Math.pow(energyInput/energyBase,-(gamma-1.0));
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
		sum_count += count;
		if(powerLawWeighted) count = count*Math.pow(energyInput/energyBase,-(gamma-1.0));
		double logCount = -15.0;
		if(count>0.0) logCount = Math.log(count)/ln10;
		System.out.println("data " + logEnergyInput + " 0.0 " +
				   logCount + " 0.0");
		logEnergyInput += Particle.getDeltaLogEnergy();
	    }
	}

	System.out.format("mssg sum(%6.3e)\n",sum_count);
	System.out.println("join");
	System.out.println("disp");
	System.out.println("endg");

    }
}
