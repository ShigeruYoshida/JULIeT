package  iceCube.uhe.propagation;

import iceCube.uhe.particles.*;
import iceCube.uhe.propagation.*;
import java.io.*;

public class ComparePropMtxWithGlashow {

    static final double ln10 = Math.log(10.0);

    // Directory path for the precalculated PropagationMatrix files.
    private static String[] pathName = {
	"/disk0/data/IceCube/JULIeT/propMtxWithoutGlashow/BB/",
	"/disk0/data/IceCube/JULIeT/propMtxWithGlashow/BB/"
    };

    public static void main(String[] args) throws IOException {

        String fileName = null;
	String matrixName = null;
	int inputFlavor = 1;
	int inputDoublet = 0;
	double logEnergyInput = 6.0;
	int outputFlavor = 1;
	int outputDoublet = 1;


        if(args.length<6){
            System.out.println(
   "Usage: ComparePropMtxWithGlashow file-name inFlavor inDoublet logEin outFlavor outDoublet");
	    System.exit(0);
        }else{
            fileName = args[0];
	    inputFlavor = Integer.valueOf(args[1]).intValue();
	    inputDoublet = Integer.valueOf(args[2]).intValue();
	    logEnergyInput = Double.valueOf(args[3]).doubleValue();
	    outputFlavor = Integer.valueOf(args[4]).intValue();
	    outputDoublet = Integer.valueOf(args[5]).intValue();
        }

	// Input particle object
	double energyInput = Math.pow(10.0,logEnergyInput);
	Particle inputParticle = 
	    new Particle(inputFlavor,inputDoublet,energyInput);
	System.err.println("Primary particle is " + 
			   inputParticle.particleName(inputFlavor,inputDoublet));

	// Output particle object
	Particle outputParticle = 
	    new Particle(outputFlavor,outputDoublet);
	System.err.println("Post-propagation particle is " + 
			   outputParticle.particleName(outputFlavor,outputDoublet));


	// Instance of PropagationMatrixFactory	
	PropagationMatrixFactory matrixWithNoGlashow = 
	    new PropagationMatrixFactory();
	matrixWithNoGlashow.whetherPropagationMatrixWithGlashowResonance(false);
	
	PropagationMatrixFactory  matrixWithGlashow =
	    new PropagationMatrixFactory();

	// Read the serialized object of the Neutrino Charged Interaction Matrix
	String matrixDataFileName = pathName[0].concat(fileName);
	DataInputStream in = 
	    new DataInputStream(new FileInputStream(matrixDataFileName));
	matrixWithNoGlashow.readMatrix(in);
	in.close( );

	matrixDataFileName = pathName[1].concat(fileName);
	in = new DataInputStream(new FileInputStream(matrixDataFileName));
	matrixWithGlashow.readMatrix(in);
	in.close( );


	System.out.println("zone 2 1");

	System.out.println("titx Log E[GeV]");
	System.out.println("tity Log Count");
	System.out.println("gwin 0.2 0.9 0.2 0.9");
	System.out.println("scal 6.0 12.0 -12.0 1.0");

	// Matrix with No Glashow
	System.out.println("lnth 4");
	System.out.println("lncl 3");
	double logEnergyOutput = Particle.getLogEnergyMinimum();
	while(logEnergyOutput <= 12.0){ // E < 10^12 GeV
	    double energyOutput = Math.pow(10.0,logEnergyOutput);
	    outputParticle.putLogEnergy(logEnergyOutput);
	    outputParticle.putEnergy(logEnergyOutput);
	    double count = matrixWithNoGlashow.getDF(inputParticle,outputParticle);
	    double logCount = -15.0;
	    if(count>0.0) logCount = Math.log(count)/ln10;
	    System.out.println("data " + logEnergyOutput + " 0.0 " +
			       logCount + " 0.0");
	    logEnergyOutput += Particle.getDeltaLogEnergy();
	}
	System.out.println("join");
	System.out.println("disp");
	System.out.println("cont");

	System.out.println("lnth 1");
	System.out.println("lncl 4");
	// Matrix with  Glashow
	logEnergyOutput = Particle.getLogEnergyMinimum();
	while(logEnergyOutput <= 12.0){ // E < 10^12 GeV
	    double energyOutput = Math.pow(10.0,logEnergyOutput);
	    outputParticle.putLogEnergy(logEnergyOutput);
	    outputParticle.putEnergy(logEnergyOutput);
	    double count = matrixWithGlashow.getDF(inputParticle,outputParticle);
	    double logCount = -15.0;
	    if(count>0.0) logCount = Math.log(count)/ln10;
	    System.out.println("data " + logEnergyOutput + " 0.0 " +
			       logCount + " 0.0");
	    logEnergyOutput += Particle.getDeltaLogEnergy();
	}

	System.out.println("join");
	System.out.println("disp");
	System.out.println("endg");

	// ratio
	System.out.println("titx Log E[GeV]");
	System.out.println("tity Log Ratio (Glashow/no Glashow)");
	System.out.println("gwin 0.2 0.9 0.2 0.9");
	System.out.println("scal 6.0 12.0 -1.0 1.0");

	logEnergyOutput = Particle.getLogEnergyMinimum();
	while(logEnergyOutput <= 12.0){ // E < 10^12 GeV
	    double energyOutput = Math.pow(10.0,logEnergyOutput);
	    outputParticle.putLogEnergy(logEnergyOutput);
	    outputParticle.putEnergy(logEnergyOutput);
	    double count = matrixWithGlashow.getDF(inputParticle,outputParticle);
	    double logCountWithGlashow = -15.0;
	    if(count>0.0) logCountWithGlashow = Math.log(count)/ln10;
	    count = matrixWithNoGlashow.getDF(inputParticle,outputParticle);
	    double logCountWithNoGlashow = -15.0;
	    if(count>0.0) logCountWithNoGlashow = Math.log(count)/ln10;
	    double ratio = logCountWithGlashow-logCountWithNoGlashow;
	    System.out.println("data " + logEnergyOutput + " 0.0 " +
			       ratio + " 0.0");
	    logEnergyOutput += Particle.getDeltaLogEnergy();
	}

	System.out.println("join");
	System.out.println("disp");
	System.out.println("endg");
    }
}
