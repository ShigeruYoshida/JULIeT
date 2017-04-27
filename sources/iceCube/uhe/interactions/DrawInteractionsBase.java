package iceCube.uhe.interactions;

import iceCube.uhe.particles.*;
import iceCube.uhe.interactions.*;
import java.io.*;

/** Draw the total cross section and the energy loss("beta term") by reading
    the pre-calculated and serialized IntertactionMatrix object stored in the file. */
public class DrawInteractionsBase {

    public static void main(String[] args) throws IOException {

        String fileName = null;
        if(args.length!=1){
            System.out.println("Usage: DrawInteractionsBase IntMatrix-file-name");
	    System.exit(0);
        }else{
            fileName = args[0];
        }

	// Read the serialized object of the Interaction Matrix
        FileInputStream in = new FileInputStream(fileName);
	InteractionsMatrix intMtx = 
	    InteractionsMatrixInput.inputInteractionsMatrix(in);
	in.close( );

	// Build the InteractionsBase object
	InteractionsBase intBase = new InteractionsBase(intMtx);


	// Drawing
	System.out.println("titx Log Energy [GeV]");
	System.out.println("tity Cummulative Probability");
	System.out.println("gwin 0.2 0.9 0.2 0.9");
	System.out.println("scal 1.0 12.0 0.0 1.0");

	// Primary Particle Energy Loop
	for(int iLogE=0;iLogE<= Particle.getDimensionOfLogEnergyMatrix();iLogE+=50){
	    if(iLogE == Particle.getDimensionOfLogEnergyMatrix()) iLogE--;

	    double logE = Particle.getLogEnergyMinimum( ) + 
		Particle.getDeltaLogEnergy( )*(double )iLogE;

	    // Produced Particle Energy loop
	    int expandedDim =  Particle.getDimensionOfLogEnergyMatrix() +
		(int )((Particle.getLogEnergyMinimum()-
			intBase.getLogEnergyProducedMinimum())/Particle.getDeltaLogEnergy());
	    for(int jLogE=0;jLogE<expandedDim;jLogE++){
		double logProducedEnergy = intBase.getLogEnergyProducedMinimum( )+
		    Particle.getDeltaLogEnergy( )*(double )jLogE;
		double prob = intBase.getCumulativeProbability(logE, logProducedEnergy);
		System.out.println("data " + logProducedEnergy + " 0.0 " +
				   prob + " 0.0");
	    }
	    System.out.println("join");
	    System.out.println("disp");
	    System.out.println("cont");
	}
	System.out.println("endg");

    }


}
