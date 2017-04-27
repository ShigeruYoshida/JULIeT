package iceCube.uhe.interactions;

import iceCube.uhe.particles.*;
import iceCube.uhe.interactions.*;
import iceCube.uhe.points.*;
import java.io.*;

/** Check the generated transfer matrix of the  Neutrino Charged Interactions */
public class InfoFromMatrixExample {

    public static void main(String[] args) throws IOException {

        String fileName = null;
        if(args.length!=1){
            System.out.println("Usage:InfoFromMatrixExample  file-name");
	    System.exit(0);
        }else{
            fileName = args[0];
        }

	// Read the serialized object of the Interaction Matrix
        FileInputStream in = new FileInputStream(fileName);
	InteractionsMatrix intMtx = 
	    InteractionsMatrixInput.inputInteractionsMatrix(in);
	in.close( );

        // Reference the ParticlePoint class in the Interaction Matrix
        ParticlePoint s = intMtx.interactions.s;
	System.out.println("The medium's material number is " +
			   s.getMaterialNumber());

	// Primary Input Particle
	Particle primaryParticle = intMtx.interactions.p;
	System.out.println("The input particle is " +
		   primaryParticle.particleName(primaryParticle.getFlavor(),
						primaryParticle.getDoublet()
						));

	// Produced particle 
	int producedFlavor = intMtx.interactions.producedFlavor;
	System.out.println("The produced particle is " +
			   Particle.particleName(producedFlavor,1));

        // Showing the total cross section as a function of energy
	int iLogE;
	for(iLogE=0;iLogE<= Particle.getDimensionOfLogEnergyMatrix();iLogE+=50){

	    if(iLogE== Particle.getDimensionOfLogEnergyMatrix()) iLogE--;
	    // Total Cross Section
	    double sigma = intMtx.getSigmaMatrix(iLogE);

	    System.out.println("sigma " + sigma);
	}

    }
}
