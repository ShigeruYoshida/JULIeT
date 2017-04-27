package iceCube.uhe.interactions;

import iceCube.uhe.particles.*;
import iceCube.uhe.interactions.*;
import iceCube.uhe.points.*;
import java.io.*;

/** Check the generated transfer matrix of the  Neutrino Charged Interactions */
public class CheckNeutrinoChargeMtx {

    public static void main(String[] args) throws IOException {

        String fileName = null;
        if(args.length!=1){
            System.out.println("Usage: CheckNeutrinoChargeMtx file-name");
	    System.exit(0);
        }else{
            fileName = args[0];
        }

        // Generate the ParticlePoint class.
        ParticlePoint s = new ParticlePoint(0.0, 5.0*Math.PI/180.0,1);//Rock
	double massNumber = 0.0;
        for(int i=0;i<s.NumberOfSpecies[s.getMaterialNumber( )];i++){
	    massNumber += s.getNumberOfAtoms(i)*s.getAtomicNumber(i);
	}

	// Read the serialized object of the Neutrino Charged Interaction Matrix
        FileInputStream in = new FileInputStream(fileName);
	InteractionsMatrix nuCCMtx = 
	    InteractionsMatrixInput.inputInteractionsMatrix(in);
	in.close( );

	int iLogE;
	for(iLogE=0;iLogE<= Particle.getDimensionOfLogEnergyMatrix();iLogE+=50){

	    if(iLogE== Particle.getDimensionOfLogEnergyMatrix()) iLogE--;
	    // Total Cross Section
	    double sigma = nuCCMtx.getSigmaMatrix(iLogE);

	    // Average Energy loss term
	    double beta = s.NA/massNumber*nuCCMtx.getInelasticityMatrix(iLogE);

	    // Transfer Matrix Summation
	    int jLogE;
	    double sigmaY = 0.0; double sigmaZ = 0.0; 
	    for(jLogE=0;jLogE<= iLogE;jLogE++){
		sigmaY += nuCCMtx.getTransferMatrix(iLogE,jLogE);
		sigmaZ += nuCCMtx.getLeptonTransferMatrix(iLogE,jLogE);
	    }


	    System.out.println(iLogE);
	    System.out.println("sigma " + sigma + " sigmaY " + sigmaY + 
			       " sigmaZ " + sigmaZ);
	    System.out.println("beta [x e-6 g^-1 cm^2] " + beta*1.0e6);
	}
	System.out.println(nuCCMtx.interactions.s.getMaterialNumber());
	System.out.println(nuCCMtx.interactions.p.getDoublet());
	System.out.println(nuCCMtx.interactions.p.getFlavor());

    }
}
