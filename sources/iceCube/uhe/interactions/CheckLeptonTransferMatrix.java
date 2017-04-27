package iceCube.uhe.interactions;

import iceCube.uhe.particles.*;
import iceCube.uhe.interactions.*;
import java.io.*;

/** Check Lepton Transfer matrix when inelasiticity is extremele small
such as tau's pair creation and calculation of
the differential cross section for z = 1- y ~ 1 could not avoid
numerical error.

Reading the original interaction matrix and output the calibrated matrix.
*/

public class CheckLeptonTransferMatrix {

    public static void main(String[] args) throws IOException {

        String fileName = null;
        if(args.length!=1){
            System.out.println("Usage: CheckLeptonTransferMatrix inputMatrixFileName");
	    System.exit(0);
        }else{
            fileName = args[0];
        }

	// Read the serialized object of the Neutrino Charged Interaction Matrix
        FileInputStream in = new FileInputStream(fileName);
	InteractionsMatrix intMtx = 
	    InteractionsMatrixInput.inputInteractionsMatrix(in);
	in.close( );

	int iLogE;
	for(iLogE=0;iLogE< Particle.getDimensionOfLogEnergyMatrix();iLogE++){

	    // Total Cross Section
	    double sigma = intMtx.getSigmaMatrix(iLogE);

	    double sumZ = 0.0; 
	    for(int jLogE=0;jLogE<=iLogE;jLogE++){
		sumZ += intMtx.getLeptonTransferMatrix(iLogE,jLogE);
	    }
	    double ratioZ = sumZ/sigma;

	    double sum = 0.0; 
	    for(int jLogE=0;jLogE<= iLogE;jLogE++){
		sum += intMtx.getTransferMatrix(iLogE,jLogE);
	    }
	    double ratio = sum/sigma;

	    System.out.println(iLogE + " Ratio(Z) = " + ratioZ +
			       " Ratio(Y) = " + ratio);
	}


    }
}
