package iceCube.uhe.interactions;

import iceCube.uhe.particles.*;
import iceCube.uhe.interactions.*;
import java.io.*;

/** Calibrate Lepton Transfer matrix when inelasiticity is extremele small
such as tau's pair creation and calculation of
the differential cross section for z = 1- y ~ 1 could not avoid
numerical error.

Reading the original interaction matrix and output the calibrated matrix.
*/

public class CalibrateLeptonTransferMatrix {

    public static void main(String[] args) throws IOException {

        String fileName = null;
	String outputFileName = null;
	boolean calibrated = false;
        if(args.length!=2){
            System.out.println("Usage: CalibrateLeptonTransferMatrix inputMatrixFileName outPutMatrixFileName");
	    System.exit(0);
        }else{
            fileName = args[0];
            outputFileName = args[1];
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
	    if(sigma==0.0){
		for(int jLogE=0;jLogE<=iLogE;jLogE++){
		    intMtx.transferAMtx[iLogE][jLogE] = 0.0;
		    intMtx.transferMtx[iLogE][jLogE] = 0.0;
		}
		calibrated = true;
		continue;
	    }
		

	    // Calibrate LeptionTransferMatrix
	    double sumZ = 0.0; 
	    for(int jLogE=0;jLogE<=iLogE;jLogE++){
		sumZ += intMtx.getLeptonTransferMatrix(iLogE,jLogE);
	    }
	    double ratioZ = sumZ/sigma;
	    if(ratioZ>1.0){
		double sum = 0.0; 
		calibrated = true;
		for(int jLogE=0;jLogE< iLogE;jLogE++){
		    sum += intMtx.getLeptonTransferMatrix(iLogE,jLogE);
		}
		double sigmaZ = sigma-sum;
		System.out.println("calibrated! " + iLogE + " Ratio(" + ratioZ +
				   ") leptonTransferMtx=" +
				   intMtx.getLeptonTransferMatrix(iLogE,iLogE) +
				   " sigmaZ=" + sigmaZ);
		intMtx.transferAMtx[iLogE][iLogE] = sigmaZ;
	    }

	    // Calibrate TransferMatrix
	    double sumY = 0.0; 
	    for(int jLogE=0;jLogE<=iLogE;jLogE++){
		sumY += intMtx.getTransferMatrix(iLogE,jLogE);
	    }
	    double ratio = sumY/sigma;
	    if(ratio>1.0){
		double sum = 0.0; 
		calibrated = true;
		for(int jLogE=0;jLogE< iLogE;jLogE++){
		    sum += intMtx.getTransferMatrix(iLogE,jLogE);
		}
		double sigmaY = sigma-sum;
		System.out.println("calibrated! " + iLogE + " Ratio(" + ratio +
				   ") transferMtx=" +
				   intMtx.getTransferMatrix(iLogE,iLogE) +
				   " sigmaY=" + sigmaY);
		intMtx.transferMtx[iLogE][iLogE] = sigmaY;
		if(sigmaY<0.0){
		    int kLogE =1;
		    do{
			intMtx.transferMtx[iLogE][iLogE-kLogE+1] = 0.0;
			sigmaY += intMtx.getTransferMatrix(iLogE,iLogE-kLogE);
			System.out.println("  #### further calibration " +
				   intMtx.getTransferMatrix(iLogE,iLogE-kLogE) +
				   " sigmaY=" + sigmaY);
			intMtx.transferMtx[iLogE][iLogE-kLogE] = sigmaY;
			kLogE++;
		    }while(sigmaY<0.0 && ((iLogE-kLogE) >=0) );
		}

	    }

	}

	// Output the calibrated interaction matrix
	if(calibrated){
	    FileOutputStream out = new FileOutputStream(outputFileName);
	    InteractionsMatrixOutput.outputInteractionsMatrix(intMtx, out);
	    out.close( );
	}


    }
}
