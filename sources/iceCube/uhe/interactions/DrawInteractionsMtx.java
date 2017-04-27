package iceCube.uhe.interactions;

import iceCube.uhe.particles.*;
import iceCube.uhe.interactions.*;
import iceCube.uhe.points.*;
import java.io.*;

/** Draw the total cross section and the energy loss("beta term") by reading
    the pre-calculated and serialized IntertactionMatrix object stored in the file. */
public class DrawInteractionsMtx {

    public static void main(String[] args) throws IOException {

        String fileName = null;
	boolean isExtended = false;
        if(args.length<1){
            System.out.println("Usage: CheckNeutrinoChargeMtx file-name (extended?)");
	    System.exit(0);
        }else{
            fileName = args[0];
	    if(args.length == 2) isExtended = true;
        }

        // Generate the ParticlePoint class.
        ParticlePoint s = new ParticlePoint(0.0, 5.0*Math.PI/180.0,1);//Rock
	double massNumber = 0.0;
        for(int i=0;i<s.NumberOfSpecies[s.getMaterialNumber( )];i++){
	    massNumber += s.getNumberOfAtoms(i)*s.getAtomicNumber(i);
	}

	// Read the serialized object of the Interaction Matrix
        FileInputStream in = new FileInputStream(fileName);
	InteractionsMatrix nuCCMtx = 
	    InteractionsMatrixInput.inputInteractionsMatrix(in);
	in.close( );

	Particle lepton = nuCCMtx.interactions.p;
	if(isExtended){
	    lepton.putLogEnergyMinimum(0.0); // 1 GeV threshold
	    lepton.putDimensionOfLogEnergyMatrix(1200); // up to 10^12 GeV
	}


        System.out.println("zone 2 2");


	System.out.println("titx Energy [GeV]");
	System.out.println("tity Cross Section [cm^2]");
	System.out.println("gwin 0.2 0.9 0.2 0.9");
	if(!isExtended) System.out.println("scal 5.0 12.0 -34.0 -24.0");
	else System.out.println("scal 0.0 12.0 -34.0 -24.0");

	int iLogE;
	for(iLogE=0;iLogE< lepton.getDimensionOfLogEnergyMatrix();iLogE++){

	    // Total Cross Section
	    double sigma = nuCCMtx.getSigmaMatrix(iLogE);

	    double logE = lepton.getLogEnergyMinimum( ) + 
		lepton.getDeltaLogEnergy( )*(double )iLogE;
	    System.out.println("data " + logE + " 0.0 " + Math.log(sigma)/Math.log(10.0) +
			       " 0.0");
	}
	System.out.println("join");
	System.out.println("disp");
	System.out.println("endg");

	System.out.println("titx Energy [GeV]");
	System.out.println("tity beta [g^-1 cm^2]");
	System.out.println("gwin 0.2 0.9 0.2 0.9");
	if(!isExtended) System.out.println("scal 5.0 12.0 -12.0 -5.0");
	else System.out.println("scal 0.0 12.0 -12.0 -5.0");

	for(iLogE=0;iLogE< lepton.getDimensionOfLogEnergyMatrix();iLogE++){

	    // Total Cross Section
	    double sigma = nuCCMtx.getSigmaMatrix(iLogE);

	    double logE = lepton.getLogEnergyMinimum( ) + 
		lepton.getDeltaLogEnergy( )*(double )iLogE;

	    // Average Energy loss term
	    double beta = s.NA/massNumber*nuCCMtx.getInelasticityMatrix(iLogE);

	    System.out.println("data " + logE + " 0.0 " + Math.log(beta)/Math.log(10.0) +
			       " 0.0");
	}

	System.out.println("join");
	System.out.println("disp");
	System.out.println("endg");

	System.out.println("titx Energy [GeV]");
	System.out.println("tity Differential Cross Section [GeV]");
	System.out.println("gwin 0.2 0.9 0.2 0.9");
	if(!isExtended) System.out.println("scal 5.0 12.0 -40.0 -24.0");
	else System.out.println("scal 0.0 12.0 -40.0 -24.0");

	int jLogE;
	for(iLogE=0;iLogE< lepton.getDimensionOfLogEnergyMatrix();iLogE+=50){
	    for(jLogE=0;jLogE< lepton.getDimensionOfLogEnergyMatrix();jLogE++){

		// Differential Cross Section
		double sigma = nuCCMtx.getLeptonTransferMatrix(iLogE,jLogE);

		double logE = lepton.getLogEnergyMinimum( ) + 
		    lepton.getDeltaLogEnergy( )*(double )jLogE;
		if(sigma>0.0)
		System.out.println("data " + logE + " 0.0 " + 
				   Math.log(sigma)/Math.log(10.0) + " 0.0");
	    }
	    System.out.println("join");
	    System.out.println("disp");
	    System.out.println("cont");
	}

	System.out.println("endg");
    

	System.out.println("titx Energy [GeV]");
	System.out.println("tity Differential Cross Section [GeV]");
	System.out.println("gwin 0.2 0.9 0.2 0.9");
	if(!isExtended) System.out.println("scal 5.0 12.0 -40.0 -24.0");
	else System.out.println("scal 0.0 12.0 -40.0 -24.0");

	for(iLogE=0;iLogE< lepton.getDimensionOfLogEnergyMatrix();iLogE+=50){
	    for(jLogE=0;jLogE< lepton.getDimensionOfLogEnergyMatrix();jLogE++){

		// Differential Cross Section
		double sigma = nuCCMtx.getTransferMatrix(iLogE,jLogE);

		double logE = lepton.getLogEnergyMinimum( ) + 
		    lepton.getDeltaLogEnergy( )*(double )jLogE;
		if(sigma>0.0)
		System.out.println("data " + logE + " 0.0 " + 
				   Math.log(sigma)/Math.log(10.0) + " 0.0");
	    }
	    System.out.println("join");
	    System.out.println("disp");
	    System.out.println("cont");
	}

	System.out.println("endg");
    
    }

}
