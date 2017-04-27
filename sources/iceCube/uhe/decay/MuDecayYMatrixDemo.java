package iceCube.uhe.decay;

import iceCube.uhe.particles.*;
import iceCube.uhe.decay.*;
import java.io.*;

/** Make the mu decay matrix and calculate some numbers derived by it.*/
public class MuDecayYMatrixDemo{

    public static void main(String[] args) {

	int flavor = 1;        // Mu Flavor
        int doublet =1;        // Charged Lepton
	int material = 1;      // Rock
	double energy = 1.0e9; // 1EeV= 10^9 GeV

        // Generate the particle class.
        Particle mu = 
	    new Particle(flavor, doublet, energy); // Mu


        System.err.println("The Particle Name is " + 
        mu.particleName(mu.getFlavor(), mu.getDoublet()));

	//Generate object of the mu decay Y-Matrix
	MuDecayYMatrix muDecayYMtx = new MuDecayYMatrix(mu);


	int iLogE;
	for(iLogE=0;iLogE<mu.getDimensionOfLogEnergyMatrix();iLogE++){
	    double logEnergy = mu.getLogEnergyMinimum()+
		mu.getDeltaLogEnergy()*(double )iLogE;
	    energy = Math.pow(10.0,logEnergy);
	    if(iLogE%100==0) System.err.println("The Incident energy " + energy + " GeV");

	    // LifeTime [sec]
	    muDecayYMtx.setLifeTimeMatrix(iLogE);
	}

	int iLogY;
	for(iLogY=0;iLogY<mu.getDimensionOfLogEnergyMatrix();iLogY++){
	    // Decay Matrix
	    muDecayYMtx.setMuDecayMatrix(iLogY);
	}



	System.err.println("Matrix calculation done");

	for(iLogE=0;iLogE<mu.getDimensionOfLogEnergyMatrix();iLogE+=100){
	    double logEnergy = mu.getLogEnergyMinimum()+
		mu.getDeltaLogEnergy()*(double )iLogE;
	    energy = Math.pow(10.0,logEnergy);
	    System.out.println(energy + " GeV " + "lifeTime " + 
			       muDecayYMtx.getLifeTimeMatrix(iLogE)*1.0e9 + " nsec");
	}


	double sumNuMu = 0.0;	double sumNu = 0.0;
	double sumLeptons = 0.0;
	for(iLogY=0;iLogY<mu.getDimensionOfLogEnergyMatrix();iLogY++){
	    sumNuMu += muDecayYMtx.getMuToNuMuDecayMatrix(iLogY);
	    sumNu += muDecayYMtx.getMuToNuEDecayMatrix(iLogY);
	    sumLeptons += muDecayYMtx.getMuToEDecayMatrix(iLogY);
	}
	System.out.println("NuMu " + sumNuMu + " NuE " + sumNu);
	System.out.println("Lepton " + sumLeptons);

    }
}
