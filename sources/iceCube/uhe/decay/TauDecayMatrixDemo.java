package iceCube.uhe.decay;

import iceCube.uhe.particles.*;
import iceCube.uhe.decay.*;
import java.io.*;

/** Make the tau decay matrix and calculate some numbers derived by it.*/
public class TauDecayMatrixDemo{

    public static void main(String[] args) {

	int flavor = 2;        // Tau Flavor
        int doublet =1;        // Charged Lepton
	int material = 1;      // Rock
	double energy = 1.0e9; // 1EeV= 10^9 GeV

        // Generate the particle class.
        Particle tau = 
	    new Particle(flavor, doublet, energy); // Tau


        System.err.println("The Particle Name is " + 
        tau.particleName(tau.getFlavor(), tau.getDoublet()));

	//Generate object of the tau decay Matrix
	TauDecayMatrix tauDecayMtx = new TauDecayMatrix(tau);


	int iLogE;
	for(iLogE=0;iLogE<tau.getDimensionOfLogEnergyMatrix();iLogE++){
	    double logEnergy = tau.getLogEnergyMinimum()+
		tau.getDeltaLogEnergy()*(double )iLogE;
	    energy = Math.pow(10.0,logEnergy);
	    if(iLogE%100==0) System.err.println("The Incident energy " + energy + " GeV");

	    // LifeTime [sec]
	    tauDecayMtx.setLifeTimeMatrix(iLogE);


	    // Decay Matrix
	    int jLogE;
	    for(jLogE=0;jLogE<tau.getDimensionOfLogEnergyMatrix();jLogE++){
		tauDecayMtx.setTauDecayMatrix(iLogE,jLogE);
	    }

	}

	System.err.println("Matrix calculation done");

	for(iLogE=0;iLogE<tau.getDimensionOfLogEnergyMatrix();iLogE+=100){
	    double logEnergy = tau.getLogEnergyMinimum()+
		tau.getDeltaLogEnergy()*(double )iLogE;
	    energy = Math.pow(10.0,logEnergy);
	    System.out.println(energy + " GeV " + "lifeTime " + 
			       tauDecayMtx.getLifeTimeMatrix(iLogE)*1.0e9 + " nsec");
	    int jLogE;
	    double sumNuTau = 0.0;	double sumNu = 0.0;
	    double sumLeptons = 0.0;double sumHadrons = 0.0;
	    for(jLogE=0;jLogE<tau.getDimensionOfLogEnergyMatrix();jLogE++){
		sumNuTau += tauDecayMtx.getTauToNuTauDecayMatrix(iLogE,jLogE);
		sumNu += tauDecayMtx.getTauToNuDecayMatrix(iLogE,jLogE);
		sumLeptons += tauDecayMtx.getTauToChargedLeptonDecayMatrix(iLogE,jLogE);
		sumHadrons += tauDecayMtx.getTauToHadronDecayMatrix(iLogE,jLogE);
	    }
	    System.out.println("NuTau " + sumNuTau + " NuE/NuMu " + sumNu);
	    System.out.println("Lepton " + sumLeptons + " Hadron " + sumHadrons);

	    double sum = sumNu+sumLeptons;
	    System.out.println("NuE/NuMu+Leptons " + sum);
	}

    }
}
