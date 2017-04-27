package iceCube.uhe.interactions;

import iceCube.uhe.particles.*;
import iceCube.uhe.interactions.*;
import iceCube.uhe.points.*;
import java.io.*;

/** Make the energy transfer matrix of PairCreation Distribution Interactions */
public class MakePairCreationMtx {

    public static void main(String[] args) throws IOException {

        String fileName = null;
	int flavor = 1;
	int producedFlavor = 0;     // Recoiled particle flavor
        int doublet =1;        // Charged Lepton
	int material = 1;      // Rock
	double energy = 1.0e9; // 1EeV= 10^9 GeV
	double epsilon = 1.0e-1;
	double energyCut = 1.0e4; // 10000 GeV
	double ln10 = Math.log(10.0);


        if(args.length!=3){
            System.out.println("Usage: MakePairCreationMtx file-name flavor producedFlavor");
	    System.exit(0);
        }else{
            fileName = args[0];
            flavor = Integer.valueOf(args[1]).intValue();
            producedFlavor = Integer.valueOf(args[2]).intValue();
        }


        // Generate the particle class.
        Particle lepton = 
	    new Particle(flavor, doublet, energy); // Charged Leptons


        System.err.println("The Particle Name is " + 
        lepton.particleName(lepton.getFlavor(), lepton.getDoublet()));

	// Generate the ParticlePoint class.
	ParticlePoint s = new ParticlePoint(0.0, 5.0*Math.PI/180.0,material);

	for(int i=0;i<s.NumberOfSpecies[material];i++){
	    System.out.println("Charge " + s.getCharge(i));
	    System.out.println("Atomic Number " + s.getAtomicNumber(i));
	}

	//Generate object of the Pair Creation interaction.
	PairCreation pairC = new PairCreation(lepton, s, producedFlavor);
	System.err.println(pairC.interactionName( ));

	//Generate object of the Interaction Matrix
	InteractionsMatrix pairCMtx = new InteractionsMatrix(pairC);


	numRecipes.Integration.setRelativeAccuracy(epsilon);
	// integration accuracy for save CPUtime

	int iLogE;
	for(iLogE=0;iLogE<lepton.getDimensionOfLogEnergyMatrix();iLogE++){
	    pairC.setIncidentParticleEnergy(iLogE);
	    energy = pairC.getIncidentParticleEnergy();
	    System.err.println("The Incident energy " + energy + " GeV");

	    double logEnergy = Math.log(energy)/ln10;
	    if(logEnergy >(lepton.getLogEnergyMinimum( )+2.0)) 
		logEnergy = lepton.getLogEnergyMinimum( )+2.0;

	    energyCut = Math.pow(10.0,logEnergy-4.0); //GeV
	    pairC.setEnergyCut(energyCut);


	    // Total Cross Section
	    pairCMtx.setSigmaMatrix(iLogE);
	    System.err.println("  Total cross section done");


	    // Transfer Matrix
	    int jLogE;
	    for(jLogE=0;jLogE<lepton.getDimensionOfLogEnergyMatrix();jLogE++){
		pairCMtx.setTransferMatrix(iLogE,jLogE);
		//if(jLogE%100==0) System.err.println("  Transfer Matrix " + jLogE);
	    }

	}

        FileOutputStream out = new FileOutputStream(fileName);
        InteractionsMatrixOutput.outputInteractionsMatrix(pairCMtx, out);
        out.close( );

    }
}
