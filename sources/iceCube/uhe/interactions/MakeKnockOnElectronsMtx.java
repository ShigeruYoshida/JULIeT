package iceCube.uhe.interactions;

import iceCube.uhe.particles.*;
import iceCube.uhe.interactions.*;
import iceCube.uhe.points.*;
import java.io.*;

/** Make the energy transfer matrix of KnockOnElectrons Distribution Interactions */
public class MakeKnockOnElectronsMtx {

    public static void main(String[] args) throws IOException {

        String fileName = null;
	int flavor = 1;
        int doublet =1;        // Charged Lepton
	int material = 0;      // Ice->0, Rock->1
	double energy = 1.0e9; // 1EeV= 10^9 GeV
	double epsilon = 5.0e-3;
	double energyCut = 1.0e4; // 10000 GeV


        if(args.length!=2){
            System.out.println("Usage: MakeKnockOnElectronsMtx file-name flavor");
	    System.exit(0);
        }else{
            fileName = args[0];
            flavor = Integer.valueOf(args[1]).intValue();
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

	//Generate object of the Brremstrahlung interaction.
	KnockOnElectrons knockOn = new KnockOnElectrons(lepton, s);
	System.err.println(knockOn.interactionName( ));

	energyCut = Math.pow(10.0,lepton.getLogEnergyMinimum( )-2.0); //GeV
	knockOn.setEnergyCut(energyCut);

	//Generate object of the Interaction Matrix
	InteractionsMatrix knockOnMtx = new InteractionsMatrix(knockOn);


	numRecipes.Integration.setRelativeAccuracy(epsilon);
	// integration accuracy for save CPUtime

	int iLogE;
	for(iLogE=0;iLogE<lepton.getDimensionOfLogEnergyMatrix();iLogE++){
	    knockOn.setIncidentParticleEnergy(iLogE);
	    if(iLogE%50==0) System.err.println("The Incident energy " + 
			   knockOn.getIncidentParticleEnergy() + " GeV");

	    // Total Cross Section
	    knockOnMtx.setSigmaMatrix(iLogE);
	    //System.err.println("  Total cross section done");


	    // Transfer Matrix
	    int jLogE;
	    for(jLogE=0;jLogE<lepton.getDimensionOfLogEnergyMatrix();jLogE++){
		knockOnMtx.setTransferMatrix(iLogE,jLogE);
		//if(jLogE%100==0) System.err.println("  Transfer Matrix " + jLogE);
	    }

	}

        FileOutputStream out = new FileOutputStream(fileName);
        InteractionsMatrixOutput.outputInteractionsMatrix(knockOnMtx, out);
        out.close( );

    }
}
