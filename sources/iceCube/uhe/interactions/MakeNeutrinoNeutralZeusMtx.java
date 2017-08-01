package iceCube.uhe.interactions;

import iceCube.uhe.particles.*;
import iceCube.uhe.interactions.*;
import iceCube.uhe.points.*;
import java.io.*;

/** Make the energy transfer matrix of Neutrino Charged Interactions */
public class MakeNeutrinoNeutralZeusMtx {

    public static void main(String[] args) throws IOException {

        String fileName = null;
        int flavor =2;         // e->0, mu->1, tau->2
        int doublet =0;        // Neutrinos
	int material = 0;      // Ice->0, Rock->1
	double energy = 1.0e9; // 1EeV= 10^9 GeV
	double epsilon = 1.0e-4;


        if(args.length!=3){
            System.out.println("Usage: MakeNeutrinoNeutralZeusMtx file-name flavor material");
	    System.exit(0);
        }else{
            fileName = args[0];
	    flavor = Integer.valueOf(args[1]).intValue();
	    material = Integer.valueOf(args[2]).intValue();
        }


        // Generate the particle class.
        Particle nu = 
	    new Particle(flavor, doublet, energy); 

        System.err.println("The Particle Name is " + 
        nu.particleName(nu.getFlavor(), nu.getDoublet()));

	// Generate the ParticlePoint class.
	ParticlePoint s = new ParticlePoint(0.0, 5.0*Math.PI/180.0,material);

	for(int i=0;i<s.NumberOfSpecies[material];i++){
	    System.out.println("Charge " + s.getCharge(i));
	    System.out.println("Atomic Number " + s.getAtomicNumber(i));
	}

	//Generate object of the Charged Current Neutrino interaction.
	NeutrinoNeutralZeusNewPDF nuNC = new NeutrinoNeutralZeusNewPDF(nu, s);
	System.err.println(nuNC.interactionName( ));

	//Generate object of the Interaction Matrix
	InteractionsMatrix nuNCMtx = new InteractionsMatrix(nuNC);


	numRecipes.Integration.setRelativeAccuracy(epsilon);
	// integration accuracy for save CPUtime

	int iLogE;
	for(iLogE=0;iLogE<nu.getDimensionOfLogEnergyMatrix();iLogE++){
	    nuNC.setIncidentParticleEnergy(iLogE);
	    if(iLogE%100==0) System.err.println("The Incident energy " + 
			   nuNC.getIncidentParticleEnergy() + " GeV");

	    // Total Cross Section
	    nuNCMtx.setSigmaMatrix(iLogE);
	    //System.err.println("  Total cross section done");


	    // Transfer Matrix
	    int jLogE;
	    for(jLogE=0;jLogE<nu.getDimensionOfLogEnergyMatrix();jLogE++){
		nuNCMtx.setTransferMatrix(iLogE,jLogE);
		//if(jLogE%100==0) System.err.println("  Transfer Matrix " + jLogE);
	    }

	}

        FileOutputStream out = new FileOutputStream(fileName);
        InteractionsMatrixOutput.outputInteractionsMatrix(nuNCMtx, out);
        out.close( );

    }
}
