package iceCube.uhe.interactions;

import iceCube.uhe.particles.*;
import iceCube.uhe.interactions.*;
import iceCube.uhe.points.*;
import java.io.*;

public class DrawIonization {

    public static void main(String[] args) throws IOException {

        int materialNumber = 0; //ice
        if(args.length!=1){
            System.out.println("Usage: DrawIonization material-numner");
	    System.exit(0);
        }else{
            materialNumber = Integer.valueOf(args[0]).intValue();
        }

        // Generate the ParticlePoint class.
        ParticlePoint s = new ParticlePoint(0.0, 5.0*Math.PI/180.0,materialNumber);//Rock
	double massNumber = 0.0;
        for(int i=0;i<s.NumberOfSpecies[s.getMaterialNumber( )];i++){
	    massNumber += s.getNumberOfAtoms(i)*s.getAtomicNumber(i);
	}

	/** Directory path for the dumped InteractionsMatrix objects */
	String interactionsMatrixDirectoryInIce  = "iceCube/uhe/interactions/ice/";
	String interactionsMatrixDirectoryInRock = "iceCube/uhe/interactions/rock/";
	String interactionsMatrixDirectory;

        if(materialNumber==0)  interactionsMatrixDirectory = interactionsMatrixDirectoryInIce;
        else interactionsMatrixDirectory = interactionsMatrixDirectoryInRock;

	// Muon Bremss Matrix for comparison
	InputStream in = ClassLoader.getSystemResourceAsStream(interactionsMatrixDirectory.concat("muBremsstrahlungMtx"));

	InteractionsMatrix muonBrems = InteractionsMatrixInput.inputInteractionsMatrix(in);

	Ionization ion = new Ionization(new Particle(1,1),s); // muon
	System.err.println(ion.interactionName());

	//
	// Drawing
	//
	System.out.println("zone 2 1");

	System.out.println("titx Energy [GeV]");
	System.out.println("tity beta [g^-1 cm^2]");
	System.out.println("gwin 0.2 0.9 0.2 0.9");
	System.out.println("scal 0.0 12.0 -17.0 -5.0");

	// Bremss
	for(int iLogE=0;iLogE< Particle.getDimensionOfLogEnergyMatrix();iLogE++){

            double logE = Particle.getLogEnergyMinimum( ) + 
                Particle.getDeltaLogEnergy( )*(double )iLogE;

	    // Average Energy loss term
	    double beta = s.NA/massNumber*muonBrems.getInelasticityMatrix(iLogE);

	    System.out.println("data " + logE + " 0.0 " + Math.log(beta)/Math.log(10.0) +
			       " 0.0");
	}

	System.out.println("join");
	System.out.println("disp");
	System.out.println("cont");

	//Ionization
	double logE = 0.0;
	System.out.println("lncl 4");
	while(logE <=12.0){
	    double energy = Math.pow(10.0,logE);
	    ion.setIncidentParticleEnergy(energy);
	    double beta = ion.getDEDX()/energy;

	    System.out.println("data " + logE + " 0.0 " + Math.log(beta)/Math.log(10.0) +
			       " 0.0");
	    logE += 0.1;
	}

	System.out.println("join");
	System.out.println("disp");
	System.out.println("endg");


	System.out.println("titx Energy [GeV]");
	System.out.println("tity Energy loss [GeV g^-1 cm^2]");
	System.out.println("gwin 0.2 0.9 0.2 0.9");
	System.out.println("scal 1.0 1e12 1e-3 1e5");

	// Bremss
	for(int iLogE=0;iLogE< Particle.getDimensionOfLogEnergyMatrix();iLogE++){

            logE = Particle.getLogEnergyMinimum( ) + 
                Particle.getDeltaLogEnergy( )*(double )iLogE;
	    double energy = Math.pow(10.0,logE);

	    // Average Energy loss term
	    double beta = s.NA/massNumber*muonBrems.getInelasticityMatrix(iLogE)*energy;

	    System.out.println("data " + energy + " 0.0 " + beta +
			       " 0.0");
	}

	System.out.println("logx");
	System.out.println("logy");
	System.out.println("join");
	System.out.println("disp");
	System.out.println("cont");

	//Ionization
	logE = 0.0;
	System.out.println("lncl 4");
	while(logE <=12.0){
	    double energy = Math.pow(10.0,logE);
	    ion.setIncidentParticleEnergy(energy);
	    double beta = ion.getDEDX();

	    System.out.println("data " + energy + " 0.0 " + beta +
			       " 0.0");
	    logE += 0.1;
	}

	System.out.println("logx");
	System.out.println("logy");
	System.out.println("join");
	System.out.println("disp");
	System.out.println("cont");

	logE = 0.0;
	System.out.println("lncl 6");
	while(logE <=12.0){
	    double energy = Math.pow(10.0,logE);
	    ion.setIncidentParticleEnergy(energy);
	    double beta = s.NA/massNumber*ion.getYDSigmaDy(ion.getYmin(),ion.getYmax());

	    System.out.println("data " + energy + " 0.0 " + beta*energy +
			       " 0.0");
	    logE += 0.2;
	}

	System.out.println("logx");
	System.out.println("logy");
	System.out.println("join");
	System.out.println("disp");
	System.out.println("endg");

    
    }

}
