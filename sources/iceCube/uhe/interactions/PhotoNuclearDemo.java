package iceCube.uhe.interactions;

import iceCube.uhe.particles.*;
import iceCube.uhe.interactions.*;
import iceCube.uhe.points.*;

/** PhotoNuclear.class Demo program */
public class PhotoNuclearDemo {

    public static void main(String[] args){

        int flavor =2;         // Tuon
        int doublet =1;
	int material = 1;      // Rock
	double energy = 1.0e9; // 1EeV= 10^9 GeV
	double y = 0.5;
	double yLogLow = 0.0;
	double yLogUp = 0.0;
	double yUp = 0.0;
	double yLow = 0.0;
	double epsilon = 1.0e-2;
	double energyCut;
	double sigma = 1.0;
	double yAverage = 0.5;
	double zAverage = 0.5;


        if(args.length!=2){
            System.out.println("Usage: PhotoNuclearDemo y epsilon");
	    System.exit(0);
        }else{
            y = Double.valueOf(args[0]).doubleValue();
            epsilon = Double.valueOf(args[1]).doubleValue();
        }


        // Generate the particle class.
        Particle mu = 
	    new Particle(flavor, doublet, energy); // Muon

        System.out.println("The Particle Name is " + 
        mu.particleName(mu.getFlavor(), mu.getDoublet()));

	// Generate the ParticlePoint class.
	ParticlePoint s = new ParticlePoint(0.0, 5.0*Math.PI/180.0,material);

	for(int i=0;i<s.NumberOfSpecies[material];i++){
	    System.out.println("Charge " + s.getCharge(i));
	    System.out.println("Atomic Number " + s.getAtomicNumber(i));
	}

	//Generate the PhotoNuclear classes
	PhotoNuclear muPhoto = new PhotoNuclear(mu, s);
	System.out.println(muPhoto.interactionName( ));

	numRecipes.Integration.setRelativeAccuracy(epsilon);
	// integration accuracy for save CPUtime

	// Cross Section 
	energyCut = Math.pow(10.0,mu.getLogEnergyMinimum( )-1.0); //GeV
	muPhoto.setEnergyCut(energyCut);
	System.out.println("epsilon " + epsilon + " energyCut " + 
			   energyCut + " GeV");

	int iLogE;
	for(iLogE=0;iLogE<=500;iLogE+=100){
	    muPhoto.setIncidentParticleEnergy(iLogE);
	    System.out.println("The Incident energy " + 
		   muPhoto.getIncidentParticleEnergy() + " GeV");

	    // Differntial Cross Section
	    System.out.println("  dSigmaB/dY(" + y + ") " +
			       muPhoto.getDSigmaDy(y));



	    // Partial Cross Section
	    yLow = muPhoto.getYmin( );
	    yLogLow = Math.log(yLow)/Math.log(10.0);
	    yLogUp = yLogLow + mu.getDeltaLogEnergy( );
	    yUp = Math.pow(10.0,yLogUp);
	    System.out.println("  YdSigmaB/dY(" + yLow + ") " +
			       muPhoto.getYDSigmaDy(yLow,yUp));

	    yUp = muPhoto.getYmax( );
	    yLogUp = Math.log(yUp)/Math.log(10.0);
	    yLogLow = yLogUp - mu.getDeltaLogEnergy( );
	    yLow = Math.pow(10.0,yLogLow);
	    System.out.println("  YdSigmaB/dY(" + yLow + ") " +
			       muPhoto.getYDSigmaDy(yLow,yUp));

	    // Total Cross Section
	    sigma = muPhoto.getSigma();
	    System.out.println("  Total cross section B " + sigma);

	    // Inelasticity
	    yAverage = muPhoto.getYDSigmaDy(muPhoto.getYmin(),
					      muPhoto.getYmax())/sigma;
	    System.out.println("  Inelasticity B <y> " + yAverage);
	    zAverage =  muPhoto.getZDSigmaDZ(1.0-muPhoto.getYmax(),
					       1.0-muPhoto.getYmin())/sigma;
	    System.out.println("  Inelasticity B <1-y> " + zAverage);
	    System.out.println("  sum " + (yAverage+zAverage));

	}



    }
}
