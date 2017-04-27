package iceCube.uhe.interactions;

import iceCube.uhe.particles.*;
import iceCube.uhe.interactions.*;
import iceCube.uhe.points.*;

/** Bremsstrahlung.class and KnockOnElectrons.class Demo program */
public class BremsstrahlungDemo {

    public static void main(String[] args){

        int flavor =1;         // Muon
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
	double sigmaB = 1.0;
	double sigmaK = 1.0;
	double yAverage = 0.5;
	double zAverage = 0.5;


        if(args.length!=2){
            System.out.println("Usage: BremsstrahlungDemo y epsilon");
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

	//Generate the Bremss/Knockon electron classes
	Bremsstrahlung muToPhotonB = new Bremsstrahlung(mu, s);
	System.out.println(muToPhotonB.interactionName( ));
	KnockOnElectrons muToPhotonK = new KnockOnElectrons(mu, s);
	System.out.println(muToPhotonK.interactionName( ));


	numRecipes.Integration.setRelativeAccuracy(epsilon);
	// integration accuracy for save CPUtime

	// Cross Section 
	energyCut = Math.pow(10.0,mu.getLogEnergyMinimum( )-1.0); //GeV
	muToPhotonB.setEnergyCut(energyCut);
	muToPhotonK.setEnergyCut(energyCut);
	System.out.println("epsilon " + epsilon + " energyCut " + 
			   energyCut + " GeV");

	int iLogE;
	for(iLogE=0;iLogE<=500;iLogE+=100){
	    muToPhotonB.setIncidentParticleEnergy(iLogE);
	    muToPhotonK.setIncidentParticleEnergy(iLogE);
	    System.out.println("The Incident energy " + 
			   muToPhotonK.getIncidentParticleEnergy() + " GeV");

	    // Differntial Cross Section
	    System.out.println("  dSigmaB/dY(" + y + ") " +
			       muToPhotonB.getDSigmaDy(y));
	    System.out.println("  dSigmaK/dY(" + y + ") " +
			       muToPhotonK.getDSigmaDy(y));



	    // Partial Cross Section
	    yLow = muToPhotonB.getYmin( );
	    yLogLow = Math.log(yLow)/Math.log(10.0);
	    yLogUp = yLogLow + mu.getDeltaLogEnergy( );
	    yUp = Math.pow(10.0,yLogUp);
	    System.out.println("  YdSigmaB/dY(" + yLow + ") " +
			       muToPhotonB.getYDSigmaDy(yLow,yUp));

	    yLow = muToPhotonK.getYmin( );
	    yLogLow = Math.log(yLow)/Math.log(10.0);
	    yLogUp = yLogLow + mu.getDeltaLogEnergy( );
	    yUp = Math.pow(10.0,yLogUp);
	    System.out.println("  YdSigmaK/dY(" + yLow + ") " +
			       muToPhotonK.getYDSigmaDy(yLow,yUp));

	    yUp = muToPhotonB.getYmax( );
	    yLogUp = Math.log(yUp)/Math.log(10.0);
	    yLogLow = yLogUp - mu.getDeltaLogEnergy( );
	    yLow = Math.pow(10.0,yLogLow);
	    System.out.println("  YdSigmaB/dY(" + yLow + ") " +
			       muToPhotonB.getYDSigmaDy(yLow,yUp));
	    yUp = muToPhotonK.getYmax( );
	    yLogUp = Math.log(yUp)/Math.log(10.0);
	    yLogLow = yLogUp - mu.getDeltaLogEnergy( );
	    yLow = Math.pow(10.0,yLogLow);
	    System.out.println("  YdSigmaK/dY(" + yLow + ") " +
			       muToPhotonK.getYDSigmaDy(yLow,yUp));


	    // Total Cross Section
	    sigmaB = muToPhotonB.getSigma();
	    System.out.println("  Total cross section B " + sigmaB);

	    sigmaK = muToPhotonK.getSigma();
	    System.out.println("  Total cross section K " + sigmaK);


	    // Inelasticity
	    yAverage = muToPhotonB.getYDSigmaDy(muToPhotonB.getYmin(),
					      muToPhotonB.getYmax())/sigmaB;
	    System.out.println("  Inelasticity B <y> " + yAverage);
	    zAverage =  muToPhotonB.getZDSigmaDZ(1.0-muToPhotonB.getYmax(),
					       1.0-muToPhotonB.getYmin())/sigmaB;
	    System.out.println("  Inelasticity B <1-y> " + zAverage);
	    System.out.println("  sum " + (yAverage+zAverage));

	    yAverage = muToPhotonK.getYDSigmaDy(muToPhotonK.getYmin(),
					      muToPhotonK.getYmax())/sigmaK;
	    System.out.println("  Inelasticity K <y> " + yAverage);
	    zAverage =  muToPhotonK.getZDSigmaDZ(1.0-muToPhotonK.getYmax(),
					       1.0-muToPhotonK.getYmin())/sigmaK;
	    System.out.println("  Inelasticity K <1-y> " + zAverage);
	    System.out.println("  sum " + (yAverage+zAverage));

	}



    }
}
