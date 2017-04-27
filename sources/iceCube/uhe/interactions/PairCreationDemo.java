package iceCube.uhe.interactions;

import iceCube.uhe.particles.*;
import iceCube.uhe.interactions.*;
import iceCube.uhe.points.*;

/** PairCreation.class Demo program */
public class PairCreationDemo {

    public static void main(String[] args){

        int flavor =1;         // Muon
        int doublet =1;
	int recoilFlavor = 0;  // Electron
	int material = 1;      // Rock
	double energy = 1.0e9; // 1EeV= 10^9 GeV
	double y = 0.5;
	double yLogLow = 0.0;
	double yLogUp = 0.0;
	double yUp = 0.0;
	double yLow = 0.0;
	double yPlus = 0.1;
	double epsilon = 1.0e-2;
	double energyCut;
	double sigma = 1.0;
	double sigmaPair = 1.0;
	double yAverage = 0.5;
	double zAverage = 0.5;
	double yPlusAverage = 0.5;


        if(args.length!=3){
            System.out.println("Usage: PairCreationDemo y yPlus epsilon");
	    System.exit(0);
        }else{
            y = Double.valueOf(args[0]).doubleValue();
            yPlus = Double.valueOf(args[1]).doubleValue();
            epsilon = Double.valueOf(args[2]).doubleValue();
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

	//Generate the pair creation class
	PairCreation muToEpair = new PairCreation(mu, s, recoilFlavor);
	System.out.println(muToEpair.interactionName( ));


	numRecipes.Integration.setRelativeAccuracy(epsilon);
	// integration accuracy for save CPUtime

	// Cross Section 
	energyCut = Math.pow(10.0,mu.getLogEnergyMinimum( )-1.0); //GeV
	muToEpair.setEnergyCut(energyCut);
	System.out.println("epsilon " + epsilon + " energyCut " + 
			   energyCut + " GeV");
	int iLogE;
	for(iLogE=0;iLogE<=500;iLogE+=100){
	    muToEpair.setIncidentParticleEnergy(iLogE);
	    System.out.println("The Incident energy " + 
			   muToEpair.getIncidentParticleEnergy() + " GeV");

	    // Differntial Cross Section
	    System.out.println("  dSigma/dY(" + y + ") " +
			       muToEpair.getDSigmaDy(y));
	    // Differential Cross Section over electron pair
	    System.out.println("  dSigma/dYplus(" +  yPlus + ") " +
			       muToEpair.getDSigmaDyPlus(yPlus));
	    // Partial Cross Section
	    yLow = muToEpair.getYmin( );
	    yLogLow = Math.log(yLow)/Math.log(10.0);
	    yLogUp = yLogLow + mu.getDeltaLogEnergy( );
	    yUp = Math.pow(10.0,yLogUp);
	    System.out.println("  YdSigma/dY(" + yLow + ") " +
			       muToEpair.getYDSigmaDy(yLow,yUp));
	    yUp = muToEpair.getYmax( );
	    yLogUp = Math.log(yUp)/Math.log(10.0);
	    yLogLow = yLogUp - mu.getDeltaLogEnergy( );
	    yLow = Math.pow(10.0,yLogLow);
	    System.out.println("  YdSigma/dY(" + yLow + ") " +
			       muToEpair.getYDSigmaDy(yLow,yUp));
	    // Total Cross Section
	    sigma = muToEpair.getSigma();
	    System.out.println("  Total cross section " + sigma);
	    // Total Cross Section by integration over energy of electron
	    sigmaPair = muToEpair.integralDSigmaDyPlus(muToEpair.getYPlusMin(),
						 muToEpair.getYPlusMax());
	    System.out.println("  Total cross section Plus " + sigmaPair);
	    System.out.println("  Ratio " + sigmaPair/sigma);
	    // Inelasticity
	    yAverage = muToEpair.getYDSigmaDy(muToEpair.getYmin(),
					      muToEpair.getYmax())/sigma;
	    System.out.println("  Inelasticity <y> " + yAverage);
	    zAverage =  muToEpair.getZDSigmaDZ(1.0-muToEpair.getYmax(),
					       1.0-muToEpair.getYmin())/sigma;
	    System.out.println("  Inelasticity <1-y> " + zAverage);
	    System.out.println("  sum " + (yAverage+zAverage));
	    yPlusAverage = muToEpair.getYPlusDSigmaDyPlus(muToEpair.getYPlusMin(),
					      muToEpair.getYPlusMax())/sigmaPair;
	    System.out.println("  Inelasticity <yPlus> " + yPlusAverage);
	}



    }
}
