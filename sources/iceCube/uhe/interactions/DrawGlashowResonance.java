package iceCube.uhe.interactions;

import iceCube.uhe.particles.*;
import iceCube.uhe.interactions.*;
import iceCube.uhe.points.*;
import java.io.*;
import java.util.*;

public class DrawGlashowResonance {

    public static void main(String[] args) throws IOException {

        // Generate the ParticlePoint class.
        ParticlePoint s = new ParticlePoint(0.0, 5.0*Math.PI/180.0,0); //ice
	double massNumber = 0.0;
        for(int i=0;i<s.NumberOfSpecies[s.getMaterialNumber( )];i++){
	    massNumber += s.getNumberOfAtoms(i)*s.getAtomicNumber(i);
	}

	// Read the serialized object of the Charged Current cross section matrix
	// for the reference
	String ccMatrixFileName = "iceCube/uhe/interactions/ice/ENeutrinoChargeMtx";
	InputStream in = ClassLoader.getSystemResourceAsStream(ccMatrixFileName);
	InteractionsMatrix nuCCMtx = 
	    InteractionsMatrixInput.inputInteractionsMatrix(in);
	in.close( );

	// Generate the Glashow Resonance object
	GlashowResonanceLeptonic grLepton =  new GlashowResonanceLeptonic(s,1);
	               // muonic decay -- produced flavor =1
	grLepton.calculateCrossSectionAsPerElectron();
	GlashowResonanceHadronic grHadron =  new GlashowResonanceHadronic(s);
	grHadron.calculateCrossSectionAsPerElectron();


	System.out.println("titx Energy [GeV]");
	System.out.println("tity Cross Section [cm^2!]");
	System.out.println("gwin 0.2 0.9 0.2 0.9");
	System.out.println("scal 1.0e5 1.0e12 1.0e-36 1.0e-30");

	// CC reaction
	int iLogE;
	for(iLogE=0;iLogE< Particle.getDimensionOfLogEnergyMatrix();iLogE++){

	    // Total Cross Section
	    double sigma = nuCCMtx.getSigmaMatrix(iLogE);

	    double logE = Particle.getLogEnergyMinimum( ) + 
		Particle.getDeltaLogEnergy( )*(double )iLogE;
	    double energy = Math.pow(10.0,logE);
	    System.out.println("data " + energy + " 0.0 " + sigma + " 0.0");
	}
	System.out.println("logx");
	System.out.println("logy");
	System.out.println("join");
	System.out.println("disp");
	System.out.println("cont");

	// Glashow Resonance reaction
	for(iLogE=0;iLogE< Particle.getDimensionOfLogEnergyMatrix();iLogE++){

	    // Total Cross Section
	    grLepton.setIncidentParticleEnergy(iLogE);

	    double sigma = grLepton.getSigma();

	    double logE = Particle.getLogEnergyMinimum( ) + 
		Particle.getDeltaLogEnergy( )*(double )iLogE;
	    double energy = Math.pow(10.0,logE);
	    System.out.println("data " + energy + " 0.0 " + sigma + " 0.0");
	}
	System.out.println("logx");
	System.out.println("logy");
	System.out.println("join");
	System.out.println("disp");
	System.out.println("cont");
	// W Resonance behavior checks
	double wResEnergy = grLepton.massW*grLepton.massW/(2.0*grLepton.massE);
	grLepton.setIncidentParticleEnergy(wResEnergy);
	grHadron.setIncidentParticleEnergy(wResEnergy);
	double sigmaRes = grLepton.getSigma()+grHadron.getSigma();
	double zRes1PeV = 1.0e6/wResEnergy;
	double sigmaRes1PeV = grLepton.integralDSigmaDz(zRes1PeV*0.65,zRes1PeV*1.35);
	System.out.println(" Res Energy " + wResEnergy + " sigma = " + sigmaRes + 
			   " 1PeV/ResE = " + zRes1PeV + " part.sigma = " + sigmaRes1PeV);

	for(iLogE=0;iLogE< Particle.getDimensionOfLogEnergyMatrix();iLogE++){

	    // Total Cross Section
	    grHadron.setIncidentParticleEnergy(iLogE);

	    double sigma = grHadron.getSigma();

	    double logE = Particle.getLogEnergyMinimum( ) + 
		Particle.getDeltaLogEnergy( )*(double )iLogE;
	    double energy = Math.pow(10.0,logE);
	    System.out.println("data " + energy + " 0.0 " + sigma + " 0.0");
	}
	System.out.println("logx");
	System.out.println("logy");
	System.out.println("join");
	System.out.println("disp");
	System.out.println("endg");


    }

}
