package  iceCube.uhe.event;

import iceCube.uhe.particles.*;
import iceCube.uhe.interactions.*;
import iceCube.uhe.neutrinoModel.*;
import java.io.*;

/**

This class calculates differential flux dF/dLogE [/cm^2 sec sr]
of muons/taus as a function of emg/hadron cascade energy deposited
in the detector volume expected from a neutrino flux with monochromatic energy
E^2dF/dE = 10^-9 GeV/cm^2 sec sr. 
The differential flux of NuMu/NuTau's which enter into
the detector volume after propagating in the earth is given apriori
by PropagationMatrix.java in the propagation package and
the results are readout from the Data Input Stream which
is given when the method readMatrix(in) is called.
The cascade energy distribution is given by EventMatrix.class
in the Event package. This class is used for neutrino aperture
estimation. 

A difference from EventMonoEnergyFlux.java is that
this class calculates the case when mu/tau are 
generated INSIDE the volume. So The intensity
of NEUTRINOs before entering the volume is
initial flux and the cascade distribution for initial neutrinos
is used to evaluate the resultant differential flux
of mu/tau.

*/

public class NeutEventMonoEnergyFlux {

    static int dimension = Particle.getDimensionOfLogEnergyMatrix();
    int inputParticle, outputParticle;
    int inLogE = -1;
    double logYmin,logYmax;
    double[][] FnuEToNuE,FnuMuToNuE,FnuTauToNuE,FmuToNuE,FtauToNuE;
    double[][]           FnuMuToNuMu,FnuTauToNuMu,FmuToNuMu,FtauToNuMu;
    double[][]           FnuMuToNuTau,FnuTauToNuTau,FmuToNuTau,FtauToNuTau;
    double[][] FnuEToE,FnuMuToE,FnuTauToE,FmuToE,FtauToE;
    double[][]           FnuMuToMu,FnuTauToMu,FmuToMu,FtauToMu;
    double[][]           FnuMuToTau,FnuTauToTau,FmuToTau,FtauToTau;
    double[][] FnuEToHadron,FnuMuToHadron,FnuTauToHadron,FmuToHadron,FtauToHadron;
    private static final double ln10 = Math.log(10.0);
    static double logEcascadeMin = InteractionsBase.getLogEnergyProducedMinimum();
    static int offset = (int )(( Particle.getLogEnergyMinimum()-logEcascadeMin)/
			       Particle.getDeltaLogEnergy());
    public static final double E2dFdE = 1.0e-9; // Intensity of Flux. [GeV/cm sec sr]
    public EventMatrix eventMatrix;

    /** Constructor. Array for propagationMatrix and  EventMatrix classes are
     generated. Note th objects still needs to read the matix
    file to perform any futher calculation. This can be done by calling the methods
    readMatrix(DataInputStream in).
    */
    public NeutEventMonoEnergyFlux(DataInputStream eventIn) throws IOException {

	// Generate PropagatingNeurinoFlux class
        FnuEToNuE= new double[dimension][dimension];
        FnuEToE= new double[dimension][dimension];
        FnuEToHadron= new double[dimension][dimension];
        FnuMuToNuE= new double[dimension][dimension];
        FnuMuToNuMu= new double[dimension][dimension];
        FnuMuToNuTau= new double[dimension][dimension];
        FnuMuToE= new double[dimension][dimension];
        FnuMuToMu= new double[dimension][dimension];
        FnuMuToTau= new double[dimension][dimension];
        FnuMuToHadron= new double[dimension][dimension];
        FnuTauToNuE= new double[dimension][dimension];
        FnuTauToNuMu= new double[dimension][dimension];
        FnuTauToNuTau= new double[dimension][dimension];
        FnuTauToE= new double[dimension][dimension];
        FnuTauToMu= new double[dimension][dimension];
        FnuTauToTau= new double[dimension][dimension];
        FnuTauToHadron= new double[dimension][dimension];
        FmuToNuE= new double[dimension][dimension];
        FmuToNuMu= new double[dimension][dimension];
        FmuToNuTau= new double[dimension][dimension];
        FmuToE= new double[dimension][dimension];
        FmuToMu= new double[dimension][dimension];
        FmuToTau= new double[dimension][dimension];
        FmuToHadron= new double[dimension][dimension];
        FtauToNuE= new double[dimension][dimension];
        FtauToNuMu= new double[dimension][dimension];
        FtauToNuTau= new double[dimension][dimension];
        FtauToE= new double[dimension][dimension];
        FtauToMu= new double[dimension][dimension];
        FtauToTau= new double[dimension][dimension];
        FtauToHadron= new double[dimension][dimension];

	eventMatrix = new EventMatrix(eventIn);
    }


    /** Read the calculated propagation matrix */
    public void readMatrix(DataInputStream in) throws IOException {
        int iLogE;
        for(iLogE=0;iLogE<dimension;iLogE++){
            int jLogE;
            for(jLogE=0;jLogE<=iLogE;jLogE++){

                FnuEToNuE[iLogE][jLogE]= in.readDouble( );
                FnuEToE[iLogE][jLogE]= in.readDouble( );
                FnuEToHadron[iLogE][jLogE]= in.readDouble( );
                FnuMuToNuE[iLogE][jLogE]= in.readDouble( );
                FnuMuToNuMu[iLogE][jLogE]= in.readDouble( );
                FnuMuToNuTau[iLogE][jLogE]= in.readDouble( );
                FnuMuToE[iLogE][jLogE]= in.readDouble( );
                FnuMuToMu[iLogE][jLogE]= in.readDouble( );
                FnuMuToTau[iLogE][jLogE]= in.readDouble( );
                FnuMuToHadron[iLogE][jLogE]= in.readDouble( );
                FnuTauToNuE[iLogE][jLogE]= in.readDouble( );
                FnuTauToNuMu[iLogE][jLogE]= in.readDouble( );
                FnuTauToNuTau[iLogE][jLogE]= in.readDouble( );
                FnuTauToE[iLogE][jLogE]= in.readDouble( );
                FnuTauToMu[iLogE][jLogE]= in.readDouble( );
                FnuTauToTau[iLogE][jLogE]= in.readDouble( );
                FnuTauToHadron[iLogE][jLogE]= in.readDouble( );
                FmuToNuE[iLogE][jLogE]= in.readDouble( );
                FmuToNuMu[iLogE][jLogE]= in.readDouble( );
                FmuToNuTau[iLogE][jLogE]= in.readDouble( );
                FmuToE[iLogE][jLogE]= in.readDouble( );
                FmuToMu[iLogE][jLogE]= in.readDouble( );
                FmuToTau[iLogE][jLogE]= in.readDouble( );
                FmuToHadron[iLogE][jLogE]= in.readDouble( );
                FtauToNuE[iLogE][jLogE]= in.readDouble( );
                FtauToNuMu[iLogE][jLogE]= in.readDouble( );
                FtauToNuTau[iLogE][jLogE]= in.readDouble( );
                FtauToE[iLogE][jLogE]= in.readDouble( );
                FtauToMu[iLogE][jLogE]= in.readDouble( );
                FtauToTau[iLogE][jLogE]= in.readDouble( );
                FtauToHadron[iLogE][jLogE]= in.readDouble( );

            }
        }

        in.close( );

    }


    /** Calculate dF/dLogE [/cm^2 sec sr] for EMG cascade 
	neutrnoFlavor is defined by the Particle class. i.e. 0 for nu-e 1 for nu-mu etc.
     */
    public double getDFEmgCascadeDLogE(int neutrinoFlavor, double logEneutrino, double logEcascade){

	double count = 0.0;
	int jLogE = (int)((logEcascade - logEcascadeMin)/Particle.getDeltaLogEnergy());
	int iLogEmin = jLogE - offset;
	if(iLogEmin<0) iLogEmin = 0; //starts from Particle.getLogEnergyMinimum()
	int iLogEmax = (int)((logEneutrino - Particle.getLogEnergyMinimum())/
			     Particle.getDeltaLogEnergy());


	for(int iLogE=iLogEmin;iLogE<=iLogEmax;iLogE++){

	    double logEprimary = Particle.getLogEnergyMinimum( ) + 
		Particle.getDeltaLogEnergy( )*(double )iLogE;
	    if(eventMatrix.getPropFlavor()==1 && neutrinoFlavor == 1){ // nu-mu to Muon
		count += FnuMuToNuMu[iLogEmax][iLogE]*
		    eventMatrix.getEmgCascadeFlux(logEprimary, logEcascade);
	    }if(eventMatrix.getPropFlavor()==1 && neutrinoFlavor == 2){ // nu-tau to Muon
		count += FnuTauToNuMu[iLogEmax][iLogE]*
		    eventMatrix.getEmgCascadeFlux(logEprimary, logEcascade);
	    }if(eventMatrix.getPropFlavor()==2 && neutrinoFlavor == 1){ // nu-mu to Tau
		count += FnuMuToNuTau[iLogEmax][iLogE]*
		    eventMatrix.getEmgCascadeFlux(logEprimary, logEcascade);
	    }if(eventMatrix.getPropFlavor()==2 && neutrinoFlavor == 2){ // nu-tau to Tau
		count += FnuTauToNuTau[iLogEmax][iLogE]*
		    eventMatrix.getEmgCascadeFlux(logEprimary, logEcascade);
	    }
	}
        double energy = Math.pow(10.0,logEneutrino); // [GeV]
	return(E2dFdE*count*ln10/energy);
    }

    /** Calculate dF/dLogE [/cm^2 sec sr] for Hadron cascade 
	neutrnoFlavor is defined by the Particle class. i.e. 0 for nu-e 1 for nu-mu etc.
     */
    public double getDFHadronCascadeDLogE(int neutrinoFlavor, double logEneutrino, double logEcascade){

	double count = 0.0;
	int jLogE = (int)((logEcascade - logEcascadeMin)/Particle.getDeltaLogEnergy());
	int iLogEmin = jLogE - offset;
	if(iLogEmin<0) iLogEmin = 0; //starts from Particle.getLogEnergyMinimum()
	int iLogEmax = (int)((logEneutrino - Particle.getLogEnergyMinimum())/
			     Particle.getDeltaLogEnergy());


	for(int iLogE=iLogEmin;iLogE<=iLogEmax;iLogE++){

	    double logEprimary = Particle.getLogEnergyMinimum( ) + 
		Particle.getDeltaLogEnergy( )*(double )iLogE;
	    if(eventMatrix.getPropFlavor()==1 && neutrinoFlavor == 1){ // nu-mu to Muon
		count += FnuMuToNuMu[iLogEmax][iLogE]*
		    eventMatrix.getHadronCascadeFlux(logEprimary, logEcascade);
	    }if(eventMatrix.getPropFlavor()==1 && neutrinoFlavor == 2){ // nu-tau to Muon
		count += FnuTauToNuMu[iLogEmax][iLogE]*
		    eventMatrix.getHadronCascadeFlux(logEprimary, logEcascade);
	    }if(eventMatrix.getPropFlavor()==2 && neutrinoFlavor == 1){ // nu-mu to Tau
		count += FnuMuToNuTau[iLogEmax][iLogE]*
		    eventMatrix.getHadronCascadeFlux(logEprimary, logEcascade);
	    }if(eventMatrix.getPropFlavor()==2 && neutrinoFlavor == 2){ // nu-tau to Tau
		count += FnuTauToNuTau[iLogEmax][iLogE]*
		    eventMatrix.getHadronCascadeFlux(logEprimary, logEcascade);
	    }
	}
        double energy = Math.pow(10.0,logEneutrino); // [GeV]
	return(E2dFdE*count*ln10/energy);
    }

    /** Calculate dF/dLogE [/cm^2 sec sr] for total cascade 
	neutrnoFlavor is defined by the Particle class. i.e. 0 for nu-e 1 for nu-mu etc.
     */
    public double getDFTotalCascadeDLogE(int neutrinoFlavor, double logEneutrino, double logEcascade){

	double count = 0.0;
	int jLogE = (int)((logEcascade - logEcascadeMin)/Particle.getDeltaLogEnergy());
	int iLogEmin = jLogE - offset;
	if(iLogEmin<0) iLogEmin = 0; //starts from Particle.getLogEnergyMinimum()
	int iLogEmax = (int)((logEneutrino - Particle.getLogEnergyMinimum())/
			     Particle.getDeltaLogEnergy());


	for(int iLogE=iLogEmin;iLogE<=iLogEmax;iLogE++){

	    double logEprimary = Particle.getLogEnergyMinimum( ) + 
		Particle.getDeltaLogEnergy( )*(double )iLogE;
	    if(eventMatrix.getPropFlavor()==1 && neutrinoFlavor == 1){ // nu-mu to Muon
		count += FnuMuToNuMu[iLogEmax][iLogE]*
		    eventMatrix.getTotalCascadeFlux(logEprimary, logEcascade);
	    }if(eventMatrix.getPropFlavor()==1 && neutrinoFlavor == 2){ // nu-tau to Muon
		count += FnuTauToNuMu[iLogEmax][iLogE]*
		    eventMatrix.getTotalCascadeFlux(logEprimary, logEcascade);
	    }if(eventMatrix.getPropFlavor()==2 && neutrinoFlavor == 1){ // nu-mu to Tau
		count += FnuMuToNuTau[iLogEmax][iLogE]*
		    eventMatrix.getTotalCascadeFlux(logEprimary, logEcascade);
	    }if(eventMatrix.getPropFlavor()==2 && neutrinoFlavor == 2){ // nu-tau to Tau
		count += FnuTauToNuTau[iLogEmax][iLogE]*
		    eventMatrix.getTotalCascadeFlux(logEprimary, logEcascade);
	    }
	}
        double energy = Math.pow(10.0,logEneutrino); // [GeV]
	return(E2dFdE*count*ln10/energy);
    }



}
