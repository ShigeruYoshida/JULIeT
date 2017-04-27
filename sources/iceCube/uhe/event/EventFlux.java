package  iceCube.uhe.event;

import iceCube.uhe.particles.*;
import iceCube.uhe.interactions.*;
import iceCube.uhe.neutrinoModel.*;
import java.io.*;

/**

This class calculates differential flux dF/dLogE [/cm^2 sec sr]
of muons/taus as a function of emg/hadron cascade energy deposited
in the detector volume expected irom a given UHE neutrino model
such as GZK. The differential flux of mu/tau's which enter into
the detector volume after propagating in the earth is given
by PropagatingNeutrinoFlux.class in the neutrinoModel package.
The cascade energy distribution is given by EventMatrix.class
in the Event package.

*/

public class EventFlux {

    static int dimension = Particle.getDimensionOfLogEnergyMatrix();
    static double logEcascadeMin = InteractionsBase.getLogEnergyProducedMinimum();
    static int offset = (int )(( Particle.getLogEnergyMinimum()-logEcascadeMin)/
			       Particle.getDeltaLogEnergy());
    public PropagatingNeutrinoFlux propLeptonFlux;
    public EventMatrix eventMatrix;
    private static final double ln10 = Math.log(10.0);

    /** Constructor. PropagatingNeutrinoFlux and EventMatrix classes are
     generated. Note that both of these objects still needs to read the matix
    file to perform any futher calculation. This can be done by calling the methods
    propLeptonFlux.readMatrix(DataInputStream in) and
    eventMatrix.readMatrix(DataInputStream in). You can run these method
    through this object because both propLeptonFlux and eventMatrix
    are public class variables.
    */
    public EventFlux(int model) throws IOException {

	// Generate PropagatingNeurinoFlux class
        propLeptonFlux = new PropagatingNeutrinoFlux(model);
	eventMatrix = new EventMatrix();
    }

    /** Same as the constructor EventFlux(int model) except
	it reads the event matrix from the DataInputStream eventIn. 
    */
    public EventFlux(int model, DataInputStream eventIn) throws IOException {

	// Generate PropagatingNeurinoFlux class
        propLeptonFlux = new PropagatingNeutrinoFlux(model);
	eventMatrix = new EventMatrix(eventIn);
    }


    /** Calculate dF/dLogE [/cm^2 sec sr] for EMG cascade */
    public double getDFEmgCascadeDLogE(double logEcascade) throws IOException {

	double count = 0.0;
	int jLogE = (int)((logEcascade - logEcascadeMin)/Particle.getDeltaLogEnergy());
	int iLogEmin = jLogE - offset;
	if(iLogEmin<0) iLogEmin = 0; //starts from Particle.getLogEnergyMinimum()


	for(int iLogE=iLogEmin;iLogE<dimension;iLogE++){

	    double logEprimary = Particle.getLogEnergyMinimum( ) + 
		Particle.getDeltaLogEnergy( )*(double )iLogE;
	    if(eventMatrix.getPropFlavor()==1){ //Muon
		count += propLeptonFlux.getDFMuDLogE(iLogE)*
		    eventMatrix.getEmgCascadeFlux(logEprimary, logEcascade);
	    }
	    if(eventMatrix.getPropFlavor()==2){ //Tau
		count += propLeptonFlux.getDFTauDLogE(iLogE)*
		    eventMatrix.getEmgCascadeFlux(logEprimary, logEcascade);
	    }
	}
	return(count);
    }

    /** Calculate dF/dLogE [/cm^2 sec sr] for Hadron cascade */
    public double getDFHadronCascadeDLogE(double logEcascade)  throws IOException {

	double count = 0.0;
	int jLogE = (int)((logEcascade - logEcascadeMin)/Particle.getDeltaLogEnergy());
	int iLogEmin = jLogE - offset;
	if(iLogEmin<0) iLogEmin = 0; //starts from Particle.getLogEnergyMinimum()


	for(int iLogE=iLogEmin;iLogE<dimension;iLogE++){

	    double logEprimary = Particle.getLogEnergyMinimum( ) + 
		Particle.getDeltaLogEnergy( )*(double )iLogE;
	    if(eventMatrix.getPropFlavor()==1){ //Muon
		count += propLeptonFlux.getDFMuDLogE(iLogE)*
		    eventMatrix.getHadronCascadeFlux(logEprimary, logEcascade);
	    }
	    if(eventMatrix.getPropFlavor()==2){ //Tau
		count += propLeptonFlux.getDFTauDLogE(iLogE)*
		    eventMatrix.getHadronCascadeFlux(logEprimary, logEcascade);
	    }
	}
	return(count);
    }


    /** Calculate dF/dLogE [/cm^2 sec sr] for total cascade */
    public double getDFTotalCascadeDLogE(double logEcascade) throws IOException {

	double count = 0.0;
	int jLogE = (int)((logEcascade - logEcascadeMin)/Particle.getDeltaLogEnergy());
	int iLogEmin = jLogE - offset;
	if(iLogEmin<0) iLogEmin = 0; //starts from Particle.getLogEnergyMinimum()


	for(int iLogE=iLogEmin;iLogE<dimension;iLogE++){

	    double logEprimary = Particle.getLogEnergyMinimum( ) + 
		Particle.getDeltaLogEnergy( )*(double )iLogE;
	    if(eventMatrix.getPropFlavor()==1){ //Muon
		count += propLeptonFlux.getDFMuDLogE(iLogE)*
		    eventMatrix.getTotalCascadeFlux(logEprimary, logEcascade);
	    }
	    if(eventMatrix.getPropFlavor()==2){ //Tau
		count += propLeptonFlux.getDFTauDLogE(iLogE)*
		    eventMatrix.getTotalCascadeFlux(logEprimary, logEcascade);
	    }
	}
	return(count);
    }


}

