package  iceCube.uhe.event;

import iceCube.uhe.particles.*;
import iceCube.uhe.interactions.*;
import iceCube.uhe.muonModel.*;
import java.io.*;

/**

This class calculates differential flux dF/dLogE [/cm^2 sec sr]
of atm muons as a function of emg/hadron cascade energy deposited
in the detector volume.
The differential flux of mu/tau's which enter into
the detector volume after propagating in the earth is given
by PropagatingAtmMuonFlux.class in the muonModel package.
The cascade energy distribution is given by EventMatrix.class
in the Event package.

*/

public class EventAtmMuonFlux {

    static int dimension = Particle.getDimensionOfLogEnergyMatrix();
    static double logEcascadeMin = InteractionsBase.getLogEnergyProducedMinimum();
    static int offset = (int )(( Particle.getLogEnergyMinimum()-logEcascadeMin)/
			       Particle.getDeltaLogEnergy());
    public PropagatingAtmMuonFlux propMuonFlux;
    public EventMatrix eventMatrix;
    private static final double ln10 = Math.log(10.0);

    /** Constructor. PropagatingAtmMuonFlux and EventMatrix classes are
     generated. Note that both of these objects still needs to read the matix
    file to perform any futher calculation. This can be done by calling the methods
    propMuonFlux.readMatrix(DataInputStream in) and
    eventMatrix.readMatrix(DataInputStream in). You can run these method
    through this object because both propMuonFlux and eventMatrix
    are public class variables.
    */
    public EventAtmMuonFlux() throws IOException {

	// Generate PropagatingNeurinoFlux class
        propMuonFlux = new PropagatingAtmMuonFlux();
	eventMatrix = new EventMatrix();
    }

    /** Same as the constructor EventFlux(int model) except
	it reads the event matrix from the DataInputStream eventIn. 
    */
    public EventAtmMuonFlux(DataInputStream eventIn) throws IOException {

	// Generate PropagatingNeurinoFlux class
        propMuonFlux = new PropagatingAtmMuonFlux();
	eventMatrix = new EventMatrix(eventIn);
    }


    /** Calculate dF/dLogE [/cm^2 sec sr] for EMG cascade */
    public double getDFEmgCascadeDLogE(double logEcascade, double cosTheta){

	double count = 0.0;
	int jLogE = (int)((logEcascade - logEcascadeMin)/Particle.getDeltaLogEnergy());
	int iLogEmin = jLogE - offset;
	if(iLogEmin<0) iLogEmin = 0; //starts from Particle.getLogEnergyMinimum()


	for(int iLogE=iLogEmin;iLogE<dimension;iLogE++){

	    double logEprimary = Particle.getLogEnergyMinimum( ) + 
		Particle.getDeltaLogEnergy( )*(double )iLogE;
	    if(eventMatrix.getPropFlavor()==1){ //Muon
		count += propMuonFlux.getDFMuDLogE(logEprimary,cosTheta)*
		    eventMatrix.getEmgCascadeFlux(logEprimary, logEcascade);
	    }
	}
	return(count);
    }

    /** Calculate dF/dLogE [/cm^2 sec sr] for Hadron cascade */
    public double getDFHadronCascadeDLogE(double logEcascade, double cosTheta){

	double count = 0.0;
	int jLogE = (int)((logEcascade - logEcascadeMin)/Particle.getDeltaLogEnergy());
	int iLogEmin = jLogE - offset;
	if(iLogEmin<0) iLogEmin = 0; //starts from Particle.getLogEnergyMinimum()


	for(int iLogE=iLogEmin;iLogE<dimension;iLogE++){

	    double logEprimary = Particle.getLogEnergyMinimum( ) + 
		Particle.getDeltaLogEnergy( )*(double )iLogE;
	    if(eventMatrix.getPropFlavor()==1){ //Muon
		count += propMuonFlux.getDFMuDLogE(logEprimary,cosTheta)*
		    eventMatrix.getHadronCascadeFlux(logEprimary, logEcascade);
	    }
	}
	return(count);
    }


    /** Calculate dF/dLogE [/cm^2 sec sr] for total cascade */
    public double getDFTotalCascadeDLogE(double logEcascade, double cosTheta){

	double count = 0.0;
	int jLogE = (int)((logEcascade - logEcascadeMin)/Particle.getDeltaLogEnergy());
	int iLogEmin = jLogE - offset;
	if(iLogEmin<0) iLogEmin = 0; //starts from Particle.getLogEnergyMinimum()


	for(int iLogE=iLogEmin;iLogE<dimension;iLogE++){

	    double logEprimary = Particle.getLogEnergyMinimum( ) + 
		Particle.getDeltaLogEnergy( )*(double )iLogE;
	    if(eventMatrix.getPropFlavor()==1){ //Muon
		count += propMuonFlux.getDFMuDLogE(logEprimary,cosTheta)*
		    eventMatrix.getTotalCascadeFlux(logEprimary, logEcascade);
	    }
	}
	return(count);
    }


}

