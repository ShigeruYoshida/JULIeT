package  iceCube.uhe.neutrinoModel;

import iceCube.uhe.particles.*;
import iceCube.uhe.interactions.*;
import iceCube.uhe.points.*;
import iceCube.uhe.propagation.*;
import java.io.*;
import java.util.*;

/**
   This class calculates differential flux dF/dLogE [/cm^2 sec sr]
   of neutrinos and charged leptons after propagation in the earth
   for a given model of primary cosmic neutrino production in the Universe.
   The primary flux of UHE cosmic neutrinos is given by the NeutrinoFlux
   class with taking into account the neutrino oscillation effect. 
   The particle propagation calculation is done by protected member 
   QuickNeutrinoPropagator propagator instead of the propagation matrix
   as the super class PropagatingNeutrinoFlux does.

   The argument for the constructor, "model" is for the NeutrinoFlux class.
   Consult the detais to the API document of the NeutrinoFlux.java
   in this package.
*/

public class QuickPropagatingNeutrinoFlux extends PropagatingNeutrinoFlux {

    /** Neutrino Quick Propagator */
    protected NeutrinoQuickPropagator propagator = null;
    private boolean hasPropagated = false;

    /** Constructor */
    public QuickPropagatingNeutrinoFlux(int model) throws IOException {
	super(model,false);
    }


    /** Calculate dF/dLogE [/cm^2 sec sr] for particle (flavor doublet).
	The neutrino flux is calculated for the model
        you specified in the argument.
    */
    private double getDFDLogE(int model, int jLogE, int flavor, int doublet) throws IOException {

	NeutrinoFlux neutFlux = getNeutrinoFluxFromList(model);
	if(neutFlux == null){ // This is a new model. Needs to generate the object
	    generateNeutrinoFluxObject(model);
	    neutFlux = getNeutrinoFluxFromList(model);
	    System.err.println("JULIeT: new NeutrinoFlux object generated with model "
			       + model);
	}


	double logE;
	double count = 0.0;
	int iLogE;
	double epsilon = 1.0e-4;
	double logOutputE = Particle.getLogEnergyMinimum( ) + 
		Particle.getDeltaLogEnergy( )*(double )jLogE+epsilon;

	for(iLogE=jLogE;iLogE<dimension;iLogE++){

	    logE = Particle.getLogEnergyMinimum( ) + 
		Particle.getDeltaLogEnergy( )*(double )iLogE+epsilon;


	    // get nuetrino flux at the earth surface after propagation in space
	    // with neutrino oscillation
            double[] nuflux_osci = neutFlux.getDFDLogEwzOsci(logE);

	    if(propagator.getDF(0,logE,flavor,doublet,logOutputE)>0.0){ // from nuE
		count += nuflux_osci[0]*propagator.getDF(0,logE,flavor,doublet,logOutputE);
	    }
	    if(propagator.getDF(1,logE,flavor,doublet,logOutputE)>0.0){ // from nuMu
		count += nuflux_osci[1]*propagator.getDF(1,logE,flavor,doublet,logOutputE);
	    }
	    if(propagator.getDF(2,logE,flavor,doublet,logOutputE)>0.0){ // from nuTau
		count += nuflux_osci[2]*propagator.getDF(2,logE,flavor,doublet,logOutputE);
	    }

	}
	return(count);
    }

    /** Calculate dF/dLogE [/cm^2 sec sr] for nu-e 
	The neutrino flux is calculated for the model
        you specified in the argument.
     */
    public double getDFNuEDLogE(int model, int jLogE) throws IOException {

	return(getDFDLogE(model,jLogE,0,0));
    }


    /** Calculate dF/dLogE [/cm^2 sec sr] for nu-mu 
	The neutrino flux is calculated for the model
        you specified in the argument.
     */
    public double getDFNuMuDLogE(int model, int jLogE) throws IOException {

	return(getDFDLogE(model,jLogE,1,0));
    }
 


    /** Calculate dF/dLogE [/cm^2 sec sr] for nu-tau 
	The neutrino flux is calculated for the model
        you specified in calling the constructor.
     */
    public double getDFNuTauDLogE(int model, int jLogE) throws IOException {

	return(getDFDLogE(model,jLogE,2,0));
    }
 

    /** Calculate dF/dLogE [/cm^2 sec sr] for mu 
	The neutrino flux is calculated for the model
        you specified in calling the constructor.
     */
    public double getDFMuDLogE(int model, int jLogE) throws IOException {

	return(getDFDLogE(model,jLogE,1,1));
    }
 

    /** Calculate dF/dLogE [/cm^2 sec sr] for tau 
	The neutrino flux is calculated for the model
        you specified in the argument.
     */
    public double getDFTauDLogE(int model, int jLogE) throws IOException {

	return(getDFDLogE(model,jLogE,2,1));
    }
 
}
