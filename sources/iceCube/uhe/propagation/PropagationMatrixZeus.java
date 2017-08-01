package iceCube.uhe.propagation;

import iceCube.uhe.particles.*;
import iceCube.uhe.interactions.*;
import iceCube.uhe.decay.*;
import iceCube.uhe.points.*;
import java.io.*;

/**                                                                                                                                    
PropagationMatrix using the Zeus based neutrino-nucleon cross section

*/

public class PropagationMatrixZeus extends PropagationMatrix {

    protected static String nuCCMtxObjectFile = "ENeutrinoChargeZeusNewMtx";
    protected static String nuNCMtxObjectFile = "ENeutrinoNeutralZeusNewMtx";

    public PropagationMatrixZeus(Particle nuE, Particle nuMu, Particle nuTau,
                             Particle e,   Particle mu,   Particle tau,
                             Particle pi,  ParticlePoint s,
                             int interactionsSwitch, int decaySwitch,
                             double neutrinoFactor)
	throws IOException{

	super(nuE, nuMu, nuTau, e, mu, tau, pi, s, interactionsSwitch, decaySwitch, neutrinoFactor,
	      nuCCMtxObjectFile,nuNCMtxObjectFile);
    }

    public PropagationMatrixZeus(Particle nuE, Particle nuMu, Particle nuTau,
                             Particle e,   Particle mu,   Particle tau,
                             Particle pi,  ParticlePoint s,
			     int interactionsSwitch, int decaySwitch) throws IOException{

	super(nuE, nuMu, nuTau, e, mu, tau, pi, s, interactionsSwitch, decaySwitch,
	      nuCCMtxObjectFile,nuNCMtxObjectFile);
    }

}


