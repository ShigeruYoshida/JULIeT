package iceCube.uhe.propagation;

import iceCube.uhe.particles.*;
import iceCube.uhe.interactions.*;
import iceCube.uhe.points.*;
import java.io.*;

public class RunPropagationMatrixZeus extends RunPropagationMatrix{

    /** 
	<pre>
	Constructor. Generate the relevant particles objects and 
	the propagation matrix object. The medium is assumed to be
        the standard rock.

	double nadirAngle:   Nadir angle [deg] of trajectory of the incoming particles.
	int inSwitch:        Interaction Switch to turn on/off the individual interaction
	                     channel.
	int decaySwitch:     Decay Switch to turn on/off the individual decay channel.
	int mediumNumber:    Medium number 0 ice 1 rock
	</pre>
    */
    public RunPropagationMatrixZeus(double nadirAngle,
				int intSwitch,int decaySwitch, int mediumNumber) throws IOException {

	/** For Glashow Resonance **/
	//if((0<intSwitch && intSwitch <256) && (0<decaySwitch && decaySwitch <256)){
	if((0<intSwitch && intSwitch <512) && (0<decaySwitch && decaySwitch <512)){

	    // Generate the ParticlePoint class.
	    s = new ParticlePoint(0.0, nadirAngle*Math.PI/180.0,mediumNumber);
	    System.err.println("Axis length of propagation trajectory " + s.getAxisLength( )
			       + " [cm]");

	    // Generate Particles involved
	    nuE = new Particle(0,0);    //Electron Neutrinos
	    nuMu = new Particle(1,0);   //Muon Neutrinos
	    nuTau = new Particle(2,0);  //Tau  Neutrinos
	    e = new Particle(0,1);      //Electron
	    mu = new Particle(1,1);     //Muon
	    tau = new Particle(2,1);    //Tauon
	    pi =  new Particle(3,1);    // hadron.. pi+

	    // Generate the Propagation Matrix
	    propMtx = new PropagationMatrixZeus(nuE,nuMu,nuTau,e,mu,tau,pi,
					    s,intSwitch,decaySwitch);
	}
    }

    public RunPropagationMatrixZeus(double nadirAngle,
				int intSwitch,int decaySwitch, int mediumNumber, 
				double neutrinoFactor) throws IOException {

	/** For Glashow Resonance **/
	//if((0<intSwitch && intSwitch <256) && (0<decaySwitch && decaySwitch <256)){
	if((0<intSwitch && intSwitch <512) && (0<decaySwitch && decaySwitch <512)){

	    // Generate the ParticlePoint class.
	    s = new ParticlePoint(0.0, nadirAngle*Math.PI/180.0,mediumNumber);
	    System.err.println("Axis length of propagation trajectory " + s.getAxisLength( )
			       + " [cm]");

	    // Generate Particles involved
	    nuE = new Particle(0,0);    //Electron Neutrinos
	    nuMu = new Particle(1,0);   //Muon Neutrinos
	    nuTau = new Particle(2,0);  //Tau  Neutrinos
	    e = new Particle(0,1);      //Electron
	    mu = new Particle(1,1);     //Muon
	    tau = new Particle(2,1);    //Tauon
	    pi =  new Particle(3,1);    // hadron.. pi+

	    // Generate the Propagation Matrix
	    System.err.println("Setting the neutrino factor " + neutrinoFactor);
	    propMtx = new PropagationMatrixZeus(nuE,nuMu,nuTau,e,mu,tau,pi,
					    s,intSwitch,decaySwitch,neutrinoFactor);
	}
    }

}
