package iceCube.uhe.propagation;

import iceCube.uhe.particles.*;
import iceCube.uhe.interactions.*;
import iceCube.uhe.points.*;
import java.io.*;

/****** 
<pre>
	Run the Propagation matrix to calculate the energy distribution
	and the flux of particles after propagation in the Earth. 

	The Incoming Nadir angle, the switches to control
	the interaction and decay channels involved, and the file name
	to record the calculated matrix are read thorugh the arguments.

interactionsSwitch
bit8    bit7     bit6      bit5     bit4     bit3     bit2     bit1     bit0
GlaRes  LepWeak  PhotoNucl Bremss  KnockOn  PairCHeavy PairC  Neutral  Charged

decaySwitch
bit7-2   bit 1     bit0
Reserv.  TauDecay  MuDecay
</pre>
*/

public class RunPropagationMatrix {

    /** Electron-Neutrino. */
    Particle nuE = null;   
    /** Muon-Neutrino. */
    Particle nuMu = null;
    /** Tau-Neutrino. */
    Particle nuTau = null;
    /** Electrons. */
    Particle e = null;
    /** Muons. */
    Particle mu = null;
    /** Taus. */
    Particle tau = null;
    /** pions. */
    Particle pi = null;

    /** ParticlePoint objects to define the particle location. */
    ParticlePoint s = null;

    /** Propagation Matrix object. */
    PropagationMatrix propMtx = null; 

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
    public RunPropagationMatrix(double nadirAngle,
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
	    propMtx = new PropagationMatrix(nuE,nuMu,nuTau,e,mu,tau,pi,
					    s,intSwitch,decaySwitch);
	}
    }

    /** Constructor to do nothing - for subclass modules */
    public RunPropagationMatrix(){};

    public RunPropagationMatrix(double nadirAngle,
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
	    propMtx = new PropagationMatrix(nuE,nuMu,nuTau,e,mu,tau,pi,
					    s,intSwitch,decaySwitch,neutrinoFactor);
	}
    }



    /** 
	<pre>
	Trace particles running over a finite distance 
	from the current location l to the length of the trajectory.
	This method is actually doing the following processes:

	(1) Calculate how many times the infinitesimal propagation
	    (dX [g/cm^2]) has to be invoked for the entire full
	    propagation of the particles from the start point to
	    the end point.
	(2) Propagate particles 10 x dX[g/cm^2] by calling propagateDX( )
	    in the PropagationMatrix object 10 times.
	(3) Propagate the particles involved over a 2^k 10dX [g/cm^2] 
	    by calling propagateDXpowered( ) in the PropagationMatrix object 
	    k times.
	(4) During (3), store the propagation matrix when particles reaches
	    3 % of the total path in their journey. This stored matrix
	    is used for the final propagation steps in (5) next.
        (5) Propagate the particles over 0.03 x total path by mulitplication
	    of the progation matrix stored in (4). Finish the calculation
	    when the particles reaches the final point.

	(6) Repeat (1)-(5) until the particles finishes their travel.
	    By doing so, you can reasonably take into account
	    the r-dependence of the local medium density which decides
	    relative constributions of decay over interactions.
	</pre>
    */
    public void traceParticles(double trajectoryLength){

	System.err.println("Entering traceParticles."); 

	double lNow = s.getParticleLocation( ); // The current particle location
	                                        // along the trajectory
	double deltaL;
	double l;
	double Xsum = 1.0;
	int section;




	for(section=1;section<=5;section++){
	    // Building the elementary infinitesimal propagation matrix
	    propMtx.calculateTransferMatrix( );

	    // Trace the propagation step until the particles reaches the end point
	    // and obtain how many times (int n) step you need.
	    int n = 0;
	    l = lNow;
	    while(l<(0.2*trajectoryLength*(double )(section))){
		deltaL = propMtx.dX/s.getMediumDensity( ); 
		l += deltaL;n++;
		s.setParticleLocation(l);
	    }
	    s.setParticleLocation(lNow);


	    // Propagate the particles over 10 x dX [g/cm^2].
	    l = lNow;
	    for(int i=0;i<10;i++){
		deltaL = propMtx.dX/s.getMediumDensity( );// Propagation distance [cm]
		propMtx.propagateDX( );                   // Propagate particles over DX[g/cm]
		l += deltaL;
		s.setParticleLocation(l);

		System.err.println("Location " + l/100.0 + " [m]" 
				   + (trajectoryLength-l)/100.0 + " [m] to go.");
		System.err.println("nuMu to nuMu " + propMtx.getFnuMuToNuMu(300,300) + 
				   " nuMu to nuTau " + propMtx.getFnuMuToNuTau(300,300));
	    }


	    // Propagate the particles over 2^k times.
	    int k = (int )(Math.log((double )n/(double )10)/Math.log(2.0));
	    int times = 1;
	    double X = propMtx.dX*(double )10;
	    Xsum = 0.0;
	    int index = 1;
	    System.err.println("Doubling the finite matrix by " + k + " times");
	    for(int i=0;i<k;i++){
		propMtx.propagateDXpowered( );
		for(int j=0;j<times;j++){
		    deltaL = X/s.getMediumDensity( ); 
		    l += deltaL;
		    if(index == 1) Xsum += X;
		    s.setParticleLocation(l);
		}

		if(l>(0.2*(3.0e-2+(double )(section-1))*trajectoryLength) 
		   && index ==1){  
		    propMtx.copyTransferMatrix( );
		    index = 0;
		    System.err.println("Store the matrix at " + 
				   X +" g/cm^2 (" + l/100.0 + ") [m]");
		}


		times = times*2;
		System.err.println("Location " + l/100.0 + " [m]" 
				   + (trajectoryLength-l)/100.0 + " [m] to go.");
		System.err.println("nuMu to nuMu " + propMtx.getFnuMuToNuMu(300,300) + 
				   " nuMu to nuTau " + propMtx.getFnuMuToNuTau(300,300));
	    }

	    // Store the results
	    propMtx.storePropagateMatrix( );
	    // Initialize the propagation matrix
	    propMtx.init( );
	    // The current particle location
	    lNow = s.getParticleLocation( ); 
	}

	//Copy the stored matrix back to the main propagation matrix.
	propMtx.copyTransferMatrixFromStore( );


	// Final steps. A Step propagate the particles over Xsum [g/cm].
	while((l = s.getParticleLocation( ))<trajectoryLength){
	    deltaL = Xsum/s.getMediumDensity( ); // propagation distance [cm]
	    propMtx.propagateX( );
	    l += deltaL;
	    s.setParticleLocation(l);
	    System.err.println("Location " + l/100.0 + " [m]" 
			       + (trajectoryLength-l)/100.0 + " [m] to go.");
	    System.err.println("nuMu to nuMu " + propMtx.getFnuMuToNuMu(300,300) + 
			       " nuMu to nuTau " + propMtx.getFnuMuToNuTau(300,300));
	}
    }



    /**
       <pre>
       Same as traceParticles(double trajectoryLength)
       but the end point is decided to be where the particles
       emerged from the underground.
       </pre>
    */
    public void traceParticles( ){
	double trajectoryLength = s.getAxisLength( );
	traceParticles(trajectoryLength);
    }


    /** Save the calculated propagatin matrix */
    public void saveMatrix(DataOutputStream out) throws IOException {
	int dimension = Particle.getDimensionOfLogEnergyMatrix();
	int iLogE;
	for(iLogE=0;iLogE<dimension;iLogE++){
	    int jLogE;
	    for(jLogE=0;jLogE<=iLogE;jLogE++){
		/** For Glashow Resonance -begin **/
		out.writeDouble(propMtx.getFnuEToNuE(iLogE,jLogE));
		out.writeDouble(propMtx.getFnuEToNuMu(iLogE,jLogE));
		out.writeDouble(propMtx.getFnuEToNuTau(iLogE,jLogE));
		out.writeDouble(propMtx.getFnuEToE(iLogE,jLogE));
		out.writeDouble(propMtx.getFnuEToMu(iLogE,jLogE));
		out.writeDouble(propMtx.getFnuEToTau(iLogE,jLogE));
		out.writeDouble(propMtx.getFnuEToHadron(iLogE,jLogE));
		/** For Glashow Resonance -begin **/
		out.writeDouble(propMtx.getFnuMuToNuE(iLogE,jLogE));
		out.writeDouble(propMtx.getFnuMuToNuMu(iLogE,jLogE));
		out.writeDouble(propMtx.getFnuMuToNuTau(iLogE,jLogE));
		out.writeDouble(propMtx.getFnuMuToE(iLogE,jLogE));
		out.writeDouble(propMtx.getFnuMuToMu(iLogE,jLogE));
		out.writeDouble(propMtx.getFnuMuToTau(iLogE,jLogE));
		out.writeDouble(propMtx.getFnuMuToHadron(iLogE,jLogE));
		out.writeDouble(propMtx.getFnuTauToNuE(iLogE,jLogE));
		out.writeDouble(propMtx.getFnuTauToNuMu(iLogE,jLogE));
		out.writeDouble(propMtx.getFnuTauToNuTau(iLogE,jLogE));
		out.writeDouble(propMtx.getFnuTauToE(iLogE,jLogE));
		out.writeDouble(propMtx.getFnuTauToMu(iLogE,jLogE));
		out.writeDouble(propMtx.getFnuTauToTau(iLogE,jLogE));
		out.writeDouble(propMtx.getFnuTauToHadron(iLogE,jLogE));
		out.writeDouble(propMtx.getFmuToNuE(iLogE,jLogE));
		out.writeDouble(propMtx.getFmuToNuMu(iLogE,jLogE));
		out.writeDouble(propMtx.getFmuToNuTau(iLogE,jLogE));
		out.writeDouble(propMtx.getFmuToE(iLogE,jLogE));
		out.writeDouble(propMtx.getFmuToMu(iLogE,jLogE));
		out.writeDouble(propMtx.getFmuToTau(iLogE,jLogE));
		out.writeDouble(propMtx.getFmuToHadron(iLogE,jLogE));
		out.writeDouble(propMtx.getFtauToNuE(iLogE,jLogE));
		out.writeDouble(propMtx.getFtauToNuMu(iLogE,jLogE));
		out.writeDouble(propMtx.getFtauToNuTau(iLogE,jLogE));
		out.writeDouble(propMtx.getFtauToE(iLogE,jLogE));
		out.writeDouble(propMtx.getFtauToMu(iLogE,jLogE));
		out.writeDouble(propMtx.getFtauToTau(iLogE,jLogE));
		out.writeDouble(propMtx.getFtauToHadron(iLogE,jLogE));
	    }
	}

	out.close( );

    }



}
