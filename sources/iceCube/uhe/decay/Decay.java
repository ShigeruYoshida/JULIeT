package iceCube.uhe.decay;

import java.io.*;
import java.util.*;

import iceCube.uhe.particles.*;

/** Particle decay such as tau -> mu -> e is handled in this class.
    The particle livetimes themselves are provided by the Particl class
    (Particle.java), but the decay matrix to give the energy distributions
    of the secondary particles, and the terms due to the Laorentz duration
    effects are supplied by the methods here. 
*/
public class Decay{

    /** The squared mass raio (m_mu/m_pi)^2. This term
	is frequently used in the differential energy
	distribution in the 2-body decay determined
	by the simple kinetics.
    */
    public static final double rMuPI = 5.731085e-1;


    /** The squared mass raio (m_pi/m_tau)^2. This term
	is frequently used in the differential energy
	distribution in the 2-body decay determined
	by the simple kinetics. In this paticular case,
	hadronic decay of taus producing pions.
    */
    public static final double rPiTau = 6.1197e-3;


    /** The squared mass raio (m_rho/m_tau)^2. This term
	is frequently used in the differential energy
	distribution in the 2-body decay determined
	by the simple kinetics. In this paticular case,
	hadronic decay of taus producing rhos.
    */
    public static final double rRhoTau = 1.8535e-1;


    /** The squared mass raio (m_A1/m_tau)^2. This term
	is frequently used in the differential energy
	distribution in the 2-body decay determined
	by the simple kinetics. 
    */
    public static final double rA1Tau = 4.9877e-1;


    /** The squared mass raio (m_X/m_tau)^2. This term
	is frequently used in the differential energy
	distribution in the 2-body decay determined
	by the simple kinetics. m_X is the "averaged mass"
	of produced hadrons in the decay mode other than
	pi,rho,and a1.
    */
    public static final double rXTau = 0.7;



    /** The Branching ratio of Leptonic tau decay. */
    public static final double BRatioTau2Leptons = 0.36;
    /** The Branching ratio of Hadronic tau decay puroducing pi's. */
    public static final double BRatioTau2Pi = 0.12;
    /** The Branching ratio of Hadronic tau decay puroducing rho's. */
    public static final double BRatioTau2Rho = 0.26;
    /** The Branching ratio of Hadronic tau decay puroducing a1's. */
    public static final double BRatioTau2A1 = 0.13;
    /** The Branching ratio of Hadronic tau decay puroducing others. */
    public static final double BRatioTau2X = 0.13;

    /** The Particle class involved with the decay */
    Particle p;     
    /** The decay particle energy [GeV]. */
    double energy;
    /** lifetime of the decay particle [sec].*/
    double tau;



    /** Constructor: Register the Particle involved.
	<pre>
	Particle P: The Particle class involved with the decay.
                    This should be either muon or tauons.
       </pre>
    */
    public Decay(Particle p){
        if(isValidParticle(p)){
            this.p = p;
	    this.energy = p.getEnergy( );
	    tau = p.getLifeTime( );
	}else{
            System.err.println("This particle " + 
                   p.particleName(p.getFlavor(), p.getDoublet()) + 
                   " should not be involved in the decay");
            System.exit(0);
        }
    }

    /** Checking the particle kind involved with
        a given decay channel. Only muons, taus, pi's
        but not electrons and neutrinos are allowed to be involved
	for the obvious reasons.
    */
    public boolean isValidParticle(Particle p){
        if(p.getDoublet()==1 && p.getFlavor()>0)
        return true;
        else return false;
    }


    /** Calculate the Weak Decay probability per inelasticity "y"
	that the decay product see its energy of y*Energy of 
	the decaying particle.
        The product originated in the charged current from
        the decaying particle is taken care of.

	<pre>
        y       Inelasticity 
        parity  Spin polarization vector. 0 for no-polarization 
	</pre>

    */
    public static double getWeakDecayProbToW(double y, double parity){
	double nonParityTerm = 5.0/3.0 - 3.0*y*y + 4.0/3.0*y*y*y;
	if(parity == 0.0){
	    return nonParityTerm;
	}else{
	    double parityTerm = 1.0/3.0 - 3.0*y*y + 8.0/3.0*y*y*y;
	    return nonParityTerm + parity*parityTerm;
	}
    }

    /** Calculate the Weak Decay probability per inelasticity "y"
	that the decay product see its energy of y*Energy of 
	the decaying particle.
        The produced neutrino originated in the charged current on 
	the another side  is taken care of.

	<pre>
        y       Inelasticity 
        parity  Spin polarization vector. 0 for no-polarization 
	</pre>
    */
    public static double getWeakDecayProbFromW(double y, double parity){
	double nonParityTerm = 2.0 - 6.0*y*y + 4.0*y*y*y;
	if(parity == 0.0){
	    return nonParityTerm;
	}else{
	    double parityTerm = -2.0 +12.0*y - 18.0*y*y + 8.0*y*y*y;
	    return nonParityTerm + parity*parityTerm;
	}
    }


    /** Integral of the weak decay probability from lowerY to upperY.
        The integration is analyticaly made in this method */
    public static double integralWeakDecayProbToW(double lowerY,
					  double upperY, double parity){
	double yLow;
	if(lowerY<getYmin()) yLow = getYmin();
	else yLow = lowerY;
	double yUp;
	if(getYmax()<upperY) yUp = getYmax();
	else yUp = upperY;

	double nonParityTermUp = 5.0/3.0*yUp - yUp*yUp*yUp + 1.0/3.0*yUp*yUp*yUp*yUp;
	double nonParityTermLow = 
	    5.0/3.0*yLow - yLow*yLow*yLow + 1.0/3.0*yLow*yLow*yLow*yLow;
	if(parity == 0.0){
	    return (nonParityTermUp-nonParityTermLow);
	}else{
	    double parityTermUp = 1.0/3.0*yUp - yUp*yUp*yUp + 2.0/3.0*yUp*yUp*yUp*yUp;
	    double parityTermLow = 
		1.0/3.0*yLow - yLow*yLow*yLow + 2.0/3.0*yLow*yLow*yLow*yLow;
	    return (nonParityTermUp-nonParityTermLow+ parity*(parityTermUp-parityTermLow));
	}
    }

    /** Same as integralWeakDecayProbToW( ), but handling the particles
	of the charged current on the another side.*/
    public static double integralWeakDecayProbFromW(double lowerY,
						  double upperY, double parity){
	double yLow;
	if(lowerY<getYmin()) yLow = getYmin();
	else yLow = lowerY;
	double yUp;
	if(getYmax()<upperY) yUp = getYmax();
	else yUp = upperY;

	double nonParityTermUp = 2.0*yUp - 2.0*yUp*yUp*yUp + yUp*yUp*yUp*yUp;
	double nonParityTermLow = 2.0*yLow - 2.0*yLow*yLow*yLow + yLow*yLow*yLow*yLow;
	if(parity == 0.0){
	    return (nonParityTermUp-nonParityTermLow);
	}else{
	    double parityTermUp = 
		-2.0*yUp + 6.0*yUp*yUp - 6.0*yUp*yUp*yUp + 2.0*yUp*yUp*yUp*yUp;
	    double parityTermLow = 
		-2.0*yLow + 6.0*yLow*yLow - 6.0*yLow*yLow*yLow + 2.0*yLow*yLow*yLow*yLow;
	    return (nonParityTermUp-nonParityTermLow+ parity*(parityTermUp-parityTermLow));
	}
    }

    public static double getYmin( ){
	return 0.0;
    }

    public static double getYmax( ){
	return 1.0;
    }


    /** Calculate the Tau hadron Decay probability per inelasticity "y"
	that the decay product (in this case nu-tau) see 
	its energy of y*Energy of tau.
    */
    public static double getTauHadronDecayProbToW(double y, double massRatio){
	return 1.0/(1.0-massRatio);
    }
    /** Calculate the Tau hadron Decay probability per inelasticity "y"
	that the decay product (in this case hadron) see 
	its energy of y*Energy of tau.
    */
    public static double getTauHadronDecayProbFromW(double y, double massRatio){
	return 1.0/(1.0-massRatio);
    }

    public static double getYmax(double massRatio){
	return (1.0-massRatio);
    }


    /** Integral of the tau hadron decay probability from lowerY to upperY.
        The integration is analyticaly made in this method */
    public static double integralTauHadronDecayProbToW(double lowerY,
					  double upperY, double massRatio){
	double yLow;
	if(lowerY<getYmin()) yLow = getYmin();
	else yLow = lowerY;
	double yUp;
	if(getYmax(massRatio)<upperY) yUp = getYmax(massRatio);
	else yUp = upperY;

	double y = yUp-yLow;
	return y/(1.0-massRatio);
    }


    public static double integralTauHadronDecayProbFromW(double lowerY,
					  double upperY, double massRatio){

	return integralTauHadronDecayProbToW(1.0-upperY,
					     1.0-lowerY, massRatio);
    }
}

