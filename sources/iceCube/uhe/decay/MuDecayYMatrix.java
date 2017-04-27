package iceCube.uhe.decay;

import java.io.*;

import iceCube.uhe.interactions.*;
import iceCube.uhe.particles.*;


/** 
<pre>
Matrix of the Energy Transfer dN/dLogY by the mu decays
    The matrix elements are calculated by the methods supplied
    by the Decay class.

The bin width is defined in the Particle class. 

Actually each matrix element dN/dLogYdecay is calculated by the integral,
\int dN/dY dY from logY - 0.5xbinWidth to logY + 0.5*binWidth where
logY = logEdecay - logEmu. This is more accurate way than simply
calculating dN/dLogEdecay.

The transfer matrix to the nu-mu
is acquired by the method getMuToNuMuDecayMatrix( ),
that to the electron neutrinos by getMuToNuEMatrix( )
that to the other charged leptons such as e by getMuToEDecayMatrix( )
</pre>
*/

public class MuDecayYMatrix {

    /** Array for the energy tansfer probability from Mu to nu-Mu or electrons. */
    double[] nuMuMtx;
    /** Array for the energy tansfer probability from Mu to nu-E. */
    double[] nuMtx;
    /** Array for the lifetime of muon considering the Lorentz duration. */
    double[] lifetimeMtx;
    /** Array Dimension */
    int dimension = Particle.getDimensionOfLogEnergyMatrix();
    /** The bin width of the matrix element.*/
    double delta = Particle.getDeltaLogEnergy();
    /** The Particle class. Should be muons. This is checked by the constructor.*/
    Particle p;

    /** Constructor: Generate the matrix array */
    public MuDecayYMatrix(Particle p){
	if(p.getDoublet()==1 && p.getFlavor()==1){ // Requires mu's
	    this.p = p;
	    nuMuMtx = new double[dimension];
	    nuMtx = new double[dimension];
	    lifetimeMtx = new double[dimension];
	}else{
	    System.err.println("This particle " +
			       p.particleName(p.getFlavor(), p.getDoublet()) +
			       " is not MUONs!!");
	    System.exit(0);
	}
    }	


    /** Calculate the decay matrix from mu to neutrinos and and electrons*/
    public void setMuDecayMatrix(int iLogY){


	if(isValidIndex(iLogY)){
	    double logY = -delta*(double )(iLogY);
	    double logYLow = logY-0.5*delta;
	    double logYUp = logY+0.5*delta;


	    double yLow = Math.pow(10.0,logYLow);
	    double yLowRange = yLow;
	    if(yLow<=Decay.getYmin())  yLowRange = Decay.getYmin( );

	    double yUp = Math.pow(10.0,logYUp);
	    double yUpRange = yUp;
	    if(yUp>=Decay.getYmax())   yUpRange = Decay.getYmax( );



	    // Weak lepton decay 
	    double decayElement = 0.0;
	    if(yLowRange<yUpRange){
		decayElement = 
		    Decay.integralWeakDecayProbToW(yLowRange,yUpRange, -1.0);
	    }
	    nuMuMtx[iLogY] = decayElement;
		              // to nu-mu
	    decayElement = 0.0;
	    if(yLowRange<yUpRange){
		decayElement=
		    Decay.integralWeakDecayProbFromW(yLowRange,yUpRange, -1.0);
	    }
	    nuMtx[iLogY] = decayElement;
	                     // to nu-e

	}
    }

    /** Calculate the life time matrix considering the Lorentz duration */
    public void setLifeTimeMatrix(int iLogE){
        double logEnergy = p.getLogEnergyMinimum()+
            p.getDeltaLogEnergy()*(double )iLogE;
        double energy = Math.pow(10.0,logEnergy);
        p.putEnergy(energy);
        p.putLogEnergy(logEnergy);
	if(isValidIndex(iLogE)){
	    lifetimeMtx[iLogE] = p.getLifeTime()*p.getEnergy( )/p.getMass( );
	}
    }

    /** get the element of the decay matrix  of m->nuMu */
    public double getMuToNuMuDecayMatrix(int iLogY){
	if(isValidIndex(iLogY)){
	    return nuMuMtx[iLogY];
	}
	else{
	    return 0.0;
	}
    }

    /** get the element of the decay matrix  of mu->nuE */
    public double getMuToNuEDecayMatrix(int iLogY){
	if(isValidIndex(iLogY)){
	    return nuMtx[iLogY];
	}
	else{
	    return 0.0;
	}
    }

    /** get the element of the decay matrix  of mu-> e */
    public double getMuToEDecayMatrix(int iLogY){
	if(isValidIndex(iLogY)){
	    return nuMuMtx[iLogY];
	}
	else{
	    return 0.0;
	}
    }

    /** get the element of the LifeTime matrix */
    public double getLifeTimeMatrix(int iLogE){
	if(isValidIndex(iLogE)){
	    return lifetimeMtx[iLogE];
	}
	else{
	    return 0.0;
	}
    }


    /** Checking the energy index */
    public boolean isValidIndex(int iLogY){
	if(0<=iLogY && iLogY<dimension) return true;
	else return false;
    }

}
