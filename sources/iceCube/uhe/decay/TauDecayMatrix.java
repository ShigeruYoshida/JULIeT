package iceCube.uhe.decay;

import java.io.*;

import iceCube.uhe.interactions.*;
import iceCube.uhe.particles.*;


/** 
<pre>
Matrix of the Energy Transfer by the tau decays
    The matrix elements are calculated by the methods supplied
    by the Decay class.

        /------------------------------------------------------------\
logEmin | 0   0  ...................................   EmindN/dEmin  |
logE1   | 0   0  ............................E1dN/dE1  EmindN/dEmin  |
  .     |                                                            |x binWidth x ln(10)
  .     |                                                            |
logEmax |EmaxdN/dEmax................................  EmindN/dEmin  |
        \------------------------------------------------------------/

The bin width and logEmin are defined in the Particle class.

Actually each matrix element dN/dLogEdecay is calculated by the integral,
\int dN/dY dY from logY - 0.5xbinWidth to logY + 0.5*binWidth where
logY = logEdecay - logEtau. This is more accurate way than simply
calculating dN/dLogEdecay.

The transfer matrix to the nu-tau
is acquired by the method getTauToNuTauDecayMatrix( ),
that to the other neutrinos such as nu-e by getTauToNuDecayMatrix( )
that to the other charged leptons such as e by getTauToChargedLeptonDecayMatrix( )
while the decay into hadrons
is obtained by getTauToHadronsDecayMatrix( ).
</pre>
*/
public class TauDecayMatrix {

    /** Array for the energy tansfer probability from Tau to nu-Tau.*/
    double[][] nuTauMtx;
    /** Array for the energy tansfer probability from Tau to nu-Mu or nuE.*/
    double[][] nuMtx;
    /** Array for the energy tansfer probability from Tau to Electon or Muon.*/
    double[][] leptonsMtx;
    /** Array for the energy tansfer probability from Tau to Hadrons.
     Adds the contributions from all hadrocnic decay modes. */
    double[][] hadronsMtx;
    /** Array for the lifetime of tauon considering the Lorentz duration. */
    double[] lifetimeMtx;
    /** Array Dimension */
    int dimension = Particle.getDimensionOfLogEnergyMatrix();
    /** The bin width of the matrix element.*/
    double delta = Particle.getDeltaLogEnergy();
    /** The Particle class. Should be tau. This is checked by the constructor.*/
    Particle p;

    /** Constructor: Generate the matrix array */
    public TauDecayMatrix(Particle p){
        if(p.getDoublet()==1 && p.getFlavor()==2){ // Requires taus
	    this.p = p;
	    nuTauMtx = new double[dimension][dimension];
	    nuMtx = new double[dimension][dimension];
	    leptonsMtx = new double[dimension][dimension];
	    hadronsMtx = new double[dimension][dimension];
	    lifetimeMtx = new double[dimension];
	}else{
            System.err.println("This particle " + 
                   p.particleName(p.getFlavor(), p.getDoublet()) + 
                   " is not TAUs!!");
            System.exit(0);
	}
    }



    /** Calculate the decay matrix  from tau to neutrinos, leptons, and hadrons*/
    public void setTauDecayMatrix(int iLogE, int jLogE){

        double logEnergy = p.getLogEnergyMinimum()+
            p.getDeltaLogEnergy()*(double )iLogE;
        double energy = Math.pow(10.0,logEnergy);
        p.putEnergy(energy);
        p.putLogEnergy(logEnergy);

	if(isValidIndex(iLogE) && isValidIndex(jLogE)){
	    if(jLogE>iLogE){
		nuTauMtx[iLogE][jLogE]=0.0;
		nuMtx[iLogE][jLogE]=0.0;
		leptonsMtx[iLogE][jLogE]=0.0;
		hadronsMtx[iLogE][jLogE]=0.0;
	    }else{
		double logY = delta*(double )(jLogE-iLogE);
		double logYLow = logY-0.5*delta;
		double logYUp = logY+0.5*delta;
		//System.err.println("logYLow " + logYLow + " logYUp " + logYUp);


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
		nuTauMtx[iLogE][jLogE] = Decay.BRatioTau2Leptons*decayElement;
		              // to nu-tau
		leptonsMtx[iLogE][jLogE] = 0.5*Decay.BRatioTau2Leptons*decayElement;
		              // to electrons and muons, respectively
		decayElement = 0.0;
		if(yLowRange<yUpRange){
		    decayElement=Decay.integralWeakDecayProbFromW(yLowRange,yUpRange, -1.0);
		}
		nuMtx[iLogE][jLogE] = 0.5*Decay.BRatioTau2Leptons*decayElement;
		              // to nu-e and nu-mu, respectively



		// Hadron decay
		yUpRange = yUp;
		yLowRange = yLow;
		if(yUp>=Decay.getYmax(Decay.rPiTau)) 
		    yUpRange = Decay.getYmax(Decay.rPiTau);
                if(yLowRange<yUpRange){
		    decayElement = 
		    Decay.integralTauHadronDecayProbToW(yLowRange,yUpRange,Decay.rPiTau);
		    nuTauMtx[iLogE][jLogE] += Decay.BRatioTau2Pi*decayElement;
		             // to nuTau via the tau -> pi
		}
		yUpRange = yUp;
		yLowRange = yLow;
                if(yLow<=(1.0-Decay.getYmax(Decay.rPiTau))) 
                    yLowRange = 1.0-Decay.getYmax(Decay.rPiTau);
                if(yUp>=(1.0-Decay.getYmin()))  
                    yUpRange = 1.0-Decay.getYmin( );
                if(yLowRange<yUpRange){
		    decayElement = 
		    Decay.integralTauHadronDecayProbFromW(yLowRange,yUpRange,Decay.rPiTau);
		    hadronsMtx[iLogE][jLogE] = Decay.BRatioTau2Pi*decayElement;
		}else{
		    hadronsMtx[iLogE][jLogE] = 0.0;
		}
		             // to pions via the tau -> pi

		yUpRange = yUp;
		yLowRange = yLow;
		if(yUp>=Decay.getYmax(Decay.rRhoTau))   
		    yUpRange = Decay.getYmax(Decay.rRhoTau);
                if(yLowRange<yUpRange){
		    decayElement = 
		    Decay.integralTauHadronDecayProbToW(yLowRange,yUpRange,Decay.rRhoTau);
		    nuTauMtx[iLogE][jLogE] += Decay.BRatioTau2Rho*decayElement;
		             // to nuTau via the tau -> rho
		}
		yUpRange = yUp;
		yLowRange = yLow;
                if(yLow<=(1.0-Decay.getYmax(Decay.rRhoTau))) 
                    yLowRange = 1.0-Decay.getYmax(Decay.rRhoTau);
                if(yUp>=(1.0-Decay.getYmin()))  
                    yUpRange = 1.0-Decay.getYmin( );
                if(yLowRange<yUpRange){
		    decayElement = 
		    Decay.integralTauHadronDecayProbFromW(yLowRange,yUpRange,Decay.rRhoTau);
		    hadronsMtx[iLogE][jLogE] += Decay.BRatioTau2Rho*decayElement;
		}
		             // to hadrons via the tau -> rho

		yUpRange = yUp;
		yLowRange = yLow;
		if(yUp>=Decay.getYmax(Decay.rA1Tau))   
		    yUpRange = Decay.getYmax(Decay.rA1Tau);
                if(yLowRange<yUpRange){
		    decayElement = 
		    Decay.integralTauHadronDecayProbToW(yLowRange,yUpRange,Decay.rA1Tau);
		    nuTauMtx[iLogE][jLogE] += Decay.BRatioTau2A1*decayElement;
		             // to nuTau via the tau -> a1
		}
		yUpRange = yUp;
		yLowRange = yLow;
                if(yLow<=(1.0-Decay.getYmax(Decay.rA1Tau))) 
                    yLowRange = 1.0-Decay.getYmax(Decay.rA1Tau);
                if(yUp>=(1.0-Decay.getYmin()))  
                    yUpRange = 1.0-Decay.getYmin( );
                if(yLowRange<yUpRange){
		    decayElement = 
	            Decay.integralTauHadronDecayProbFromW(yLowRange,yUpRange,Decay.rA1Tau);
		    hadronsMtx[iLogE][jLogE] += Decay.BRatioTau2A1*decayElement;
		}
		             // to hadrons via the tau -> a1

		yUpRange = yUp;
		yLowRange = yLow;
		if(yUp>=Decay.getYmax(Decay.rXTau))   
		    yUpRange = Decay.getYmax(Decay.rXTau);
                if(yLowRange<yUpRange){
		    decayElement = 
		    Decay.integralTauHadronDecayProbToW(yLowRange,yUpRange,Decay.rXTau);
		    nuTauMtx[iLogE][jLogE] += Decay.BRatioTau2X*decayElement;
		             // to nuTau via the tau -> X
		}
		yUpRange = yUp;
		yLowRange = yLow;
                if(yLow<=(1.0-Decay.getYmax(Decay.rXTau))) 
                    yLowRange = 1.0-Decay.getYmax(Decay.rXTau);
                if(yUp>=(1.0-Decay.getYmin()))  
                    yUpRange = 1.0-Decay.getYmin( );
                if(yLowRange<yUpRange){
		    decayElement = 
		    Decay.integralTauHadronDecayProbFromW(yLowRange,yUpRange,Decay.rXTau);
		    hadronsMtx[iLogE][jLogE] += Decay.BRatioTau2X*decayElement;
		}
		             // to hadrons via the tau -> X

	    }
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

    /** get the element of the decay matrix  of tau->nuTau */
    public double getTauToNuTauDecayMatrix(int iLogE, int jLogE){
	if(isValidIndex(iLogE) && isValidIndex(jLogE)){
	    return nuTauMtx[iLogE][jLogE];
	}
	else{
	    return 0.0;
	}
    }

    /** get the element of the decay matrix  of tau->nuE or nuMu */
    public double getTauToNuDecayMatrix(int iLogE, int jLogE){
	if(isValidIndex(iLogE) && isValidIndex(jLogE)){
	    return nuMtx[iLogE][jLogE];
	}
	else{
	    return 0.0;
	}
    }

    /** get the element of the decay matrix  of tau-> l^+- */
    public double getTauToChargedLeptonDecayMatrix(int iLogE, int jLogE){
	if(isValidIndex(iLogE) && isValidIndex(jLogE)){
	    return leptonsMtx[iLogE][jLogE];
	}
	else{
	    return 0.0;
	}
    }

    /** get the element of the decay matrix  of tau-> hadrons */
    public double getTauToHadronDecayMatrix(int iLogE, int jLogE){
	if(isValidIndex(iLogE) && isValidIndex(jLogE)){
	    return hadronsMtx[iLogE][jLogE];
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
    public boolean isValidIndex(int iLogE){
	if(0<=iLogE && iLogE<dimension) return true;
	else return false;
    }

}
