package iceCube.uhe.decay;

import java.io.*;

import iceCube.uhe.interactions.*;
import iceCube.uhe.particles.*;


/** 
<pre>
Matrix of the Energy Transfer dN/dLogY  by the tau decays
    The matrix elements are calculated by the methods supplied
    by the Decay class.

The bin width is defined in the Particle class.

Actually each matrix element dN/dLogYdecay is calculated by the integral,
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
public class TauDecayYMatrix {

    /** Array for the energy tansfer probability from Tau to nu-Tau.*/
    double[] nuTauMtx;
    /** Array for the energy tansfer probability from Tau to nu-Mu or nuE.*/
    double[] nuMtx;
    /** Array for the energy tansfer probability from Tau to Electon or Muon.*/
    double[] leptonsMtx;
    /** Array for the energy tansfer probability from Tau to Hadrons.
     Adds the contributions from all hadrocnic decay modes. */
    double[] hadronsMtx;
    /** Array for the lifetime of tauon considering the Lorentz duration. */
    double[] lifetimeMtx;
    /** Array Dimension */
    int dimension = Particle.getDimensionOfLogEnergyMatrix();
    /** The bin width of the matrix element.*/
    double delta = Particle.getDeltaLogEnergy();
    /** The Particle class. Should be tau. This is checked by the constructor.*/
    Particle p;

    /** Constructor: Generate the matrix array */
    public TauDecayYMatrix(Particle p){
        if(p.getDoublet()==1 && p.getFlavor()==2){ // Requires taus
	    this.p = p;
	    nuTauMtx = new double[dimension];
	    nuMtx = new double[dimension];
	    leptonsMtx = new double[dimension];
	    hadronsMtx = new double[dimension];
	    lifetimeMtx = new double[dimension];
	}else{
            System.err.println("This particle " + 
			       p.particleName(p.getFlavor(), p.getDoublet()) + 
			       " is not TAUs!!");
            System.exit(0);
	}
    }



    /** Calculate the decay matrix  from tau to neutrinos, leptons, and hadrons*/
    public void setTauDecayMatrix(int iLogY){

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
	    nuTauMtx[iLogY] = Decay.BRatioTau2Leptons*decayElement;
		              // to nu-tau
	    leptonsMtx[iLogY] = 0.5*Decay.BRatioTau2Leptons*decayElement;
		              // to electrons and muons, respectively
	    decayElement = 0.0;
	    if(yLowRange<yUpRange){
		decayElement=Decay.integralWeakDecayProbFromW(yLowRange,yUpRange, -1.0);
	    }
	    nuMtx[iLogY] = 0.5*Decay.BRatioTau2Leptons*decayElement;
		              // to nu-e and nu-mu, respectively



	    // Hadron decay
	    yUpRange = yUp;
	    yLowRange = yLow;
	    if(yUp>=Decay.getYmax(Decay.rPiTau)) 
		yUpRange = Decay.getYmax(Decay.rPiTau);
	    if(yLowRange<yUpRange){
		decayElement = 
		    Decay.integralTauHadronDecayProbToW(yLowRange,yUpRange,Decay.rPiTau);
		nuTauMtx[iLogY] += Decay.BRatioTau2Pi*decayElement;
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
		hadronsMtx[iLogY] = Decay.BRatioTau2Pi*decayElement;
	    }else{
		hadronsMtx[iLogY] = 0.0;
	    }
		             // to pions via the tau -> pi

	    yUpRange = yUp;
	    yLowRange = yLow;
	    if(yUp>=Decay.getYmax(Decay.rRhoTau))   
		yUpRange = Decay.getYmax(Decay.rRhoTau);
	    if(yLowRange<yUpRange){
		decayElement = 
		    Decay.integralTauHadronDecayProbToW(yLowRange,yUpRange,Decay.rRhoTau);
		nuTauMtx[iLogY] += Decay.BRatioTau2Rho*decayElement;
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
		hadronsMtx[iLogY] += Decay.BRatioTau2Rho*decayElement;
	    }
		             // to hadrons via the tau -> rho

	    yUpRange = yUp;
	    yLowRange = yLow;
	    if(yUp>=Decay.getYmax(Decay.rA1Tau))   
		yUpRange = Decay.getYmax(Decay.rA1Tau);
	    if(yLowRange<yUpRange){
		decayElement = 
		    Decay.integralTauHadronDecayProbToW(yLowRange,yUpRange,Decay.rA1Tau);
		nuTauMtx[iLogY] += Decay.BRatioTau2A1*decayElement;
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
		hadronsMtx[iLogY] += Decay.BRatioTau2A1*decayElement;
	    }
		             // to hadrons via the tau -> a1

	    yUpRange = yUp;
	    yLowRange = yLow;
	    if(yUp>=Decay.getYmax(Decay.rXTau))   
		yUpRange = Decay.getYmax(Decay.rXTau);
	    if(yLowRange<yUpRange){
		decayElement = 
		    Decay.integralTauHadronDecayProbToW(yLowRange,yUpRange,Decay.rXTau);
		nuTauMtx[iLogY] += Decay.BRatioTau2X*decayElement;
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
		hadronsMtx[iLogY] += Decay.BRatioTau2X*decayElement;
	    }
		             // to hadrons via the tau -> X

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
    public double getTauToNuTauDecayMatrix(int iLogY){
	if(isValidIndex(iLogY)){
	    return nuTauMtx[iLogY];
	}
	else{
	    return 0.0;
	}
    }

    /** get the element of the decay matrix  of tau->nuE or nuMu */
    public double getTauToNuDecayMatrix(int iLogY){
	if(isValidIndex(iLogY)){
	    return nuMtx[iLogY];
	}
	else{
	    return 0.0;
	}
    }

    /** get the element of the decay matrix  of tau-> l^+- */
    public double getTauToChargedLeptonDecayMatrix(int iLogY){
	if(isValidIndex(iLogY)){
	    return leptonsMtx[iLogY];
	}
	else{
	    return 0.0;
	}
    }

    /** get the element of the decay matrix  of tau-> hadrons */
    public double getTauToHadronDecayMatrix(int iLogY){
	if(isValidIndex(iLogY)){
	    return hadronsMtx[iLogY];
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
