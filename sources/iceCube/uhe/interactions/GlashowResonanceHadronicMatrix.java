package iceCube.uhe.interactions;

import java.io.*;

import iceCube.uhe.particles.*;
import iceCube.uhe.points.*;


/** 
<pre>
    Matrix of the Energy Transfer by particle Interactions.
    The matrix elements are calculated by the methods supplied
    by the Interaction class.

        /------------------------------------------------------------\
logEmin | 0   0  ...................................   EmindS/dEmin  |
logE1   | 0   0  ............................E1dS/dE1  EmindS/dEmin  |
  .     |                                                            |x binWidth x ln(10)
  .     |                                                            |
logEmax |EmaxdS/dEmax................................  EmindS/dEmin  |
        \------------------------------------------------------------/

The bin width and logEmin are defined in the Particle class.

Actually each matrix element dSigma/dLogE is calculated by the integral,
\int dSigma/dY dY from logY - 0.5xbinWidth to logY + 0.5*binWidth where
logY = logEproduced - logEincoming. This is more accurate way than simply
calculating dSigma/dLogE.

The transfer matrix of the recolied lepton energy (1-y)*Eincoming 
is acquired by the method getLeptonTransferMatrix( )
while the energy of the transfered matrix for the counter current 
y*Eincoming is obtained by getTransferMatrix( ).
</pre>
*/

//public class GlashowResonanceHadronicMatrix extends InteractionsMatrix implements Serializable {
public class GlashowResonanceHadronicMatrix extends InteractionsMatrix {

    /** Constructor: Generate the matrix array */
    public GlashowResonanceHadronicMatrix(Interactions interactions){

	super(interactions);

	// Interaction must be "GlashowResonanceHadronic"
	if(!interactions.interactionName().startsWith("Glashow Resonance Hadronic")){
            System.err.println(interactions.interactionName() + "is not GlashowResonanceHadronic interaction");
            System.exit(0);
	}
	System.err.println("GlashowResonanceHadronicMatrix has been constructed. ");
    }

    /** Calculate the total cross section matrix */
    public void setSigmaMatrix(int iLogE){

        interactions.setIncidentParticleEnergy(iLogE);/** Setup the primary E*/
        if(isValidIndex(iLogE)){
            sigmaMtx[iLogE] = interactions.getSigma( );
            inelasticityMtx[iLogE] = sigmaMtx[iLogE];
        }
    }

    /** Calculate the transfer matrix

        This matrix should be a diagonal matrix with just "Sigma",
        because all the incident energy is deposited as a hadronic cascade.
        So, you don't need integral differencial cross-section,
        and do need itneractions.getSigma( ) method.
    */
    public void setTransferMatrix(int iLogE, int jLogE){

	interactions.setIncidentParticleEnergy(iLogE);/** Setup the primary E*/

	if(isValidIndex(iLogE) && isValidIndex(jLogE)){
            if(jLogE==iLogE){
                transferMtx[iLogE][jLogE]=interactions.getSigma( );
                transferAMtx[iLogE][jLogE]=0.0;
            }else{
                transferMtx[iLogE][jLogE]=0.0;
                transferAMtx[iLogE][jLogE]=0.0;
            }
	}
    }
    
}

