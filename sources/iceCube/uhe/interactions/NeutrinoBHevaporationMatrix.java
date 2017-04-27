package iceCube.uhe.interactions;

import java.io.*;

import iceCube.uhe.particles.*;
import iceCube.uhe.points.*;


/** 
*/

public class NeutrinoBHevaporationMatrix extends InteractionsMatrix {


    /** Constructor: Generate the matrix array */
    public NeutrinoBHevaporationMatrix(NeutrinoBHevaporation interactions){

	super(interactions);

	// Interaction must be "Neutrino BH evaporation"
	if(!interactions.interactionName().startsWith("Neutrino Black Hole")){
            System.err.println(interactions.interactionName() + "is not BH-evaporation interaction");
            System.exit(0);
	}
    }

 
    /** Calculate the total cross section matrix */
    public void setSigmaMatrix(int iLogE){
	    
	interactions.setIncidentParticleEnergy(iLogE);/** Setup the primary E*/
	NeutrinoBHevaporation nuBH = (NeutrinoBHevaporation)interactions;
	nuBH.switchToShower();
	if(isValidIndex(iLogE)){
	    sigmaMtx[iLogE] = interactions.getSigma( );
	    inelasticityMtx[iLogE] = 
	    interactions.getYDSigmaDy(interactions.getYmin()+interactions.roundOffError,
				      interactions.getYmax()-interactions.roundOffError);
	}
    }



    /** Calculate the transfer matrix */
    public void setTransferMatrix(int iLogE, int jLogE){

	interactions.setIncidentParticleEnergy(iLogE);/** Setup the primary E*/

	if(isValidIndex(iLogE) && isValidIndex(jLogE)){
	    if(jLogE>iLogE){
		transferMtx[iLogE][jLogE]=0.0;
		transferAMtx[iLogE][jLogE]=0.0;
	    }
	    else{
		double logY = delta*(double )(jLogE-iLogE);
		double logYLow = logY-0.5*delta;
		double logYUp = logY+0.5*delta;
		//System.err.println("logYLow " + logYLow + " logYUp " + logYUp);

		double yLow = Math.pow(10.0,logYLow);
		double yLowRange = yLow;
		NeutrinoBHevaporation nuBH = (NeutrinoBHevaporation)interactions;
		nuBH.switchToShower();
		if(yLow<=interactions.getYmin()+2.0*interactions.roundOffError)  
		    yLowRange = interactions.getYmin( ) + 2.0*interactions.roundOffError;

		double yUp = Math.pow(10.0,logYUp);
		double yUpRange = yUp;
		if(yUp>=interactions.getYmax()-2.0*interactions.roundOffError)   
		    yUpRange = interactions.getYmax( ) - 2.0*interactions.roundOffError;


		if(yLowRange<yUpRange){
		    transferMtx[iLogE][jLogE]=
			interactions.integralDSigmaDy(yLowRange,yUpRange);
		}else{
		    transferMtx[iLogE][jLogE]=0.0;
		}

		yLowRange = yLow; yUpRange = yUp;
		nuBH.switchToChargedLepton();
		if(yLow<=(1.0-interactions.getYmax()+2.0*interactions.roundOffError)) 
		    yLowRange = 1.0-interactions.getYmax( )+2.0*interactions.roundOffError;
		if(yUp>=(1.0-interactions.getYmin()-2.0*interactions.roundOffError))  
		    yUpRange = 1.0-interactions.getYmin( )-2.0*interactions.roundOffError;

		if(yLowRange<yUpRange){
		    transferAMtx[iLogE][jLogE]=
			interactions.integralDSigmaDz(yLowRange,yUpRange);
		}else{
		    transferAMtx[iLogE][jLogE]=0.0;
		}

	    }
	}
    }


}


