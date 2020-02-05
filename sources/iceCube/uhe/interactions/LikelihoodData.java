package iceCube.uhe.interactions;

/**
   a class to store the likelihood value and the propagating particle energy 
   after this interaction occures. Essentially this is a data class.
*/
public class LikelihoodData {
    double interactionLikelihood = 0.0;
    double trackllh = 0.0;
    double particleEnergyAfterThisInteraction = 0.0;

    /** constructor */
    public LikelihoodData(){};

    /** Return the likelihood value */
    public double getLikelihoodValue(){
	return interactionLikelihood;
    }
    
    /** Return the propagating particle energy [GeV]  */
    public double getEnergyAfterThisInteraction(){
	return particleEnergyAfterThisInteraction;
    }

    /** add the likelihood of this single interaction to
	the overall interactions log-likelihood along a track */
    public void addThisLikelihoodValueToTrackInteractionLikelihood(double likelihoodValue){
	trackllh += Math.log(likelihoodValue);
    }

    public void initializeTrackLikelihood(){
	trackllh = 0.0;
    }

    /** return the present log-likelihood of track */
    public double getLogTrackInteractionLikelihodd(){
	return trackllh;
    }
}
	
    

