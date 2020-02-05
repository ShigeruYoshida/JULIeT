package iceCube.uhe.interactions;

import numRecipes.*;
import geometry.*;
import iceCube.uhe.geometry.*;
import iceCube.uhe.particles.*;
import iceCube.uhe.points.*;
import iceCube.uhe.event.*;

import java.io.*;
import java.util.*;

public class InteractionsLikelihoodBuilder {

    private final static double ln10 = Math.log(10.0);

    /* Calculate interactions likelihood of a given lepton track stored in the list containers
       created by the JulietEventGenerator object.
       <pre>
       InteractionsLikelihood intLikelihood : InteractionsLikelihood object to calculate likelihood 
       ListIterator trackIterator           : list iterator to access the track particle list
       ListIterator particleIterator        : list iterator to access the secondary cascade particle list
       ListIterator particleLocationIterator: list iterator to access the cascade location list
       J3Vector startPosition_ice3          : location of the in-ice track start position in IceCube coordinate
       double localDensity                  : density of medium (ice) [g/cm2]
       </pre>
       */
    public static double 
	calculateInteractionsLikelihoodFromJulietLists(InteractionsLikelihood intLikelihood,
						       ListIterator trackIterator,
						       ListIterator particleIterator,
						       ListIterator particleLocationIterator,
						       J3Vector startPosition_ice3,
						       double localDensity){

	double logEthreshold = intLikelihood.logThresholdEnergy;
	double logLikelihoodValue = 0.0;
	J3Vector r_start = new J3Vector();
	r_start.putVector(startPosition_ice3);
	double inIce_initial_energy = ((Particle )(trackIterator.next())).getEnergy();
	//System.out.format("initial inice energy =%e log(threshold E)=%5.2f\n",inIce_initial_energy,logEthreshold);
	double iniceEnergy = inIce_initial_energy;
	double stochastic_int_pathLength = 0.0;
	double iniceEnergy_stochastic = inIce_initial_energy;
	int ndf = 0;
	while(particleLocationIterator.hasNext()){
	    Particle particle = (Particle )(particleIterator.next());
	    String particleName = particle.particleName(particle.getFlavor(),
							particle.getDoublet());
	    double cascade_energy = particle.getEnergy();
	    double logProducedEnergy = Math.log(cascade_energy)/ln10;
		
	    J3Vector r_end = (J3Vector )(particleLocationIterator.next());
	    J3Vector distanceBetweenVertex = J3Vector.subtract(r_start,r_end);
	    double pathLength = distanceBetweenVertex.getLength()*localDensity;
	    stochastic_int_pathLength += pathLength;
	    iniceEnergy -= cascade_energy;
	    r_start.putVector(r_end);

	    if(logProducedEnergy> logEthreshold){ // Stochastic energy loss
		double logInIceEnergy = Math.log(iniceEnergy_stochastic)/ln10;
		LikelihoodData likelihoodData =
		    intLikelihood.getInteractionLikelihood(logInIceEnergy,logProducedEnergy,
							   stochastic_int_pathLength);
		double prob = likelihoodData.getLikelihoodValue();
		if(prob>0.0){
		    logLikelihoodValue += -Math.log(prob);
		    ndf++;
		}
		//System.out.println(" secondary " + particleName + " " + 
		//		   cascade_energy + " [GeV] " +
		//		   stochastic_int_pathLength + " prob = " + prob);
		iniceEnergy_stochastic = iniceEnergy;
		stochastic_int_pathLength = 0.0;
	    }
	}
	if(ndf>0){
	    return logLikelihoodValue/(double)ndf;
	}else{
	    return Double.POSITIVE_INFINITY;
	}

    }

}
