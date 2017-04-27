package iceCube.uhe.event;

import iceCube.uhe.particles.*;
import iceCube.uhe.decay.*;
import iceCube.uhe.interactions.*;

import numRecipes.*;

import java.util.*;

/**
<pre>

This is the abstract class to define the methods for both intereactions and decay
that determines the pathlength and produced energy with the Monte Carlo method.
The acutual implementaion of the mothods are made by InteractionsBase.java
(in the intereactions package) and MuDecayBase/TauDecayBase.java (in the decay package).

This class also defines the method to tell you particle flavor and doublet involved
so that the Event class can decide which interaction (or decay) should be
taken for a given primary particle.
</pre>
*/


public abstract class MonteCarloBase {
    
    /** Get the pathlength */                                                        
    public abstract double getPathLength(int iLogE, RandomGenerator rand);
    public abstract double getPathLength(double logEnergy, RandomGenerator rand);
    public abstract double getNeutrinoPathLength(int iLogE, RandomGenerator rand);
    public abstract double getNeutrinoPathLength(double logEnergy, RandomGenerator rand);

    /** Get the logEnergy of the produced particle */
    public abstract double getProducedEnergy(int iLogE, RandomGenerator rand);
    public abstract double getProducedEnergy(double logEnergy, RandomGenerator rand);

    /** Get information of the particle propagating*/
    public abstract int getPropFlavor();
    public abstract int getPropDoublet();

    /** Get information of the produced particle */
    public abstract int getProducedFlavor();

    /** Get name of the interaction */
    public abstract String getInteractionName();

    /** Get type of the interaction (Interaction->0; Decay->1) */
    public abstract int getTypeOfInteraction();

}

