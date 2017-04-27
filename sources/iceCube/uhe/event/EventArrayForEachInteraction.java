package  iceCube.uhe.event;

import iceCube.uhe.particles.*;
import iceCube.uhe.interactions.*;
import java.io.*;

/**
   This object handles a matrix table generated by 
   RunManager.runEventOnMatrix(DataOutputStream out, int numEvent).
   It gives a dN/dLogEcascade(mu or tau) where Ecascade is a total
   energy deposit in form of electrons or hadrons initiated by
   UHE muons or taus running in the detector voluce.
   Propagation of muons and tauons are treated by Event.java.

   The matrix table this object reads out is similar with 
   PropagationMatrix.getFtauToHadron(int iLogE,int jLogE) or
   PropagationMatrix.getFmuToHadron(int iLogE,int jLogE),
   for example, 
   in the propagation package. A main difference is 
   the PropagationMatrix is calculated by numerically
   solving the transport equation which is suitable
   for particle propagation in the Earth over long range
   while THIS matrix table is calculated by the Monte Carlo
   method by Event.java which works especially for particle
   propagation inside the detector volume where the fluctuation
   of the energy loss profile and total energy deposit from an event
   do matter.
  
*/

public class EventArrayForEachInteraction {

    static final int mc  = 7; // length of mcBases[]  muStandard=7, tauStandard=15  
    static int dim = Particle.getDimensionOfLogEnergyMatrix();
    static int expandedDim = dim + 
                       (int )(( Particle.getLogEnergyMinimum()-InteractionsBase.getLogEnergyProducedMinimum())/
			      Particle.getDeltaLogEnergy());
    private double[][] cascadeMtx;
    private int primaryFlavor = 0; //Flavor of the primary propagating particles.
    private int primaryDoublet = 0; //Doublet of the primary propagating particles.
    private String interactionName = null;

    /** Constructor. Generate matrix array.*/
    public EventArrayForEachInteraction() {
        cascadeMtx    = new double[mc][expandedDim];
    }

    public EventArrayForEachInteraction(DataInputStream in ) throws IOException {
        cascadeMtx    = new double[mc][expandedDim];
	readArray(in);
    }

    /** Read the calculated event matrix by RunManager.runEventOnMatrix */
    public void readArray(DataInputStream in) throws IOException {
	
        for(int n=0; n<mc; n++){
            for(int jLogE=0;jLogE<expandedDim;jLogE++){
                cascadeMtx[n][jLogE]= in.readDouble();
		System.err.println(cascadeMtx[n][jLogE]);
            }
        }

        in.close();
    }

    /** Obtain primarily propagating particle flavor.
	As defined in the Particle class, it should be either
	1 (muon) or 2 (tau). This value is read from the Matrix file
	by the method readArray() in this class.
    */

    public int getPrimaryFlavor(){
	return(primaryFlavor);
    }

    /** Obtain primarily propagating particle doublet.
	As defined in the Particle class, it should be either
	0 (neutrino) or 1 (charged lepton). This value is read from the Matrix file
	by the method readArray() in this class.
    */

    public int getPrimaryDoublet(){
	return(primaryDoublet);
    }

    /** Obtain the name of interaction.
	This is read from the Matrix file
	by the method readArray() in this class.
    */

    public String getInteractionName(){
	return(interactionName);
    }



    /** Obtain dN/dLogE */
    public double getCascadeFlux(int n, double logE){

	int jLogE = (int)((logE - InteractionsBase.getLogEnergyProducedMinimum())
			  /Particle.getDeltaLogEnergy());

	double count = 0.0;
	if((0<=n && n<mc) && (0<=jLogE && jLogE<expandedDim)){
	    count = cascadeMtx[n][jLogE];
	}
	return(count);
   }

}

