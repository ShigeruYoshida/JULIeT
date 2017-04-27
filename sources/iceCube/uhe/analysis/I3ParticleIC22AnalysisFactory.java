package iceCube.uhe.analysis;

import iceCube.uhe.particles.*;
import iceCube.uhe.analysis.*;

import hep.aida.*;

import java.io.*;
import java.util.*;

/** I3Particles Analysis Factory for the IC22 data.

    Written originally by S. Yoshida for the IceCube EHE analysis.
    2008/6/13
*/
public class I3ParticleIC22AnalysisFactory extends I3ParticleAnalysisFactory {

    //protected double observationTime = 1.01088e7; // 117 days in [sec]
    //protected double observationTime = 2.083968e7; //  241.2 days in [sec]
    protected double observationTime = 2.091744e7; //  242.1 days in [sec]
    //public final static int[] badRunIC22 = { // including short runs
    //108895, 108896, 108911, 108912, 109049, 109097, 109117, 109135,
    //109136, 109137, 109480, 109495, 109512};
    public final static int[] badRunIC22 = {
	108895, 108896, 108911, 108912,109049,109097,110303,110527,110648,110655,110715};
    //public static final int startRunIC22 = 108874;  // 1st Run ID to be analized
    public static final int startRunIC22 = 107933;  // 1st Run ID to be analized
    //public static final int endRunIC22 = 109792;  // end Run ID to be analized
    public static final int endRunIC22 = 110773;  // end Run ID to be analized
    public static final int[] excludedRunIC22 = {108715,108826}; // June16-30

    public I3ParticleIC22AnalysisFactory(InputStream in ) throws IOException{
	super(in);
	setObservationTime(observationTime);
    }

    public I3ParticleIC22AnalysisFactory(InputStream in,  
					 boolean filterOutBadRunData) throws IOException{
	super(in,filterOutBadRunData);
	setObservationTime(observationTime);
    }

    /** Judge if this event has to be excluded because of
        the bad run. You have to set filterOutBadRunData(true)
        to get this method effective. Implemeneted fowllowing the IC22 bad run list
    */
    protected boolean isBadRunData(I3Particle iceParticle){
	if(!filterOutBadRunData) return false;
	int runID = iceParticle.getIceCubeData().getEventNumber();
	if(startRunIC22 <= runID && runID <= endRunIC22){
	    if(excludedRunIC22[0] <= runID && runID <= excludedRunIC22[1]){
		System.err.println("Bad Run! " + runID);
		return true;
	    }
	    for(int i=0; i<badRunIC22.length; i++){
		if(badRunIC22[i] == runID){
		    System.err.println("Bad Run! " + runID);
		    return true;
		}
	    }
	    return false;
	}else{
	    System.err.println("This Run is too early or too late! " + runID);
	    return true;
	}
    }



}
