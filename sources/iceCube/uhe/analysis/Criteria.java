package iceCube.uhe.analysis;

import iceCube.uhe.particles.*;
import iceCube.uhe.analysis.*;
import geometry.*;

/**
   Criteria class defines requrements of the IceCube EHE events
   that are subject to futher analysis. This will be "signal cut".

   <pre>
   The criteria you can set are the followings :
   (a) log Npe > certain value
          AND
   (a') (cosZenith - cosZenMin)*(logNpeMax-logNpeMin)
               -(cosZenMax-cosZenMin)(logNpe-logNpeMin)<0
   (a'')  OR
         cosZenith < highestOfCosZenith && logNpe > lowestOfLogNpe
         AND
   (a''') zenith angle simple cut
         cosZenith < cosZenCut downCut = true
         cosZenith > cosZenCut downCut = false
	 AND
   (a'''') Npe Simple cut
          logNpe > logNpeBound npeMinimumCut = true
          logNpe < logNpeBound npeMinimumCut = false

         AND
   (b)  Number of DOMs > certain value
         AND
   (c)  first guess Quality > certain value
         AND
   (d) COB/Cascade vertex - origin you set < max threshold of distance

   Theshold values in each of these can be set by calling the relevant
   method.
   </pre>

   When the method setEHESuperCut() is called,,however,
   the criterua set for the EHE super CUT determines
   if a given event passes the criteria. The present super CUT
   developed by Aya Ishihara (aya@hepburn.s.chiba-u.jp) categorize
   events into the four groups as
   <pre>

    Npe-ATWD > Npe-FADC  &&  | Npe-ATWD < Npe-FADC  &&
    cobZ =[min max]          | cobZ =[min max]  
    ----------------------------------------------------
    Npe-ATWD > Npe-FADC  &&  | Npe-ATWD < Npe-FADC  &&
    cobZ out of [min max]    | cobZ  out of [min max]  


   Written originally by S. Yoshida for the IceCube EHE anaysis
   2007/1/26

   Rewritten by S. Yoshida to include the EHE super Cut 2007/12/01
   </pre>

*/
public class Criteria {
    private final double ln10 = Math.log(10.0);

    private double minRangeOfLogNpe =  I3ParticleAnalysisFactory.minLogNpe;
    private double maxRangeOfLogNpe =  I3ParticleAnalysisFactory.minLogNpe;
    private double thresholdOfLogNpe = I3ParticleAnalysisFactory.minLogNpe;
    private double lowestOfLogNpe = I3ParticleAnalysisFactory.maxLogNpe;
    private double highestOfCosZenith = -1.0;
    private double minRangeOfCosZenith = -1.0;
    private double maxRangeOfCosZenith = 1.0;
    private double cosZenCut = 1.0; private boolean downCut = true;
    private double logNpeBound = 20.0; private boolean npeMinimumCut = false;
    private int minNDOMs = 0;
    private double fgThreshold = 0.0;
    private double maxDistance = 1.0e6; // 1,0e6 cm = 10 km
    private double logNpeScaleFactor = 0.0;
    private J3Vector origin = null;

    //The variables concerning the EHE super cut
    protected boolean isEHESuperCut = false;
    protected boolean isCOBZCut = false;
    protected double minCobZ = -4.0e4; // -4.0e4 cm = -400 m
    protected double maxCobZ =  4.0e4; //  4.0e4 cm =  400 m
    protected double min2CobZ = 4.0e6; // 4.0e6 cm = 40 km - not used.
    protected double max2CobZ = 5.0e6; // 5.0e6 cm = 50 km - not used
    protected int minNDOMsInEHESuperCut = 80;
    protected J3Vector[][] vertexLocation = null; // [category][vertex] (lonNpe, cosTheta)
    /**  number of the event category */
    static protected int numberOfEventCategory = 4;
    /**  maximum number of the vertex points on the logNpe-cosTheta plane
	 to define the EHE super cut */
    protected int maxNumberOfVertex = 5;
    /**  number of the vertex points on the logNpe-cosTheta plane
	 set by the method setEHESuperCut() */
    protected int[] numberOfVertex;
    /** The default vertex X location on the logNpe-CosZenith plane
	in the EHE super cut. Used in the setEHESuperCut()
    */
    protected static double[][] vertexDefaultLocationX = 
    {
	{4.6, 4.6, 5.6, 6.2, 6.5},
	{5.2, 5.2, 5.4, 6.0, 6.0},
	{4.3, 4.3, 4.7, 5.45, 5.45},
	{4.1, 4.1, 4.7, 5.45, 5.45}

    };
    /** The default vertex y location on the logNpe-CosZenith plane
	in the EHE super cut. Used in the setEHESuperCut()
    */
    protected static double[][] vertexDefaultLocationY = 
    {
	{-1.0, 0.1, 0.2, 0.4, 1.0},
	{-1.0, 0.0, 0.0, 0.4, 1.0},
	{-1.0, 0.1, 0.1, 0.5, 1.0},
	{-1.0, -0.5, -0.5, 0.3, 1.0}
    };


    /** constructor. Nothing is determined */
    public Criteria(){
    }

    /** Set the threshold of Log10(Npe) 
	no matter where an event comes from.
	This requirement corresponds
	to (a) in the description you find at top of this class.
    */
    public void setThresholdOfLogNpe(double logNpe){
	thresholdOfLogNpe = logNpe;
    }

    /** Sets the range of the line on cosZenith -logNpe plain.
	This requirement corresponds
	to (a') in the description you find at top of this class.
    */
    public void setRangeOfCosineOfZenith(double minCosTheta, double maxCosTheta,
					 double minLogNpe, double maxLogNpe){
	minRangeOfCosZenith = minCosTheta;
	maxRangeOfCosZenith = maxCosTheta;
	minRangeOfLogNpe = minLogNpe;
	maxRangeOfLogNpe = maxLogNpe;
    }

    /** Set the highest cosZenith and the lowest LogNpe -- Criterion (a'')*/
    public void setMinimumBound(double highestOfCosZenith, double lowestOfLogNpe){
	this.highestOfCosZenith = highestOfCosZenith;
	this.lowestOfLogNpe = lowestOfLogNpe;
    }


    /** Set the simple zenith Cut  -- Criterion (a''')/

     */
    public void setSimpleCosZenithCut(double cosZenCut, boolean downCut){
	this.cosZenCut = cosZenCut;
	this.downCut = downCut;
    }

    /** Set the simple Npe Cut  -- Criterion (a'''')/

     */
    public void setSimpleNpeCut(double logNpeBound, boolean npeMinimumCut){
	this.logNpeBound = logNpeBound;
	this.npeMinimumCut = npeMinimumCut;
    }


    /** Set the threshold of Number Of DOMs. Condition (b) */
    public void setThresholdOfNDOMs(int nDOMs){
	minNDOMs = nDOMs;
    }

    /** Sets the threshold of the First Guess fitting quality. Condition (c)*/
    public void setThresholdOfFirstGuessQuality(double fgThreshold){
	this.fgThreshold = fgThreshold;
    }

    /** Sets the maximum distance : The COB/Cascade vertex - origin
        must be shorter than this value. Condition (d)
    */
    public void setMaxDistance(double distance, J3Vector origin){
	maxDistance = distance;
	this.origin = origin;
    }

    /** Sets the COBZ cut */
    public void setCOBZcut(double cobZmin, double cobZmax){
	isCOBZCut = true;
	minCobZ = cobZmin;
	maxCobZ = cobZmax;
    }
    /** Sets the COBZ cut with the Default COB-Z range */
    public void setCOBZcut(){
	isCOBZCut = true;
    }
    /** Unsets the COBZ cut */
    public void unsetCOBZcut(){
	isCOBZCut = false;
    }

    /** Set the NPE scaling factor (default:1.0). This factor
	accounts a possible rescaling of Energy-NPE due to 
        the systematic errors. NPE of I3Particle is multiplied
        by this factor whenever the event cut is applied.
    */
    public void setNPEScalingFactor(double factor){
	if(factor>0.0) logNpeScaleFactor = Math.log(factor)/ln10;
    }

    /** Determine if the iceParticle satisfies the criteria */
    public boolean doesThisEventPass(I3Particle iceParticle){

	if(!isEHESuperCut) return doesThisEventPassSimpleCut(iceParticle);
	else return doesThisEventPassEHESuperCut(iceParticle);

    }


    /** 
	Tell Whether this event passes the "simple cut" .
        Called inside the method doesThisEventPass() when
        the EHE super cut does not apply.
     */
    protected boolean doesThisEventPassSimpleCut(I3Particle iceParticle){
	// Number Of PEs and hitted DOMs of this event
	double logNpe = iceParticle.getIceCubeData().getLogBestNpe()
	    + logNpeScaleFactor;
	int nDOMs = iceParticle.getIceCubeData().getNDOMsLaunch();

	// energy
	double logEnergy = iceParticle.getLogEnergy();
	double logRecoEnergy = iceParticle.getLogRecoEnergy();

	// direction
	J3UnitVector n = iceParticle.getDirectionInIceCubeCoordinate();
	double cosZenith = -n.getZ(); // Reversed vector

	// First Guess fitting Quality index
	double fgQ = iceParticle.getFirstGuessQuality();

	// Now Check if this event passes the criteria
	boolean passThisEvent = false;
	if((logNpe >= thresholdOfLogNpe) && 
	   ((cosZenith <= cosZenCut && downCut)|| (cosZenCut <= cosZenith && !downCut)) &&
	   ((logNpeBound <= logNpe && npeMinimumCut) || (logNpe < logNpeBound && !npeMinimumCut)) &&
	   ((((cosZenith-minRangeOfCosZenith)*(maxRangeOfLogNpe-minRangeOfLogNpe)
	      -(maxRangeOfCosZenith-minRangeOfCosZenith)*(logNpe-minRangeOfLogNpe))<0.0)||
	    (cosZenith< highestOfCosZenith && lowestOfLogNpe < logNpe))){

	    if(fgQ>=fgThreshold && nDOMs >= minNDOMs){
		if(origin!= null){ // Then check the COB/Cascade vertex position
		    J3Vector r = iceParticle.getR0InIceCubeCoordinate();
		    J3Vector cobIC9 = J3Vector.subtract(r,origin);
		    double cobDistance = cobIC9.getLength();
		    if(maxDistance>=cobDistance) passThisEvent = true;
		    else passThisEvent = false;
		}else if(isCOBZCut){ // COB Z cut
		    return isCOBZwithinRange(iceParticle);
		}else{
		    passThisEvent = true;
		}
	    }else{
		passThisEvent = false;
	    }

	}

	return passThisEvent;
    }


    /** 
	Tell Whether this event passes the "EHE super cut" .
        Called inside the method doesThisEventPass() when
        the EHE super cut applies. The criteria is defined
        for each of the four following event categories :
       <pre>

       <category 0>               <category 2>
       Npe-ATWD > Npe-FADC  &&  | Npe-ATWD < Npe-FADC  &&
       cobZ =[min max]          | cobZ =[min max]  
       ----------------------------------------------------
       <category 1>               <category 3>
       Npe-ATWD > Npe-FADC  &&  | Npe-ATWD < Npe-FADC  &&
       cobZ out of [min max]    | cobZ  out of [min max]  

       </pre>

     */
    protected boolean doesThisEventPassEHESuperCut(I3Particle iceParticle){

	int eventCategory;
	if(isATWDNpeLarger(iceParticle)){
	    if(isCOBZwithinRange(iceParticle)) eventCategory = 0;
	    else eventCategory = 1;
	}else{
	    if(isCOBZwithinRange(iceParticle)) eventCategory = 2;
	    else eventCategory = 3;
	}

	if(isEHESuperCut){
	    return isThisEventInTheGZKRegion(iceParticle,eventCategory);
	}else{ // The relevant parameters not setting!
	    System.err.println("Must call setEHESuperCut() to set the parameters");
	    System.exit(0);
	    return false;
	}
    }

    /**
       Tell if this event is within the GZK boundary
       set by the method setEHESuperCUT(category, vertex)
       in the corresponding event category. Note that
       in the "EHE Super Cut", the GZK boundary in the logNPE-cosZenith
       plane is defined individually in each of the category.
    */
    protected boolean isThisEventInTheGZKRegion(I3Particle iceParticle,
						int category){

	int nDOMs = iceParticle.getIceCubeData().getNDOMsLaunch();
	double logNpe = iceParticle.getIceCubeData().getLogBestNpe()
	    + logNpeScaleFactor;
	J3UnitVector n = iceParticle.getDirectionInIceCubeCoordinate();
	double cosZenith = -n.getZ(); // Reversed vector

	boolean isOnTheRightSide = false;
	if(nDOMs<minNDOMsInEHESuperCut) return isOnTheRightSide;

	for(int number = 1; number<numberOfVertex[category];number++){

	    double logNpeLeft = vertexLocation[category][number-1].getX();
	    double logNpeRight = vertexLocation[category][number].getX();
	    double cosZenithLeft = vertexLocation[category][number-1].getY();
	    double cosZenithRight = vertexLocation[category][number].getY();

	    if(cosZenithLeft!=cosZenithRight){
		// See if this event is within the range of cosZenith in the line
		if(((cosZenith-cosZenithLeft)*(cosZenith-cosZenithRight))<0.0){
		    double gradiant = 
			(logNpeRight-logNpeLeft)/(cosZenithRight-cosZenithLeft);
		    double logNpeBoundary = 
			logNpeLeft + gradiant*(cosZenith-cosZenithLeft);
		    if(logNpe >= logNpeBoundary) isOnTheRightSide = true;
		    else{
			isOnTheRightSide = false;
			break;
		    }
		}
	    }

	}

	return isOnTheRightSide;
    }

    /** set all the GZK boundary values of Npe and cosZenith
        in the "EHE Super Cut" for each of the four event categories:
       <pre>

       <category 0>               <category 2>
       Npe-ATWD > Npe-FADC  &&  | Npe-ATWD < Npe-FADC  &&
       cobZ =[min max]          | cobZ =[min max]  
       ----------------------------------------------------
       <category 1>               <category 3>
       Npe-ATWD > Npe-FADC  &&  | Npe-ATWD < Npe-FADC  &&
       cobZ out of [min max]    | cobZ  out of [min max]  

       </pre>
       The GZK boundary are defined by connecting the lines
       from the i-th point to the i+1 th point. When you call
       this method for the first time, the first point is set
       by J3Vector vertex. The 2nd point is set in the next call
       of this method. Just keep going.
       By doing so you can set the successive points
       to form the GZK boundary for a given category.
       The vector Vertex should contain (x,y) = (logNpe,cosZenith).
    */
    protected void setEHESuperCut(int category, J3Vector vertex){

	if(!isEHESuperCut){
	    vertexLocation = new J3Vector[numberOfEventCategory][maxNumberOfVertex];
	    numberOfVertex = new int[numberOfEventCategory];
	    for(int i_category = 0;i_category<numberOfEventCategory;i_category++)
		numberOfVertex[i_category] = 0;
	}
	isEHESuperCut = true;

	if(category<maxNumberOfVertex && 
	   (numberOfVertex[category]<maxNumberOfVertex)){
	    vertexLocation[category][numberOfVertex[category]] = vertex;
	    numberOfVertex[category]++;
	}else{
	    System.err.println("either category number(" +
			       category + ") is out of range");
	    System.err.println(" or number of the vertex you set (" +
			      numberOfVertex[category] + ") is out of range");
	    System.exit(0);
	}
    }


    /** set all the GZK boundary values of Npe and cosZenith
        in the "EHE Super Cut" for a given event categories
        with the default settings.
    */
    protected void setEHESuperCut(int category){
	for(int number = 0; number<maxNumberOfVertex;number++){
	    J3Vector vertex = new J3Vector(vertexDefaultLocationX[category][number],
					   vertexDefaultLocationY[category][number],
					   0.0);
	    setEHESuperCut(category,vertex);
	}

    }


    /** set all the GZK boundary values of Npe and cosZenith
        in the "EHE Super Cut" for all of the event categories
        with the default settings.
    */
    protected void setEHESuperCut(){
	for(int i_category = 0;i_category<numberOfEventCategory;i_category++)
	    setEHESuperCut(i_category);
    }


    /*** Tell whether this event has ATWD-based NPEs larger than 
	 ones based on FADC. Used in the EHE super Cut */
    protected boolean isATWDNpeLarger(I3Particle iceParticle){

	double atwdNPE = iceParticle.getIceCubeData().getNpeATWD();
	double fadcNPE = iceParticle.getIceCubeData().getNpeFADC();
	if(atwdNPE>=fadcNPE) return true;
	else return false;

    }

    /*** Tell whether this event has z depth of the Center of Brightness
         vertex within a given range set by the method setEHESuperCut().
	 Used in the EHE super Cut 
    */
    protected boolean isCOBZwithinRange(I3Particle iceParticle){

	boolean wasMCTruth = true;
	if(!iceParticle.isMCTruth()) wasMCTruth = false;

	iceParticle.switchToReco(); // now the first-guess geometry
	double cobZ = iceParticle.getR0InIceCubeCoordinate().getZ();

	// back to the original status
	if(wasMCTruth) iceParticle.switchToMCTruth();
	else iceParticle.switchToReco();

	if(((minCobZ <= cobZ ) && (cobZ <= maxCobZ)) || 
	   ((min2CobZ <= cobZ ) && (cobZ <= max2CobZ))){
		return true;
	}else{
	    return false;
	}
    }





}


