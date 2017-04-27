package iceCube.uhe.analysis;

import iceCube.uhe.particles.*;
import iceCube.uhe.propagation.*;
import iceCube.uhe.analysis.*;
import iceCube.uhe.muonModel.*;
import geometry.*;
import numRecipes.*;

import hep.aida.*;
import hep.aida.ext.*;
import hep.aida.util.*;
import hep.aida.util.comparison.*;

import java.io.*;
import java.util.*;

/** I3Particles Analysis Factory. Utility methods to help your analysis
    using I3Particles. A main function is making a histogram that 
    consists of events satisfying the criteria defined by Criteria class. 

    Written originally by S. Yoshida for the IceCube EHE analysis.
    2007/1/21
*/
public class I3ParticleAnalysisFactory {

    protected static final double ln10 = Math.log(10.0);
    protected static double xCenterOfIC9 = 3.10e4; // Center of IC9 array - X[cm]
    protected static double yCenterOfIC9 = 1.89e4; // Center of IC9 array - Y[cm]

    // Parameters on  histogram
    protected double minLogNumberOfEvents = -2.0;

    private double deltaLogNpe = 0.1;
    protected static double maxLogNpe = 7.0;
    protected static double minLogNpe = 3.0;
    private double deltaLogEnergy = 0.1;
    protected static double maxLogEnergy = 11.0; // 10^11 [GeV]
    protected static double minLogEnergy = 5.0;  // 10^5  [GeV]

    private double deltaCosZenith = 0.5;

    private double maxNDOMsOnDrawing = 601.0;
    private double deltaNDOMs = 50.0;

    private double deltaFgQuality = 0.1*PropagationMatrix.c;
    private double maxFgQuality = 2.0*PropagationMatrix.c;

    /** Criteria of events that are subject to analysis */
    private Criteria criteria = null;
    /** The lowest level Requirement of NDOMs. For saving memory.
	Only I3Particles to pass this condition will be stored
        in the container and subject to further analysis with
        higher level cuts defined by the Criteral class
     */
    protected static int minNDOMsToAnalize = 32;
    /** The lowest level Requirement of NPEs. For saving memory.
	Only I3Particles to pass this condition will be stored
        in the container and subject to further analysis with
        higher level cuts defined by the Criteral class
     */
    protected static double minLogNPEToAnalize = 0.0;

    /** MCTruth Flag for the I3Particles.*/
    private boolean isMCTruth = false;

    /** Flag for weighting/unweighting */
    private boolean isWeighted = false;
    private boolean isAtmMuon = true;

    // Model name of GZK/Atm Muon fluxes stored in I3Particle 
    private String modelName = null;

    // Options for the "grafig"-based Xfig drawing
    private boolean plotHistogram = false; // if false, plot data points 
                                           //with error bars
    private boolean drawByAutoScale = true;

    private int numberOfEvents = 0;

    private RandomGenerator randomGen = null;

    // Container of the I3Particles
    private List i3particleList = null;
    private ListIterator i3particleIterator = null;

    private List  numberOfAllI3ParticlesList;
    private List  numberOfFilledI3ParticlesList;
    protected int numberOfFilledI3Particles = 0;

    /** Area where juliet particles are thrown in the MC simulation .
	Unit should be [cm^2].
     */
    protected double mcArea = 880.0*880.0*Math.PI*1.0e4; // 880*800*pi [m^2]
    /** Solid angle where juliet particles are thrown in 
	the MC simulation */
    protected double mcOmega = 4.0*Math.PI;
    /** Area of the Corsika MMC particles are thrown */
    protected double mcCorsikaArea = Math.PI*880.0*(880.0+1760.0)*1.0e4; //[cm^2]
    /** Solid angle of the Corsika mmc particles */
    protected double mcCorsikaOmega = Math.PI;
    /** Observation Time [sec]*/
    protected double observationTime = 1.0726387e7; // 124.148 days (IC9 data)

    /** AtmMuonBundleFlux object for calculating the primary cosmic ray energy */
    protected AtmMuonBundleFlux muonBundle = null;

    /**
       Bad Run Numbers.. for the real data analsis only
    */
    public final static int[][] badRunID = {
	{89658, 89712} // FY 2006 data sample.
    };

    /** Flag for filtering out the bad run data */
    protected boolean filterOutBadRunData = false;

    private double epsilon = 1.0e-2;

    private int dimensionLogNpe = (int )((maxLogNpe-minLogNpe)/deltaLogNpe + epsilon);
    private int dimensionLogEnergy = 
	(int )((maxLogEnergy-minLogEnergy)/deltaLogEnergy + epsilon);
    private int dimensionNDOMs = (int )((maxNDOMsOnDrawing + epsilon)/deltaNDOMs);
    private int dimensionCosZenith = (int )(2.0/deltaCosZenith + epsilon);
    private int dimensionFgQuality = (int)(maxFgQuality/deltaFgQuality+epsilon);

    protected double[][] histLogNpeCosZenith;
    protected double[][] histLogNpeNDOM;
    protected double[][] histLogNpeFgQuality;
    protected double[][] histCosZenithFgQuality;
    protected double[][] histLogNpeLogEnergy;
    protected double[][] histLogEnergyCosZenith;

    /** Default Constructor 
	<pre>
	InputStream in  : Stream to readout a series of I3Particles that is subject
                          to your analysis
        </pre>
        It reads all I3Particles from InputStream and hold them in form
        of List.
     */
    public I3ParticleAnalysisFactory(InputStream in) throws IOException{
	i3particleList          = new LinkedList();
	numberOfAllI3ParticlesList = new LinkedList();
	numberOfFilledI3ParticlesList = new LinkedList();
	numberOfFilledI3Particles = 0;
	readI3Particles(in);
    }

    /** Constructor 
	<pre>
	InputStream in  : Stream to readout a series of I3Particles that is subject
                          to your analysis
        filterOutBadRunData  :  do not use events in the bad run. The real data analysis only.
        </pre>
        It reads all I3Particles from InputStream and hold them in form
        of List.
     */
    public I3ParticleAnalysisFactory(InputStream in, boolean filterOutBadRunData) throws IOException{
	this.filterOutBadRunData = filterOutBadRunData;
	i3particleList          = new LinkedList();
	numberOfAllI3ParticlesList = new LinkedList();
	numberOfFilledI3ParticlesList = new LinkedList();
	numberOfFilledI3Particles = 0;
	readI3Particles(in);
    }


    /** reading the I3Particle objects */
    protected void readI3Particles(InputStream in) throws IOException{

	I3Particle iceParticle = null; 
	int numberOfAllI3Particles = 0;
	while((iceParticle = I3ParticleInputStream.inputI3Particle(in))
	      !=null){
	    if(!isBadRunData(iceParticle)) numberOfAllI3Particles++;
	    if(minNDOMsToAnalize<=iceParticle.getIceCubeData().getNDOMsLaunch() && 
	       minLogNPEToAnalize<=iceParticle.getIceCubeData().getLogBestNpe() && 
	       !isBadRunData(iceParticle)){
		I3Particle iceLite = getI3ParticleLite(iceParticle);
		i3particleList.add(iceLite);
		numberOfFilledI3Particles++;
	    }
	}
	numberOfAllI3ParticlesList.add(new Integer(numberOfAllI3Particles));
	numberOfFilledI3ParticlesList.add(new Integer(numberOfFilledI3Particles));

	System.err.println("Number of I3Particles before applying any criteria" +
			   numberOfAllI3Particles);
	System.err.println("Number of Filled I3Particles" + 
			   numberOfFilledI3Particles);
    }

    /** make I3Particle as a temporal dataclass for the analysis.
	This method is introduced for saving memory
    */
    protected I3Particle getI3ParticleLite(I3Particle iceParticle){
	int flavor = iceParticle.getFlavor();
	int doublet = iceParticle.getDoublet();
	// Generate an I3Particle object
	I3Particle iceLite = new I3Particle(flavor,doublet);

	// sets geometry
	if(!filterOutBadRunData){
	    iceParticle.switchToMCTruth();
	    iceLite.switchToMCTruth();
	    J3UnitVector n = iceParticle.getDirectionInIceCubeCoordinate();
	    J3Vector cog = iceParticle.getR0InIceCubeCoordinate();
	    J3Line axisInIce3 = new J3Line(cog,n);
	    iceLite.setParticleAxisInIceCubeCoordinate(axisInIce3);
	    iceLite.transformParticleAxisToEarthCenterCoordinate();
	}
	iceParticle.switchToReco();
	iceLite.switchToReco();
	J3UnitVector n = iceParticle.getDirectionInIceCubeCoordinate();
	J3Vector cog = iceParticle.getR0InIceCubeCoordinate();
	J3Line axisInIce3 = new J3Line(cog,n);
	iceLite.setParticleAxisInIceCubeCoordinate(axisInIce3);
	iceLite.transformParticleAxisToEarthCenterCoordinate();

	// sets energy
	iceLite.putEnergy(iceParticle.getEnergy());
	iceLite.putLogEnergy(iceParticle.getLogEnergy());
	iceLite.putRecoEnergy(iceParticle.getRecoEnergy());

	// Sets the distance from the Earth surface
	double distance = iceParticle.getDistanceFromEarthSurfaceToIceCube();
	iceLite.setDistanceFromEarthSurfaceToIceCube(distance);

	//Sets the first guess quality
	double fgQ = iceParticle.getFirstGuessQuality();
	iceLite.setFirstGuessQuality(fgQ);

	// Fills the IceCube data
	int eventNumber = iceParticle.getIceCubeData().getEventNumber();
	iceLite.getIceCubeData().setEventNumber(eventNumber);
	double npeATWD = iceParticle.getIceCubeData().getNpeATWD();
	iceLite.getIceCubeData().setNpeATWD(npeATWD);
	double npeFADC = iceParticle.getIceCubeData().getNpeFADC();
	iceLite.getIceCubeData().setNpeFADC(npeFADC);
	double npeBest = iceParticle.getIceCubeData().getBestNpe();
	iceLite.getIceCubeData().setBestNpe(npeBest);
	int nDOMsATWD = iceParticle.getIceCubeData().getNDOMsATWD();
	iceLite.getIceCubeData().setNDOMsATWD(nDOMsATWD);
	int nDOMsFADC = iceParticle.getIceCubeData().getNDOMsFADC();
	iceLite.getIceCubeData().setNDOMsFADC(nDOMsFADC);
	int nDOMsLaunch = iceParticle.getIceCubeData().getNDOMsLaunch();
	iceLite.getIceCubeData().setNDOMsLaunch(nDOMsLaunch);
	//System.err.println(" logNPE= " + iceLite.getIceCubeData().getLogBestNpe());
	//System.err.println(" NDOM= " +  iceLite.getIceCubeData().getNDOMsLaunch());

	// set the MC primary spectrum weight
	if(!filterOutBadRunData){
	    iceLite.switchToMCTruth();
	    iceLite.setMCPrimarySpectrumWeight(iceParticle.getMCPrimarySpectrumWeight());
	}

	// flux weight
	Iterator iteGZK = iceParticle.iteratorOfGZKNeutrinoFlux();
	while(iteGZK.hasNext()){
	    Map.Entry entry = (Map.Entry )(iteGZK.next());
	    Double flux = (Double )(entry.getValue());
	    String name = (String )(entry.getKey());
	    iceLite.setGZKNeutrinoFlux(flux.doubleValue(),name);
	}

	Iterator iteAtm = iceParticle.iteratorOfAtmosphericMuonFlux();
	while(iteAtm.hasNext()){
	    Map.Entry entry = (Map.Entry )(iteAtm.next());
	    Double flux = (Double )(entry.getValue());
	    String name = (String )(entry.getKey());
	    iceLite.setAtmosphericMuonFlux(flux.doubleValue(),name);
	    //System.err.println(" atm Flux(" + name + ")= "+
	    //		       iceLite.getAtmosphericMuonFlux(name) +
	    //		       " distance= " + 
	    //		       iceLite.getDistanceFromEarthSurfaceToIceCube() + 
	    //		       " log(energy)= " + iceLite.getLogEnergy());
	}

	return iceLite;
    }

    /** Return the Particle Iterator loping over I3Particles*/
    public ListIterator getParticleIterator(){
	    return i3particleList.listIterator();
    }

    /** Set the MC area [cm^2] */
    public void setMCArea(double area){
	mcArea = area;
    }

    /** Set the MC solid angle [sr] */
    public void setMCSolidAngle(double omega){
	mcOmega = omega;
    }

    /** Change the MC area and solid angle to the Corska configuration */
    public void changeAOmega(){
	mcOmega = mcCorsikaOmega;
	mcArea = mcCorsikaArea;
    }


    /** Set the MC solid angle [sec] */
    public void setObservationTime(double time){
	observationTime = time;
    }

    /** Judge if this event has to be excluded because of
        the bad run. You have to set filterOutBadRunData(true)
        to get this method effective.
    */
    protected boolean isBadRunData(I3Particle iceParticle){
	if(!filterOutBadRunData) return false;
	int runID = iceParticle.getIceCubeData().getEventNumber();
	for(int i=0; i<badRunID.length; i++){
	    if(badRunID[i][0]<= runID && runID <= badRunID[i][1]) return true;
	}
	return false;
    }

    public void autoScale(boolean flag){
	drawByAutoScale = flag;
    }

    /** Calculate the MC Primary Spectrum dN/dLogE 
	except normalization constatnt
        determined by the total numner of MC events.
	<pre>
	double powerLaw    :   dN/dE = E**(-powerLaw)
	double energy      :   Energy of this I3Particle [GeV]
	double logEnergyMinimum :  log(Energy Minimum [GeV]) in the spectral range
	double logEnergyMaximum :  log(Energy Maximum [GeV]) in the spectral range
	</pre>
     */
    protected static double getDNDLogE(double powerLaw, double energy,
				       double logEnergyMinimum,
				       double logEnergyMaximum){
	double flux;
	if(powerLaw == 1.0){
	    flux = 1.0/(logEnergyMaximum - logEnergyMinimum);
	}else {
	    double minTerm = Math.pow(10.0,logEnergyMinimum*(1.0-powerLaw));
	    double maxTerm = Math.pow(10.0,logEnergyMaximum*(1.0-powerLaw));
	    double energyTerm = Math.pow(energy,(1.0-powerLaw));
	    flux = (powerLaw-1.0)*energyTerm/(minTerm - maxTerm)*ln10;
	}
	return flux;
    }


    /** Set historgram bin size 
	<pre>
	double deltaLogEnergy  bin size of logEnergy
	double deltaLogNpe     bin size of logNpe
	double deltaCosZenith  bin size of cos(zenith angle)
	double deltaFgQuality  bin size of "First Guess Quality index" in unit of beta
	</pre>
     */
    public void setBinSize(double deltaLogEnergy,double deltaLogNpe, 
			   double deltaCosZenith,
			   double deltaFgQuality){
	setBinSize(deltaLogNpe,deltaCosZenith,deltaFgQuality);
	if(deltaLogEnergy>0.0) this.deltaLogEnergy = deltaLogEnergy;
	dimensionLogEnergy = 
	(int )((maxLogEnergy-minLogEnergy)/deltaLogEnergy + epsilon);
    }

    public void setBinSize(double deltaLogNpe, double deltaCosZenith,
			   double deltaFgQuality){
	setBinSize(deltaLogNpe,deltaCosZenith);
	if(deltaFgQuality>0) this.deltaFgQuality = deltaFgQuality*PropagationMatrix.c;
	dimensionFgQuality = (int)(maxFgQuality/this.deltaFgQuality+epsilon);
    }

    public void setBinSize(double deltaLogNpe, double deltaCosZenith){
	if(deltaLogNpe>0.0) this.deltaLogNpe = deltaLogNpe;
	if(deltaCosZenith>0.0) this.deltaCosZenith = deltaCosZenith;
	dimensionLogNpe = (int )((maxLogNpe-minLogNpe)/deltaLogNpe + epsilon);
	dimensionCosZenith = (int )(2.0/deltaCosZenith + epsilon);
    }

    /** Set Max Number Of DOMs in Drawing */
    public void setMaxRangeOfNDOMsInDrawing(int nDOMs){
	maxNDOMsOnDrawing = (double )nDOMs;
	dimensionNDOMs = (int )((maxNDOMsOnDrawing + epsilon)/deltaNDOMs);
    }


    /** Switch to parameters concerned with MCtrue.*/
    public void switchToMCTruth(){
	isMCTruth = true;
    }

    /** Switch to parameters concerned with Reco results.*/
    public void switchToReco(){
	isMCTruth = false;
    }

    /** Plotting histogram by jointed-line */
    public void plotByJointLine(){
	plotHistogram = true;
    }


    /** Plotting histogram by points with error bars */
    public void plotByPointsWithErrorBars(){
	plotHistogram = false;
    }


    /** draw atmospheric muon weighted-events */
    public void drawEventsWithAtmMuonWeights(String modelName){
	isWeighted = true;
	isAtmMuon = true;
	this.modelName = modelName;
    }

    /** draw GZK flux weighted-events */
    public void drawEventsWithGZKWeights(String modelName){
	isWeighted = true;
	isAtmMuon = false;
	this.modelName = modelName;
    }

    /** draw unweighted-events */
    public void drawEventsWithNoWeights(){
	isWeighted = false;
    }

    /** Sets the criteria on making Histogram */
    public void setCriteria(Criteria criteria){
	this.criteria = criteria;
    }


    /** 
	Making historgram for drawing by "grafig" graphics package on Xfig.
	You have to call this method before drawing any plots on xfig.
     */
    public void makeHistogram(){

	// initialization of histograms
	histLogNpeCosZenith =  new double[dimensionLogNpe][dimensionCosZenith];
	histLogEnergyCosZenith =  new double[dimensionLogEnergy][dimensionCosZenith];
	histLogNpeNDOM =  new double[dimensionLogNpe][dimensionNDOMs];
	histLogNpeLogEnergy =  new double[dimensionLogNpe][dimensionLogEnergy];
	histLogNpeFgQuality =  new double[dimensionLogNpe][dimensionFgQuality];
	histCosZenithFgQuality =  new double[dimensionCosZenith][dimensionFgQuality];
	for(int iNpe = 0; iNpe<dimensionLogNpe ; iNpe++){
	    for(int iZenith = 0; iZenith<dimensionCosZenith; iZenith++){
		histLogNpeCosZenith[iNpe][iZenith] = 0.0;
	    }
	    for(int nDOM= 0; nDOM<dimensionNDOMs; nDOM++){
		histLogNpeNDOM[iNpe][nDOM] = 0.0;
	    }
	    for(int iFg= 0; iFg<dimensionFgQuality; iFg++){
		histLogNpeFgQuality[iNpe][iFg] = 0.0;
	    }
	    for(int iE= 0; iE<dimensionLogEnergy; iE++){
		histLogNpeLogEnergy[iNpe][iE] = 0.0;
	    }
	}
	for(int iZenith = 0; iZenith<dimensionCosZenith; iZenith++){
	    for(int iFg= 0; iFg<dimensionFgQuality; iFg++){
		histCosZenithFgQuality[iZenith][iFg] = 0.0;
	    }
	    for(int iE= 0; iE<dimensionLogEnergy; iE++){
		histLogEnergyCosZenith[iE][iZenith] = 0.0;
	    }
	}

	// loop over I3Particle Objects
	i3particleIterator = i3particleList.listIterator();
	numberOfEvents = 0;
	while(i3particleIterator.hasNext()){

	    I3Particle iceParticle = (I3Particle)i3particleIterator.next();

	    if(isMCTruth) iceParticle.switchToMCTruth();
   	    else iceParticle.switchToReco();

	    numberOfEvents++;

	    if(criteria == null || criteria.doesThisEventPass(iceParticle)){
		System.err.print("- EventNumber is " + 
				 iceParticle.getIceCubeData().getEventNumber());
		System.err.println(" The I3Particle Name is " + 
			   iceParticle.particleName(iceParticle.getFlavor(), 
						    iceParticle.getDoublet()));

		J3UnitVector n = iceParticle.getDirectionInIceCubeCoordinate();
		double cosZenith = -n.getZ(); // Reversed vector

		double logEnergy = iceParticle.getLogEnergy();
		if(logEnergy<minLogEnergy) logEnergy = minLogEnergy;

		// Number Of PEs and hitted DOMs of this event
		double logNpe = iceParticle.getIceCubeData().getLogBestNpe();
		int nDOMs = iceParticle.getIceCubeData().getNDOMsLaunch();

		double fgQ = iceParticle.getFirstGuessQuality();

		double w = 1.0;
		if(isWeighted){
		    double primaryFluxWeight = iceParticle.getMCPrimarySpectrumWeight();
		    double modelFluxWeight;
		    if(isAtmMuon){
			modelFluxWeight = iceParticle.getAtmosphericMuonFlux(modelName);
		    }else{
			modelFluxWeight = iceParticle.getGZKNeutrinoFlux(modelName);
		    }
		    ListIterator numberIterator = numberOfAllI3ParticlesList.listIterator();
		    ListIterator numFilled = numberOfFilledI3ParticlesList.listIterator();
		    int numberOfAllI3Particles = 0;
		    while(numberIterator.hasNext()){
			numberOfAllI3Particles = ((Integer )numberIterator.next()).intValue();
			if(numberOfEvents <= 
			   ((Integer )numFilled.next()).intValue()) break;
		    }
		    w = (modelFluxWeight/primaryFluxWeight)*
			mcArea*mcOmega*observationTime/(double )numberOfAllI3Particles;
		    //System.err.println(" fluxWeight= " + modelFluxWeight + 
		    //		       " MCPrimaryWeight " + primaryFluxWeight +
		    //		       " Total Number of Evnets " + numberOfAllI3Particles +
		    //		       " w= " + w);
		    //System.err.println(" logE=" + logEnergy + " cos(zenith)= " + cosZenith);
		}

		histLogNpeCosZenith[(int )((logNpe-minLogNpe)/deltaLogNpe)][(int )((cosZenith+1.0)/deltaCosZenith)] += w;
		histLogEnergyCosZenith[(int )((logEnergy-minLogEnergy)/deltaLogEnergy)][(int )((cosZenith+1.0)/deltaCosZenith)] += w;
		histLogNpeLogEnergy[(int )((logNpe-minLogNpe)/deltaLogNpe)][(int )((logEnergy-minLogEnergy)/deltaLogEnergy)] += w;
		histLogNpeNDOM[(int )((logNpe-minLogNpe)/deltaLogNpe)][(int )(nDOMs/deltaNDOMs)] += w;
		histLogNpeFgQuality[(int )((logNpe-minLogNpe)/deltaLogNpe)][(int )(fgQ/deltaFgQuality)] += w;
		histCosZenithFgQuality[(int )((cosZenith+1.0)/deltaCosZenith)][(int )(fgQ/deltaFgQuality)] += w;

	    }
	}

    }

    private double drawOnXfig(double[][] histogram, 
			      double minX, double maxX, int dimensionX, double deltaX,
			      double minY, double maxY, int dimensionY, double deltaY,
			      String xLabel, String yLabel,
			      double markerSize, int sliceStep,
			      boolean plotXwithSliceOfY){

	double logNumberOfEvents = Math.log(numberOfEvents)/ln10;
	double maxLogNumberOfEvents = Math.ceil(logNumberOfEvents);

	System.out.println("titx " + xLabel);
	System.out.println("tity " + yLabel);
	System.out.println("gwin 0.2 0.9 0.2 0.9");
	if(drawByAutoScale){
	    if(plotXwithSliceOfY)
		System.out.println("scal " + minX + " " + maxX + " " + 
				   minLogNumberOfEvents + " " + 
				   maxLogNumberOfEvents);
	    else
		System.out.println("scal " + minY + " " + maxY + " " + 
				   minLogNumberOfEvents + " " + 
				   maxLogNumberOfEvents);
	}
	System.out.println("mssg Number of Events (" + numberOfEvents + ")");

	int dimensionSlice; int dimensionPlot; double deltaPlot;
	if(plotXwithSliceOfY){
	    dimensionSlice = dimensionY;
	    dimensionPlot = dimensionX;
	    deltaPlot = deltaX;
	}else{
	    dimensionSlice = dimensionX;
	    dimensionPlot = dimensionY;
	    deltaPlot = deltaY;
	}
	for(int iSlice = 0; iSlice<dimensionSlice; iSlice+=sliceStep){
	    int colorMark = dimensionSlice-iSlice-sliceStep;
	    if(plotHistogram){
		System.out.println("hicl " + colorMark);
		System.out.println("hfil -1");
		System.out.println("lnth 2");
	    }else{
		System.out.println("mksz " + markerSize);
		System.out.println("mkcl " + colorMark);
		System.out.println("lncl " + colorMark);
	    }
	    for(int iPlot = 0; iPlot<dimensionPlot ; iPlot++){
		double x; double nEvents;
		if(plotXwithSliceOfY){
		    x = minX + deltaPlot*(double )iPlot;
		    nEvents = histogram[iPlot][iSlice];
		}else{
		    x = minY + deltaPlot*(double )iPlot;
		    nEvents = histogram[iSlice][iPlot];
		}
		double logNEvents; 
		double logMaxNEvents;double logMinNEvents;
		if(nEvents>0.0){
		    logNEvents = Math.log(nEvents)/ln10;
		    logMaxNEvents = Math.log(nEvents+Math.sqrt(nEvents))/ln10;
		    double minEvents= nEvents-Math.sqrt(nEvents);
		    if(minEvents>0.0){
			logMinNEvents = Math.log(minEvents)/ln10;
		    }else{
			logMinNEvents = minLogNumberOfEvents;
		    }
		}else{ 
		    logNEvents = minLogNumberOfEvents; 
		    logMaxNEvents = logNEvents;logMinNEvents = logNEvents;
		}

		if(!plotHistogram && nEvents>0.0){
		    System.out.println("data " + x + " " + 
				       0.5*deltaPlot + " " + logNEvents + " 0.0");
		    System.out.println("line " + x + " " + logMinNEvents +
				       " " + x + " " + logMaxNEvents);
		}else if(plotHistogram){
		    System.out.println("data " + x + " 0.0 " + 
				       logNEvents + " 0.0");
		}
	    }

	    if(plotHistogram){
		System.out.println("hist");
	    }else{
		System.out.println("plot");
	    }
	    System.out.println("disp");
	    System.out.println("cont");
	}

	return maxLogNumberOfEvents;
    }



    private double drawOnXfig(double[][] histogram, 
			      double minX, double maxX, int dimensionX, double deltaX,
			      double minY, double maxY, int dimensionY, double deltaY,
			      String xLabel, String yLabel,
			      double markerSize, 
			      boolean plotXwithSliceOfY){

	return drawOnXfig(histogram, 
			  minX, maxX, dimensionX, deltaX,
			  minY, maxY, dimensionY, deltaY,
			  xLabel, yLabel,
			  markerSize, 1, plotXwithSliceOfY);
    }


    /** Draw Npe distribution with slices of cos(zenith) with its bin width */
    public double drawNpeDistributionOnXfig(){

	return drawOnXfig(histLogNpeCosZenith,
			  minLogNpe, maxLogNpe, dimensionLogNpe, deltaLogNpe,
			  -1.0, 1.0, dimensionCosZenith, deltaCosZenith,
			  "Log(Npe)","Log(Number of Events)",
			   0.7, true);
    }

    /** Draw Npe distribution with slices of Log(Energy) with its bin width */
    public double drawNpeDistributionWithSliceOfEnergyOnXfig(){

	return drawOnXfig(histLogNpeLogEnergy,
			  minLogNpe, maxLogNpe, dimensionLogNpe, deltaLogNpe,
			  minLogEnergy, maxLogEnergy, dimensionLogEnergy, deltaLogEnergy,
			  "Log(Npe)","Log(Number of Events)",
			  0.7, 5,true);
    }

    /** Draw Energy distribution with slices of cos(zenith) with its bin width */
    public double drawEnergyDistributionOnXfig(){

	return drawOnXfig(histLogEnergyCosZenith,
			  minLogEnergy, maxLogEnergy, dimensionLogEnergy, deltaLogEnergy,
			  -1.0, 1.0, dimensionCosZenith, deltaCosZenith,
			  "Log(Energy)","Log(Number of Events)",
			   0.7, true);
    }

    /** Draw cos(zenith) distribution with slices of Npe with its bin width */
    public double drawZenithAngleDistributionOnXfig(){

	return drawOnXfig(histLogNpeCosZenith,
			  minLogNpe, maxLogNpe, dimensionLogNpe, deltaLogNpe,
			  -1.0, 1.0, dimensionCosZenith, deltaCosZenith,
			  "cos(Zenith Angle)","Log(Number of Events)",
			   1.0, false);
    }


    /** Draw cos(zenith) distribution with slices of 
	First Guess quality with its bin width */
    public double drawZenithAngleDistributionWithFGsliceOnXfig(){

	return drawOnXfig(histCosZenithFgQuality,
			  -1.0, 1.0, dimensionCosZenith, deltaCosZenith,
			  0.0, maxFgQuality/PropagationMatrix.c, dimensionFgQuality, 
			  deltaFgQuality/PropagationMatrix.c,
			  "cos(Zenith Angle)","Log(Number of Events)",
			   1.0, true);
    }

    /** Draw NPE distribution with slices of NDOM with its bin width */
    public double drawNpeDistributionWithNDOMsliceOnXfig(){

	return drawOnXfig(histLogNpeNDOM,
			  minLogNpe, maxLogNpe, dimensionLogNpe, deltaLogNpe,
			  0.0, maxNDOMsOnDrawing, dimensionNDOMs, deltaNDOMs,
			  "Log(Npe)","Log(Number of Events)",
			   0.7, true);
    }


    /** Draw First Guess Quality distribution with slices of Npe */
    public double drawFirstGuessQualityDistributionOnXfig(){

	return drawOnXfig(histLogNpeFgQuality,
			  minLogNpe, maxLogNpe, dimensionLogNpe, deltaLogNpe,
			  0.0, maxFgQuality/PropagationMatrix.c, dimensionFgQuality, 
			  deltaFgQuality/PropagationMatrix.c,
			  "velocity [beta]","Log(Number of Events)",
			   1.0, false);
    }

    /** Make and Fill IHistogram1D by reading out the variables of I3Particles
	<pre>
	String option
	"logEnergy"      :  plot logE
	"logRecoEnergy"  :  plot logRecoEnergy
	"logNpe"         :  plot logNpe
	"cosZenith"      :  plot cos(Zenith Angle)
	"nDOMs"          :  plot numberOfDOMs
	"firstGuess"     :  plot beta - first guess quality
	boolean bootstrap 
	false : default. No bootstraping
	true  : resampling the (weighted) I3Particles 
                following the poisson distribution.
	</pre>
    */
    protected IHistogram1D makeJaida1DHistogram(String option, String histName,
						boolean bootstrap,
						IHistogramFactory jaidaHistoFactory){
	double minX = 0.0; double maxX = 0.0;
	int dimensionX = 0;
	if(option.startsWith("logE") || option.startsWith("logRecoE")){
	    minX = minLogEnergy;
	    maxX = maxLogEnergy;
	    dimensionX = dimensionLogEnergy;
	}else if(option.startsWith("logNp")){
	    minX = minLogNpe;
	    maxX = maxLogNpe;
	    dimensionX = dimensionLogNpe;
	}else if(option.startsWith("cosZe")){
	    minX = -1.0;
	    maxX = 1.0;
	    dimensionX = dimensionCosZenith;
	}else if(option.startsWith("nDOM")){
	    minX = 0.0;
	    maxX = maxNDOMsOnDrawing;
	    dimensionX = dimensionNDOMs;
	}else if(option.startsWith("firstGu")){
	    minX = 0.0;
	    maxX = maxFgQuality/PropagationMatrix.c;
	    dimensionX = dimensionFgQuality;
	}else{
	    System.err.println("Option (" + option + 
			       ") is wrong. Cannot make histogram!");
	    System.exit(0);
	}

	IHistogram1D h1 = 
	    jaidaHistoFactory.createHistogram1D(histName,dimensionX,minX,maxX);


	// Reading I3Particles
	// loop over I3Particle Objects
	i3particleIterator = i3particleList.listIterator();
	numberOfEvents = 0;
	while(i3particleIterator.hasNext()){

	    I3Particle iceParticle = (I3Particle)i3particleIterator.next();

	    if(isMCTruth) iceParticle.switchToMCTruth();
   	    else iceParticle.switchToReco();

	    numberOfEvents++;

	    if(criteria == null || criteria.doesThisEventPass(iceParticle)){

		double xval = 0.0;
		if(option.startsWith("logEnerg")){
		    xval = iceParticle.getLogEnergy();

		}else if(option.startsWith("logECR")){
		    if(muonBundle == null) muonBundle = new AtmMuonBundleFlux();
		    iceParticle.switchToMCTruth();// must be MC truth zenith
		    J3UnitVector n = iceParticle.getDirectionInIceCubeCoordinate();
		    double cosZenith = -n.getZ(); // Reversed vector
		    if(isMCTruth) iceParticle.switchToMCTruth();
		    else iceParticle.switchToReco();
		    double logInIceMuonEnergy = iceParticle.getLogEnergy();
		    double distance = 
			iceParticle.getDistanceFromEarthSurfaceToIceCube();
		    double slantDepth = distance*0.917; // [g/cm^2] for ice
		    double beta_loss = 4.4619776009127244E-6;
		    if(logInIceMuonEnergy >= 7.0) 
			beta_loss = iceCube.uhe.interactions.CELbeta.getBeta(logInIceMuonEnergy);
		    double cosmicRayEnergy = 
			muonBundle.getEffectiveEnergyOfCRs(logInIceMuonEnergy,
							   cosZenith,beta_loss,slantDepth);
		    double logCosmicRayEnergy = Math.log(cosmicRayEnergy)/ln10;
		    //System.err.println("cos(Zenith)=" + cosZenith + 
		    //       " distance=" + distance +
		    //       " [cm] logMuE=" + logInIceMuonEnergy + 
		    //       " logECR=" + logCosmicRayEnergy);
		    xval = logCosmicRayEnergy;

		}else if(option.startsWith("logRecoE")){
		    xval = iceParticle.getLogRecoEnergy();

		}else if(option.startsWith("logNp")){
		    double logNpe = iceParticle.getIceCubeData().getLogBestNpe();
		    if(logNpe>=maxLogNpe) logNpe = maxLogNpe;
		    xval = logNpe;

		}else if(option.startsWith("cosZe")){
		    J3UnitVector n = iceParticle.getDirectionInIceCubeCoordinate();
		    double cosZenith = -n.getZ(); // Reversed vector
		    xval = cosZenith;

		}else if(option.startsWith("nDOM")){
		    int nDOMs = iceParticle.getIceCubeData().getNDOMsLaunch();
		    xval = (double )nDOMs;

		}else if(option.startsWith("firstGu")){
		    double fgQ = iceParticle.getFirstGuessQuality()/PropagationMatrix.c;
		    xval = fgQ;
		}

		double w = 1.0;
		if(isWeighted){
		    double primaryFluxWeight = iceParticle.getMCPrimarySpectrumWeight();
		    double modelFluxWeight;
		    if(isAtmMuon){
			modelFluxWeight = iceParticle.getAtmosphericMuonFlux(modelName);
		    }else{
			modelFluxWeight = iceParticle.getGZKNeutrinoFlux(modelName);
		    }
		    ListIterator numberIterator = numberOfAllI3ParticlesList.listIterator();
		    ListIterator numFilled = numberOfFilledI3ParticlesList.listIterator();
		    int numberOfAllI3Particles = 0;
		    while(numberIterator.hasNext()){
			numberOfAllI3Particles = ((Integer )numberIterator.next()).intValue();
			if(numberOfEvents <= 
			   ((Integer )numFilled.next()).intValue()) break;
		    }
		    w = (modelFluxWeight/primaryFluxWeight)*
			mcArea*mcOmega*observationTime/(double )numberOfAllI3Particles;
		}

		// now bootsrappiing
		if(bootstrap){
		    if(randomGen == null) randomGen = new RandomGenerator();
		    long statW = randomGen.GetPoissonian(1.0);
		    w = w*(double )statW;
		}

		// fill 1D histogram
		h1.fill(xval,w);
	    } // loop for passing Criteria ends

	}// loop over I3Particles ends


	return h1;

    }

    /** Histogram name is assinged by the default  - same as "option" */
    protected IHistogram1D makeJaida1DHistogram(String option,
						IHistogramFactory jaidaHistoFactory){
	return makeJaida1DHistogram(option,option,false,jaidaHistoFactory);
    }

    /** Make and Fill IHistogram2D by reading out the variables of I3Particles
	<pre>
	String option
	"logE-Npe"        :  plot logE-logNpe
	"logECR-Npe"      :  plot log(Primary Cosmic Ray E)-logNpe
	"logECR-logE"     :  plot log(Primary Cosmic Ray E)-logE
	"logRecoE-Npe"    :  plot logRecoE-logNpe
	"logE-cosZenith"  :  plot logE-cosZenith
	"logNpe-cosZenith":  plot logNpe-cosZenith
	"logNpe-CobZ"     :  plot logNpe-ConterOfBrightnessZ
	"CobR-CobZ"       :  plot CenterOfBrightness(sqrt(x^2+y^2))-ConterOfBrightnessZ
	boolean bootstrap 
	false : default. No bootsraping
	true  : resampling the (weighted) I3Particles 
                following the poisson distribution.
	</pre>
    */
    protected IHistogram2D makeJaida2DHistogram(String option, String histName,
						boolean bootstrap,
						IHistogramFactory jaidaHistoFactory){
	double minX = 0.0; double maxX = 0.0;
	int dimensionX = 0;
	double minY = 0.0; double maxY = 0.0;
    	int dimensionY = 0;
	if(option.startsWith("logE-Npe")){
	    minX = minLogEnergy;
	    maxX = maxLogEnergy;
	    dimensionX = dimensionLogEnergy;
	    minY = minLogNpe;
	    maxY = maxLogNpe;
	    dimensionY = dimensionLogNpe;
	}else if(option.startsWith("logECR-Npe")){
	    minX = minLogEnergy;
	    maxX = maxLogEnergy;
	    dimensionX = dimensionLogEnergy;
	    minY = minLogNpe;
	    maxY = maxLogNpe;
	    dimensionY = dimensionLogNpe;
	    if(muonBundle == null) muonBundle = new AtmMuonBundleFlux();
	}else if(option.startsWith("logECR-logE")){
	    minX = minLogEnergy;
	    maxX = maxLogEnergy;
	    dimensionX = dimensionLogEnergy;
	    minY = minLogEnergy;
	    maxY = maxLogEnergy;
	    dimensionY = dimensionLogEnergy;
	    if(muonBundle == null) muonBundle = new AtmMuonBundleFlux();
	}else if(option.startsWith("logRecoE-Npe")){
	    minX = minLogEnergy;
	    maxX = maxLogEnergy;
	    dimensionX = dimensionLogEnergy;
	    minY = minLogNpe;
	    maxY = maxLogNpe;
	    dimensionY = dimensionLogNpe;
	}else if(option.startsWith("logE-CosZ")){
	    minX = minLogEnergy;
	    maxX = maxLogEnergy;
	    dimensionX = dimensionLogEnergy;
	    minY = -1.0;
	    maxY = 1.0;
	    dimensionY = dimensionCosZenith;
	}else if(option.startsWith("logNpe-CosZ")){
	    minX = minLogNpe;
	    maxX = maxLogNpe;
	    dimensionX = dimensionLogNpe;
	    minY = -1.0;
	    maxY = 1.0;
	    dimensionY = dimensionCosZenith;
	}else if(option.startsWith("logNpe-Cob")){
	    minX = minLogNpe;
	    maxX = maxLogNpe;
	    dimensionX = dimensionLogNpe;
	    minY = -7.5e4; // [cm]
	    maxY = 7.5e4;  // [cm]
	    dimensionY = 50;
	}else if(option.startsWith("CobR-CobZ")){
	    minX = 0.0; // [cm]
	    maxX = 1.0e5;  // [cm]
	    dimensionX = 100;
	    minY = -1.0e5; // [cm]
	    maxY = 1.0e5;  // [cm]
	    dimensionY = 200;
	}else{
	    System.err.println("Option (" + option + 
			       ") is wrong. Cannot make histogram!");
	    System.exit(0);
	}

	IHistogram2D h2 = 
	    jaidaHistoFactory.createHistogram2D(histName,dimensionX,minX,maxX,
						dimensionY,minY,maxY);


	// Reading I3Particles
	// loop over I3Particle Objects
	i3particleIterator = i3particleList.listIterator();
	numberOfEvents = 0;
	while(i3particleIterator.hasNext()){

	    I3Particle iceParticle = (I3Particle)i3particleIterator.next();

	    if(isMCTruth) iceParticle.switchToMCTruth();
   	    else iceParticle.switchToReco();

	    numberOfEvents++;

	    if(criteria == null || criteria.doesThisEventPass(iceParticle)){

		double logNpe = iceParticle.getIceCubeData().getLogBestNpe();
		if(logNpe>=maxLogNpe) logNpe = maxLogNpe;

		J3UnitVector n = iceParticle.getDirectionInIceCubeCoordinate();
		double cosZenith = -n.getZ(); // Reversed vector

		double xval = 0.0; double yval = 0.0;
		if(option.startsWith("logE-Np")){
		    xval = iceParticle.getLogEnergy();
		    yval = logNpe;
		}else if(option.startsWith("logECR")){
		    iceParticle.switchToMCTruth();// must be MC truth zenith
		    J3UnitVector nTruth = iceParticle.getDirectionInIceCubeCoordinate();
		    double cosZenithTruth = -nTruth.getZ(); // Reversed vector
		    if(isMCTruth) iceParticle.switchToMCTruth();
		    else iceParticle.switchToReco();
		    double logInIceMuonEnergy = iceParticle.getLogEnergy();
		    double distance = 
			iceParticle.getDistanceFromEarthSurfaceToIceCube();
		    double slantDepth = distance*0.917; // [g/cm^2] for ice
		    double beta_loss = 4.4619776009127244E-6;
		    if(logInIceMuonEnergy >= 7.0) 
			beta_loss = iceCube.uhe.interactions.CELbeta.getBeta(logInIceMuonEnergy);
		    double cosmicRayEnergy = 
			muonBundle.getEffectiveEnergyOfCRs(logInIceMuonEnergy,
							   cosZenithTruth,beta_loss,
							   slantDepth);
		    double logCosmicRayEnergy = Math.log(cosmicRayEnergy)/ln10;
		    xval = logCosmicRayEnergy;
		    if(option.startsWith("logECR-Npe")) yval = logNpe;
		    else yval = iceParticle.getLogEnergy();
		}else if(option.startsWith("logRecoE-Npe")){
		    xval = iceParticle.getLogRecoEnergy();
		    yval = logNpe;

		}else if(option.startsWith("logE-Cos")){
		    xval = iceParticle.getLogRecoEnergy();
		    yval = cosZenith;

		}else if(option.startsWith("logNpe-Cos")){
		    xval = logNpe;
		    yval = cosZenith;

		}else if(option.startsWith("logNpe-CobX")){
		    xval = logNpe;
		    yval = iceParticle.getR0InIceCubeCoordinate().getX();

		}else if(option.startsWith("logNpe-CobY")){
		    xval = logNpe;
		    yval = iceParticle.getR0InIceCubeCoordinate().getY();

		}else if(option.startsWith("logNpe-CobZ")){
		    xval = logNpe;
		    yval = iceParticle.getR0InIceCubeCoordinate().getZ();

		}else if(option.startsWith("CobR-CobZ")){
		    double cob_x = iceParticle.getR0InIceCubeCoordinate().getX();
		    double cob_y = iceParticle.getR0InIceCubeCoordinate().getY();
		    xval = Math.sqrt(cob_x*cob_x + cob_y*cob_y);
		    yval = iceParticle.getR0InIceCubeCoordinate().getZ();

		}

		double w = 1.0;
		if(isWeighted){
		    double primaryFluxWeight = iceParticle.getMCPrimarySpectrumWeight();
		    double modelFluxWeight;
		    if(isAtmMuon){
			modelFluxWeight = iceParticle.getAtmosphericMuonFlux(modelName);
		    }else{
			modelFluxWeight = iceParticle.getGZKNeutrinoFlux(modelName);
		    }
		    ListIterator numberIterator = numberOfAllI3ParticlesList.listIterator();
		    ListIterator numFilled = numberOfFilledI3ParticlesList.listIterator();
		    int numberOfAllI3Particles = 0;
		    while(numberIterator.hasNext()){
			numberOfAllI3Particles = ((Integer )numberIterator.next()).intValue();
			if(numberOfEvents <= 
			   ((Integer )numFilled.next()).intValue()) break;
		    }
		    w = (modelFluxWeight/primaryFluxWeight)*
			mcArea*mcOmega*observationTime/(double )numberOfAllI3Particles;
		}

		// now bootsrappiing
		if(bootstrap){
		    if(randomGen == null) randomGen = new RandomGenerator();
		    long statW = randomGen.GetPoissonian(1.0);
		    w = w*(double )statW;
		}

		// fill 1D histogram
		h2.fill(xval,yval,w);
	    } // loop for passing Criteria ends

	}// loop over I3Particles ends


	return h2;

    }

    /** Histogram name is assinged by the default  - same as "option" */
    protected IHistogram2D makeJaida2DHistogram(String option,
						IHistogramFactory jaidaHistoFactory){
	return makeJaida2DHistogram(option,option,false,jaidaHistoFactory);
    }


    /** Chi^2 comparison - analysis utility method */
    public static IComparisonResult calcChi2(IHistogram1D mcHistogram,
					     IHistogram1D dataHistogram,
					     int binFirst, int binLast){


	if(StatisticalComparison.canCompare(mcHistogram,dataHistogram,"chi2")){

	    int nbin = mcHistogram.axis().bins();
	    if(binFirst<0 || binLast>nbin){
		System.err.println("bin out of range");
		System.exit(0);
	    }
	    double chi2 = 0.0;
	    int nDof = 0;
	    for(int i= binFirst;i<binLast;i++){
		double nEventsData = dataHistogram.binHeight(i);
		double nEntriesData = dataHistogram.binEntries(i);
		double nEventsMC = mcHistogram.binHeight(i);
		double nEntriesMC = mcHistogram.binEntries(i);
		double nErrorMC = mcHistogram.binError(i);
		//double weight = nEventsData + nEntriesMC;
		double weight = nEventsData + nErrorMC*nErrorMC;
		//if(weight>0.0 && nEventsData > 0.0){
 	        if(weight>0.0){
		    double dw = nEventsData - nEventsMC;
		    chi2 += dw*dw/weight;
		    nDof++;
		    //System.err.println("x=" + dataHistogram.axis().binLowerEdge(i) +
		    //	       " NeventsData =" + nEventsData +
		    //	       " NentryData = " + nEntriesData +
		    //	       " NeventsMC =" + nEventsMC +
		    //	       " NentryMC = " + nEntriesMC);

		}
	    }

	    //IComparisonResult result =
	    //StatisticalComparison.compare(mcHistogram,dataHistogram,"chi2",
	    //			      "rejectionLevel=0.95");

	    ComparisonResult result = new ComparisonResult();
	    result.setQuality(chi2);
	    result.setnDof(nDof);
	    return result;
	}else{
	    System.out.println("Cannot compare the results");
	    return null;
	}
    }

    /** Chi^2 comparison - analysis utility method.
     [xFirst xMast] in the histograms participates the chi^2 calculation*/
    public static IComparisonResult calcChi2(IHistogram1D mcHistogram,
					     IHistogram1D dataHistogram,
					     double xFirst, double xLast){

	int binFirst = mcHistogram.coordToIndex(xFirst);
	int binLast = mcHistogram.coordToIndex(xLast);
	//System.err.println("binF= " + binFirst + " binL = " + binLast);

	return calcChi2(mcHistogram,dataHistogram,binFirst, binLast);
    }
    /** Chi^2 comparison - analysis utility method.
     All bins in the histograms participates the chi^2 calculation*/
    public static IComparisonResult calcChi2(IHistogram1D mcHistogram,
					     IHistogram1D dataHistogram ){

	return calcChi2(mcHistogram,dataHistogram,0,mcHistogram.axis().bins());
    }
 

}
