package iceCube.uhe.particles;

import iceCube.uhe.particles.*;
import iceCube.uhe.geometry.*;
import geometry.*;

import java.io.*;
import java.util.*;

/**
   This particle class inherited from Particle.java describes
   a digested IceCube event (both MC and real) in form of
   JULIeT particle. It has a geometry, flux weights (if MC),
   a propgation matix data, and digested data of
   the detector signals such as NPEs.
<P>
   If you want to get/set valuables on MC truth, then 
   run the methonds of this class object, p, as following.

   <pre>
   % p.switchToMCTruth();
   % J3UnitVector n = p.getDirectionInIceCubeCoordinate();
   </pre>

   for the "reconstruction" results, then

   <pre>
   % p.switchToReco();
   % J3UnitVector n = p.getDirectionInIceCubeCoordinate();
   </pre>

<P>
   This class can be serialized and, thus, can be stored
   in disk once you fill it with the data.

<P>
   Written originaly for the IceCube EHE analysis
   by S.Yoshida 2006/12/26

*/

public class I3Particle extends Particle implements Serializable {

    private static final long serialVersionUID = 974651822323863237L;
    private final double ln10 = Math.log(10.0);

    /** Geometry of propagating particle represended by EarthCenterCoordinate
	 - MC Truth -
     */
    J3Line particleAxis_J3Line_center_MC = null;
    /** Geometry of propagating particle represended by IceCubeCoordinate
	 - MC Truth -
     */
    J3Line particleAxis_J3Line_ice3_MC = null;

    /** Geometry of propagating particle represended by EarthCenterCoordinate
	 - Reconstruction track -
     */
    J3Line particleAxis_J3Line_center_reco = null;
    /** Geometry of propagating particle represended by IceCubeCoordinate
	 - Reconstruction track -
     */
    J3Line particleAxis_J3Line_ice3_reco = null;

    /** Energy estimated by reconstruction */
    private double energyReco;

    /** log(Energy) estimated by reconstruction */
    private double logEnergyReco;

    /** First guess quality index */
    private double fgQuality;

    /** propagation distance from the earth surface */
    private double distanceFromEarthSurface;

    /** I3 Data class */
    private I3Data ice3Data;

    /** Flux of GZK Neutrinos. Calculated by PropagatingNeutrinoFlux class
	in the neutrinoModel package
    */
    private Map dFdLogE_GZKMap;

    /** Flux of GZK Neutrinos. Calculated by PropagatingAtmMuonFlux class
	in the muonModel package
    */
    private Map dFdLogE_AtmMuonMap;

    /** MC primary spectrum weight */
    private double primarySpecWeight; 


    /** This variable gives two cases -
	<pre>
	true - sets/returns the parameters on MC truth
        false - sets/returns the parameters on reco results
        </pre>
     */
    private boolean isMCTruth = false;


    /** Constructor.
	   <pre>
	   flabor ... flavor valuable
	   doblet ... doublet valuable
	   energy ... initial Energy [GeV]
	   </pre>
    */
    public I3Particle(int flavor, int doublet, double energy){
	super(flavor,doublet,energy);
        dFdLogE_GZKMap = new LinkedHashMap();
        dFdLogE_AtmMuonMap = new LinkedHashMap();
	energyReco = -1.0; logEnergyReco = Double.NEGATIVE_INFINITY;
	fgQuality = -1.0;
	distanceFromEarthSurface = -1.0;
	ice3Data = new I3Data();
    }

    /** Constructor.
	   <pre>
	   flabor ... flavor valuable
	   doblet ... doublet valuable
	   </pre>
    */
    public I3Particle(int flavor, int doublet){
	super(flavor,doublet);
        dFdLogE_GZKMap = new LinkedHashMap();
        dFdLogE_AtmMuonMap = new LinkedHashMap();
	energyReco = -1.0; logEnergyReco = Double.NEGATIVE_INFINITY;
	fgQuality = -1.0;
	distanceFromEarthSurface = -1.0;
	ice3Data = new I3Data();
    }


    /** Switch to parameters concerned with MCtrue.*/
    public void switchToMCTruth(){
	isMCTruth = true;
    }

    /** Switch to parameters concerned with Reco results.*/
    public void switchToReco(){
	isMCTruth = false;
    }

    /** Tells if this object returns valuables on MC Truth.
	If not, it returns valuables  resulted from reconstruction */
    public boolean isMCTruth(){
	return isMCTruth;
    }

    /** put the energy estimated by a reconstruction program */
    public void putRecoEnergy(double energy){
	energyReco =  energy;
	if(energyReco>0.0) logEnergyReco = Math.log(energy)/ln10;
	else logEnergyReco = Double.NEGATIVE_INFINITY;
    }

    /** get the energy estimated by a reconstruction program */
    public double getRecoEnergy(){ return energyReco;}

    /** get the energy estimated by a reconstruction program */
    public double getLogRecoEnergy(){ return logEnergyReco;}

    /** put the First Guess fit quality parameter such as the velocity*/
    public void setFirstGuessQuality(double quality){
	fgQuality = quality;
    }

    /** Set the distance from the Earth Surface that this particle
        has propagated before reching to the IceCube volume
    */
    public void setDistanceFromEarthSurfaceToIceCube(double distance){
	distanceFromEarthSurface = distance;
    }

    /** Get the distance from the Earth Surface that this particle
        has propagated before reching to the IceCube volume
    */
    public double getDistanceFromEarthSurfaceToIceCube(){
	return distanceFromEarthSurface;
    }

    /** Returns the First Guess fit quality parameter such as the velocity*/
    public double getFirstGuessQuality(){ return fgQuality;}

    /** Sets the particle axis defined in the IceCube coordinate */
    public void setParticleAxisInIceCubeCoordinate(J3Line axis){
	if(isMCTruth) particleAxis_J3Line_ice3_MC = axis;
	else particleAxis_J3Line_ice3_reco = axis;
    }

    /** Sets the particle axis defined in the Earth Center coordinate */
    public void setParticleAxisInEarthCenterCoordinate(J3Line axis){
	if(isMCTruth) particleAxis_J3Line_center_MC = axis;
	else particleAxis_J3Line_center_reco = axis;
    }


    /** Returns the unit vector of the axis direction in 
	the IceCube coordinate */
    public J3UnitVector getDirectionInIceCubeCoordinate(){
	if(isMCTruth && particleAxis_J3Line_ice3_MC!=null){
	    return particleAxis_J3Line_ice3_MC.getDirection();
	}else if((!isMCTruth) && particleAxis_J3Line_ice3_reco!=null){
	    return particleAxis_J3Line_ice3_reco.getDirection();
	}else{
	    System.err.println("I3Particle: particle axis vector is null");
	    System.exit(0);
	    return null;
	}
    }

    /** Returns the unit vector of the axis direction in 
	the Earth Center coordinate */
    public J3UnitVector getDirectionInEarthCenterCoordinate(){
	if(isMCTruth && particleAxis_J3Line_center_MC!=null){
	    return particleAxis_J3Line_center_MC.getDirection();
	}else if((!isMCTruth) && particleAxis_J3Line_center_reco!=null){
	    return particleAxis_J3Line_center_reco.getDirection();
	}else{
	    System.err.println("I3Particle: particle axis vector is null");
	    System.exit(0);
	    return null;
	}
    }


    /** Returns the axis location (possibly, center of brightness) in 
	the Earth Center coordinate */
    public J3Vector getR0InIceCubeCoordinate(){
	if(isMCTruth && particleAxis_J3Line_ice3_MC!=null){
	    return particleAxis_J3Line_ice3_MC.getR0();
	}else if((!isMCTruth) && particleAxis_J3Line_ice3_reco!=null){
	    return particleAxis_J3Line_ice3_reco.getR0();
	}else{
	    System.err.println("I3Particle: particle axis vector is null");
	    System.exit(0);
	    return null;
	}
    }


    /** Returns the axis location (possibly, center of brightness) in 
	the IceCube coordinate */
    public J3Vector getR0InEarthCenterCoordinate(){
	if(isMCTruth && particleAxis_J3Line_center_MC!=null){
	    return particleAxis_J3Line_center_MC.getR0();
	}else if((!isMCTruth) && particleAxis_J3Line_center_reco!=null){
	    return particleAxis_J3Line_center_reco.getR0();
	}else{
	    System.err.println("I3Particle: particle axis vector is null");
	    System.exit(0);
	    return null;
	}
    }

    /** Return the trajectory of the particle in form of J3Line */
    public J3Line getParticleAxis(){
	if(isMCTruth && particleAxis_J3Line_center_MC!=null){
	    return particleAxis_J3Line_center_MC;
	}else if((!isMCTruth) && particleAxis_J3Line_center_reco!=null){
	    return particleAxis_J3Line_center_reco;
	}else{
	    System.err.println("I3Particle: particle axis vector is null");
	    System.exit(0);
	    return null;
	}
    }


    /** Transofrm the particle axis vector from the IceCube coordinate
	to the Earth Cencer coordinate */
    public void transformParticleAxisToEarthCenterCoordinate(){
	J3Line particleAxis_J3Line_ice3;
	J3Line particleAxis_J3Line_center = null;
	if(isMCTruth){
	    particleAxis_J3Line_ice3 = particleAxis_J3Line_ice3_MC;
	}else{
	    particleAxis_J3Line_ice3 = particleAxis_J3Line_ice3_reco;
	}

	if(particleAxis_J3Line_ice3 != null){
	    IceCubeCoordinate ice3Coordinate = new IceCubeCoordinate();
	    J3Vector passingPoint_J3Vector_center = 
		ice3Coordinate.transformVectorToEarthCenter(
                  particleAxis_J3Line_ice3.getR0());
	    J3UnitVector direction_J3UnitVector_center =
		ice3Coordinate.transformUnitVectorToEarthCenter(
                  particleAxis_J3Line_ice3.getDirection());
	    particleAxis_J3Line_center = new J3Line(passingPoint_J3Vector_center,
						    direction_J3UnitVector_center);
	}else{
	    System.err.println("I3Particle: particle axis vector is null");
	}

	if(isMCTruth){
	    particleAxis_J3Line_center_MC = particleAxis_J3Line_center;
	}else{
	    particleAxis_J3Line_center_reco = particleAxis_J3Line_center;
	}
    }

    /** Transofrm the particle axis vector from the Earth Center coordinate
	to the IceCube coordinate */
    public void transformParticleAxisToIceCubeCoordinate(){
	J3Line particleAxis_J3Line_ice3 = null;
	J3Line particleAxis_J3Line_center;
	if(isMCTruth){
	    particleAxis_J3Line_center = particleAxis_J3Line_center_MC;
	}else{
	    particleAxis_J3Line_center = particleAxis_J3Line_center_reco;
	}

	if(particleAxis_J3Line_center != null){
	    IceCubeCoordinate ice3Coordinate = new IceCubeCoordinate();
	    EarthCenterCoordinate earthCoordinate = new EarthCenterCoordinate();
	    J3Vector particleLocation_J3Vector_ice3 =
		ice3Coordinate.transformVectorToThisCoordinate(
                  particleAxis_J3Line_center.getR0(),earthCoordinate);

	    J3Vector direction_ice3 =
		ice3Coordinate.transformVectorToThisCoordinate(
                  particleAxis_J3Line_center.getDirection(),earthCoordinate);

	    J3UnitVector directionUnitVector_ice3 =
		new J3UnitVector(direction_ice3.getX(),
				 direction_ice3.getY(),direction_ice3.getZ());
	    particleAxis_J3Line_ice3 =
		new J3Line(particleLocation_J3Vector_ice3,
			   directionUnitVector_ice3);
	}else{
	    System.err.println("I3Particle: particle axis vector is null");
	}

	if(isMCTruth){
	    particleAxis_J3Line_ice3_MC = particleAxis_J3Line_ice3;
	}else{
	    particleAxis_J3Line_ice3_reco = particleAxis_J3Line_ice3;
	}
    }

    /** 
	Sets the GZK neutrino flux dF/dLogE [/cm^2 sec sr] for weight.
	Calculated by PropagatingNeutrinoFlux class
	in the neutrinoModel package.
        <pre>
	flux    :   dF/dLogE [/cm^2 sec sr]
        fluxName:   model name like "YT m=4 Zmax=4"
	</pre>
    */
    public void setGZKNeutrinoFlux(double flux, String fluxName){
	dFdLogE_GZKMap.put(fluxName,new Double(flux));
    }

    /** 
	Sets the Atmospheric Muon flux dF/dLogE [/cm^2 sec sr] for weight.
	Calculated by PropagatingAtmMuonFlux class
	in the muonModel package.
        <pre>
	flux    :   dF/dLogE [/cm^2 sec sr]
        fluxName:   model name like "Elbert bundle model"
	</pre>
    */
    public void setAtmosphericMuonFlux(double flux, String fluxName){
	dFdLogE_AtmMuonMap.put(fluxName,new Double(flux));
    }

    /** 
	Return the GZK neutrino flux dF/dLogE [/cm^2 sec sr]
	stored by String "fluxName".
        If it finds no data stored, return -1.0 with the warning message.
    */
    public double getGZKNeutrinoFlux(String fluxName){
	if(dFdLogE_GZKMap.containsKey(fluxName)){
	    Double flux_obj = (Double )(dFdLogE_GZKMap.get(fluxName));
	    double flux = flux_obj.doubleValue();
	    return flux;
	}else{
	    System.err.println(" No GZK neutrino flux data stored by " +
			       fluxName);
	    return -1.0;
	}
    }

    /** Remove the GZK neutrino flux data stored by String "fluxName" */
    public void removeGZKNeutrinoFlux(String fluxName){
	if(dFdLogE_GZKMap.containsKey(fluxName)){
	    dFdLogE_GZKMap.remove(fluxName);
	}
    }

    /** 
	Return the Atmospheric Muon flux dF/dLogE [/cm^2 sec sr]
	stored by String "fluxName".
        If it finds no data stored, return -1.0 with the warning message.
    */
    public double getAtmosphericMuonFlux(String fluxName){
	if(dFdLogE_AtmMuonMap.containsKey(fluxName)){
	    Double flux_obj = (Double )(dFdLogE_AtmMuonMap.get(fluxName));
	    double flux = flux_obj.doubleValue();
	    return flux;
	}else{
	    System.err.println(" No Atmospheric Muon flux data stored by " +
			       fluxName);
	    return -1.0;
	}
    }

    /** Remove the Atmospheric Muon flux data stored by String "fluxName" */
    public void removeAtmosphericMuonFlux(String fluxName){
	if(dFdLogE_AtmMuonMap.containsKey(fluxName)){
	    dFdLogE_AtmMuonMap.remove(fluxName);
	}
    }

    /**
       Return iterator to the map of GZK neutrino fluxes.
       You can use this iterator to loop over all
       fluxes data stored by setGZKNeutrinoFLux(flux, fluxname).
       <pre>
       A example :
	Iterator gzkIterator = iceParticle.iteratorOfGZKNeutrinoFlux();
	while(gzkIterator.hasNext()){
	    Map.Entry entry = (Map.Entry )(gzkIterator.next());
	    Double flux = (Double )(entry.getValue());
	    String name = (String )(entry.getKey());
        }
	</pre>
    */
    public Iterator iteratorOfGZKNeutrinoFlux(){
	return dFdLogE_GZKMap.entrySet().iterator();
    }

    /**
       Return iterator to the map of GZK neutrino fluxes.
       You can use this iterator to loop over all
       fluxes data stored by setAtmosphericMuonFLux(flux, fluxname)
    */
    public Iterator iteratorOfAtmosphericMuonFlux(){
	return dFdLogE_AtmMuonMap.entrySet().iterator();
    }


    /** 
	Copy the given propagation matrix to the energy distribution matrix
        dN/dLogEin.
        The propagation matrix calculated 
	by PropagationMatrix.java in the propagation package is in form
	of something like FnuEToNuE(logEin,logEout). logEout here should
        be equal to THIS particle energy while logEin is that of parimary
        particle entering into the earth surface. This method automatically
        builds dN/dLogEin from the propagation matrix data and store it
       in the particle class's energy distribution matrix.
       <pre>
       double[][] fNuToThis  :   A Propagation Matrix
       </pre>
    */
    public void copyPropagationMatrixData(double[][] fNuToThis){
	int dimension = getDimensionOfLogEnergyMatrix();
	int jLogE = (int)((getLogEnergy() - getLogEnergyMinimum())
			  /getDeltaLogEnergy());
	for(int iLogE=jLogE;iLogE<dimension;iLogE++){
	    double logE = getLogEnergyMinimum() + 
		getDeltaLogEnergy()*(double )iLogE;
	    putLogEnergyMatrix(logE,fNuToThis[iLogE][jLogE]);
	}
    }

    /** Set the MC primary spectrum weight
	You have to switch to MCTruth (by calling switchToMCTruth())
	to put this weight, due to the safety reasons.
    */
    public void setMCPrimarySpectrumWeight(double weight){
	if(isMCTruth) primarySpecWeight = weight;
    }

    /** Return the MC primary spectrum weight */
    public double getMCPrimarySpectrumWeight() {return primarySpecWeight;}


    /** The method to access inner IceCube data class */
    public I3Data getIceCubeData() { return ice3Data;}

    /**
       IceCube data class.
    */
    public class I3Data implements Serializable {
	private static final long serialVersionUID = 8129368670320780033L;
	int eventNumber = 0;
	double npeFADC = 0.0; double npeATWD = 0.0; double npeBest = 0.0;
	int nDOMsFADC = 0; int nDOMsATWD = 0; int nDOMsLaunch = 0;
	double logNpeFADC = Double.NEGATIVE_INFINITY; 
	double logNpeATWD = Double.NEGATIVE_INFINITY; 
	double logNpeBest = Double.NEGATIVE_INFINITY; 

	public int getEventNumber() {return eventNumber;}
	public void setEventNumber(int number) {eventNumber = number;}
	public double getNpeATWD() {return npeATWD;}
	public double getLogNpeATWD() {return logNpeATWD;}
	public void setNpeATWD(double npe) {
	    npeATWD = npe;
	    if(npeATWD>0.0) logNpeATWD = Math.log(npeATWD)/ln10;
	}
	public double getNpeFADC() {return npeFADC;}
	public double getLogNpeFADC() {return logNpeFADC;}
	public void setNpeFADC(double npe) {
	    npeFADC = npe;
	    if(npeFADC>0.0) logNpeFADC = Math.log(npeFADC)/ln10;
	}
	public double getBestNpe() {return npeBest;}
	public double getLogBestNpe() {return logNpeBest;}
	public void setBestNpe(double npe) {
	    npeBest = npe;
	    if(npeBest>0.0) logNpeBest = Math.log(npeBest)/ln10;
	}
	public int getNDOMsATWD() { return nDOMsATWD;}
	public void setNDOMsATWD(int nDOMs) { nDOMsATWD = nDOMs;}
	public int getNDOMsFADC() { return nDOMsFADC;}
	public void setNDOMsFADC(int nDOMs) { nDOMsFADC = nDOMs;}
	public int getNDOMsLaunch() { return nDOMsLaunch;}
	public void setNDOMsLaunch(int nDOMs) { nDOMsLaunch = nDOMs;}
    }

}
