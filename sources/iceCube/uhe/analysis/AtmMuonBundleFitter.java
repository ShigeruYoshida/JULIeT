package iceCube.uhe.analysis;

import iceCube.uhe.particles.*;
import iceCube.uhe.points.*;
import iceCube.uhe.interactions.*;
import iceCube.uhe.muonModel.*;
import iceCube.uhe.geometry.*;
import iceCube.uhe.analysis.*;
import geometry.*;

import hep.aida.*;
import hep.aida.ext.*;
import hep.aida.util.*;
import hep.aida.util.comparison.*;

import java.io.*;
import java.util.*;

/** 
    This class provides methods to fit I3Particles with 
    the CEL approximated AtmMuonBundle flux or
    precalculated fluxes by PropagatingAtmMuonFlux.

    Use AtmMuonBundleFlux.java in the muonModel package
    to calculate the energy spectrum of muon bundles
    at the IceCube depth based on the CEL approximation, i.e.,
    AtmMuonBundleFlux.getDFDLogE(logEnergy,cosTheta,beta_loss,slantDepth).

    An andvanced method 
    <pre>
    public static void fillAtmMuonBundleFluxWeight(I3Particle iceParticle, 
					   ElbertFluxTableFactory muonFluxTable)
    </pre>
    is also available for using caluclation based upon PropagatingAtmMuonFlux.

    See the API document of the JULIeT's muonModel package for the details.

    Written by S. Yoshida 2007 February 16
    Revised by S. Yoshida 2008 June 2nd
*/

public class AtmMuonBundleFitter {

    protected static String fluxWeightName = "Elbert CEL Model";
    protected static double muonBeta = 4.4619776009127244E-6; 
    // Muon Inelasticity [/g/cm^2] E=10^7 GeV


    /** 
	Fill I3Particle with the Atmospheric Muon Bundle Flux weights.

	The flux at the IceCube depth is given by the CEL approximated
        calculation in getDFDLogE(logE,cos(theta),beta,depth) in AtmMuonBundleFlux.
    */
    public static void fillAtmMuonBundleFluxWeight(I3Particle iceParticle, 
						   ParticlePoint s,
						   AtmMuonBundleFlux muonFlux){
	// All the parameters should be MC truth 
	iceParticle.switchToMCTruth();

	J3UnitVector n_ice = iceParticle.getDirectionInIceCubeCoordinate();
	double cosZenith = -n_ice.getZ(); // Reversed vector

	// Get the muon bundle flux at the IceCube depth
	double slantDepth = 
	    iceParticle.getDistanceFromEarthSurfaceToIceCube()*s.getMediumDensity();
	double logEnergy = iceParticle.getLogEnergy();

	double beta = muonBeta;
	if(logEnergy >= 7.0) beta = CELbeta.getBeta(logEnergy);
	// energy dependence of beta must be took into account for logE > 7.0
	// for the better matching with the propagation matrix-based calculation.

	double flux = muonFlux.getDFDLogE(logEnergy,cosZenith,
					  beta,slantDepth);

	// The primary cosmic ray energy responsible for this I3Particle
	double cosmicRayEnergy = 
	    muonFlux.getEffectiveEnergyOfCRs(logEnergy,cosZenith,
					     beta,slantDepth);
	// Fill the flux weight
	iceParticle.putRecoEnergy(cosmicRayEnergy);
	iceParticle.removeAtmosphericMuonFlux(fluxWeightName);
	                             // Should be overwritten
	iceParticle.setAtmosphericMuonFlux(flux,fluxWeightName);
    }


    /**
       Set the Atmospheric Muon Bundle flux to all the I3Particles
       stored in the I3ParticleAnalysisFactory.

       The set flux is subject to be compared with the real data distribution.
     */
    public static void setAtmMuonBundleFlux(I3ParticleAnalysisFactory analizer,
					    ParticlePoint s,
					    AtmMuonBundleFlux muonFlux){

	ListIterator i3particleIterator =analizer.getParticleIterator();
	while(i3particleIterator.hasNext()){
	    I3Particle iceParticle = (I3Particle)i3particleIterator.next();
	    fillAtmMuonBundleFluxWeight(iceParticle,s,muonFlux);
	}
    }


    /** 
	Fill I3Particle with the Atmospheric Muon Bundle Flux weights.

	The flux at the IceCube depth is given by the calculation
	based on PropagatingAtmMuonBundleFlux, where muon propagation
        is processed by the full numerical propagation matrix 
	without any CEL approximation. Because PropagatingAtmMuonBundleFlux
	must read our pre-calculataed matrix data which takes time,
	this method uses ElbertFluxTableFactory object in the muonModel package
	which reads compact table data files containing the calculatated fluxes
	values by PropagatingAtmMuonBundleFlux.
	
    */
    public static void fillAtmMuonBundleFluxWeight(I3Particle iceParticle, 
					   ElbertFluxTableFactory muonFluxTable){
	// All the parameters should be MC truth 
	iceParticle.switchToMCTruth();
	// Propagation Distance [cm]
	double distance = iceParticle.getDistanceFromEarthSurfaceToIceCube();
	// log(Muon Energy [GeV])
	double logEnergy = iceParticle.getLogEnergy();

	double flux = muonFluxTable.getDFMuDLogE(distance,logEnergy);

	// Fill the flux weight
	iceParticle.removeAtmosphericMuonFlux(fluxWeightName);
	                             // Should be overwritten
	iceParticle.setAtmosphericMuonFlux(flux,fluxWeightName);
    }

    /**
       Set the Atmospheric Muon Bundle flux to all the I3Particles
       stored in the I3ParticleAnalysisFactory.

       The set flux is subject to be compared with the real data distribution.
     */
    public static void setAtmMuonBundleFlux(I3ParticleAnalysisFactory analizer,
					   ElbertFluxTableFactory muonFluxTable){

	ListIterator i3particleIterator =analizer.getParticleIterator();
	while(i3particleIterator.hasNext()){
	    I3Particle iceParticle = (I3Particle)i3particleIterator.next();
	    fillAtmMuonBundleFluxWeight(iceParticle,muonFluxTable);
	}
    }

    public static void fillAtmMuonBundleFluxWeight(I3Particle iceParticle, 
						   String modelName,
						   RelativeElbertFluxTableMaker muonFluxTable) throws IOException {
	// All the parameters should be MC truth 
	iceParticle.switchToMCTruth();

	// log(Muon Energy [GeV])
	double logEnergy = iceParticle.getLogEnergy();

	J3UnitVector n_ice = iceParticle.getDirectionInIceCubeCoordinate();
	double cosZenith = -n_ice.getZ(); // Reversed vector

	// The Elbert Flux alerady weighted by the default parameters
	double muonFlux = iceParticle.getAtmosphericMuonFlux(modelName);

	double flux = muonFluxTable.getRelativeFlux(logEnergy,cosZenith)*muonFlux;

	// Fill the flux weight
	iceParticle.removeAtmosphericMuonFlux(fluxWeightName);
	                             // Should be overwritten
	iceParticle.setAtmosphericMuonFlux(flux,fluxWeightName);
    }

    public static void setAtmMuonBundleFlux(I3ParticleAnalysisFactory analizer,
					    String modelName,
					    RelativeElbertFluxTableMaker muonFluxTable) throws IOException{

	ListIterator i3particleIterator =analizer.getParticleIterator();
	while(i3particleIterator.hasNext()){
	    I3Particle iceParticle = (I3Particle)i3particleIterator.next();
	    fillAtmMuonBundleFluxWeight(iceParticle,modelName,muonFluxTable);
	}
    }

    /** This is for the debugging. */
    protected static 
	void fillAtmMuonBundleFluxRatioWeight(I3Particle iceParticle, 
					      String modelName,
					      ElbertFluxTableFactory muonFluxTable,
					      RelativeElbertFluxTableMaker relativeFluxTable) 
	throws IOException {

	// All the parameters should be MC truth 
	iceParticle.switchToMCTruth();

	// log(Muon Energy [GeV])
	double logEnergy = iceParticle.getLogEnergy();

	J3UnitVector n_ice = iceParticle.getDirectionInIceCubeCoordinate();
	double cosZenith = -n_ice.getZ(); // Reversed vector

	// The Elbert Flux alerady weighted by the default parameters
	double muonFlux = iceParticle.getAtmosphericMuonFlux(modelName);

	double flux = relativeFluxTable.getRelativeFlux(logEnergy,cosZenith)*muonFlux;

	// Propagation Distance [cm]
	double distance = iceParticle.getDistanceFromEarthSurfaceToIceCube();
	double elbertFlux = muonFluxTable.getDFMuDLogE(distance,logEnergy);

	// take a ratio
	double ratio = -1.0;
	if(elbertFlux>0.0) ratio = flux/elbertFlux;

	// Fill the flux weight
	iceParticle.putRecoEnergy(ratio);
    }

    /** This is for the debugging. */
    protected static void setAtmMuonBundleFluxRatio(I3ParticleAnalysisFactory analizer,
						    String modelName,
						    ElbertFluxTableFactory muonFluxTable,
						    RelativeElbertFluxTableMaker relativeFluxTable) throws IOException{

	ListIterator i3particleIterator =analizer.getParticleIterator();
	while(i3particleIterator.hasNext()){
	    I3Particle iceParticle = (I3Particle)i3particleIterator.next();
	    fillAtmMuonBundleFluxRatioWeight(iceParticle,modelName,muonFluxTable,relativeFluxTable);
	}
    }
    /**
       Make a chi2-based statistical comparison between
       the MC with the muon bundle flux and the real data
     */
    public static IComparisonResult compare(I3ParticleAnalysisFactory mcAnalizer,
					    I3ParticleAnalysisFactory dataAnalizer,
					    IHistogramFactory histoFactory){

	mcAnalizer.switchToMCTruth();
	mcAnalizer.drawEventsWithAtmMuonWeights(fluxWeightName);
	IHistogram1D h1LogNpeMC = 
	    mcAnalizer.makeJaida1DHistogram("logNpe",histoFactory);
	IHistogram1D h1LogNpeData = 
	    dataAnalizer.makeJaida1DHistogram("logNpe",histoFactory);


	if(StatisticalComparison.canCompare(h1LogNpeMC,h1LogNpeData,"chi2")){

	    int nbin = h1LogNpeMC.axis().bins();
	    double chi2 = 0.0;
	    int nDof = 0;
	    for(int i= 0;i<nbin;i++){
		double nEventsData = h1LogNpeData.binHeight(i);
		double nEntriesData = h1LogNpeData.binEntries(i);
		double nEventsMC = h1LogNpeMC.binHeight(i);
		double nEntriesMC = h1LogNpeMC.binEntries(i);
		double nErrorMC = h1LogNpeMC.binError(i);
		//double weight = nEventsData + nEntriesMC;
		double weight = nEventsData + nErrorMC*nErrorMC;
		if(weight>0.0 && nEventsData > 0.0){
		    double dw = nEventsData - nEventsMC;
		    chi2 += dw*dw/weight;
		    nDof++;
		    //System.err.println("x=" + h1LogNpeData.axis().binLowerEdge(i) +
		    //	       " NeventsData =" + nEventsData +
		    //	       " NentryData = " + nEntriesData +
		    //	       " NeventsMC =" + nEventsMC +
		    //	       " NentryMC = " + nEntriesMC);

		}
	    }

	    //IComparisonResult result =
	    //StatisticalComparison.compare(h1LogNpeMC,h1LogNpeData,"chi2",
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

    /**
       Make a chi2-based statistical comparison between
       the two MCs (one with E**-1, another with E**-2, for example) with the muon bundle flux 
       and the real data
     */
    public static IComparisonResult compare(I3ParticleAnalysisFactory mcAnalizer1,
					    I3ParticleAnalysisFactory mcAnalizer2,
					    I3ParticleAnalysisFactory dataAnalizer,
					    IHistogramFactory histoFactory){

	mcAnalizer1.switchToMCTruth();
	mcAnalizer1.drawEventsWithAtmMuonWeights(fluxWeightName);

	mcAnalizer2.switchToMCTruth();
	mcAnalizer2.drawEventsWithAtmMuonWeights(fluxWeightName);

	IHistogram1D h1LogNpeMC1 = 
	    mcAnalizer1.makeJaida1DHistogram("logNpe",histoFactory);
	IHistogram1D h1LogNpeMC2 = 
	    mcAnalizer2.makeJaida1DHistogram("logNpe",histoFactory);
	IHistogram1D h1LogNpeData = 
	    dataAnalizer.makeJaida1DHistogram("logNpe",histoFactory);


	if(StatisticalComparison.canCompare(h1LogNpeMC1,h1LogNpeData,"chi2") &&
	   StatisticalComparison.canCompare(h1LogNpeMC2,h1LogNpeData,"chi2")){

	    int nbin = h1LogNpeMC1.axis().bins();
	    double chi2 = 0.0;
	    int nDof = 0;
	    for(int i= 0;i<nbin;i++){
		double nEventsData = h1LogNpeData.binHeight(i);
		double nEntriesData = h1LogNpeData.binEntries(i);
		double nEventsMC = h1LogNpeMC1.binHeight(i) + h1LogNpeMC2.binHeight(i);
		double nEntriesMC = h1LogNpeMC1.binEntries(i) + h1LogNpeMC2.binEntries(i);
		double nErrorMC = h1LogNpeMC1.binError(i) + h1LogNpeMC2.binError(i);
		//double weight = nEventsData + nEntriesMC;
		double weight = nEventsData + nErrorMC*nErrorMC;
		if(weight>0.0 && nEventsData > 0.0){
		    double dw = nEventsData - nEventsMC;
		    chi2 += dw*dw/weight;
		    nDof++;
		}
	    }

	    ComparisonResult result = new ComparisonResult();
	    result.setQuality(chi2);
	    result.setnDof(nDof);
	    return result;
	}else{
	    System.out.println("Cannot compare the results");
	    return null;
	}
    }

    /**
       Make a chi2-based statistical comparison between
       the two MCs (one with E**-1, another with E**-2, for example) 
       with the muon bundle flux  and the real data.

       The chi2 is built with mcAnalizer2 for logNpe < logNpeBoundary
       and with mcAnalizer1 otherwise.
    */
    public static IComparisonResult compare(I3ParticleAnalysisFactory mcAnalizer1,
					    I3ParticleAnalysisFactory mcAnalizer2,
					    I3ParticleAnalysisFactory dataAnalizer,
					    double logNpeBoundary,
					    IHistogramFactory histoFactory){

	mcAnalizer1.switchToMCTruth();
	mcAnalizer1.drawEventsWithAtmMuonWeights(fluxWeightName);

	mcAnalizer2.switchToMCTruth();
	mcAnalizer2.drawEventsWithAtmMuonWeights(fluxWeightName);

	IHistogram1D h1LogNpeMC1 = 
	    mcAnalizer1.makeJaida1DHistogram("logNpe",histoFactory);
	IHistogram1D h1LogNpeMC2 = 
	    mcAnalizer2.makeJaida1DHistogram("logNpe",histoFactory);
	IHistogram1D h1LogNpeData = 
	    dataAnalizer.makeJaida1DHistogram("logNpe",histoFactory);


	if(StatisticalComparison.canCompare(h1LogNpeMC1,h1LogNpeData,"chi2") &&
	   StatisticalComparison.canCompare(h1LogNpeMC2,h1LogNpeData,"chi2")){

	    int nbin = h1LogNpeMC1.axis().bins();
	    double chi2 = 0.0;
	    int nDof = 0;
	    for(int i= 0;i<nbin;i++){
		double nEventsData = h1LogNpeData.binHeight(i);
		double nEntriesData = h1LogNpeData.binEntries(i);
		double nEventsMC;
		double nEntriesMC;
		double nErrorMC;
		if(h1LogNpeData.axis().binLowerEdge(i)<logNpeBoundary){
		    nEventsMC = h1LogNpeMC2.binHeight(i);
		    nEntriesMC = h1LogNpeMC2.binEntries(i);
		    nErrorMC = h1LogNpeMC2.binError(i);
		}else{
		    nEventsMC = h1LogNpeMC1.binHeight(i);
		    nEntriesMC = h1LogNpeMC1.binEntries(i);
		    nErrorMC = h1LogNpeMC1.binError(i);
		}
		//double weight = nEventsData + nEntriesMC;
		double weight = nEventsData + nErrorMC*nErrorMC;
		if(weight>0.0 && nEventsData > 0.0){
		    double dw = nEventsData - nEventsMC;
		    chi2 += dw*dw/weight;
		    nDof++;
		}
	    }

	    ComparisonResult result = new ComparisonResult();
	    result.setQuality(chi2);
	    result.setnDof(nDof);
	    return result;
	}else{
	    System.out.println("Cannot compare the results");
	    return null;
	}
    }
}
