package iceCube.uhe.analysis;

import iceCube.uhe.analysis.*;
import iceCube.uhe.muonModel.*;
import geometry.*;

import hep.aida.*;
import hep.aida.ext.*;
import hep.aida.util.*;
import hep.aida.util.comparison.*;

import java.io.*;
import java.util.*;

/**
   Fit the IC9 data with the Elbert Formula using
   AtmMuonBundleFetter class. The muon bundle intensity
   at the IceCube depth expected by the Elbert model
   is given by ElbertFluxTableFactory in the muon model
   package.

   Written by S. Yoshida 2008 June 2nd

 */

public class FitElbertModelFixedEthWithIC9Data {

    public static void main(String[] args) throws IOException{

	boolean doPlot = true; // plot the Npe and Cos(Zenith) spectrum

	// ElbertFlux calculator
	ElbertFluxTableFactory.tablePath="../data/muonBundleFluxTable3TeVthreshold/";
	ElbertFluxTableFactory.alphaMin= 1.8;
	ElbertFluxTableFactory.numberOfAlphaSteps = 56;
	ElbertFluxTableFactory.maxNumberOfMuEThSteps = 1;
	ElbertFluxTableFactory muonFluxTable = new ElbertFluxTableFactory();

	String realDataFileName = "iceCube/uhe/analysis/EHERealI3Particles";
	String mcSlope2DataFileName = 
	    "iceCube/uhe/analysis/EHEMCReducedIC9Slope2I3Particles";
	String mcSlope1DataFileName = 
	    "iceCube/uhe/analysis/EHEMCReducedIC9Slope1I3Particles";

	I3ParticleAnalysisFactory.minNDOMsToAnalize = 60;

	//MC data Analyzer (E**-2)
	InputStream in = ClassLoader.getSystemResourceAsStream(mcSlope2DataFileName);
	I3ParticleAnalysisFactory mcSlope2Analizer = new I3ParticleAnalysisFactory(in);
	in.close();

	//MC data Analyzer (E**-1)
	in = ClassLoader.getSystemResourceAsStream(mcSlope1DataFileName);
	I3ParticleAnalysisFactory mcSlope1Analizer = new I3ParticleAnalysisFactory(in);
	in.close();

	// Real data
	in = ClassLoader.getSystemResourceAsStream(realDataFileName);
	I3ParticleAnalysisFactory dataAnalizer = new I3ParticleAnalysisFactory(in,true);
	                 // Filter out all the data in the bad run. This is real data analysis.
	in.close();

	J3Vector ic9Center = new J3Vector(mcSlope1Analizer.xCenterOfIC9,
					  mcSlope1Analizer.yCenterOfIC9, 
					  0.0);
	// IceCube 9 string array center


	// Sets the Criteria
	Criteria criteria = new Criteria();
	criteria.setThresholdOfLogNpe(3.8);
	criteria.setThresholdOfNDOMs(80);
	dataAnalizer.setCriteria(criteria);

	Criteria criteriaSlope2 = new Criteria();
	criteriaSlope2.setThresholdOfLogNpe(3.8);
	criteriaSlope2.setThresholdOfNDOMs(80);
	criteriaSlope2.setSimpleNpeCut(4.2,false); // logNpe < 4.2
	mcSlope2Analizer.setCriteria(criteriaSlope2);

	Criteria criteriaSlope1 = new Criteria();
	criteriaSlope1.setThresholdOfLogNpe(3.8);
	criteriaSlope1.setThresholdOfNDOMs(80);
	criteriaSlope1.setSimpleNpeCut(4.2,true); // logNpe >= 4.2
	mcSlope1Analizer.setCriteria(criteriaSlope1);

	// Set the Anaizer's bin size and range
	mcSlope2Analizer.minLogNpe = 4.0;
	mcSlope2Analizer.maxLogNpe = 5.0;
	mcSlope2Analizer.setBinSize(0.05,0.05,0.05,0.5); 
	                // Delta (LogE) Delta(LogNpe), Delta cos(Zenith), DeltaFGquality

	mcSlope1Analizer.minLogNpe = 4.0;
	mcSlope1Analizer.maxLogNpe = 5.0;
	mcSlope1Analizer.setBinSize(0.05,0.05,0.05,0.5); 
	                // Delta (LogE) Delta(LogNpe), Delta cos(Zenith), DeltaFGquality

                	// a half of the observation time must be set for averaing with E**-1 data
	dataAnalizer.minLogNpe = 4.0;
	dataAnalizer.maxLogNpe = 5.0;
	dataAnalizer.setBinSize(0.05,0.05,0.05,0.5); 
	                // Delta (LogE) Delta(LogNpe), Delta cos(Zenith), DeltaFGquality

	/// Jaida FreeHep objects
	IAnalysisFactory jaidaFactory = IAnalysisFactory.create();
	ITree jaidaTree = jaidaFactory.createTreeFactory().create();
	IHistogramFactory jaidaHistoFactory = 
	    jaidaFactory.createHistogramFactory(jaidaTree);
	IFunctionFactory jaidaFuncFactory = null;
	IDataPointSetFactory jaidaDpsFactory = 
	    jaidaFactory.createDataPointSetFactory(jaidaTree);

	//
	//
	// Put the Atmuospheric Muon Bundle Fluxes and fit
	//
	//

	// initial values
	double chi2min = 1.0e10; 
	double alphaMin = 1.75;
	int alphabinMin = 0; 
	double[] chi2minAlpha = new double[muonFluxTable.numberOfAlphaSteps];

	for(int ialpha=0;ialpha<muonFluxTable.numberOfAlphaSteps;ialpha++){
	    chi2minAlpha[ialpha] = 1.0e10;

	    muonFluxTable.setElbertParameters(ialpha);
	    double alpha = muonFluxTable.getAlpha();
	    double muonThresholdEnergy = muonFluxTable.getMuETh();

	    // Set the flux weight to I3Particle's
	    AtmMuonBundleFitter.setAtmMuonBundleFlux(mcSlope2Analizer,
						     muonFluxTable);
	    AtmMuonBundleFitter.setAtmMuonBundleFlux(mcSlope1Analizer,
						     muonFluxTable);

	    // Comparison with the real data
	    IComparisonResult result = 
		AtmMuonBundleFitter.compare(mcSlope1Analizer,mcSlope2Analizer,
					    dataAnalizer,jaidaHistoFactory);
	    //double chi2 =  result.quality()/result.nDof();
	    double chi2 =  result.quality();
	    chi2minAlpha[ialpha] = chi2;
	    if(chi2<chi2min){
		chi2min = chi2;
		alphaMin = alpha;
		alphabinMin = ialpha;
	    }

	    System.err.println(alpha + 
			       " " + muonThresholdEnergy +
			       " " + chi2);

	}

	// Now Drawing
	if(doPlot){

	    muonFluxTable.setElbertParameters(alphabinMin);
	    AtmMuonBundleFitter.setAtmMuonBundleFlux(mcSlope2Analizer,
						     muonFluxTable);
	    AtmMuonBundleFitter.setAtmMuonBundleFlux(mcSlope1Analizer,
						     muonFluxTable);

	    mcSlope2Analizer.drawEventsWithAtmMuonWeights(AtmMuonBundleFitter.fluxWeightName);
	    mcSlope2Analizer.switchToReco();
	    mcSlope1Analizer.drawEventsWithAtmMuonWeights(AtmMuonBundleFitter.fluxWeightName);
	    mcSlope1Analizer.switchToReco();

	    dataAnalizer.switchToReco();

	    IHistogram1D h1LogNpeMC1 = 
		mcSlope1Analizer.makeJaida1DHistogram("logNpe",jaidaHistoFactory);
	    IHistogram1D h1LogNpeMC2 = 
		mcSlope2Analizer.makeJaida1DHistogram("logNpe",jaidaHistoFactory);

	    IHistogram1D h1LogNpeMC = jaidaHistoFactory.add("Slope1+Slope2",
							    h1LogNpeMC1,h1LogNpeMC2);

	    IHistogram1D h1LogNpeData = 
		dataAnalizer.makeJaida1DHistogram("logNpe",jaidaHistoFactory);
	    IDataPointSet d1LogNpeData = jaidaDpsFactory.create("dps1DFromHist",
								h1LogNpeData);

	    IPlotterFactory plotterFactory = jaidaFactory.createPlotterFactory();
	    IPlotter plotter = plotterFactory.create("Event Distribution");
	    plotter.destroyRegions();
	    plotter.createRegion(0,0,0.66,1);
	    plotter.createRegions(2,1);

	    IPlotterStyle npeStyle = plotter.region(0).style();
	    JaidaPlotStyleSetter.setPlotStyle(npeStyle,
					      "Log(Npe)","Number of Events");
	    plotter.region(0).plot(h1LogNpeMC);
	    plotter.region(0).plot(d1LogNpeData);

	    IPlotterStyle cosStyle = plotter.region(1).style();
	    JaidaPlotStyleSetter.setPlotStyle(cosStyle,
				      "cos(Zenith) ","Number of Events");
	    IHistogram1D h1cosZenithMC1 = 
		mcSlope1Analizer.makeJaida1DHistogram("cosZenith",jaidaHistoFactory);
	    IHistogram1D h1cosZenithMC2 = 
		mcSlope2Analizer.makeJaida1DHistogram("cosZenith",jaidaHistoFactory);

	    IHistogram1D h1cosZenithMC = jaidaHistoFactory.add("Zenith Slope1+Slope2",
							       h1cosZenithMC1,h1cosZenithMC2);

	    IHistogram1D h1cosZenithData = 
		dataAnalizer.makeJaida1DHistogram("cosZenith",jaidaHistoFactory);
	    IDataPointSet d1cosZenithData = jaidaDpsFactory.create("dps1DFromHist",
							   h1cosZenithData);
	    plotter.region(1).plot(h1cosZenithMC);
	    plotter.region(1).plot(d1cosZenithData);
	    plotter.show();

	    int nbin = chi2minAlpha.length;
	    IDataPointSet chi2ps = jaidaDpsFactory.create("chi2D","Two dimensional chi plot",2);
	    int numberOfPoints = 0;
	    for(int ialpha = 0; ialpha< nbin; ialpha++){

		muonFluxTable.setElbertParameters(ialpha);
		double alpha = muonFluxTable.getAlpha();
		chi2ps.addPoint();
		chi2ps.point(numberOfPoints).coordinate(0).setValue(alpha);
		chi2ps.point(numberOfPoints).coordinate(1).setValue(chi2minAlpha[ialpha]);
		numberOfPoints++;
	    }

	    IPlotter plotterChi2 = plotterFactory.create("chi2");
	    JaidaPlotStyleSetter.setPlotStyle(plotterChi2.region(0).style(),
				      "alpha ","chi squared");
	    plotterChi2.region(0).plot(chi2ps);
	    plotterChi2.show();


	    System.err.println(" alpha(" + alphaMin + 
			       ") prob = " + chi2min);
	    System.err.println(" alphabin(" + alphabinMin);
	}

	if(args.length>0){ // f2k out
	    int nbin = chi2minAlpha.length;
	    muonFluxTable.setElbertParameters(0);
	    double alphaLeft = muonFluxTable.getAlpha();
	    muonFluxTable.setElbertParameters(nbin-1);
	    double alphaRight = muonFluxTable.getAlpha();
	    System.out.println(nbin + " " + alphaLeft + " " + alphaRight);

	    for(int ialpha = 0; ialpha< nbin; ialpha++){

		muonFluxTable.setElbertParameters(ialpha);
		double alpha = muonFluxTable.getAlpha();
		double muonThresholdEnergy = muonFluxTable.getMuETh();

		System.out.println(alpha + " " + chi2minAlpha[ialpha] + " " + 
				   muonThresholdEnergy);
	    }

	}

    }

}

