package iceCube.uhe.analysis;

import iceCube.uhe.analysis.*;
import iceCube.uhe.points.*;
import iceCube.uhe.muonModel.*;
import geometry.*;

import hep.aida.*;
import hep.aida.ext.*;
import hep.aida.util.*;
import hep.aida.util.comparison.*;

import java.io.*;
import java.util.*;

public class RunAtmMuonBundleFitter {

    public static void main(String[] args) throws IOException{

	boolean doPlot = true;
	String realDataFileName = "iceCube/uhe/analysis/EHERealI3Particles";
	//String mcDataFileName = 
	//  "iceCube/uhe/analysis/EHEMCReducedIC9WeightFilledI3Particles";
	String mcDataFileName = 
	    "iceCube/uhe/analysis/EHEMCReducedIC9Slope2I3Particles";

	//MC data Analyzer
	InputStream in = ClassLoader.getSystemResourceAsStream(mcDataFileName);
	I3ParticleAnalysisFactory mcAnalizer = new I3ParticleAnalysisFactory(in);
	in.close();

	in = ClassLoader.getSystemResourceAsStream(realDataFileName);
	I3ParticleAnalysisFactory dataAnalizer = new I3ParticleAnalysisFactory(in,true);
	                 // Filter out all the data in the bad run. This is real data analysis.
	in.close();

	J3Vector ic9Center = new J3Vector(mcAnalizer.xCenterOfIC9,
					  mcAnalizer.yCenterOfIC9, 
					  0.0);
	// IceCube 9 string array center


	// Sets the Criteria
	Criteria criteria = new Criteria();
	criteria.setThresholdOfLogNpe(3.8);
	criteria.setThresholdOfNDOMs(80);
	mcAnalizer.setCriteria(criteria);
	dataAnalizer.setCriteria(criteria);

	// Set the Anaizer's bin size and range
	mcAnalizer.minLogNpe = 3.8;
	mcAnalizer.maxLogNpe = 5.0;
	mcAnalizer.setBinSize(0.025,0.025,0.05,0.5); 
	                // Delta (LogE) Delta(LogNpe), Delta cos(Zenith), DeltaFGquality
	dataAnalizer.minLogNpe = 3.8;
	dataAnalizer.maxLogNpe = 5.0;
	dataAnalizer.setBinSize(0.025,0.025,0.05,0.5); 
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
	AtmMuonBundleFlux muonFlux = new AtmMuonBundleFlux();
	ParticlePoint s = new ParticlePoint(0.0,0.0,0); // ice

	// initial values
	double chi2min = 1.0e10; 
	double alphaMin = 1.75; double muEmin = 0.0;
	double alpha = 1.75;
	double alphaMax = 2.152;
	double muonEnergyLowest = 1.0e1;// 10 GeV
	while(alpha<=alphaMax){
	    double muonThresholdEnergy = muonEnergyLowest;
	    double muonEnergyMax = 2.0e3; // 2 TeV
	    double chi2minAlpha = 1.0e10;
	    while(muonThresholdEnergy <= muonEnergyMax){

		muonFlux.setMuonThresholdEnergy(muonThresholdEnergy);
		muonFlux.setAlpha(alpha);

		// Set the flux
		System.err.println("--- now Setting the bundle fluxes to MC I3Particles");
		AtmMuonBundleFitter.setAtmMuonBundleFlux(mcAnalizer,s,muonFlux);

		// Comparison with the real data
		IComparisonResult result = 
		    AtmMuonBundleFitter.compare(mcAnalizer,dataAnalizer,jaidaHistoFactory);
		double chi2 =  result.quality()/result.nDof();
		if(chi2<chi2min){
		    chi2min = chi2;
		    muEmin = muonThresholdEnergy;
		    alphaMin = alpha;
		}
		if(chi2<chi2minAlpha){
		    chi2minAlpha = chi2;
		    muonEnergyLowest = muonThresholdEnergy;
		}
		System.out.println(" alpha(" + alpha + 
				   ") Ethreshold(" + muonThresholdEnergy +
				   ") prob = " + chi2);

		muonThresholdEnergy += 5.0; //  add 5 GeV
		if((chi2minAlpha+1.0)<chi2) break;
	    }
	    alpha += 0.02;
	    muonEnergyLowest -= 10.0; // subtract 10 GeV for mergin
	}

	// Now Drawing
	if(doPlot){
	    muonFlux.setMuonThresholdEnergy(muEmin);
	    muonFlux.setAlpha(alphaMin);
		AtmMuonBundleFitter.setAtmMuonBundleFlux(mcAnalizer,s,
							 muonFlux);
	    mcAnalizer.drawEventsWithAtmMuonWeights(AtmMuonBundleFitter.fluxWeightName);
	    mcAnalizer.switchToReco();
	    dataAnalizer.switchToReco();

	    IHistogram1D h1LogNpeMC = 
		mcAnalizer.makeJaida1DHistogram("logNpe",jaidaHistoFactory);
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
				      "cos(ZenithAngle)","Number of Events");
	    IHistogram1D h1cosZenithMC = 
		mcAnalizer.makeJaida1DHistogram("cosZenith",jaidaHistoFactory);
	    IHistogram1D h1cosZenithData = 
		dataAnalizer.makeJaida1DHistogram("cosZenith",jaidaHistoFactory);
	    IDataPointSet d1cosZenithData = jaidaDpsFactory.create("dps1DFromHist",
							   h1cosZenithData);
	    plotter.region(1).plot(h1cosZenithMC);
	    plotter.region(1).plot(d1cosZenithData);
	    plotter.show();
	    System.out.println(" alpha(" + alphaMin + 
			       ") Ethreshold(" + muEmin +
			       ") prob = " + chi2min);
	}
    }


}
