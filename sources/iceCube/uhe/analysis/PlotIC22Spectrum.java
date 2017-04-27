package iceCube.uhe.analysis;

import iceCube.uhe.analysis.*;
import geometry.*;

import hep.aida.*;

import java.io.*;
import java.util.*;

public class PlotIC22Spectrum {

    public static void main(String[] args) throws IOException{

	String realDataFileName = "iceCube/uhe/analysis/EHEIC22Real_NPE3k_I3Particles";
	String mcNuMuDataFileName = "iceCube/uhe/analysis/EHEMCIC22Slope1NuMu_mtx_flux_I3Particles";
	String mcNuEDataFileName = "iceCube/uhe/analysis/EHEMCIC22Slope1NuE_mtx_flux_I3Particles";
	String mcNuTauDataFileName = "iceCube/uhe/analysis/EHEMCIC22Slope1NuTau_mtx_flux_I3Particles";
        String mcMuonDataFileName = "iceCube/uhe/analysis/EHEMCIC22Slope1Mu_mtx_allfluxes_I3Particles";
        String mcTauDataFileName = "iceCube/uhe/analysis/EHEMCIC22Slope1TAU_mtx_flux_I3Particles";
        String mcMuonSlope1DataFileName = 
	    "iceCube/uhe/analysis/EHEMCIC22Slope1Mu_mtx_allfluxes_I3Particles";
        String mcMuonSlope2DataFileName = 
	    "iceCube/uhe/analysis/EHEMCIC22Slope2Mu_atmFlux_I3Particles";

	double logNpeBoundary = 4.7; // slope1 > 4.7 : 4.7 > slope2
	double cosZenLow = -1.0;
	double cosZenUp = 1.0;

	I3ParticleAnalysisFactory.minNDOMsToAnalize = 80; // 80 DOMs cut
	I3ParticleAnalysisFactory.minLogNpe = 4.0; // NPE cut
	InputStream in = ClassLoader.getSystemResourceAsStream(realDataFileName);
	I3ParticleIC22AnalysisFactory analizerData = new I3ParticleIC22AnalysisFactory(in,true);

	// Atm Muon MC slope1
	in = ClassLoader.getSystemResourceAsStream(mcMuonSlope1DataFileName);
	I3ParticleIC22AnalysisFactory analizerBGSlope1 = 
	    new I3ParticleIC22AnalysisFactory(in);
	in.close();

	// Atm Muon MC slope2
	in = ClassLoader.getSystemResourceAsStream(mcMuonSlope2DataFileName);
	I3ParticleIC22AnalysisFactory analizerBGSlope2 = 
	    new I3ParticleIC22AnalysisFactory(in);
	in.close();

	// GZK signal
	in = ClassLoader.getSystemResourceAsStream(mcMuonDataFileName);
	I3ParticleIC22AnalysisFactory analizerSIG = new I3ParticleIC22AnalysisFactory(in);
	in = ClassLoader.getSystemResourceAsStream(mcTauDataFileName);
	analizerSIG.readI3Particles(in);
	in.close();
	in = ClassLoader.getSystemResourceAsStream(mcNuEDataFileName);
	analizerSIG.readI3Particles(in);
	in.close();
	in = ClassLoader.getSystemResourceAsStream(mcNuMuDataFileName);
	analizerSIG.readI3Particles(in);
	in.close();
	in = ClassLoader.getSystemResourceAsStream(mcNuTauDataFileName);
	analizerSIG.readI3Particles(in);
	in.close();

	// Sets the Criteria for Data
	CriteriaIC22 criteriaData = new CriteriaIC22();
	//criteriaData.setThresholdOfLogNpe(4.0);
	//criteriaData.setThresholdOfNDOMs(80);
	//criteriaData.setRangeOfCosineOfZenith(1.0,-1.0,1.0,1.0);// block
	//criteriaData.setMinimumBound(cosZenUp,1.0); // cosZen < cosUp logNpe>1
	//criteriaData.setSimpleCosZenithCut(cosZenLow,false);// cosZen >cosLow
	criteriaData.setCOBZcut();
	criteriaData.setEHESuperCut();

	// Sets the Criteria for MC slope 1
	CriteriaIC22 criteriaBGSlope1 = new CriteriaIC22();
	//criteriaBGSlope1.setThresholdOfLogNpe(logNpeBoundary);
	//criteriaBGSlope1.setThresholdOfNDOMs(80);
	//criteriaBGSlope1.setRangeOfCosineOfZenith(1.0,-1.0,1.0,1.0);// block
	//criteriaBGSlope1.setMinimumBound(cosZenUp,1.0); // cosZen < cosUp logNpe>1
	//criteriaBGSlope1.setSimpleCosZenithCut(cosZenLow,false);// cosZen >cosLow
	criteriaBGSlope1.setCOBZcut();
	criteriaBGSlope1.setEHESuperCut();

	// Sets the Criteria for MC slope 2
	CriteriaIC22 criteriaBGSlope2 = new CriteriaIC22();
	// logNpeBoundary> logNPE > 4.0
	//criteriaBGSlope2.setThresholdOfLogNpe(4.0);
	//criteriaBGSlope2.setSimpleNpeCut(logNpeBoundary,false);
	//criteriaBGSlope2.setThresholdOfNDOMs(80);
	//criteriaBGSlope2.setRangeOfCosineOfZenith(1.0,-1.0,1.0,1.0);// block
	//criteriaBGSlope2.setMinimumBound(cosZenUp,1.0); // cosZen < cosUp logNpe>1
	//criteriaBGSlope2.setSimpleCosZenithCut(cosZenLow,false);// cosZen >cosLow
	criteriaBGSlope2.setCOBZcut();
	criteriaBGSlope2.setThresholdOfLogNpe(10.0); // for blocking slope2 MC

	// Sets the Criteria for Signal
	CriteriaIC22 criteriaSIG = new CriteriaIC22();
	//criteriaSIG.setThresholdOfLogNpe(4.0);
	//criteriaSIG.setThresholdOfNDOMs(80);
	//criteriaSIG.setRangeOfCosineOfZenith(1.0,-1.0,1.0,1.0);// block
	//criteriaSIG.setMinimumBound(cosZenUp,1.0); // cosZen < cosUp logNpe>1
	//criteriaSIG.setSimpleCosZenithCut(cosZenLow,false);// cosZen >cosLow
	criteriaSIG.setCOBZcut();
	criteriaSIG.setEHESuperCut();

	analizerData.setCriteria(criteriaData);
	analizerBGSlope1.setCriteria(criteriaBGSlope1);
	analizerBGSlope2.setCriteria(criteriaBGSlope2);
	analizerSIG.setCriteria(criteriaSIG);

	// Set the Anaizer's bin size and range
	analizerData.maxLogNpe = 7.0;
	analizerBGSlope1.maxLogNpe = 7.0;
	analizerBGSlope2.maxLogNpe = 7.0;
	analizerSIG.maxLogNpe = 7.0;

	//analizerData.maxLogEnergy = 12.0;
	analizerBGSlope1.maxLogEnergy = 12.0; 
	analizerBGSlope2.maxLogEnergy = 12.0; 
	analizerSIG.maxLogEnergy = 12.0;
	analizerData.setBinSize(0.05,0.05,0.05,0.1); 
	analizerBGSlope1.setBinSize(0.05,0.05,0.05,0.1); 
	analizerBGSlope2.setBinSize(0.05,0.05,0.05,0.1); 
	analizerSIG.setBinSize(0.05,0.05,0.05,0.1); 
	                // Delta (LogE) Delta(LogNpe), Delta cos(Zenith), DeltaFGquality


	// Jaida FreeHep objects
	IAnalysisFactory jaidaFactory = IAnalysisFactory.create();
	ITree jaidaTree = jaidaFactory.createTreeFactory().
	    createTree("IC221DSpectrum_AtmMu_GZK_YT_4_4.aida","xml",ITreeFactory.RECREATE);
	//ITree jaidaTree = jaidaFactory.createTreeFactory().create();
	IHistogramFactory jaidaHistoFactory = 
	    jaidaFactory.createHistogramFactory(jaidaTree);
	IFunctionFactory jaidaFuncFactory = null;

	//
	//
	// Drawing
	//
	//
	analizerData.switchToReco();

	analizerBGSlope1.switchToReco();
	analizerBGSlope1.drawEventsWithAtmMuonWeights("elbert_1_97_1505");
	analizerBGSlope2.switchToReco();
	analizerBGSlope2.drawEventsWithAtmMuonWeights("elbert_1_97_1505");

	analizerSIG.switchToReco();
	analizerSIG.drawEventsWithGZKWeights("GZK_YT_4_4");
	//analizerSIG.drawEventsWithGZKWeights("TD_SUSY");
	
	//
	// 1D Histograms
	//
	jaidaTree.mkdir("/logNpe");
	jaidaTree.cd("/logNpe");
	IHistogram1D h1LogNpeData = 
	    analizerData.makeJaida1DHistogram("logNpe","Data", false, jaidaHistoFactory);
	IHistogram1D h1LogNpeBGSlope1 = 
	  analizerBGSlope1.makeJaida1DHistogram("logNpe","BG Slope1", false, jaidaHistoFactory);
	IHistogram1D h1LogNpeBGSlope2 = 
	  analizerBGSlope2.makeJaida1DHistogram("logNpe","BG Slope2", false, jaidaHistoFactory);
	IHistogram1D h1LogNpeBGCombined = 
	    jaidaHistoFactory.add("Atmospheric Muon",h1LogNpeBGSlope1,h1LogNpeBGSlope2);
	
	IHistogram1D h1LogNpeSIG = 
	   analizerSIG.makeJaida1DHistogram("logNpe","GZK", false, jaidaHistoFactory);


	jaidaTree.mkdir("/cosZenith");
	jaidaTree.cd("/cosZenith");
	IHistogram1D h1cosZenithData = 
	    analizerData.makeJaida1DHistogram("cosZenith","Data", false, jaidaHistoFactory);
	IHistogram1D h1cosZenithBGSlope1 = 
	  analizerBGSlope1.makeJaida1DHistogram("cosZenith","BG Slope1", false, jaidaHistoFactory);
	IHistogram1D h1cosZenithBGSlope2 = 
	  analizerBGSlope2.makeJaida1DHistogram("cosZenith","BG Slope2", false, jaidaHistoFactory);
	IHistogram1D h1cosZenithBGCombined = 
	    jaidaHistoFactory.add("Atmospheric Muon",h1cosZenithBGSlope1,h1cosZenithBGSlope2);
	IHistogram1D h1cosZenithSIG = 
	  analizerSIG.makeJaida1DHistogram("cosZenith","GZK", false, jaidaHistoFactory);

	jaidaTree.commit();

	//
	// plotting
	//
	IPlotterFactory plotterFactory = jaidaFactory.createPlotterFactory();
	IPlotter plotter = plotterFactory.create("Npe and Cos(Zenith) distribution");
	plotter.destroyRegions();
	plotter.createRegion(0,0,0.66,1);
	plotter.createRegions(2,1);

	IPlotterStyle npeStyle = plotter.region(0).style();
	JaidaPlotStyleSetter.setPlotStyle(npeStyle,
					  "Log(Npe)","Number of Events");
	plotter.region(0).plot(h1LogNpeData);
	plotter.region(0).plot(h1LogNpeBGCombined);
	plotter.region(0).plot(h1LogNpeSIG);

	IPlotterStyle cosStyle = plotter.region(1).style();
	JaidaPlotStyleSetter.setPlotStyle(cosStyle,
					      "cos(ZenithAngle)","Number of Events");
	plotter.region(1).plot(h1cosZenithData);
	plotter.region(1).plot(h1cosZenithBGCombined);
	plotter.region(1).plot(h1cosZenithSIG);

	plotter.show();

	int nbin = h1LogNpeSIG.axis().bins();
	double sigRate = 0.0;
	for(int ix = 0; ix< nbin; ix++){
	    sigRate += h1LogNpeSIG.binHeight(ix);
	}
	System.out.println("signal Rate = " + sigRate + " ");
	nbin = h1LogNpeSIG.axis().bins();
	double bgRate = 0.0;
	for(int ix = 0; ix< nbin; ix++){
	    bgRate += h1LogNpeBGCombined.binHeight(ix);
	}
	System.out.println("BG Rate = " + bgRate + " ");

    }


}
