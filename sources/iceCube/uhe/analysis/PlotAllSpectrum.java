package iceCube.uhe.analysis;

import iceCube.uhe.analysis.*;
import geometry.*;

import hep.aida.*;

import java.io.*;
import java.util.*;

public class PlotAllSpectrum {

    public static void main(String[] args) throws IOException{

	String realDataFileName = "iceCube/uhe/analysis/EHERealI3Particles";
	String mcNuMuDataFileName = "iceCube/uhe/analysis/EHEMCIC80Slope1NuMu_mtx_flux_I3Particles";
	String mcNuEDataFileName = "iceCube/uhe/analysis/EHEMCIC80Slope1NuE_mtx_flux_I3Particles";
	String mcNuTauDataFileName = "iceCube/uhe/analysis/EHEMCIC80Slope1NuTau_mtx_flux_I3Particles";
        String mcMuonDataFileName = "iceCube/uhe/analysis/EHEMCIC80DatSet498-594Slope1Mu_mtx_allfluxes_I3Particles";
        String mcTauDataFileName = "iceCube/uhe/analysis/EHEMCIC80DatSet499-641Slope1TAU_mtx_flux_I3Particles";

	//InputStream in = ClassLoader.getSystemResourceAsStream(realDataFileName);
	//I3ParticleAnalysisFactory analizer = new I3ParticleAnalysisFactory(in,true);
	InputStream in = ClassLoader.getSystemResourceAsStream(mcMuonDataFileName);
	I3ParticleAnalysisFactory analizer = new I3ParticleAnalysisFactory(in);
	in.close();
	in = ClassLoader.getSystemResourceAsStream(mcTauDataFileName);
	analizer.readI3Particles(in);
	in.close();
	in = ClassLoader.getSystemResourceAsStream(mcNuEDataFileName);
	analizer.readI3Particles(in);
	in.close();
	in = ClassLoader.getSystemResourceAsStream(mcNuMuDataFileName);
	analizer.readI3Particles(in);
	in.close();
	in = ClassLoader.getSystemResourceAsStream(mcNuTauDataFileName);
	analizer.readI3Particles(in);
	in.close();

	in = ClassLoader.getSystemResourceAsStream(mcMuonDataFileName);
	I3ParticleAnalysisFactory analizerBG = new I3ParticleAnalysisFactory(in);
	in.close();

	J3Vector ic9Center = new J3Vector(analizer.xCenterOfIC9,
					  analizer.yCenterOfIC9, 
					  0.0);
	// IceCube 9 string array center

	analizer.observationTime = 365.0*24.0*3600.0; // 1 year
	analizerBG.observationTime = 365.0*24.0*3600.0; // 1 year


	// Sets the Criteria
	Criteria criteria = new Criteria();
	criteria.setThresholdOfLogNpe(5.2);
	criteria.setRangeOfCosineOfZenith(0.1,1.0,6.0,6.5);
	criteria.setMinimumBound(0.1,5.2);
	criteria.setThresholdOfNDOMs(80);

	//criteria.setThresholdOfLogNpe(3.8);
	//criteria.setThresholdOfLogNpe(5.0);
	//criteria.setThresholdOfNDOMs(80);
	analizer.setCriteria(criteria);
	analizerBG.setCriteria(criteria);

	// Set the Anaizer's bin size and range
	analizer.minLogNpe = 3.8;	analizerBG.minLogNpe = 3.8;
	analizer.minLogNpe = 5.0;	analizerBG.minLogNpe = 5.0;
	analizer.maxLogNpe = 8.0;	analizerBG.maxLogNpe = 8.0;
	analizer.setBinSize(0.1,0.1,0.05,0.1); 	analizerBG.setBinSize(0.1,0.1,0.05,0.1); 
	                // Delta (LogE) Delta(LogNpe), Delta cos(Zenith), DeltaFGquality


	// Jaida FreeHep objects
	IAnalysisFactory jaidaFactory = IAnalysisFactory.create();
	//ITree jaidaTree = jaidaFactory.createTreeFactory().create("I3Particles.aida",
	//							  "",false,true);
	ITree jaidaTree = jaidaFactory.createTreeFactory().create();
	IHistogramFactory jaidaHistoFactory = 
	    jaidaFactory.createHistogramFactory(jaidaTree);
	IFunctionFactory jaidaFuncFactory = null;

	//
	//
	// Drawing
	//
	//
	//analizer.switchToMCTruth();
	analizer.switchToReco();
	//analizerBG.switchToMCTruth();
	analizerBG.switchToReco();
	
	//
	// 1D Histograms
	//
	jaidaTree.mkdir("/logNpe");
	jaidaTree.cd("/logNpe");

	//analizer.drawEventsWithGZKWeights("TD_SUSY");
	//analizer.drawEventsWithGZKWeights("GZK_sigl");
	analizer.drawEventsWithGZKWeights("GZK_YT_4_4");
	IHistogram1D h1LogNpeSig = 
	    analizer.makeJaida1DHistogram("logNpe",jaidaHistoFactory);

	analizerBG.drawEventsWithAtmMuonWeights(">Elbert_2_04_3730");
	IHistogram1D h1LogNpeBG= 
	    analizerBG.makeJaida1DHistogram("logNpe",jaidaHistoFactory);

	jaidaTree.mkdir("/cosZenith");
	jaidaTree.cd("/cosZenith");

	IHistogram1D h1cosZenithSig = 
	    analizer.makeJaida1DHistogram("cosZenith",jaidaHistoFactory);

	IHistogram1D h1cosZenithBG = 
	    analizerBG.makeJaida1DHistogram("cosZenith",jaidaHistoFactory);

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
					  "Log(Npe) ","Number of Events");
	plotter.region(0).plot(h1LogNpeSig);
	plotter.region(0).plot(h1LogNpeBG);

	IPlotterStyle cosStyle = plotter.region(1).style();
	JaidaPlotStyleSetter.setPlotStyle(cosStyle,
					  "cos(First-Guessed Zenith Angle) ","Number of Events");
	plotter.region(1).plot(h1cosZenithSig);
	plotter.region(1).plot(h1cosZenithBG);

	plotter.show();

	int nbin = h1LogNpeSig.axis().bins();
	double sigRate = 0.0;
	for(int ix = 0; ix< nbin; ix++){
	    sigRate += h1LogNpeSig.binHeight(ix);
	}
	nbin = h1LogNpeBG.axis().bins();
	double bgRate = 0.0;
	for(int ix = 0; ix< nbin; ix++){
	    bgRate += h1LogNpeBG.binHeight(ix);
	}

	System.out.println("Sig Rate = " + sigRate + " BG rate = " + bgRate);


    }

}
