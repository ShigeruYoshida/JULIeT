package iceCube.uhe.analysis;

import iceCube.uhe.analysis.*;
import geometry.*;

import hep.aida.*;

import java.io.*;
import java.util.*;

public class PlotCRSpectrum {

    public static void main(String[] args) throws IOException{

	String mcMuonDataFileName = "iceCube/uhe/analysis/EHEMCReducedIC9Slope1_mtx_bundleflux_I3Particles";
	    //String mcMuonDataFileName = "iceCube/uhe/analysis/EHEMCIC80DatSet498-594Slope1Mu_mtx_allfluxes_I3Particles";

	InputStream in = ClassLoader.getSystemResourceAsStream(mcMuonDataFileName);
	I3ParticleAnalysisFactory analizer = 
	    new I3ParticleAnalysisFactory(in);
	in.close();

	// IceCube 9 string array center
	//analizer.observationTime = 365.0*24.0*3600.0; // 1 year


	// Sets the Signal Criteria
	Criteria signalCriteria = new Criteria();
	signalCriteria.setThresholdOfLogNpe(4.6);
	signalCriteria.setRangeOfCosineOfZenith(0.1,1.0,4.7,5.8);
	signalCriteria.setMinimumBound(0.1,4.6);
	signalCriteria.setThresholdOfNDOMs(80);

	// Sets the filter Criteria
	Criteria filterCriteria = new Criteria();
	filterCriteria.setThresholdOfLogNpe(3.8);
	filterCriteria.setThresholdOfNDOMs(80);

	// Set the Anaizer's bin size and range
	analizer.maxLogEnergy = 12.0;
	analizer.setBinSize(0.1,0.1,0.05,0.1); 
	// Delta (LogE) Delta(LogNpe), Delta cos(Zenith), DeltaFGquality


	// Jaida FreeHep objects
	IAnalysisFactory jaidaFactory = IAnalysisFactory.create();
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
	
	//
	// 1D Histograms
	//
	jaidaTree.mkdir("/logEnergy");
	jaidaTree.cd("/logEnergy");

	analizer.drawEventsWithAtmMuonWeights("Elbert_2_04_3730");

	analizer.setCriteria(filterCriteria);
	IHistogram1D h1LogEnergyBefore = 
	    analizer.makeJaida1DHistogram("logECR",jaidaHistoFactory);

	analizer.setCriteria(signalCriteria);
	IHistogram1D h1LogEnergyAfter = 
	    analizer.makeJaida1DHistogram("logECR",jaidaHistoFactory);

	//
	// plotting
	//
	IPlotterFactory plotterFactory = jaidaFactory.createPlotterFactory();
	IPlotter plotter = plotterFactory.create("Cosmic Ray Energy distribution");

	IPlotterStyle logEStyle = plotter.region(0).style();
	JaidaPlotStyleSetter.setPlotStyle(logEStyle,
				  "Log(Energy [GeV]) ","Number of Events");
	plotter.region(0).plot(h1LogEnergyBefore);
	plotter.region(0).plot(h1LogEnergyAfter);

	plotter.show();

	int nbin = h1LogEnergyAfter.axis().bins();
	double sigRate = 0.0;
	for(int ix = 0; ix< nbin; ix++){
	    sigRate += h1LogEnergyAfter.binHeight(ix);
	}
	nbin = h1LogEnergyBefore.axis().bins();
	double bgRate = 0.0;
	for(int ix = 0; ix< nbin; ix++){
	    bgRate += h1LogEnergyBefore.binHeight(ix);
	}

	System.out.println("Sig Rate = " + sigRate + " BG rate = " + bgRate);


    }

}
