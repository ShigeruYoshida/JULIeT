package iceCube.uhe.analysis;

import iceCube.uhe.analysis.*;
import geometry.*;

import hep.aida.*;

import java.io.*;
import java.util.*;

public class RunJaidaAnalysis {

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

	J3Vector ic9Center = new J3Vector(analizer.xCenterOfIC9,
					  analizer.yCenterOfIC9, 
					  0.0);
	// IceCube 9 string array center

	analizer.observationTime = 365.0*24.0*3600.0; // 1 year
	//analizer.observationTime = 365.0*24.0*3600.0*1000; // 1000 year


	// Sets the Criteria
	Criteria criteria = new Criteria();
	//criteria.setThresholdOfLogNpe(4.6);
	//criteria.setRangeOfCosineOfZenith(0.1,1.0,6.0,6.5);
	//criteria.setMinimumBound(0.1,5.2);
	//criteria.setThresholdOfNDOMs(80);

	//criteria.setThresholdOfLogNpe(3.8);
	//criteria.setThresholdOfLogNpe(5.0);
	//criteria.setThresholdOfNDOMs(80);

	criteria.setEHESuperCut();

	analizer.setCriteria(criteria);

	// Set the Anaizer's bin size and range
	analizer.minLogNpe = 3.8;
	//analizer.setBinSize(0.025,0.025,0.05,0.1); 
	analizer.setBinSize(0.1,0.1,0.05,0.1); 
	//analizer.setBinSize(0.2,0.2,0.1,0.1); 
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
	boolean plot2D = true;
	//analizer.switchToMCTruth();
	analizer.switchToReco();
	//analizer.drawEventsWithAtmMuonWeights(">Elbert_2_04_3730");
	analizer.drawEventsWithGZKWeights("GZK_YT_4_4");
	//analizer.drawEventsWithGZKWeights("TD_SUSY");
	
	//
	// 1D Histograms
	//
	if(!plot2D){
	    jaidaTree.mkdir("/logNpe");
	    jaidaTree.cd("/logNpe");
	    IHistogram1D h1LogNpe = 
	    analizer.makeJaida1DHistogram("logNpe",jaidaHistoFactory);
	    //IHistogram1D h1LogNpe = 
	    //analizer.makeJaida1DHistogram("logEnergy",jaidaHistoFactory);

	    jaidaTree.mkdir("/cosZenith");
	    jaidaTree.cd("/cosZenith");
	    IHistogram1D h1cosZenith = 
		analizer.makeJaida1DHistogram("cosZenith",jaidaHistoFactory);

	    //
	    // plotting
	    //
	    IPlotterFactory plotterFactory = jaidaFactory.createPlotterFactory();
	    IPlotter plotter = plotterFactory.create("Npe and Cos(Zenith) distribution");
	    //IPlotter plotter = plotterFactory.create("Cos(Zenith) distribution");
	    plotter.destroyRegions();
	    plotter.createRegion(0,0,0.66,1);
	    plotter.createRegions(2,1);

	    IPlotterStyle npeStyle = plotter.region(0).style();
	    JaidaPlotStyleSetter.setPlotStyle(npeStyle,
					      "Log(Npe)","Number of Events");
	    plotter.region(0).plot(h1LogNpe);

	    IPlotterStyle cosStyle = plotter.region(1).style();
	    JaidaPlotStyleSetter.setPlotStyle(cosStyle,
	    				      "cos(ZenithAngle)","Number of Events");
	    plotter.region(1).plot(h1cosZenith);

	    plotter.show();
	    double logNpeDebug = 6.0;
	    int inpe = h1LogNpe.coordToIndex(logNpeDebug);
	    System.out.println("h1 (" + logNpeDebug + 
			       ") = " + h1LogNpe.binHeight(inpe));

	}
	//
	// 2D Histograms
	//
	else{
	    jaidaTree.mkdir("/logNpe-ZenithAngle");
	    jaidaTree.cd("/logNpe-ZenithAngle");
	    IHistogram2D h2 = 
		analizer.makeJaida2DHistogram("logNpe-CosZenith",jaidaHistoFactory);
		//analizer.makeJaida2DHistogram("logE-Npe",jaidaHistoFactory);
	    //analizer.makeJaida2DHistogram("logNpe-CobX",jaidaHistoFactory);
		//analizer.makeJaida2DHistogram("CobR-CobZ",jaidaHistoFactory);
	    //jaidaTree.commit();

	    //
	    // plotting
	    //
	    IPlotterFactory plotterFactory = jaidaFactory.createPlotterFactory();
	    IPlotter plotter = plotterFactory.create("Npe - Zenith");
	    //IPlotter plotter = plotterFactory.create("Energy - Npe");
	    //IPlotter plotter = plotterFactory.create("Npe - CobX");
	    //IPlotter plotter = plotterFactory.create("CobR - CobZ");
	    IPlotterStyle style = plotter.region(0).style();
	    //JaidaPlotStyleSetter.setPlotStyle(style,
	    //	      "Log(Energy at IceCube Depth [GeV]) ",
	    //      "Log(Npe)",
	    //      "hist2DStyle","colorMap");
	    JaidaPlotStyleSetter.setPlotStyle(style,
	    	      "Log(Npe) ",
	    	      "cos(First-Guessed Zenith Angle) ",
	    	      "hist2DStyle","colorMap");
	    //JaidaPlotStyleSetter.setPlotStyle(style,
	    //	      "CobR",
	    //	      "CobZ",
	    //	      "hist2DStyle","colorMap");

	    plotter.region(0).plot(h2);

	    plotter.show();
	    //plotter.writeToFile("example.ps","PS");
	    //plotter.writeToFile("example.jpg","JPEG");
	    double logNpeDebug = 6.0; double cosZenithDebug = 0.95;
	    int inpe = h2.coordToIndexX(logNpeDebug);
	    int icos = h2.coordToIndexY(cosZenithDebug);
	    System.out.println("h2 (" + logNpeDebug + " " + cosZenithDebug +
			       ") = " + h2.binHeight(inpe,icos));

	}


    }

}
