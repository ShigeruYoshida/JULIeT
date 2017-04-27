package iceCube.uhe.analysis;

import iceCube.uhe.analysis.*;
import geometry.*;

import hep.aida.*;

import java.io.*;
import java.util.*;

public class RunJaidaAnalysisIC22 {

    public static void main(String[] args) throws IOException{

	String realDataFileName = "iceCube/uhe/analysis/EHEIC22Real_NPE3k_I3Particles";
	String mcNuMuDataFileName = "iceCube/uhe/analysis/EHEMCIC22Slope1NuMu_mtx_flux_I3Particles";
	String mcNuEDataFileName = "iceCube/uhe/analysis/EHEMCIC22Slope1NuE_mtx_flux_I3Particles";
	String mcNuTauDataFileName = "iceCube/uhe/analysis/EHEMCIC22Slope1NuTau_mtx_flux_I3Particles";
        String mcMuonDataFileName = "iceCube/uhe/analysis/EHEMCIC22Slope1Mu_mtx_allfluxes_I3Particles";
        String mcTauDataFileName = "iceCube/uhe/analysis/EHEMCIC22Slope1TAU_mtx_flux_I3Particles";

	I3ParticleAnalysisFactory.minNDOMsToAnalize = 80; // 80 DOMs cut
	I3ParticleAnalysisFactory.minLogNpe = 4.0; // NPE cut

	InputStream in = ClassLoader.getSystemResourceAsStream(realDataFileName);
	I3ParticleIC22AnalysisFactory analizer = new I3ParticleIC22AnalysisFactory(in,true);
	in.close();

	//in = ClassLoader.getSystemResourceAsStream(mcMuonDataFileName);
	//I3ParticleIC22AnalysisFactory analizer = new I3ParticleIC22AnalysisFactory(in);
	//in.close();

	// Sets the Criteria
	Criteria criteria = new Criteria();
	//criteria.setThresholdOfLogNpe(4.6);
	//criteria.setRangeOfCosineOfZenith(0.1,1.0,4.7,5.8);
	//criteria.setMinimumBound(0.1,4.6);
	//criteria.setThresholdOfNDOMs(80);
	//criteria.setNPEScalingFactor(0.63)

	criteria.setThresholdOfLogNpe(4.0);
	criteria.setThresholdOfNDOMs(80);

	analizer.setCriteria(criteria);

	// Set the Anaizer's bin size and range
	analizer.minLogNpe = 4.0;
	analizer.maxLogNpe = 6.5;
	analizer.maxLogEnergy = 12.0;
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
	//analizer.drawEventsWithAtmMuonWeights("elbert_1_97_1505");
	//analizer.drawEventsWithGZKWeights("GZK_YT_4_4");
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

	    //jaidaTree.mkdir("/cosZenith");
	    //jaidaTree.cd("/cosZenith");
	    //IHistogram1D h1cosZenith = 
	    //	analizer.makeJaida1DHistogram("cosZenith",jaidaHistoFactory);
	    jaidaTree.mkdir("/logECR");
	    jaidaTree.cd("/logECR");
	    IHistogram1D h1cosZenith = 
		analizer.makeJaida1DHistogram("logECR",jaidaHistoFactory);

	    //
	    // plotting
	    //
	    IPlotterFactory plotterFactory = jaidaFactory.createPlotterFactory();
	    //IPlotter plotter = plotterFactory.create("Npe and Cos(Zenith) distribution");
	    IPlotter plotter = plotterFactory.create("Npe and E(cosmic Ray) distribution");
	    plotter.destroyRegions();
	    plotter.createRegion(0,0,0.66,1);
	    plotter.createRegions(2,1);

	    IPlotterStyle npeStyle = plotter.region(0).style();
	    JaidaPlotStyleSetter.setPlotStyle(npeStyle,
					      "Log(Npe)","Number of Events");
	    plotter.region(0).plot(h1LogNpe);

	    IPlotterStyle cosStyle = plotter.region(1).style();
	    //JaidaPlotStyleSetter.setPlotStyle(cosStyle,
	    //				      "cos(ZenithAngle)","Number of Events");
	    JaidaPlotStyleSetter.setPlotStyle(cosStyle,
	    				      "log(Energy [GeV])","Number of Events");
	    plotter.region(1).plot(h1cosZenith);

	    plotter.show();
	    double logNpeDebug = 6.0;
	    int inpe = h1LogNpe.coordToIndex(logNpeDebug);
	    System.out.println("h1 (" + logNpeDebug + 
			       ") = " + h1LogNpe.binHeight(inpe));

	    int nbin = h1LogNpe.axis().bins();
	    double sigRate = 0.0;
	    for(int ix = 0; ix< nbin; ix++){
		sigRate += h1LogNpe.binHeight(ix);
	    }
	    System.err.println("Event Rate=" + sigRate);


	}
	//
	// 2D Histograms
	//
	else{
	    jaidaTree.mkdir("/logNpe-ZenithAngle");
	    jaidaTree.cd("/logNpe-ZenithAngle");
	    IHistogram2D h2 = 
		//analizer.makeJaida2DHistogram("logECR-Npe",jaidaHistoFactory);
		//analizer.makeJaida2DHistogram("logECR-logE",jaidaHistoFactory);
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
	    //JaidaPlotStyleSetter.setPlotStyle(style,
	    //	      "Log(Primary Cosmic Ray Energy [GeV]) ",
	    //      "Log(Energy at IceCube Depth [GeV])",
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

	    if(args.length>0){ // f2k out
		int nbinx = h2.xAxis().bins();
		int nbiny = h2.yAxis().bins();
		System.out.println(nbinx + " " + analizer.minLogEnergy + " " + 
				   analizer.maxLogEnergy);
		System.out.println(nbiny + " " + analizer.minLogEnergy + " " +
				   analizer.maxLogEnergy);
		for(int iy = 0; iy< nbiny ; iy++){
		    for(int ix = 0; ix< nbinx; ix++){
			System.out.println(h2.binHeight(ix,iy));
		    }
		}
	    }

	}



    }

}
