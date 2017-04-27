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

public class DrawAtmMuonBundle {

    public static void main(String[] args) throws IOException{

	double alpha = 1.8;
	double muEThreshold = 1.0e2; // [GeV]
	if(args.length!=2){
	    System.out.println("Usage: DrawAtmMuonBundle alpha muonEth [GeV]");
	    System.exit(0);
        }else{
	    alpha = Double.valueOf(args[0]).doubleValue();
	    muEThreshold = Double.valueOf(args[1]).doubleValue();
        }

	String realDataFileName = "iceCube/uhe/analysis/EHERealI3Particles";
	String mcDataFileName = 
	    //"iceCube/uhe/analysis/EHEMCReducedIC9Slope2I3Particles";
	    "iceCube/uhe/analysis/EHEMCReducedIC9Slope1I3Particles";

	//MC data Analyzer
	InputStream in = ClassLoader.getSystemResourceAsStream(mcDataFileName);
	I3ParticleAnalysisFactory mcAnalizer = new I3ParticleAnalysisFactory(in);
	in.close();

	in = ClassLoader.getSystemResourceAsStream(realDataFileName);
	I3ParticleAnalysisFactory dataAnalizer = new I3ParticleAnalysisFactory(in,true);
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
	mcAnalizer.setBinSize(0.05,0.05,0.05,0.5); 
	                // Delta (LogE) Delta(LogNpe), Delta cos(Zenith), DeltaFGquality
	dataAnalizer.minLogNpe = 3.8;
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

	AtmMuonBundleFlux muonFlux = new AtmMuonBundleFlux();
	ParticlePoint s = new ParticlePoint(0.0,0.0,0); // ice

	// Now Drawing
	muonFlux.setMuonThresholdEnergy(muEThreshold);
	muonFlux.setAlpha(alpha);
	//muonFlux.setCutOffFeature(true);
	AtmMuonBundleFitter.setAtmMuonBundleFlux(mcAnalizer,s,
						 muonFlux);
	mcAnalizer.drawEventsWithAtmMuonWeights(AtmMuonBundleFitter.fluxWeightName);
	dataAnalizer.switchToReco();
	mcAnalizer.switchToReco();

	jaidaTree.mkdir("/logNpe");
	jaidaTree.cd("/logNpe");
	IHistogram1D h1LogNpeMC = 
	    mcAnalizer.makeJaida1DHistogram("logNpe",jaidaHistoFactory);
	IHistogram1D h1LogNpeData = 
	    dataAnalizer.makeJaida1DHistogram("logNpe",jaidaHistoFactory);
	IDataPointSet d1LogNpeData = jaidaDpsFactory.create("dps1DFromHist",
							    h1LogNpeData);

	IPlotterFactory plotterFactory = jaidaFactory.createPlotterFactory();
	IPlotter plotter = plotterFactory.create("");
	plotter.destroyRegions();
	plotter.createRegion(0,0,0.66,1);
	plotter.createRegions(2,1);

	IPlotterStyle npeStyle = plotter.region(0).style();
	JaidaPlotStyleSetter.setPlotStyle(npeStyle,
					      "Log(Npe)","Number of Events");
	plotter.region(0).plot(h1LogNpeMC);
	plotter.region(0).plot(d1LogNpeData);

	jaidaTree.mkdir("/cosZenith");
	jaidaTree.cd("/cosZenith");
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
    }


}
