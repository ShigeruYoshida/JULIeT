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

public class DrawAtmMuonBundleCombined {

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
	String mcSlope2DataFileName = 
	    "iceCube/uhe/analysis/EHEMCReducedIC9Slope2I3Particles";
	String mcSlope1DataFileName = 
	    "iceCube/uhe/analysis/EHEMCReducedIC9Slope1I3Particles";

	//MC data Analyzer (E**-2)
	InputStream in = ClassLoader.getSystemResourceAsStream(mcSlope2DataFileName);
	I3ParticleAnalysisFactory mcSlope2Analizer = new I3ParticleAnalysisFactory(in);
	in.close();

	//MC data Analyzer (E**-1)
	in = ClassLoader.getSystemResourceAsStream(mcSlope1DataFileName);
	I3ParticleAnalysisFactory mcSlope1Analizer = new I3ParticleAnalysisFactory(in);
	in.close();

	// real data
	in = ClassLoader.getSystemResourceAsStream(realDataFileName);
	I3ParticleAnalysisFactory dataAnalizer = new I3ParticleAnalysisFactory(in,true);
	in.close();

	J3Vector ic9Center = new J3Vector(mcSlope1Analizer.xCenterOfIC9,
					  mcSlope1Analizer.yCenterOfIC9, 
					  0.0);
	// IceCube 9 string array center


	// Sets the Criteria
	Criteria criteria = new Criteria();
	criteria.setThresholdOfLogNpe(3.8);
	criteria.setThresholdOfNDOMs(80);
	mcSlope2Analizer.setCriteria(criteria);
	mcSlope1Analizer.setCriteria(criteria);
	dataAnalizer.setCriteria(criteria);

	// Set the Anaizer's bin size and range
	mcSlope2Analizer.minLogNpe = 3.8;
	mcSlope2Analizer.maxLogNpe = 5.0;
	mcSlope2Analizer.setBinSize(0.025,0.025,0.05,0.5); 
	                // Delta (LogE) Delta(LogNpe), Delta cos(Zenith), DeltaFGquality
	double time = 0.5*mcSlope2Analizer.observationTime;
	mcSlope2Analizer.setObservationTime(time);
                	// a half of the observation time must be set for averaing with E**-1 data

	mcSlope1Analizer.minLogNpe = 3.8;
	mcSlope1Analizer.maxLogNpe = 5.0;
	mcSlope1Analizer.setBinSize(0.025,0.025,0.05,0.5); 
	                // Delta (LogE) Delta(LogNpe), Delta cos(Zenith), DeltaFGquality
	time = 0.5*mcSlope1Analizer.observationTime;
	mcSlope1Analizer.setObservationTime(time);
                	// a half of the observation time must be set for averaing with E**-1 data

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

	AtmMuonBundleFlux muonFlux = new AtmMuonBundleFlux();
	ParticlePoint s = new ParticlePoint(0.0,0.0,0); // ice

	// Now Drawing
	muonFlux.setMuonThresholdEnergy(muEThreshold);
	muonFlux.setAlpha(alpha);
	//muonFlux.setCutOffFeature(true);
	AtmMuonBundleFitter.setAtmMuonBundleFlux(mcSlope2Analizer,s,
						 muonFlux);
	AtmMuonBundleFitter.setAtmMuonBundleFlux(mcSlope1Analizer,s,
						 muonFlux);

	mcSlope2Analizer.drawEventsWithAtmMuonWeights(AtmMuonBundleFitter.fluxWeightName);
	mcSlope2Analizer.switchToReco();
	mcSlope1Analizer.drawEventsWithAtmMuonWeights(AtmMuonBundleFitter.fluxWeightName);
	mcSlope1Analizer.switchToReco();

	dataAnalizer.switchToReco();

	jaidaTree.mkdir("/logNpe");
	jaidaTree.cd("/logNpe");
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
    }


}
