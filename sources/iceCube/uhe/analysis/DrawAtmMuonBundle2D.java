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

public class DrawAtmMuonBundle2D {

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

        criteria.setThresholdOfLogNpe(4.6); 
        criteria.setRangeOfCosineOfZenith(0.1,1.0,4.7,5.8); 
        criteria.setMinimumBound(0.1,4.6);
        criteria.setThresholdOfNDOMs(80); 


	mcSlope2Analizer.setCriteria(criteria);
	mcSlope1Analizer.setCriteria(criteria);
	dataAnalizer.setCriteria(criteria);

	// Set the Anaizer's bin size and range
	mcSlope2Analizer.minLogNpe = 3.8;
	mcSlope2Analizer.maxLogNpe = 6.0;
	mcSlope2Analizer.maxLogEnergy = 14.0;
	mcSlope2Analizer.setBinSize(0.2,0.2,0.1,0.5); 
	                // Delta (LogE) Delta(LogNpe), Delta cos(Zenith), DeltaFGquality
	double time = 0.5*mcSlope2Analizer.observationTime;
	mcSlope2Analizer.setObservationTime(time);
                	// a half of the observation time must be set for averaing with E**-1 data

	mcSlope1Analizer.minLogNpe = 3.8;
	mcSlope1Analizer.maxLogNpe = 6.0;
	mcSlope1Analizer.maxLogEnergy = 14.0;
	mcSlope1Analizer.setBinSize(0.2,0.2,0.1,0.5); 
	                // Delta (LogE) Delta(LogNpe), Delta cos(Zenith), DeltaFGquality
	time = 0.5*mcSlope1Analizer.observationTime;
	mcSlope1Analizer.setObservationTime(time);
                	// a half of the observation time must be set for averaing with E**-1 data

	dataAnalizer.minLogNpe = 3.8;
	dataAnalizer.maxLogNpe = 6.0;
	dataAnalizer.setBinSize(0.05,0.05,0.1,0.5); 
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
	//mcSlope2Analizer.drawEventsWithNoWeights();
	mcSlope2Analizer.switchToReco();
	mcSlope1Analizer.drawEventsWithAtmMuonWeights(AtmMuonBundleFitter.fluxWeightName);
	//mcSlope1Analizer.drawEventsWithNoWeights();
	mcSlope1Analizer.switchToReco();

	dataAnalizer.switchToReco();

	jaidaTree.mkdir("/logNpe");
	jaidaTree.cd("/logNpe");
	IHistogram2D h2LogNpeMC1 = 
	    mcSlope1Analizer.makeJaida2DHistogram("logRecoE-Npe",jaidaHistoFactory);
	IHistogram2D h2LogNpeMC2 = 
	    mcSlope2Analizer.makeJaida2DHistogram("logRecoE-Npe",jaidaHistoFactory);

	IHistogram2D h2LogNpeMC = jaidaHistoFactory.add("Energy - Npe",
							h2LogNpeMC1,h2LogNpeMC2);

	IPlotterFactory plotterFactory = jaidaFactory.createPlotterFactory();
	IPlotter plotter = plotterFactory.create("CR Energy - Npe");
	IPlotterStyle npeStyle = plotter.region(0).style();
	JaidaPlotStyleSetter.setPlotStyle(npeStyle,
					  "Log(CR Energy)","Log(NPE)",
					  "hist2DStyle","colorMap"); 
	plotter.region(0).plot(h2LogNpeMC);
	plotter.show();
    }


}
