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

public class DrawAtmMuonBundleByRelativeWeight {

    public static void main(String[] args) throws IOException{

	double alpha = 1.8;
	if(args.length==0){
	    System.out.println("Usage: DrawAtmMuonBundle alpha");
	    System.exit(0);
        }else{
	    alpha = Double.valueOf(args[0]).doubleValue();
	}
	// ElbertFlux calculator
	RelativeElbertFluxTableMaker muonFluxTable = new RelativeElbertFluxTableMaker();
	muonFluxTable.setAlpha(alpha);

	String realDataFileName = "iceCube/uhe/analysis/EHERealI3Particles";

	String mcSlope1DataFileName = 
	        "iceCube/uhe/analysis/EHEMCReducedIC9Slope1_mtx_bundleflux_I3Particles";

	I3ParticleAnalysisFactory.minNDOMsToAnalize = 75;

	//MC data Analyzer (E**-1)
	InputStream in = ClassLoader.getSystemResourceAsStream(mcSlope1DataFileName);
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
	//criteria.setThresholdOfLogNpe(4.6);
	//criteria.setRangeOfCosineOfZenith(0.1,1.0,4.7,5.8);
	//criteria.setMinimumBound(0.1,4.6);
	//criteria.setThresholdOfNDOMs(80);

	criteria.setThresholdOfLogNpe(3.8);
	criteria.setThresholdOfNDOMs(80);
	dataAnalizer.setCriteria(criteria);

	Criteria criteriaSlope1 = new Criteria();
	//criteriaSlope1.setThresholdOfLogNpe(4.6);
	//criteriaSlope1.setRangeOfCosineOfZenith(0.1,1.0,4.7,5.8);
	//criteriaSlope1.setMinimumBound(0.1,4.6);
	//criteriaSlope1.setThresholdOfNDOMs(80);

	criteriaSlope1.setThresholdOfLogNpe(3.8);
	criteriaSlope1.setThresholdOfNDOMs(80);
	//criteriaSlope1.setSimpleNpeCut(4.2,true); // logNpe >= 4.2
	mcSlope1Analizer.setCriteria(criteriaSlope1);


	mcSlope1Analizer.minLogNpe = 4.0;
	mcSlope1Analizer.maxLogNpe = 5.0;
	//mcSlope1Analizer.maxLogNpe = 6.0;
	mcSlope1Analizer.setBinSize(0.05,0.05,0.05,0.5); 
	                // Delta (LogE) Delta(LogNpe), Delta cos(Zenith), DeltaFGquality

                	// a half of the observation time must be set for averaing with E**-1 data
	dataAnalizer.minLogNpe = 4.0;
	dataAnalizer.maxLogNpe = 5.0;
	//dataAnalizer.maxLogNpe = 6.0;
	dataAnalizer.setBinSize(0.05,0.05,0.05,0.5); 

	/// Jaida FreeHep objects
	IAnalysisFactory jaidaFactory = IAnalysisFactory.create();
	ITree jaidaTree = jaidaFactory.createTreeFactory().create();
	IHistogramFactory jaidaHistoFactory = 
	    jaidaFactory.createHistogramFactory(jaidaTree);
	IFunctionFactory jaidaFuncFactory = null;
	IDataPointSetFactory jaidaDpsFactory = 
	    jaidaFactory.createDataPointSetFactory(jaidaTree);

	AtmMuonBundleFitter.setAtmMuonBundleFlux(mcSlope1Analizer,"Elbert_2_04_3730",
						 muonFluxTable);

	mcSlope1Analizer.drawEventsWithAtmMuonWeights(AtmMuonBundleFitter.fluxWeightName);
	//mcSlope1Analizer.switchToMCTruth();
	mcSlope1Analizer.switchToReco();

	dataAnalizer.switchToReco();

	IHistogram1D h1LogNpeMC = 
	    mcSlope1Analizer.makeJaida1DHistogram("logNpe",jaidaHistoFactory);

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
	    mcSlope1Analizer.makeJaida1DHistogram("cosZenith",jaidaHistoFactory);

	IHistogram1D h1cosZenithData = 
		dataAnalizer.makeJaida1DHistogram("cosZenith",jaidaHistoFactory);
	IDataPointSet d1cosZenithData = jaidaDpsFactory.create("dps1DFromHist",
							       h1cosZenithData);
	plotter.region(1).plot(h1cosZenithMC);
	plotter.region(1).plot(d1cosZenithData);
	plotter.show();

	int nbin = h1LogNpeData.axis().bins();
	double numberOfData = 0.0;
	double numberOfMC = 0.0;
	for(int ix = 0; ix< nbin; ix++){
	    numberOfData += h1LogNpeData.binHeight(ix);
	    numberOfMC += h1LogNpeMC.binHeight(ix);
	}

	System.out.println("Number of data(" + numberOfData + ")" +
			   "Number of MC(" + numberOfMC + ")");


    }


}
