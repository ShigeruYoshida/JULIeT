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

public class DrawAtmMuonBundleFluxRatio {

    public static void main(String[] args) throws IOException{

	double alpha = 1.8;
	if(args.length==0){
	    System.out.println("Usage: DrawAtmMuonBundle alpha");
	    System.exit(0);
        }else{
	    alpha = Double.valueOf(args[0]).doubleValue();
	}
	// ElbertFlux calculator
	RelativeElbertFluxTableMaker relativeFluxTable = new RelativeElbertFluxTableMaker();
	relativeFluxTable.setAlpha(alpha);

	ElbertFluxTableFactory.tablePath="../data/muonBundleFluxTable2TeVthreshold/";
	ElbertFluxTableFactory.alphaMin= 1.8;
	ElbertFluxTableFactory.numberOfAlphaSteps = 56;
	ElbertFluxTableFactory.maxNumberOfMuEThSteps = 1;

	ElbertFluxTableFactory muonFluxTable = new ElbertFluxTableFactory();

	muonFluxTable.setElbertParameters(alpha);

	String realDataFileName = "iceCube/uhe/analysis/EHERealI3Particles";

	String mcSlope1DataFileName = 
	        "iceCube/uhe/analysis/EHEMCReducedIC9Slope1_mtx_bundleflux_I3Particles";

	I3ParticleAnalysisFactory.minNDOMsToAnalize = 75;

	//MC data Analyzer (E**-1)
	InputStream in = ClassLoader.getSystemResourceAsStream(mcSlope1DataFileName);
	I3ParticleAnalysisFactory mcSlope1Analizer = new I3ParticleAnalysisFactory(in);
	in.close();

	J3Vector ic9Center = new J3Vector(mcSlope1Analizer.xCenterOfIC9,
					  mcSlope1Analizer.yCenterOfIC9, 
					  0.0);
	// IceCube 9 string array center
	// Sets the Criteria

	Criteria criteriaSlope1 = new Criteria();
	//criteriaSlope1.setThresholdOfLogNpe(4.6);
	//criteriaSlope1.setRangeOfCosineOfZenith(0.1,1.0,4.7,5.8);
	//criteriaSlope1.setMinimumBound(0.1,4.6);
	//criteriaSlope1.setThresholdOfNDOMs(80);

	criteriaSlope1.setThresholdOfLogNpe(3.8);
	criteriaSlope1.setThresholdOfNDOMs(80);
	criteriaSlope1.setSimpleCosZenithCut(0.2,true);
	//criteriaSlope1.setSimpleNpeCut(4.2,true); // logNpe >= 4.2
	mcSlope1Analizer.setCriteria(criteriaSlope1);


	mcSlope1Analizer.minLogNpe = 4.0;
	//mcSlope1Analizer.maxLogNpe = 5.0;
	mcSlope1Analizer.maxLogNpe = 6.0;
	mcSlope1Analizer.minLogEnergy = -2.0;
	mcSlope1Analizer.maxLogEnergy = 2.0;
	mcSlope1Analizer.setBinSize(0.1,0.1,0.05,0.5); 
	                // Delta (LogE) Delta(LogNpe), Delta cos(Zenith), DeltaFGquality

	/// Jaida FreeHep objects
	IAnalysisFactory jaidaFactory = IAnalysisFactory.create();
	ITree jaidaTree = jaidaFactory.createTreeFactory().create();
	IHistogramFactory jaidaHistoFactory = 
	    jaidaFactory.createHistogramFactory(jaidaTree);
	IFunctionFactory jaidaFuncFactory = null;
	IDataPointSetFactory jaidaDpsFactory = 
	    jaidaFactory.createDataPointSetFactory(jaidaTree);

	AtmMuonBundleFitter.setAtmMuonBundleFluxRatio(mcSlope1Analizer,"Elbert_2_04_3730",
						      muonFluxTable,relativeFluxTable);

	mcSlope1Analizer.switchToMCTruth();
	//mcSlope1Analizer.switchToReco();


	IHistogram2D h1LogNpeMC = 
	    mcSlope1Analizer.makeJaida2DHistogram("logRecoE-Npe",jaidaHistoFactory);


	IPlotterFactory plotterFactory = jaidaFactory.createPlotterFactory();
	IPlotter plotter = plotterFactory.create("Event Distribution");
	//plotter.destroyRegions();
	//plotter.createRegion(0,0,0.66,1);
	//plotter.createRegions(2,1);

	IPlotterStyle npeStyle = plotter.region(0).style();
	JaidaPlotStyleSetter.setPlotStyle(npeStyle,
					      "Log(Ratio)", "Log(Npe)");
	plotter.region(0).plot(h1LogNpeMC);
	plotter.show();

    }


}
