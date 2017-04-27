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

public class DrawAtmMuonBundleByElbertFluxTable {

    public static void main(String[] args) throws IOException{

	double alpha = 1.8;
	double muEThreshold = 1.0e2; // [GeV]
	boolean muEset = false;
	if(args.length==0){
	    System.out.println("Usage: DrawAtmMuonBundle alpha muonEth [GeV]");
	    System.exit(0);
        }else if(args.length==1){
	    alpha = Double.valueOf(args[0]).doubleValue();
        }else{
	    alpha = Double.valueOf(args[0]).doubleValue();
	    muEThreshold = Double.valueOf(args[1]).doubleValue();
	    muEset = true;
	}
	// ElbertFlux calculator
	if(!muEset){ // Eth fixed mode
	    // ElbertFlux calculator
	    ElbertFluxTableFactory.tablePath="../data/muonBundleFluxTable2TeVthreshold/";
	    ElbertFluxTableFactory.alphaMin= 1.8;
	    ElbertFluxTableFactory.numberOfAlphaSteps = 56;
	    ElbertFluxTableFactory.maxNumberOfMuEThSteps = 1;
	}
	ElbertFluxTableFactory muonFluxTable = new ElbertFluxTableFactory();
	if(!muEset){ // Eth fixed mode
	    muonFluxTable.setElbertParameters(alpha);
	}else{
	    muonFluxTable.setElbertParameters(alpha,muEThreshold);
	}

	String realDataFileName = "iceCube/uhe/analysis/EHERealI3Particles";
	String mcSlope2DataFileName = 
	        "iceCube/uhe/analysis/EHEMCReducedIC9Slope2I3Particles";
	//"iceCube/uhe/analysis/EHEMCIC22Slope2Mu_I3Particles";
	String mcSlope1DataFileName = 
	        "iceCube/uhe/analysis/EHEMCReducedIC9Slope1I3Particles";
	//"iceCube/uhe/analysis/EHEMCIC22Slope1Mu_I3Particles";

	I3ParticleAnalysisFactory.minNDOMsToAnalize = 60;

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
	dataAnalizer.setCriteria(criteria);

	Criteria criteriaSlope2 = new Criteria();
	criteriaSlope2.setThresholdOfLogNpe(3.8);
	criteriaSlope2.setThresholdOfNDOMs(80);
	criteriaSlope2.setSimpleNpeCut(4.2,false); // logNpe < 4.2
	mcSlope2Analizer.setCriteria(criteriaSlope2);

	Criteria criteriaSlope1 = new Criteria();
	criteriaSlope1.setThresholdOfLogNpe(3.8);
	criteriaSlope1.setThresholdOfNDOMs(80);
	criteriaSlope1.setSimpleNpeCut(4.2,true); // logNpe >= 4.2
	mcSlope1Analizer.setCriteria(criteriaSlope1);

	// Set the Anaizer's bin size and range
	mcSlope2Analizer.minLogNpe = 4.0;
	mcSlope2Analizer.maxLogNpe = 5.0;
	mcSlope2Analizer.setBinSize(0.05,0.05,0.05,0.5); 
	                // Delta (LogE) Delta(LogNpe), Delta cos(Zenith), DeltaFGquality

	mcSlope1Analizer.minLogNpe = 4.0;
	mcSlope1Analizer.maxLogNpe = 5.0;
	mcSlope1Analizer.setBinSize(0.05,0.05,0.05,0.5); 
	                // Delta (LogE) Delta(LogNpe), Delta cos(Zenith), DeltaFGquality

                	// a half of the observation time must be set for averaing with E**-1 data
	dataAnalizer.minLogNpe = 4.0;
	dataAnalizer.maxLogNpe = 5.0;
	dataAnalizer.setBinSize(0.05,0.05,0.05,0.5); 

	/// Jaida FreeHep objects
	IAnalysisFactory jaidaFactory = IAnalysisFactory.create();
	ITree jaidaTree = jaidaFactory.createTreeFactory().create();
	IHistogramFactory jaidaHistoFactory = 
	    jaidaFactory.createHistogramFactory(jaidaTree);
	IFunctionFactory jaidaFuncFactory = null;
	IDataPointSetFactory jaidaDpsFactory = 
	    jaidaFactory.createDataPointSetFactory(jaidaTree);

	AtmMuonBundleFitter.setAtmMuonBundleFlux(mcSlope2Analizer,
						 muonFluxTable);
	AtmMuonBundleFitter.setAtmMuonBundleFlux(mcSlope1Analizer,
						 muonFluxTable);

	mcSlope2Analizer.drawEventsWithAtmMuonWeights(AtmMuonBundleFitter.fluxWeightName);
	mcSlope2Analizer.switchToReco();
	mcSlope1Analizer.drawEventsWithAtmMuonWeights(AtmMuonBundleFitter.fluxWeightName);
	mcSlope1Analizer.switchToReco();

	dataAnalizer.switchToReco();

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
