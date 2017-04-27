package iceCube.uhe.analysis;

import iceCube.uhe.analysis.*;
import geometry.*;

import hep.aida.*;

import java.io.*;
import java.util.*;

public class MakeBootstrapBGSpectrum {

    public static void main(String[] args) throws IOException{

	int nBootstrap = 10000; // Bootstrapping 10000 times
	double logNpeBoundary = 4.4; // slope1 > 4.4 : 4.4 > slope2
	double cosZenLow = -1.0;
	double cosZenUp = 1.0;

	if(args.length<4){
	    System.err.println(
       "MakeBootstrapBGSpectrum nBootstrap logNpeBoundary cosZenLow cosZenUp");
	    System.exit(0);
	}else{
	    nBootstrap = Integer.valueOf(args[0]).intValue();
	    logNpeBoundary = Double.valueOf(args[1]).doubleValue();
	    cosZenLow = Double.valueOf(args[2]).doubleValue();
	    cosZenUp = Double.valueOf(args[3]).doubleValue();
	}
	System.err.println("Connect slope1 and 2 at logNPE=" +
			   logNpeBoundary);
	System.err.println("Bootstrapping " + nBootstrap + " times");
	System.err.println(cosZenLow + "< cosZenith < " + cosZenUp);

	String aidaFileHeaderName = "BootstrapBGSpectrum_" + logNpeBoundary;
	String cosRange = "_from" + cosZenUp + "to" + cosZenLow;
	String aidaFileName = 
	    aidaFileHeaderName.concat(cosRange).concat(".aida");
	System.err.println("Output file " + aidaFileName);

	String realDataFileName = "iceCube/uhe/analysis/EHEIC22Real_NPE3k_I3Particles";
        String mcMuonSlope1DataFileName = 
	    "iceCube/uhe/analysis/EHEMCIC22Slope1Mu_mtx_allfluxes_I3Particles";
        String mcMuonSlope2DataFileName = 
	    "iceCube/uhe/analysis/EHEMCIC22Slope2Mu_atmFlux_I3Particles";

	I3ParticleAnalysisFactory.minNDOMsToAnalize = 80; // 80 DOMs cut
	I3ParticleAnalysisFactory.minLogNPEToAnalize = 4.0; // NPE cut

	// Real data
	InputStream in = ClassLoader.getSystemResourceAsStream(realDataFileName);
	I3ParticleIC22AnalysisFactory analizerData = 
	    new I3ParticleIC22AnalysisFactory(in,true);

	// Atm Muon MC slope1
	in = ClassLoader.getSystemResourceAsStream(mcMuonSlope1DataFileName);
	I3ParticleIC22AnalysisFactory analizerBGSlope1 = 
	    new I3ParticleIC22AnalysisFactory(in);
	in.close();

	// Atm Muon MC slope2
	in = ClassLoader.getSystemResourceAsStream(mcMuonSlope2DataFileName);
	I3ParticleIC22AnalysisFactory analizerBGSlope2 = 
	    new I3ParticleIC22AnalysisFactory(in);
	in.close();


	// Sets the Criteria for data
	Criteria criteria = new Criteria();
	criteria.setThresholdOfLogNpe(4.0);
	criteria.setThresholdOfNDOMs(80);
	criteria.setRangeOfCosineOfZenith(1.0,-1.0,1.0,1.0);// block
	criteria.setMinimumBound(cosZenUp,1.0); // cosZen < cosUp logNpe>1
	criteria.setSimpleCosZenithCut(cosZenLow,false);// cosZen >cosLow

	// Sets the Criteria for MC slope 1
	// logNPE > logNpeBoundady
	Criteria criteriaBGSlope1 = new Criteria();
	criteriaBGSlope1.setThresholdOfLogNpe(logNpeBoundary);
	criteriaBGSlope1.setThresholdOfNDOMs(80);
	criteriaBGSlope1.setRangeOfCosineOfZenith(1.0,-1.0,1.0,1.0);// block
	criteriaBGSlope1.setMinimumBound(cosZenUp,1.0); // cosZen < cosUp logNpe>1
	criteriaBGSlope1.setSimpleCosZenithCut(cosZenLow,false);// cosZen >cosLow

	// Sets the Criteria for MC slope 2
	Criteria criteriaBGSlope2 = new Criteria();
	// logNpeBoundary> logNPE > 4.0
	criteriaBGSlope2.setThresholdOfLogNpe(4.0);
	criteriaBGSlope2.setSimpleNpeCut(logNpeBoundary,false);
	criteriaBGSlope2.setThresholdOfNDOMs(80);
	criteriaBGSlope2.setRangeOfCosineOfZenith(1.0,-1.0,1.0,1.0);// block
	criteriaBGSlope2.setMinimumBound(cosZenUp,1.0); // cosZen < cosUp logNpe>1
	criteriaBGSlope2.setSimpleCosZenithCut(cosZenLow,false);// cosZen >cosLow

	analizerData.setCriteria(criteria);
	analizerBGSlope1.setCriteria(criteriaBGSlope1);
	analizerBGSlope2.setCriteria(criteriaBGSlope2);

	// Set the Anaizer's bin size and range
	analizerData.maxLogNpe = 6.0;
	analizerBGSlope1.maxLogNpe = 6.0;analizerBGSlope2.maxLogNpe = 6.0;
	analizerData.minLogNpe = 4.0;
	analizerBGSlope1.minLogNpe = 4.0;analizerBGSlope2.minLogNpe = 4.0;

	analizerData.setBinSize(0.05,0.05,0.05,0.1); 
	analizerBGSlope1.setBinSize(0.05,0.05,0.05,0.1); 
	analizerBGSlope2.setBinSize(0.05,0.05,0.05,0.1); 
	                // Delta (LogE) Delta(LogNpe), Delta cos(Zenith), DeltaFGquality


	// Jaida FreeHep objects
	IAnalysisFactory jaidaFactory = IAnalysisFactory.create();
	ITree jaidaTree = jaidaFactory.createTreeFactory().createTree(aidaFileName,
							      "xml",ITreeFactory.RECREATE);
	IHistogramFactory jaidaHistoFactory = 
	    jaidaFactory.createHistogramFactory(jaidaTree);

	//
	//
	// Atmospheric spectrum weight
	//
	//
	analizerData.switchToReco();

	analizerBGSlope1.switchToReco();
	analizerBGSlope1.drawEventsWithAtmMuonWeights("elbert_1_97_1505");
	analizerBGSlope2.switchToReco();
	analizerBGSlope2.drawEventsWithAtmMuonWeights("elbert_1_97_1505");

	//
	// Data Histograms
	//
	jaidaTree.mkdir("/logNpe");
	jaidaTree.cd("/logNpe");
	IHistogram1D h1LogNpeData = 
	    analizerData.makeJaida1DHistogram("logNpe","Npe Data",
					      false,jaidaHistoFactory);
	jaidaTree.mkdir("/cosZenith");
	jaidaTree.cd("/cosZenith");
	IHistogram1D h1cosZenithData = 
	    analizerData.makeJaida1DHistogram("cosZenith","Zenith Data",
					      false,jaidaHistoFactory);

	//
	// BG MC Histograms
	//
	jaidaTree.cd("/logNpe");
	IHistogram1D h1LogNpeBGSlope1 = 
	    analizerBGSlope1.makeJaida1DHistogram("logNpe","Npe MC Slope1",
						  false,jaidaHistoFactory);
	IHistogram1D h1LogNpeBGSlope2 = 
	    analizerBGSlope2.makeJaida1DHistogram("logNpe","Npe MC Slope2",
						  false,jaidaHistoFactory);

	IHistogram1D h1LogNpeBGCombined = 
	    jaidaHistoFactory.add("Npe MC",h1LogNpeBGSlope1,h1LogNpeBGSlope2);


	jaidaTree.cd("/cosZenith");
	IHistogram1D h1cosZenithBGSlope1 = 
	    analizerBGSlope1.makeJaida1DHistogram("cosZenith","Zenith MC Slope1",
						  false,jaidaHistoFactory);
	IHistogram1D h1cosZenithBGSlope2 = 
	    analizerBGSlope2.makeJaida1DHistogram("cosZenith","Zenith MC Slope2",
						  false,jaidaHistoFactory);
	IHistogram1D h1cosZenithBGCombined = 
	    jaidaHistoFactory.add("Zenith MC",h1cosZenithBGSlope1,
				  h1cosZenithBGSlope2);

	//
	// Bootstrapping
	//
	//

	jaidaTree.mkdir("/logNpeBootstrap");
	//jaidaTree.mkdir("/cosZenithBootstrap");
	IHistogram1D h1LogNpeBGSlope1BS;
	IHistogram1D h1LogNpeBGSlope2BS;
	//IHistogram1D h1cosZenithBGSlope1BS;
	//IHistogram1D h1cosZenithBGSlope2BS;
	for(int itrial = 0; itrial<nBootstrap; itrial++){
	    jaidaTree.cd("/logNpeBootstrap");
	    String npeHistogramName = "Npe Data " + itrial;
	    IHistogram1D h1LogNpeDataBS = 
	    analizerData.makeJaida1DHistogram("logNpe",npeHistogramName,
					      true,jaidaHistoFactory);

	    String npeHistogramName1 = "Npe MC Slope1 " + itrial;
	    h1LogNpeBGSlope1BS = 
		analizerBGSlope1.makeJaida1DHistogram("logNpe",npeHistogramName1,
						  true,jaidaHistoFactory);
	    String npeHistogramName2 = "Npe MC Slope2 " + itrial;
	    h1LogNpeBGSlope2BS = 
	    analizerBGSlope2.makeJaida1DHistogram("logNpe",npeHistogramName2,
						  true,jaidaHistoFactory);

	    String npeHistogramNameCombined = "Npe MC " + itrial;
	    IHistogram1D h1LogNpeBGBSCombined = 
	    jaidaHistoFactory.add(npeHistogramNameCombined,
				  h1LogNpeBGSlope1BS,h1LogNpeBGSlope2BS);

	    //jaidaTree.cd("/cosZenithBootstrap");
	    //String zenithHistogramName = "Zenith Data " + itrial;
	    //IHistogram1D h1cosZenithDataBS = 
	    //analizerData.makeJaida1DHistogram("cosZenith",zenithHistogramName,
	    //			      true,jaidaHistoFactory);

	    //String zenithHistogramName1 = "Zenith MC Slope1 " + itrial;
	    //h1cosZenithBGSlope1BS = 
	    //analizerBGSlope1.makeJaida1DHistogram("cosZenith",zenithHistogramName1,
	    //				  true,jaidaHistoFactory);

	    //String zenithHistogramName2 = "Zenith MC Slope2 " + itrial;
	    //h1cosZenithBGSlope2BS = 
	    //analizerBGSlope2.makeJaida1DHistogram("cosZenith",zenithHistogramName2,
	    //				  true,jaidaHistoFactory);

	    //String zenithHistogramNameCombined = "Zenith MC " + itrial;
	    //IHistogram1D h1cosZenithBGBSCombined = 
	    //jaidaHistoFactory.add(zenithHistogramNameCombined,
	    //		  h1cosZenithBGSlope1BS,h1cosZenithBGSlope2BS);

	    if(itrial%500==0) System.err.println(" bootstrap " + itrial + " times");
	}

	jaidaTree.commit();

	//
	// plotting
	//
	//IPlotterFactory plotterFactory = jaidaFactory.createPlotterFactory();
	//IPlotter plotter = plotterFactory.create("Npe and Cos(Zenith) distribution");
	//plotter.destroyRegions();
	//plotter.createRegion(0,0,0.66,1);
	//plotter.createRegions(2,1);

	//IPlotterStyle npeStyle = plotter.region(0).style();
	//JaidaPlotStyleSetter.setPlotStyle(npeStyle,
	//				  "Log(Npe)","Number of Events");
	//plotter.region(0).plot(h1LogNpeData);
	//plotter.region(0).plot(h1LogNpeBGCombined);

	//IPlotterStyle cosStyle = plotter.region(1).style();
	//JaidaPlotStyleSetter.setPlotStyle(cosStyle,
	//			  "cos(ZenithAngle)","Number of Events");
	//plotter.region(1).plot(h1cosZenithData);
	//plotter.region(1).plot(h1cosZenithBGCombined);

	//plotter.show();
	//jaidaTree.cd("/");
	//jaidaTree.ls(".", true, System.out);

    }


}
