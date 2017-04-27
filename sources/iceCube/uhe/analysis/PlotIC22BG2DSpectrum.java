package iceCube.uhe.analysis;

import iceCube.uhe.analysis.*;
import geometry.*;

import hep.aida.*;

import java.io.*;
import java.util.*;

public class PlotIC22BG2DSpectrum {

    public static void main(String[] args) throws IOException{

        String mcMuonSlope1DataFileName = 
	    "iceCube/uhe/analysis/EHEMCIC22Slope1Mu_mtx_allfluxes_I3Particles";
        String mcMuonSlope2DataFileName = 
	    "iceCube/uhe/analysis/EHEMCIC22Slope2Mu_atmFlux_I3Particles";

	double logNpeBoundary = 4.4; // slope1 > 4.4 : 4.4 > slope2
	double cosZenLow = -1.0;
	double cosZenUp = 1.0;

	if(args.length<1){
	    System.err.println(
       "PlotIC22BG2DSpectrum logNpeBoundary");
	    System.exit(0);
	}else{
	    logNpeBoundary = Double.valueOf(args[0]).doubleValue();
	}
	System.err.println("Connect slope1 and 2 at logNPE=" +
			   logNpeBoundary);
	System.err.println(cosZenLow + "< cosZenith < " + cosZenUp);

	String aidaFileHeaderName = "IC22BG2DSpectrum_" + logNpeBoundary;
	String aidaFileName = 
	    aidaFileHeaderName.concat(".aida");
	System.err.println("Output file " + aidaFileName);


	I3ParticleAnalysisFactory.minNDOMsToAnalize = 80; // 80 DOMs cut
	I3ParticleAnalysisFactory.minLogNPEToAnalize = 4.0; // NPE cut

	// Atm Muon MC slope1
	InputStream in = ClassLoader.getSystemResourceAsStream(mcMuonSlope1DataFileName);
	I3ParticleIC22AnalysisFactory analizerBGSlope1 = 
	    new I3ParticleIC22AnalysisFactory(in);
	in.close();

	// Atm Muon MC slope2
	in = ClassLoader.getSystemResourceAsStream(mcMuonSlope2DataFileName);
	I3ParticleIC22AnalysisFactory analizerBGSlope2 = 
	    new I3ParticleIC22AnalysisFactory(in);
	in.close();

	// livetime enhancement x 1000
	double livetime = analizerBGSlope1.observationTime;
	analizerBGSlope1.setObservationTime(livetime*1.0e3);
	analizerBGSlope2.setObservationTime(livetime*1.0e3);


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

	analizerBGSlope1.setCriteria(criteriaBGSlope1);
	analizerBGSlope2.setCriteria(criteriaBGSlope2);

	// Set the Anaizer's bin size and range
	analizerBGSlope1.maxLogNpe = 7.0;analizerBGSlope2.maxLogNpe = 7.0;
	analizerBGSlope1.minLogNpe = 4.0;analizerBGSlope2.minLogNpe = 4.0;

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

	analizerBGSlope1.switchToReco();
	analizerBGSlope1.drawEventsWithAtmMuonWeights("elbert_1_97_1505");
	analizerBGSlope2.switchToReco();
	analizerBGSlope2.drawEventsWithAtmMuonWeights("elbert_1_97_1505");

	//
	// Data Histograms
	//
	jaidaTree.mkdir("/NpeZenith");
	jaidaTree.cd("/NpeZenith");
	IHistogram2D h2BGSlope1 = 
	    analizerBGSlope1.makeJaida2DHistogram("logNpe-CosZenith","MC Slope1",
						  false,jaidaHistoFactory);
	IHistogram2D h2BGSlope2 = 
	    analizerBGSlope2.makeJaida2DHistogram("logNpe-CosZenith","MC Slope2",
						  false,jaidaHistoFactory);
	IHistogram2D h2BGCombined = 
	    jaidaHistoFactory.add("Background MC",h2BGSlope1,h2BGSlope2);

	jaidaTree.commit();

    }


}
