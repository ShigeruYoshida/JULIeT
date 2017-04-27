package iceCube.uhe.analysis;

import iceCube.uhe.analysis.*;
import geometry.*;

import hep.aida.*;

import java.io.*;
import java.util.*;

public class PlotIC22SIG2DSpectrum {

    public static void main(String[] args) throws IOException{

	String realDataFileName = "iceCube/uhe/analysis/EHEIC22Real_NPE3k_I3Particles";
	String mcNuMuDataFileName = "iceCube/uhe/analysis/EHEMCIC22Slope1NuMu_mtx_flux_I3Particles";
	String mcNuEDataFileName = "iceCube/uhe/analysis/EHEMCIC22Slope1NuE_mtx_flux_I3Particles";
	String mcNuTauDataFileName = "iceCube/uhe/analysis/EHEMCIC22Slope1NuTau_mtx_flux_I3Particles";
        String mcMuonDataFileName = "iceCube/uhe/analysis/EHEMCIC22Slope1Mu_mtx_allfluxes_I3Particles";
        String mcTauDataFileName = "iceCube/uhe/analysis/EHEMCIC22Slope1TAU_mtx_flux_I3Particles";

	String[] weightName = {"GZK_YT_4_4","GZK_sigl","TD_SUSY"};

	double cosZenLow = -1.0;
	double cosZenUp = 1.0;

	int fluxIndex = 0;
	if(args.length<1){
	    System.err.println("PlotIC22SIG2DSpectrum fluxIndex(0 GZK-YT, 1 GZK-SIgl, 2 TD-SUSY)");
	    System.exit(0);
	}else{
	    fluxIndex = Integer.valueOf(args[0]).intValue();
	    if(fluxIndex>weightName.length){
		System.err.println("fluxIndes out of range!");
		System.exit(0);
	    }
	}
	String aidaFileHeaderName = "IC22SIG2DSpectrum_" + weightName[fluxIndex];
	String aidaFileName = 
	    aidaFileHeaderName.concat(".aida");
	System.err.println("Output file " + aidaFileName);


	I3ParticleAnalysisFactory.minNDOMsToAnalize = 80; // 80 DOMs cut
	I3ParticleAnalysisFactory.minLogNPEToAnalize = 4.0; // NPE cut

	InputStream in = ClassLoader.getSystemResourceAsStream(mcMuonDataFileName);
	I3ParticleIC22AnalysisFactory analizerSIG = new I3ParticleIC22AnalysisFactory(in);
	in.close();
	in = ClassLoader.getSystemResourceAsStream(mcTauDataFileName);
	analizerSIG.readI3Particles(in);
	in.close();
	in = ClassLoader.getSystemResourceAsStream(mcNuEDataFileName);
	analizerSIG.readI3Particles(in);
	in.close();
	in = ClassLoader.getSystemResourceAsStream(mcNuMuDataFileName);
	analizerSIG.readI3Particles(in);
	in.close();
	in = ClassLoader.getSystemResourceAsStream(mcNuTauDataFileName);
	analizerSIG.readI3Particles(in);
	in.close();

	// livetime enhancement x 1000
	double livetime = analizerSIG.observationTime;
	analizerSIG.setObservationTime(livetime*1.0e3);


	// Sets the Criteria for MC slope 1
	// logNPE > logNpeBoundady
	//Criteria criteria = new Criteria();
	CriteriaIC22 criteria = new CriteriaIC22();
	//criteria.setThresholdOfLogNpe(4.0);
	//criteria.setThresholdOfNDOMs(80);
	//criteria.setRangeOfCosineOfZenith(1.0,-1.0,1.0,1.0);// block
	//criteria.setMinimumBound(cosZenUp,1.0); // cosZen < cosUp logNpe>1
	//criteria.setSimpleCosZenithCut(cosZenLow,false);// cosZen >cosLow
	criteria.setEHESuperCut();

	analizerSIG.setCriteria(criteria);

	// Set the Anaizer's bin size and range
	analizerSIG.maxLogNpe = 7.0;analizerSIG.minLogNpe = 4.0;

	analizerSIG.setBinSize(0.05,0.05,0.05,0.1); 
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

	analizerSIG.switchToReco();
	analizerSIG.drawEventsWithGZKWeights(weightName[fluxIndex]);

	//
	// Data Histograms
	//
	jaidaTree.mkdir("/NpeZenith");
	jaidaTree.cd("/NpeZenith");
	IHistogram2D h2SIG = 
	    analizerSIG.makeJaida2DHistogram("logNpe-CosZenith",weightName[fluxIndex],
						  false,jaidaHistoFactory);
	jaidaTree.commit();

    }


}
