package iceCube.uhe.analysis;

import iceCube.uhe.analysis.*;
import geometry.*;

import hep.aida.*;
import hep.aida.ext.*;
import hep.aida.util.*;
import hep.aida.util.comparison.*;

import java.io.*;
import java.util.*;

public class PlotChi2ByJaidaHisto {

    public static void main(String[] args) throws IOException{

	String aidaFileName = null;
	String path2BootStrapHisto1 = null;
	String path2BootStrapHisto2 = null;
	String path2Histo1 = null;
	String path2Histo2 = null;
	int itrial = 0;
	boolean plotData = false;
	if(args.length<4){
	    System.err.println(
   "PlotChi2ByJaidaHisto aida-file-name path-to-MChisto1 path-to-Datahisto2 " +
	 " trial-number (path-to-histo3) (path-to-histo4)");
	    System.exit(0);
	}else{
	    aidaFileName = args[0];
	    path2BootStrapHisto1 = args[1];
	    path2BootStrapHisto2 = args[2];
	    itrial = Integer.valueOf(args[3]).intValue();
	    if(args.length == 6){
		plotData = true;
		path2Histo1 = args[4];
		path2Histo2 = args[5];
	    }
	}
	System.err.println("Reading " + aidaFileName);
	System.err.println("Chi2 dist for the bootstrapping " +
			   path2BootStrapHisto1 + " and " +
			   path2BootStrapHisto2);
	if(plotData){
	    System.err.println("Also chi value for " + 
			   path2Histo1 + " and " +
			   path2Histo2 + " are calculated");
	}


	// Jaida FreeHep objects
	IAnalysisFactory jaidaFactory = IAnalysisFactory.create();
	ITree jaidaTree = jaidaFactory.createTreeFactory().
	    createTree(aidaFileName,"xml",ITreeFactory.READONLY);
	//jaidaTree.ls(".", true, System.out);

	IHistogramFactory jaidaHistoFactory = 
	    jaidaFactory.createHistogramFactory(jaidaTree);
	int dimensionX = 100;
	double minX = 0.0; double maxX = 10.0;
	IHistogram1D chi2Histo = 
	    jaidaHistoFactory.createHistogram1D("Bootstrap",
						dimensionX,minX,maxX);

	double logNpeMin = 4.0;
	double logNpeMax = 5.0;
	double numberChi2LessThan1 = 0.0;
	double numberChi2LessThan2 = 0.0;

	// Take the bootstrapping histos and calculate chi2
	for(int i = 0; i< itrial; i++){
	    String histo1Name = path2BootStrapHisto1 + " " + i;
	    String histo2Name = path2BootStrapHisto2 + " " + i;

	    IHistogram1D h1 = (IHistogram1D )jaidaTree.find(histo1Name);
	    IHistogram1D h2 = (IHistogram1D )jaidaTree.find(histo2Name);

	    IComparisonResult chi2 = 
		I3ParticleAnalysisFactory.calcChi2(h1,h2,logNpeMin,logNpeMax);
	    //I3ParticleAnalysisFactory.calcChi2(h1,h2);
	    if(chi2.quality()/chi2.nDof()<=1.0)numberChi2LessThan1 += 1.0;
	    if(chi2.quality()/chi2.nDof()<=2.0)numberChi2LessThan2 += 1.0;


	    chi2Histo.fill(chi2.quality()/chi2.nDof());
	    //if(i%100 ==0) System.err.println("  taking " + i + " th data");
	}

	IHistogram1D chi2DataHisto = 
	    jaidaHistoFactory.createHistogram1D("Data-MC",
						dimensionX,minX,maxX);
	if(plotData){
	    IHistogram1D h1 = (IHistogram1D )jaidaTree.find(path2Histo1);
	    IHistogram1D h2 = (IHistogram1D )jaidaTree.find(path2Histo2);

	    IComparisonResult chi2 = 
		I3ParticleAnalysisFactory.calcChi2(h1,h2,logNpeMin,logNpeMax);
	    //I3ParticleAnalysisFactory.calcChi2(h1,h2);

	    chi2DataHisto.fill(chi2.quality()/chi2.nDof(),
			       chi2Histo.maxBinHeight());
	}

	IPlotterFactory plotterFactory = jaidaFactory.createPlotterFactory();
	IPlotter plotter = plotterFactory.create("Chi2 distribution");

	IPlotterStyle chi2Style = plotter.region(0).style();
	JaidaPlotStyleSetter.setPlotStyle(chi2Style,
					  "chi2","Number of Events");
	plotter.region(0).plot(chi2Histo);
	if(plotData) plotter.region(0).plot(chi2DataHisto);


	plotter.show();

	System.out.println("trial " + itrial + 
			   " chi2<1 " + (int )(numberChi2LessThan1+0.01) + 
			   " prob=" + numberChi2LessThan1/((double ) itrial) +
			   " chi2<2 " + (int )(numberChi2LessThan2+0.01) + 
			   " prob=" + numberChi2LessThan2/((double ) itrial));


    }


}
