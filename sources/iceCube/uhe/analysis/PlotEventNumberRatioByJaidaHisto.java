package iceCube.uhe.analysis;

import iceCube.uhe.analysis.*;
import geometry.*;

import hep.aida.*;
import hep.aida.ext.*;
import hep.aida.util.*;
import hep.aida.util.comparison.*;

import java.io.*;
import java.util.*;

public class PlotEventNumberRatioByJaidaHisto {

    private static final double ln10 = Math.log(10.0);

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
	       "PlotEventNumberRatioByJaidaHisto aida-file-name " +
	       "path-to-MChisto1 path-to-Datahisto2 " +
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
	System.err.println("Event Ratio dist for the bootstrapping " +
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
	int dimensionX = 300;
	double minX = -1.0; double maxX = 1.0; // -1.0 <log(Data/MC) <1.0
	IHistogram1D ratioHisto = 
	    jaidaHistoFactory.createHistogram1D("Bootstrap",
						dimensionX,minX,maxX);

	double logNpeMin = 4.0;
	double logNpeMax = 5.0;

	// Take the bootstrapping histos and calculate chi2
	for(int i = 0; i< itrial; i++){
	    String histo1Name = path2BootStrapHisto1 + " " + i;
	    String histo2Name = path2BootStrapHisto2 + " " + i;

	    IHistogram1D h1 = (IHistogram1D )jaidaTree.find(histo1Name);
	    IHistogram1D h2 = (IHistogram1D )jaidaTree.find(histo2Name);

	    //double eventNumberInH1 = 0.0;
	    double eventNumberInH1 = h1.sumAllBinHeights();
	    //double eventNumberInH2 = 0.0;
	    double eventNumberInH2 = h2.sumAllBinHeights();
	    //int firstBin = h1.coordToIndex(logNpeMin);
	    //int lastBin = h1.coordToIndex(logNpeMax);
	    //for(int index = 0; index <= lastBin; index++){
	    //eventNumberInH1 += h1.binHeight(index);
	    //	eventNumberInH2 += h2.binHeight(index);
	    //}
	    double logRatio = Math.log(eventNumberInH2/eventNumberInH1)/ln10;
	    ratioHisto.fill(logRatio);

	    //if(i%100 ==0) System.err.println("  taking " + i + " th data");
	}

	IHistogram1D ratioDataHisto = 
	    jaidaHistoFactory.createHistogram1D("Data-MC",
						dimensionX,minX,maxX);
	double logRatioData = 0.0;
	if(plotData){
	    IHistogram1D h1 = (IHistogram1D )jaidaTree.find(path2Histo1);
	    IHistogram1D h2 = (IHistogram1D )jaidaTree.find(path2Histo2);


	    //double eventNumberInH1 = 0.0;
	    double eventNumberInH1 = h1.sumAllBinHeights();
	    //double eventNumberInH2 = 0.0;
	    double eventNumberInH2 = h2.sumAllBinHeights();
	    //int firstBin = h1.coordToIndex(logNpeMin);
	    //int lastBin = h1.coordToIndex(logNpeMax);
	    //for(int index = 0; index <= lastBin; index++){
	    //	eventNumberInH1 += h1.binHeight(index);
	    //	eventNumberInH2 += h2.binHeight(index);
	    //}
	    logRatioData = Math.log(eventNumberInH2/eventNumberInH1)/ln10;
	    ratioDataHisto.fill(logRatioData,ratioHisto.maxBinHeight());
	}

	IPlotterFactory plotterFactory = jaidaFactory.createPlotterFactory();
	IPlotter plotter = plotterFactory.create("Ratio distribution");

	IPlotterStyle ratioStyle = plotter.region(0).style();
	JaidaPlotStyleSetter.setPlotStyle(ratioStyle,
					  "log(data/MC)","Number of Events");
	plotter.region(0).plot(ratioHisto);
	if(plotData) plotter.region(0).plot(ratioDataHisto);


	plotter.show();

	//
	// Calculating the error range
	//
	double totalNumber = ratioHisto.sumAllBinHeights();
	double width = 0.025*totalNumber; //95 % CL
	// lower limit
	double lowerNumber = 0.0;
	int index = 0;
	for(index = 0; index<= ratioHisto.axis().bins();index++){
	    lowerNumber += ratioHisto.binHeight(index);
	    if(lowerNumber>= width) break;
	}
	System.err.println(" lower index(" + index + ")");
	double lowerLogRatio = ratioHisto.axis().binCenter(index);
	// upper limit
	double upperNumber = 0.0;
	for(index= ratioHisto.axis().bins(); index >=0 ;index--){
	    upperNumber += ratioHisto.binHeight(index);
	    if(upperNumber>= width) break;
	}
	double upperLogRatio = ratioHisto.axis().binCenter(index);
	System.out.println("log10(ratio)=" + logRatioData +
			   " +" + upperLogRatio + 
			   " -" + lowerLogRatio);
	System.out.println(" total(" + totalNumber + ")");

    }


}
