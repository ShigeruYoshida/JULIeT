package iceCube.uhe.analysis;

import iceCube.uhe.analysis.*;
import geometry.*;

import hep.aida.*;
import hep.aida.ext.*;
import hep.aida.util.*;
import hep.aida.util.comparison.*;

import java.io.*;
import java.util.*;

public class DrawBGSpectrumByJaidaHisto {

    public static void main(String[] args) throws IOException{

	String aidaFileName = null;
	String path2BootStrapHisto1 = null;
	String path2BootStrapHisto2 = null;
	String[] path2Histo = null;
	boolean npeSpectrum = false;
	if(args.length<3){
	    System.err.println(
   "DrawBGSpectrumByJaidaHisto aida-file-name flag(1 for NPE, 2 for zenith) " +
	 " path-to-histo1, path-to-histo2, ...");
	    System.exit(0);
	}else{
	    aidaFileName = args[0];
	    if(Integer.valueOf(args[1]).intValue()==1) npeSpectrum = true;
	    path2Histo = new String[args.length-2];
	    for(int i= 2; i<args.length ; i++){
		path2Histo[i-2] = args[i];
	    }
	}

	System.err.println("Reading " + aidaFileName);
	for(int i=0; i< path2Histo.length; i++){
	    System.err.println("Drawing histogram of " + 
			       path2Histo[i]);
	}

	// Jaida FreeHep objects
	IAnalysisFactory jaidaFactory = IAnalysisFactory.create();
	IPlotterFactory plotterFactory = jaidaFactory.createPlotterFactory();
	IPlotter plotter = plotterFactory.create("");
	IPlotterStyle plotStyle = plotter.region(0).style();
	if(npeSpectrum){
	    System.err.println(" -- Npe spectrum");
	    JaidaPlotStyleSetter.setPlotStyle(plotStyle,
					      "log(Npe) ","Number of Events");
	}else{
	    System.err.println(" -- Zenith spectrum");
	    JaidaPlotStyleSetter.setPlotStyle(plotStyle,
				      "cos(ZenithAngle) ","Number of Events");
	}
	ITree jaidaTree = jaidaFactory.createTreeFactory().
	    createTree(aidaFileName,"xml",ITreeFactory.READONLY);
	//jaidaTree.ls(".", true, System.out);

	// Take the histograms in the aida Tree
	for(int i=0; i< path2Histo.length; i++){
	    IHistogram1D h1 = (IHistogram1D )jaidaTree.find(path2Histo[i]);
	    plotter.region(0).plot(h1);
	}

	plotter.show();


    }


}
