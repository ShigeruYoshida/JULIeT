package iceCube.uhe.analysis;

import iceCube.uhe.analysis.*;
import geometry.*;

import hep.aida.*;
import hep.aida.ext.*;
import hep.aida.util.*;
import hep.aida.util.comparison.*;

import java.io.*;
import java.util.*;

public class PlotSpectrumByJaidaHisto {

    private static final double ln10 = Math.log(10.0);

    public static void main(String[] args) throws IOException{

	String aidaFileName = null;
	List path2Histo1DNPEList = null;
	List path2Histo1DZenithList = null;
	List path2Histo2D = null;
	boolean plot2D = false;
	if(args.length<1){
	    System.err.println("PlotSpectrumByJaidaHisto aida-file-name");
	    System.exit(0);
	}else{
	    aidaFileName = args[0];
	}
	System.err.println("Reading " + aidaFileName);
 
        DataInputStream input = new DataInputStream(System.in); 
        BufferedReader  d     = new BufferedReader(new InputStreamReader(input));
	String buffer; 

        System.out.print("2D histogram [yes(1)/no(0)] ->"); 
        buffer   = d.readLine(); 
        if(Integer.valueOf(buffer).intValue()==1) plot2D = true; 
 
	if(!plot2D){
	    path2Histo1DNPEList = new LinkedList();
	    path2Histo1DZenithList = new LinkedList();
	}else{
	    path2Histo2D = new LinkedList();
	}

	//
	// Reading the histogram name
	//
	if(!plot2D){ // 1D plot
	    int endMark = 0;
	    while(endMark == 0){
		System.out.print("Path to the NPE histogram ->"); 
		buffer   = d.readLine(); 
		path2Histo1DNPEList.add(new String(buffer));

		System.out.print("end [yes(1)/no(0)] ->"); 
		buffer   = d.readLine(); 
		endMark=Integer.valueOf(buffer).intValue();
	    }

	    endMark = 0;
	    while(endMark == 0){
		System.out.print("Path to the Zenith histogram ->"); 
		buffer   = d.readLine(); 
		path2Histo1DZenithList.add(new String(buffer));

		System.out.print("end [yes(1)/no(0)] ->"); 
		buffer   = d.readLine(); 
		endMark=Integer.valueOf(buffer).intValue();
	    }

	    ListIterator path2HistoIter = path2Histo1DNPEList.listIterator();
	    while(path2HistoIter.hasNext()){
		String pathName = (String )path2HistoIter.next();
		System.out.println("NPE histogram " + pathName);
	    }
	    path2HistoIter = path2Histo1DZenithList.listIterator();
	    while(path2HistoIter.hasNext()){
		String pathName = (String )path2HistoIter.next();
		System.out.println("Zenith histogram " + pathName);
	    }
	}else{ // 2D spectrum
	    int endMark = 0;
	    while(endMark == 0){
		System.out.print("Path to the 2D histogram ->"); 
		buffer   = d.readLine(); 
		path2Histo2D.add(new String(buffer));

		System.out.print("end [yes(1)/no(0)] ->"); 
		buffer   = d.readLine(); 
		endMark=Integer.valueOf(buffer).intValue();
	    }
	    ListIterator path2HistoIter = path2Histo2D.listIterator();
	    while(path2HistoIter.hasNext()){
		String pathName = (String )path2HistoIter.next();
		System.out.println("2D histogram " + pathName);
	    }
	}

	// Jaida FreeHep objects
	IAnalysisFactory jaidaFactory = IAnalysisFactory.create();
	ITree jaidaTree = jaidaFactory.createTreeFactory().
	    createTree(aidaFileName,"xml",ITreeFactory.READONLY);
	jaidaTree.ls(".", true, System.out);

	//IHistogramFactory jaidaHistoFactory = 
	//    jaidaFactory.createHistogramFactory(jaidaTree);


	//
	// plotting
	//
	if(!plot2D){ // 1D plot

	    System.err.println("now 1D plotting");
	    IPlotterFactory plotterFactory = jaidaFactory.createPlotterFactory();
	    IPlotter plotter = plotterFactory.create("Npe and Cos(Zenith) distribution");
	    plotter.destroyRegions();
	    plotter.createRegion(0,0,0.66,1);
	    plotter.createRegions(2,1);
	    IPlotterStyle npeStyle = plotter.region(0).style();
	    JaidaPlotStyleSetter.setPlotStyle(npeStyle,
					  "Log(Npe)","Number of Events");

	    ListIterator path2HistoIter = path2Histo1DNPEList.listIterator();
	    while(path2HistoIter.hasNext()){
		String pathName = (String )path2HistoIter.next();
		IHistogram1D h1 = (IHistogram1D )jaidaTree.find(pathName);
		plotter.region(0).plot(h1);
	    }

	    IPlotterStyle cosStyle = plotter.region(1).style();
	    JaidaPlotStyleSetter.setPlotStyle(cosStyle,
					      "cos(Zenith Angle) ","Number of Events");
	    path2HistoIter = path2Histo1DZenithList.listIterator();
	    while(path2HistoIter.hasNext()){
		String pathName = (String )path2HistoIter.next();
		IHistogram1D h1 = (IHistogram1D )jaidaTree.find(pathName);
		plotter.region(1).plot(h1);
	    }
	    plotter.show();

	    path2HistoIter = path2Histo1DNPEList.listIterator();
	    while(path2HistoIter.hasNext()){
		String pathName = (String )path2HistoIter.next();
		IHistogram1D h1 = (IHistogram1D )jaidaTree.find(pathName);
		int nbin = h1.axis().bins();
		double sigRate = 0.0;
		for(int ix = 0; ix< nbin; ix++){
		    sigRate += h1.binHeight(ix);
		}

		System.out.println("Rate = " + sigRate + " " + pathName);
	    }
    
	}else{ // 2D plot
	    IPlotterFactory plotterFactory = jaidaFactory.createPlotterFactory();
	    IPlotter plotter = plotterFactory.create("Npe - Zenith");
	    IPlotterStyle style = plotter.region(0).style();
	    JaidaPlotStyleSetter.setPlotStyle(style,
	    	      "Log(Npe) ",
	    	      "cos(First-Guessed Zenith Angle) ",
	    	      "hist2DStyle","colorMap");

	    ListIterator path2HistoIter = path2Histo2D.listIterator();
	    while(path2HistoIter.hasNext()){
		String pathName = (String )path2HistoIter.next();
		IHistogram2D h2 = (IHistogram2D )jaidaTree.find(pathName);
		plotter.region(0).plot(h2);
	    }

	    plotter.show();

	}

    }


}
