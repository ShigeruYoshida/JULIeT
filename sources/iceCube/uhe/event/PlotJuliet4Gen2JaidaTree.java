package iceCube.uhe.event;

import numRecipes.*;
import geometry.*;
import iceCube.uhe.geometry.*;
import iceCube.uhe.particles.*;
import iceCube.uhe.analysis.*;

import java.io.*;
import java.util.*;

import hep.aida.*;

public class PlotJuliet4Gen2JaidaTree {


    public static void main(String[] args) throws IOException{

	String jaidaTreeFileName = null;
        if(args.length<1){
            System.err.println("PlotJuliet4Gen2JaidaTree aida-file-name");
            System.exit(0);
	}else{
            jaidaTreeFileName = args[0];
        }

        // Jaida FreeHep objects 
        IAnalysisFactory jaidaFactory = IAnalysisFactory.create();
	ITree jaidaJulietTree = jaidaFactory.createTreeFactory().
            createTree(jaidaTreeFileName,"xml",ITreeFactory.READONLY);

        String path2HistoCascade = "./cascadePosition";
        String path2HistoInjection = "./injectionPosition";
        String path2HistoEnd = "./endingPosition";
	String path2HistoDistanceNadir = "./propagationDistanceNadir";
	String path2HistoDistanceEnergy = "./propagationDistanceEnergy";
	// retrive the hostograms
	IHistogram3D cascadeXYZ = (IHistogram3D )jaidaJulietTree.find(path2HistoCascade);
	IHistogram3D injectionXYZ = (IHistogram3D )jaidaJulietTree.find(path2HistoInjection);
	IHistogram3D endXYZ = (IHistogram3D )jaidaJulietTree.find(path2HistoEnd);
	IHistogram2D distanceNadir = (IHistogram2D )jaidaJulietTree.find(path2HistoDistanceNadir);
	IHistogram2D distanceEnergy = (IHistogram2D )jaidaJulietTree.find(path2HistoDistanceEnergy);

	//
	// draw
	//
	IHistogramFactory jaidaHistoFactory =
            jaidaFactory.createHistogramFactory(jaidaJulietTree);
	IPlotterFactory plotterFactory = jaidaFactory.createPlotterFactory();

        IPlotter plotterCascadeXY = plotterFactory.create("Cascade vertex XY distribution");
	IPlotterStyle styleCascadeXY = plotterCascadeXY.region(0).style();
	JaidaPlotStyleSetter.setPlotStyle(styleCascadeXY,
					  "X[m] ",
					  "Y[m] ",
					  "hist2DStyle","colorMap");

        //plotter.destroyRegions();
	//plotter.createRegion(0,0,0.66,1);
	plotterCascadeXY.region(0).plot(jaidaHistoFactory.projectionXY("Energy Deposit XY projection",cascadeXYZ));
	plotterCascadeXY.show();

        IPlotter plotterCascadeXZ = plotterFactory.create("Cascade vertex XZ distribution");
	IPlotterStyle styleCascadeXZ = plotterCascadeXZ.region(0).style();
	JaidaPlotStyleSetter.setPlotStyle(styleCascadeXZ,
					  "X[m] ",
					  "Z[m] ",
					  "hist2DStyle","colorMap");

	plotterCascadeXZ.region(0).plot(jaidaHistoFactory.projectionXZ("Energy Deposit XZ projection",cascadeXYZ));
	plotterCascadeXZ.show();


        IPlotter plotterInjectionXY = plotterFactory.create("In-Ice injection vertex XY distribution");
	IPlotterStyle styleInjectionXY = plotterInjectionXY.region(0).style();
	JaidaPlotStyleSetter.setPlotStyle(styleInjectionXY,
					  "X[m] ",
					  "Y[m] ",
					  "hist2DStyle","colorMap");
	plotterInjectionXY.region(0).plot(jaidaHistoFactory.projectionXY("In-Ice injection XY projection",injectionXYZ));
	plotterInjectionXY.show();

        IPlotter plotterInjectionXZ = plotterFactory.create("In-Ice injection vertex XZ distribution");
	IPlotterStyle styleInjectionXZ = plotterInjectionXZ.region(0).style();
	JaidaPlotStyleSetter.setPlotStyle(styleInjectionXZ,
					  "X[m] ",
					  "Z[m] ",
					  "hist2DStyle","colorMap");

	plotterInjectionXZ.region(0).plot(jaidaHistoFactory.projectionXZ("In-Ice injection XZ projection",injectionXYZ));
	plotterInjectionXZ.show();

        IPlotter plotterEndXY = plotterFactory.create("In-Ice Ending vertex XY distribution");
	IPlotterStyle styleEndXY = plotterEndXY.region(0).style();
	JaidaPlotStyleSetter.setPlotStyle(styleEndXY,
					  "X[m] ",
					  "Y[m] ",
					  "hist2DStyle","colorMap");
	plotterEndXY.region(0).plot(jaidaHistoFactory.projectionXY("In-Ice track-ending XY projection",endXYZ));
	plotterEndXY.show();


        IPlotter plotterEndXZ = plotterFactory.create("In-Ice Ending vertex XZ distribution");
	IPlotterStyle styleEndXZ = plotterEndXZ.region(0).style();
	JaidaPlotStyleSetter.setPlotStyle(styleEndXZ,
					  "X[m] ",
					  "Z[m] ",
					  "hist2DStyle","colorMap");
	plotterEndXZ.region(0).plot(jaidaHistoFactory.projectionXZ("In-Ice track-ending XZ projection",endXYZ));
	plotterEndXZ.show();

        IPlotter plotterDistanceNadir = plotterFactory.create("Distance from the surface Vs nadir angle");
	IPlotterStyle styleDistanceNadir = plotterDistanceNadir.region(0).style();
	JaidaPlotStyleSetter.setPlotStyle(styleDistanceNadir,
					  "Distance[cm] ",
					  "Nadir [deg] ",
					  "hist2DStyle","colorMap");
	plotterDistanceNadir.region(0).plot(distanceNadir);
	plotterDistanceNadir.show();

        IPlotter plotterDistanceEnergy = plotterFactory.create("Distance from the surface Vs Energy Deposit");
	IPlotterStyle styleDistanceEnergy = plotterDistanceEnergy.region(0).style();
	JaidaPlotStyleSetter.setPlotStyle(styleDistanceEnergy,
					  "Distance[cm] ",
					  "Energy Deposit [ratio] ",
					  "hist2DStyle","colorMap");
	plotterDistanceEnergy.region(0).plot(distanceEnergy);
	plotterDistanceEnergy.show();

    }

}
