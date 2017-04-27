package numRecipes;

import numRecipes.*;

import hep.aida.*;
import hep.aida.ext.*;
import hep.aida.util.*;
import hep.aida.util.comparison.*;

public class CheckKSTest {

    /** A main method for test. The Gaussian is used as a model */
    public static void main(String[] args) {

	double mean = 0.0; double sigma = 1.0;
	int n = 10; // number of data points
	if(args.length!=3){
	    System.out.println("Usage: KSTest mean sigma n");
	    System.exit(0);
	}else {
	    mean = Double.valueOf(args[0]).doubleValue();
	    sigma = Double.valueOf(args[1]).doubleValue();
	    n = Integer.valueOf(args[2]).intValue();
	}
	System.out.format("Test with the Gaussian (mean=%f sigma=%f) and %d data points\n",mean,sigma,n);

	SpecialFunctions gaussFunc = new SpecialFunctions();
	double lowerBound = mean - 30.0*sigma;
	double upperBound = mean + 30.0*sigma;
	int functionIndex = 1; // choose the Gaussian in the interface of SpecialFunctions class
	double[] parameters =  new double[2];
	parameters[0] = mean;   // for the function interface
	parameters[1] = sigma;  // for the function interface

	// Jaida FreeHep objects
	IAnalysisFactory jaidaFactory = IAnalysisFactory.create();
	ITree jaidaTree = jaidaFactory.createTreeFactory().create();
	IHistogramFactory jaidaHistoFactory = 
	    jaidaFactory.createHistogramFactory(jaidaTree);
	IFunctionFactory jaidaFuncFactory = null;
	double minX = 0.0; double maxX = 1.0;
	int dimensionX = 100; 
	IHistogram1D hTest1 = 
	    jaidaHistoFactory.createHistogram1D("Test by the gaussian with same params",dimensionX,minX,maxX);
	IHistogram1D hTest2 = 
	    jaidaHistoFactory.createHistogram1D("Test by the gaussian with different params",dimensionX,minX,maxX);

	//
	// Test with the data exactly following the model function
	//
	RandomGenerator rand = new RandomGenerator();
	double[] data = new double[n];

	for(int trial = 0; trial <2000; trial++){
	    for(int i = 0; i<n; i++){
		data[i] = rand.GetGaussianDouble(mean,sigma);
	    }
	    double ksSignificance = KSTest.getKSSignificance(data, gaussFunc, functionIndex, parameters,
				       lowerBound,upperBound,false);

	    //System.out.format(" KS significance = %f\n",ksSignificance);
	    hTest1.fill(ksSignificance);
	}

	//
	// Test with the data following the model function with DIFFERENT mean and sigma
	//
	mean += 1.0*sigma;
	sigma = 1.5*sigma;
	System.out.format("Test with the Gaussian (mean=%f sigma=%f) and %d data points\n",mean,sigma,n);
	for(int trial = 0; trial <2000; trial++){
	    for(int i = 0; i<n; i++){
		data[i] = rand.GetGaussianDouble(mean,sigma);
	    }
	    double ksSignificance = KSTest.getKSSignificance(data, gaussFunc, functionIndex, parameters,
				       lowerBound,upperBound,false);
	    //System.out.format(" KS significance = %f\n",ksSignificance);

	    //System.out.format(" KS significance = %f\n",ksSignificance);
	    hTest2.fill(ksSignificance);
	}


	//
	// draw
	//
	IPlotterFactory plotterFactory = jaidaFactory.createPlotterFactory();
	IPlotter plotter = plotterFactory.create("KS test significance");
	plotter.destroyRegions();
	plotter.createRegion(0,0,0.66,1);
	plotter.createRegions(2,1);

	IPlotterStyle plotStyle = plotter.region(0).style();
	plotStyle.xAxisStyle().setLabel("significance");
	plotStyle.yAxisStyle().setLabel("Number Of Events");
	plotStyle.xAxisStyle().tickLabelStyle().setBold(true);
	plotStyle.xAxisStyle().labelStyle().setItalic(true);
	plotStyle.xAxisStyle().tickLabelStyle().setFontSize(14);
	plotStyle.yAxisStyle().tickLabelStyle().setBold(true);
	plotStyle.yAxisStyle().labelStyle().setItalic(true);
	plotStyle.yAxisStyle().tickLabelStyle().setFontSize(14);
	plotStyle.titleStyle().textStyle().setFontSize(20);
	plotter.region(0).plot(hTest1);

	IPlotterStyle plotStyle2 = plotter.region(1).style();
	plotStyle2.xAxisStyle().setLabel("significance");
	plotStyle2.yAxisStyle().setLabel("Number Of Events");
	plotStyle2.xAxisStyle().tickLabelStyle().setBold(true);
	plotStyle2.xAxisStyle().labelStyle().setItalic(true);
	plotStyle2.xAxisStyle().tickLabelStyle().setFontSize(14);
	plotStyle2.yAxisStyle().tickLabelStyle().setBold(true);
	plotStyle2.yAxisStyle().labelStyle().setItalic(true);
	plotStyle2.yAxisStyle().tickLabelStyle().setFontSize(14);
	plotStyle2.titleStyle().textStyle().setFontSize(20);
	plotter.region(1).plot(hTest2);

	plotter.show();


    } 
}
