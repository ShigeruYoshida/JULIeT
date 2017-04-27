package iceCube.uhe.analysis;

import iceCube.uhe.analysis.*;
import geometry.*;

import java.io.*;
import java.util.*;

public class RunAnalysis {

    public static void main(String[] args) throws IOException{

	String realDataFileName = "iceCube/uhe/analysis/EHERealI3Particles";
	//String mcDataFileName = "iceCube/uhe/analysis/EHEMCI3Particles";
	String mcDataFileName = "iceCube/uhe/analysis/EHEMCReducedIC9I3Particles";
	//InputStream in = ClassLoader.getSystemResourceAsStream(realDataFileName);
	InputStream in = ClassLoader.getSystemResourceAsStream(mcDataFileName);
	I3ParticleAnalysisFactory analizer = new I3ParticleAnalysisFactory(in);
	in.close();

	J3Vector ic9Center = new J3Vector(analizer.xCenterOfIC9,
					  analizer.yCenterOfIC9, 
					  0.0);
	// IceCube 9 string array center

	// Sets the Criteria
	Criteria criteria = new Criteria();
	criteria.setThresholdOfLogNpe(3.0);
	criteria.setThresholdOfNDOMs(80);
	analizer.setCriteria(criteria);

	// Drawing
	analizer.setMaxRangeOfNDOMsInDrawing(300);

	System.out.println("zone 2 2"); // 2 x 2 panels

	//analizer.switchToMCTruth();
	analizer.makeHistogram();
	analizer.drawNpeDistributionOnXfig();
	System.out.println("endg");


	//criteria.setMaxDistance(1.0e7,ic9Center);
	analizer.switchToReco();
	//analizer.setBinSize(0.2,0.2,0.1,1.0); 
	analizer.setBinSize(0.5,0.5,0.1,1.0); 
	                // Delta (LogE) Delta(LogNpe), Delta cos(Zenith), DeltaFGquality
	analizer.makeHistogram();
	//analizer.plotByJointLine();
	//analizer.drawNpeDistributionWithSliceOfEnergyOnXfig();
	analizer.drawZenithAngleDistributionOnXfig();
	System.out.println("endg");


	//criteria.setMaxDistance(1.0e7,ic9Center);
	criteria.setThresholdOfLogNpe(4.0);	
	analizer.setBinSize(1.0,1.0,0.1,0.5); 
	                // Delta (LogE) Delta(LogNpe), Delta cos(Zenith), DeltaFGquality
	analizer.switchToReco();
	analizer.makeHistogram();
	//criteria.setThresholdOfFirstGuessQuality(1.5e10);
	analizer.plotByPointsWithErrorBars();
	analizer.drawZenithAngleDistributionWithFGsliceOnXfig();
	//analizer.drawZenithAngleDistributionOnXfig();
	System.out.println("endg");


	analizer.setBinSize(1.0,1.0,0.5,0.1); 
	                // Delta (LogE) Delta(LogNpe), Delta cos(Zenith), DeltaFGquality
	analizer.makeHistogram();
	analizer.plotByJointLine();
	analizer.drawFirstGuessQualityDistributionOnXfig();
	//analizer.drawEnergyDistributionOnXfig();
	System.out.println("endg");
    }

}
