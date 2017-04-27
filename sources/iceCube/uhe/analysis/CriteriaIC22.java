package iceCube.uhe.analysis;

import iceCube.uhe.particles.*;
import iceCube.uhe.analysis.*;
import geometry.*;

public class CriteriaIC22 extends Criteria {

    //protected double minCobZ = -2.5e4; // -2.5e4 cm = -250 m
    //protected double maxCobZ =  2.0e5; //  2.0e5 cm =  2000 m

    protected double minCobZ = -2.5e4; // -2.5e4 cm = -250 m
    protected double maxCobZ = -5.0e3; //  -5.0e3 cm =  -50 m
    protected double min2CobZ = 5.0e3; //  5.0e3 cm = 50 m



    //protected int maxNumberOfVertex = 12;
    protected int maxNumberOfVertex = 7;

    /** The default vertex X location on the logNpe-CosZenith plane
	in the EHE super cut. Used in the setEHESuperCut().
	This is for option of the cut without using the COB-z.
    */
    //protected static double[][] vertexIC22DefaultLocationXwithNoCOBcut = 
    //{
    //	{5.33, 5.38, 5.58, 5.5801, 5.78, 5.7801, 5.93, 5.98, 6.02, 6.08, 6.0801, Double.POSITIVE_INFINITY},
    //	{5.33, 5.38, 5.58, 5.5801, 5.78, 5.7801, 5.93, 5.98, 6.02, 6.08, 6.0801, Double.POSITIVE_INFINITY},
    //	{5.33, 5.38, 5.58, 5.5801, 5.78, 5.7801, 5.93, 5.98, 6.02, 6.08, 6.0801, Double.POSITIVE_INFINITY},
    //	{5.33, 5.38, 5.58, 5.5801, 5.78, 5.7801, 5.93, 5.98, 6.02, 6.08, 6.0801, Double.POSITIVE_INFINITY}

    //    };
    protected static double[][] vertexIC22DefaultLocationXwithNoCOBcut = 
    {
	{5.21, 5.21001, 5.21002, 5.21003, 6.07, 6.0701, Double.POSITIVE_INFINITY},
	{5.21, 5.21001, 5.21002, 5.21003, 6.07, 6.0701, Double.POSITIVE_INFINITY},
	{5.21, 5.21001, 5.21002, 5.21003, 6.07, 6.0701, Double.POSITIVE_INFINITY},
	{5.21, 5.21001, 5.21002, 5.21003, 6.07, 6.0701, Double.POSITIVE_INFINITY}
    };
    /** The default vertex y location on the logNpe-CosZenith plane
	in the EHE super cut. Used in the setEHESuperCut()
	This is for option of the cut without using the COB-z.
    */
    //protected static double[][] vertexIC22DefaultLocationYwithNoCOBcut = 
    //{
    //  {-1.0, -0.3, -0.1, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 1.0},
    //  {-1.0, -0.3, -0.1, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 1.0},
    //  {-1.0, -0.3, -0.1, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 1.0},
    //  {-1.0, -0.3, -0.1, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 1.0}
    //};
    protected static double[][] vertexIC22DefaultLocationYwithNoCOBcut = 
    {
	{-1.0, -0.3468, -0.34681, -0.34682, 0.65, 0.8, 1.0},
	{-1.0, -0.3468, -0.34681, -0.34682, 0.65, 0.8, 1.0},
	{-1.0, -0.3468, -0.34681, -0.34682, 0.65, 0.8, 1.0},
	{-1.0, -0.3468, -0.34681, -0.34682, 0.65, 0.8, 1.0}
    };

    /** The default vertex X location on the logNpe-CosZenith plane
	in the EHE super cut. Used in the setEHESuperCut()
    */
    //protected static double[][] vertexIC22DefaultLocationX = 
    //{
    //	{4.78, 5.23, 5.2301, 5.58, 5.63, 5.68, 5.88, 5.98, 6.03, 6.08, 6.081, Double.POSITIVE_INFINITY},
    //	{5.33, 5.38, 5.58, 5.5801, 5.78, 5.7801, 5.88, 5.8801, 5.93, 5.9301, 5.9302, Double.POSITIVE_INFINITY},
    //	{4.78, 5.23, 5.2301, 5.58, 5.63, 5.68, 5.88, 5.98, 6.03, 6.08, 6.081, Double.POSITIVE_INFINITY},
    //	{5.33, 5.38, 5.58, 5.5801, 5.78, 5.7801, 5.88, 5.8801, 5.93, 5.9301, 5.9302, Double.POSITIVE_INFINITY}
    //};
    protected static double[][] vertexIC22DefaultLocationX = 
    {
	{4.705, 4.705001, 5.195, 5.195001, 6.035, 6.035001, Double.POSITIVE_INFINITY},
	{5.385, 5.38501, 5.465, 5.859, 5.859001, 6.083, Double.POSITIVE_INFINITY},
	{4.705, 4.705001, 5.195, 5.195001, 6.035, 6.035001, Double.POSITIVE_INFINITY},
	{5.385, 5.38501, 5.465, 5.859, 5.859001, 6.083, Double.POSITIVE_INFINITY}
    };

    /** The default vertex y location on the logNpe-CosZenith plane
	in the EHE super cut. Used in the setEHESuperCut()
    */
    //protected static double[][] vertexIC22DefaultLocationY = 
    //{
    //{-1.0, -0.3, -0.1, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 1.0},
    //	{-1.0, -0.3, -0.1, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 1.0},
    //	{-1.0, -0.3, -0.1, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 1.0},
    //	{-1.0, -0.3, -0.1, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 1.0}
    //};
    protected static double[][] vertexIC22DefaultLocationY = 
    {
	{-1.0, -0.1, -0.10001, 0.1222, 0.5889, 0.8, 1.0},
	{-1.0, 0.0, 0.0001, 0.1791, 0.179101, 0.8, 1.0},
	{-1.0, -0.1, -0.10001, 0.1222, 0.5889, 0.8, 1.0},
	{-1.0, 0.0, 0.0001, 0.1791, 0.179101, 0.8, 1.0}
    };

    /** Constructor */
    public CriteriaIC22(){
	super();
	super.minCobZ = minCobZ;
	super.maxCobZ = maxCobZ;
	super.min2CobZ = min2CobZ;
	super.maxNumberOfVertex = maxNumberOfVertex;
    }

    /** set all the GZK boundary values of Npe and cosZenith
        in the "EHE Super Cut" for a given event categories
        with the default settings.
    */
    protected void setEHESuperCut(int category){
	for(int number = 0; number<maxNumberOfVertex;number++){
	    J3Vector vertex = null;
	    if(isCOBZCut){ 
		vertex = new J3Vector(vertexIC22DefaultLocationX[category][number],
				      vertexIC22DefaultLocationY[category][number],
				      0.0);
	    }else{
		vertex = new J3Vector(vertexIC22DefaultLocationXwithNoCOBcut[category][number],
				      vertexIC22DefaultLocationYwithNoCOBcut[category][number],
				      0.0);
	    }
	    setEHESuperCut(category,vertex);
	}

    }


}
