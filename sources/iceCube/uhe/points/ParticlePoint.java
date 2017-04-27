package iceCube.uhe.points;

import java.io.*;

/** 
<pre>
    This class contains the variables and methods concerining
    the particle location and the propagation medium.

    Calculations on propagation of UHE particles in the rock/ice require
    the point vector that is provided by this class.
    It also provides the parameters of the medium such as the density
    which are used for the cross section calculation made in
    the Interaction class.

    The parameters provided in this class are:

    The vector of the particle location--
        R     : Radius from the earth center.
	lAxis : Length along the propagation axis
                from the incident point to the current particle location.
        theta : Zenith angle [rad]
	phi   : Azimuth angle [rad]
	alpha : Nadir angle [rad]
	        Those angular parameters above are defined
                at the point where a particle emerges from underground.
	MaterialNumber: ice (0) rock (1)

	These are based on the Standard Earth Model
	and the related geometrical factor.

    The parameters of the material at the particle location--
        Density : Mass density [g/cm^3]
	Xslant  : Slant Depth: [g/cm^2]
            For each component..
	       Atomic Number
	       Charge
	Material name (Rock, ice, etc.)

	     Material Number  0   1
Madium                       ice rock
Number of Species             2   1
</pre>
*/

public class ParticlePoint implements Serializable {

    /** Avogadro's Number. */
    public final static double NA = 6.022e23;
    private int MaterialNumber;   // Material Number
    private final static double[] SurfaceDensity = {0.917, 2.65};//Ice, Rock
    private final static double[] IonizationPotential = {75.0, 136.4};
                                                           //Ice, Rock [eV]
    /** Number of Species of Atom. 2 for Ice (H and O) and 1 for Rock.*/
    public final static int[] NumberOfSpecies = {2,1};     //Ice, Rock
    private final static double[][] AtomicNumber = 
    {
	{15.9994, 1.00794},//Ice
	{22.0}             //Rock
    };
    private final static double[][] Charge = 
    {
	{8.0, 1.0},       //Ice
	{11.0}            //Rock
    };
    private final static double[][] RadiationConstant = 
    {
	{173.4, 202.4},       //Ice
	{189.0}               //Rock
    };
    private final static int[][] NumberOfAtoms = 
    {
	{1, 2},            //Ice
	{1}                //Rock
    };

    /** Earth radius [cm].*/
    public final static double REarth = 6.37814e8;
    private double IceRockBoundary; // Ice/Rock boundary radius [cm]
    private double R;     //Distance from the Earth Center [cm].
    private double alpha; //Nadir angle [rad]
    private double theta; //Zenith angle [rad]
    private double lAxis; //Length along the trajectory axis [cm]
    private double AxisLength;
                          //Axis length from incident point at the surface
                          // to emerging point at the earth surface.
    private double Xslant;//Slant depth along the trajectory axis [g/cm]


    /** Constructor to initialize the starting point.
	Input parameters of Length along the trajectory axis
	(0.0 for surface) and the nadir angle alpha [deg]
	is given to define the other geometrical parameter.
    */
    public ParticlePoint(double lAxis, double alpha, int MaterialNumber){
	if(isValidNadir(alpha)){
	    AxisLength = 2.0*REarth*Math.cos(alpha);
	    this.alpha = alpha;
	    theta = 2.0*Math.PI-alpha;
	    this.MaterialNumber = MaterialNumber;
	    if(lAxis <AxisLength){
		this.lAxis = lAxis;
		R = getRadiusFromEarthCenter(lAxis);
	    }else {
		System.err.println("Illegal initial parameters");
		System.exit(0);
	    }
	    IceRockBoundary = 1.01*REarth; // Default Ice/Rock boundary radius [cm]
	}else{
	    System.err.println("Illegal initial angular parameters");
	    System.exit(0);
	}
    }


    public static boolean isValidNadir(double alpha){
	if(0.0<=alpha && alpha<=2.0*Math.PI){
	    return true;
	}else {
	    return false;
	}
    }


    /** Calculate the radius from the earth center to
	the particle location.
    */
    public double getRadiusFromEarthCenter(double lAxis){
	double rSquared = REarth*REarth+(AxisLength-lAxis)*(AxisLength-lAxis)
	    -2.0*REarth*(AxisLength-lAxis)*Math.cos(alpha);
	double radius = Math.sqrt(rSquared);
	return radius;
    }

    /** Setting the particle location along the trajectory,
	i.e. lAxis.
    */
    public void setParticleLocation(double lAxis){
	if(lAxis>=0.0){
	    this.lAxis = lAxis;
	}
	if(lAxis>AxisLength){
	    //System.err.println("Warning: Your particle has now emerged from underground!");
	}

    }

    /** Obtain the particle location along the trajectory,
	i.e. lAxis.
    */
    public double getParticleLocation(){
	return lAxis;
    }
    public double getAxisLength( ){
	return  AxisLength;
    }


    /** Setting the slant depth */
    public void setSlantDepth(double Xslant){
	this.Xslant = Xslant;
    }

    /** Obtain the slant depth */
    public double getSlantDepth( ){
	return Xslant;
    }



    /** Obtain the madium density [g/cm^2] */
    public double getMediumDensity( ){
	double r = getRadiusFromEarthCenter(getParticleLocation( ));

	double density=SurfaceDensity[MaterialNumber];
	double x = r/REarth;

	if(MaterialNumber==1){ //Rock

	    if(0.0<=r && r < 1.2215e8){ //[cm]
		density = 13.0885-8.8381*x*x;
	    }else if(r < 3.48e8){
		density = 12.5815 - 1.2638*x - 3.6426*x*x -5.5281*x*x*x;
	    }else if(r < 5.701e8){
		density = 7.9565 - 6.4761*x + 5.5283*x*x - 3.0807*x*x*x;
	    }else if(r < 5.771e8){
		density = 5.3197 - 1.4836*x;
	    }else if(r < 5.971e8){
		density = 11.2494 - 8.0298*x;
	    }else if(r < 6.151e8){
		density = 7.1089 - 3.8045*x;
	    }else if(r < 6.3466e8){
		density = 2.691 + 0.6924*x;
	    }else if(r < 6.356e8){
		density = 2.9;
	    }else if(r < IceRockBoundary){
		density = SurfaceDensity[MaterialNumber];
	    }else{
		MaterialNumber = 0; //Ice
		density = SurfaceDensity[MaterialNumber];
	    }
	}

	return density;
    }


    public double getIceRockBoundaryRadius( ){
	return IceRockBoundary;
    }
    public void setIceRockBoundaryRadius(double r){
	IceRockBoundary = r;
    }


    /** Atomic-number is acquired for a given Material Number and a given species
	index of atom. */
    public double getAtomicNumber(int iThSpecies){

	double atomicNumber =0.0;

	atomicNumber = AtomicNumber[MaterialNumber][iThSpecies];
	return atomicNumber;
    }

    /** Charge-number is acquired for a given Material Number and a given species
	index of atom. */
    public double getCharge(int iThSpecies){

	double charge;

	charge = Charge[MaterialNumber][iThSpecies];
	return charge;
    }

    /** Radiation constant is acquired for a given Material Number 
        and a given species index of atom. */
    public double getRadiation(int iThSpecies){

	double r;

	r = RadiationConstant[MaterialNumber][iThSpecies];
	return r;
    }

    /** Number of Atoms is acquired for a given Material Number and a given species
	index of atom. */
    public int getNumberOfAtoms(int iThSpecies){

	int n;

	n = NumberOfAtoms[MaterialNumber][iThSpecies];
	return n;
    }

    /** Ionization Potential [eV] */
    public double getIonizationE( ){
	return IonizationPotential[MaterialNumber];
    }


    /** Sets material number to specify either Ice or standard Rock */
    public void setMaterialNumber(int MaterialNumber){
	this.MaterialNumber = MaterialNumber;
    }

    public int getMaterialNumber( ){
	return MaterialNumber;
    }

}

