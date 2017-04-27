package iceCube.uhe.analysis;

import iceCube.uhe.particles.*;
import iceCube.uhe.points.*;
import iceCube.uhe.geometry.*;
import geometry.*;

import java.io.*;
import java.util.*;

/** 
    Build IceCube Millipede results reading from the f2k data via standard input.
    By this way, You can leave the IceTray/i3 file framwork into
   the pure java world.

   The built I3Particle objects are outputed via outputStream.

   Orignally written by S. Yoshida for the IceCube EHE analysis
*/

public class I3ParticleMillipedeBuilder {

    /** flag to read out MC events with MCTruth or real events */
    private boolean isMCTruth;
    private final double ln10 = Math.log(10.0);

    private double trackLength = 0.0;      // neutrino track length for weighting
    ParticlePoint s = null;


    /** Constructor. You are forced to choose MC or real data mode. */
    public I3ParticleMillipedeBuilder(boolean isMCTruth){
	this.isMCTruth = isMCTruth;
    }

    /** the method to build I3Particle objects with data from the DataInputStream. 
	Generated I3Particle objects are serialized and written out to OutputStream.
     */
    public void process(DataInputStream in, OutputStream out) 
	throws IOException {

	// Reading data
	BufferedReader  d     = new BufferedReader(new InputStreamReader(in));
	String buffer; int sep = 0; int sepstart = 0;
        char separator = ' ';
	while((buffer = d.readLine())!=null){
	    try{

		// 1st line -- eventNumber, flavor, doublet
		sepstart = 0;
		sep = buffer.indexOf(separator,sepstart+1);
		sepstart = sep;

		sep = buffer.indexOf(separator,sepstart+1);
		int eventNumber =
		    Integer.valueOf(buffer.substring(sepstart+1,sep)).intValue( );
		sepstart = sep;

		sep = buffer.indexOf(separator,sepstart+1);
		int flavor =
		    Integer.valueOf(buffer.substring(sepstart+1,sep)).intValue( );
		sepstart = sep;

		sep = buffer.indexOf(separator,sepstart+1);
		int doublet =
		    Integer.valueOf(buffer.substring(sepstart+1,sep)).intValue( );
		sepstart = sep;

		// 2nd line energy distance (energy reco)
		buffer = d.readLine();
		sepstart = 0;
		sep = buffer.indexOf(separator,sepstart+1);
		sepstart = sep;

		sep = buffer.indexOf(separator,sepstart+1);
		double energy =
		    Double.valueOf(buffer.substring(sepstart+1,sep)).doubleValue( );
		sepstart = sep;

		sep = buffer.indexOf(separator,sepstart+1);
		double distance =
		    Double.valueOf(buffer.substring(sepstart+1,sep)).doubleValue( );
		sepstart = sep;

		sep = buffer.indexOf(separator,sepstart+1);
		double energyReco = -1.0;
		if(sep!=-1){
		    energyReco =
			Double.valueOf(buffer.substring(sepstart+1,sep)).doubleValue( );
		}
		sepstart = sep;

		// 3rd line -- direction (nx,ny,nz) (fitQuality)
		buffer = d.readLine();
		sepstart = 0;
		sep = buffer.indexOf(separator,sepstart+1);
		sepstart = sep;

		double[] n = new double[3];
		for(int jj=0;jj<3;jj++){
		    sep = buffer.indexOf(separator,sepstart+1);
		    n[jj] =
			Double.valueOf(buffer.substring(sepstart+1,sep)).doubleValue( );
		    sepstart = sep;
		}
		J3UnitVector direction = new J3UnitVector(n[0],n[1],n[2]);

		sep = buffer.indexOf(separator,sepstart+1);
		double fitQuality = -1.0;
		if(sep!=-1){
		    fitQuality =
			Double.valueOf(buffer.substring(sepstart+1,sep)).doubleValue( );
		}
		sepstart = sep;

		// 4th line -- vertex (rx,ry,rz)
		buffer = d.readLine();
		sepstart = 0;
		sep = buffer.indexOf(separator,sepstart+1);
		sepstart = sep;

		double[] r = new double[3];
		for(int jj=0;jj<3;jj++){
		    sep = buffer.indexOf(separator,sepstart+1);
		    r[jj] =
			Double.valueOf(buffer.substring(sepstart+1,sep)).doubleValue( );
		    sepstart = sep;
		}
		J3Vector vertex = new J3Vector(r[0],r[1],r[2]);
		J3Line axis = new J3Line(vertex,direction);

		// Additional line (if exists)
		J3Line axisReco = null;
		if(isMCTruth){
		    // 5th line -- direction (nx,ny,nz) (fitQuality)
		    buffer = d.readLine();
		    sepstart = 0;
		    sep = buffer.indexOf(separator,sepstart+1);
		    sepstart = sep;

		    for(int jj=0;jj<3;jj++){
			sep = buffer.indexOf(separator,sepstart+1);
			n[jj] =
			    Double.valueOf(buffer.substring(sepstart+1,sep)).doubleValue( );
			sepstart = sep;
		    }
		    direction = new J3UnitVector(n[0],n[1],n[2]);

		    sep = buffer.indexOf(separator,sepstart+1);
		    if(sep!=-1){
			fitQuality =
			Double.valueOf(buffer.substring(sepstart+1,sep)).doubleValue( );
		    }
		    sepstart = sep;


		    // 6th line -- vertex (rx,ry,rz)
		    buffer = d.readLine();
		    sepstart = 0;
		    sep = buffer.indexOf(separator,sepstart+1);
		    sepstart = sep;

		    for(int jj=0;jj<3;jj++){
			sep = buffer.indexOf(separator,sepstart+1);
			r[jj] =
			    Double.valueOf(buffer.substring(sepstart+1,sep)).doubleValue( );
			sepstart = sep;
		    }
		    vertex = new J3Vector(r[0],r[1],r[2]);
		    axisReco = new J3Line(vertex,direction);
		}

		// 5th line -- energy-deposit
		buffer = d.readLine();
		sepstart = 0;
		sep = buffer.indexOf(separator,sepstart+1);
		sepstart = sep;

		sep = buffer.indexOf(separator,sepstart+1);
		double energyDeposit =
		    Double.valueOf(buffer.substring(sepstart+1,sep)).doubleValue( );
		sepstart = sep;

		// 6th line -- track start position (rx,ry,rz)
		buffer = d.readLine();
		sepstart = 0;
		sep = buffer.indexOf(separator,sepstart+1);
		sepstart = sep;

		double[] rstart = new double[3];
		for(int jj=0;jj<3;jj++){
		    sep = buffer.indexOf(separator,sepstart+1);
		    rstart[jj] =
			Double.valueOf(buffer.substring(sepstart+1,sep)).doubleValue( );
		    sepstart = sep;
		}
		J3Vector vertex_start = new J3Vector(rstart[0],rstart[1],rstart[2]);

		// 6th line -- track end position (rx,ry,rz)
		buffer = d.readLine();
		sepstart = 0;
		sep = buffer.indexOf(separator,sepstart+1);
		sepstart = sep;

		double[] rend = new double[3];
		for(int jj=0;jj<3;jj++){
		    sep = buffer.indexOf(separator,sepstart+1);
		    rend[jj] =
			Double.valueOf(buffer.substring(sepstart+1,sep)).doubleValue( );
		    sepstart = sep;
		}
		J3Vector vertex_end = new J3Vector(rend[0],rend[1],rend[2]);

		// 7th line IceCube data - ATWDNpe ATWDNDOMs FADCNPE FADCNDOMs BestNPE BestNDOMs
		buffer = d.readLine();
		sepstart = 0;
		sep = buffer.indexOf(separator,sepstart+1);
		sepstart = sep;

		sep = buffer.indexOf(separator,sepstart+1);
		double npeATWD =
		    Double.valueOf(buffer.substring(sepstart+1,sep)).doubleValue( );
		sepstart = sep;

		sep = buffer.indexOf(separator,sepstart+1);
		int nDOMsATWD =
		    Integer.valueOf(buffer.substring(sepstart+1,sep)).intValue( );
		sepstart = sep;

		sep = buffer.indexOf(separator,sepstart+1);
		double npeFADC =
		    Double.valueOf(buffer.substring(sepstart+1,sep)).doubleValue( );
		sepstart = sep;

		sep = buffer.indexOf(separator,sepstart+1);
		int nDOMsFADC =
		    Integer.valueOf(buffer.substring(sepstart+1,sep)).intValue( );
		sepstart = sep;

		sep = buffer.indexOf(separator,sepstart+1);
		double npeBest =
		    Double.valueOf(buffer.substring(sepstart+1,sep)).doubleValue( );
		sepstart = sep;

		sep = buffer.indexOf(separator,sepstart+1);
		int nDOMsLaunch =
		    Integer.valueOf(buffer.substring(sepstart+1,sep)).intValue( );
		sepstart = sep;

		//
		// Now generate an I3Particle object
		//
		I3Particle iceParticle = null;
		if(!isMCTruth){
		    iceParticle = generateI3Particle(flavor,doublet,energy,distance,axis,
						     eventNumber,energyDeposit,
						     vertex_start,vertex_end,
						     npeFADC,npeATWD,npeBest,
						     nDOMsFADC,nDOMsATWD,
						     nDOMsLaunch);
		}else{
		    iceParticle = generateI3Particle(flavor,doublet,energy,energyReco,
						     distance,
						     axis,axisReco,
						     eventNumber,energyDeposit,
						     vertex_start,vertex_end,
						     npeFADC,npeATWD,npeBest,
						     nDOMsFADC,nDOMsATWD,
						     nDOMsLaunch);
		}
								       
		// First guess results exists
		if(fitQuality!= -1.0) iceParticle.setFirstGuessQuality(fitQuality);
		// Write out I3Particle object
		I3ParticleOutputStream.outputI3Particle(iceParticle, out);

	    }catch (EOFException e){
		buffer = null;
		break;
	    }
	}

    }

    /** Generate I3Particle from a set of the given valuables. */
    protected I3Particle generateI3Particle(int flavor, int doublet,
					    double energy, double distance,
					    J3Line axisInIce3,
					    int eventNumber,double energyDeposit,
					    J3Vector vertex_start, J3Vector vertex_end,
					    double npeFADC,double npeATWD, double npeBest,
					    int nDOMsFADC, int nDOMsATWD, 
					    int nDOMsLaunch){
	// Generate an I3Particle object
	I3Particle iceParticle = new I3Particle(flavor,doublet);

	if(isMCTruth) iceParticle.switchToMCTruth();
	else iceParticle.switchToReco();

	// Sets the Geometry
	if(!isMCTruth){
	    axisInIce3.setR0(vertex_start); // Sets R0 as the start vertex
	    double length = (vertex_end.getZ()-vertex_start.getZ())/axisInIce3.getDirection().getZ();
	    axisInIce3.setAxisLength(length); // Sets the axis length so that it points to the end vertex posiion
	}
	iceParticle.setParticleAxisInIceCubeCoordinate(axisInIce3);
	iceParticle.transformParticleAxisToEarthCenterCoordinate();

	// Sets the energy
	if(isMCTruth){
	    iceParticle.putEnergy(energy);
	    if(energy>0.0){
		double logEnergy = Math.log(energy)/ln10;
		iceParticle.putLogEnergy(logEnergy);
	    }else{
		iceParticle.putLogEnergy(Double.NEGATIVE_INFINITY);
	    }
	}else{
	    iceParticle.putRecoEnergy(energyDeposit);
	}

	// Sets the distance from the Earth surface
	iceParticle.setDistanceFromEarthSurfaceToIceCube(distance);

	// Fills the IceCube data
	iceParticle.getIceCubeData().setEventNumber(eventNumber);
	iceParticle.getIceCubeData().setNpeATWD(npeATWD);
	iceParticle.getIceCubeData().setNpeFADC(npeFADC);
	iceParticle.getIceCubeData().setBestNpe(npeBest);
	iceParticle.getIceCubeData().setNDOMsATWD(nDOMsATWD);
	iceParticle.getIceCubeData().setNDOMsFADC(nDOMsFADC);
	iceParticle.getIceCubeData().setNDOMsLaunch(nDOMsLaunch);

	return iceParticle;

    }

    /** Generate I3Particle from a set of the given valuables.
	This method will be used when you build a MC event.
     */
    protected I3Particle generateI3Particle(int flavor, int doublet,
					    double energyMCTruth, double energyReco,
					    double distance,
					    J3Line axisInIce3MCTruth,
					    J3Line axisInIce3Reco,
					    int eventNumber,double energyDeposit,
					    J3Vector vertex_start, J3Vector vertex_end,
					    double npeFADC,double npeATWD, double npeBest,
					    int nDOMsFADC, int nDOMsATWD, 
					    int nDOMsLaunch){

	// Generate an I3Particle object
	isMCTruth = true;
	I3Particle iceParticle = 
	    generateI3Particle(flavor, doublet, energyMCTruth, distance,
			       axisInIce3MCTruth,
			       eventNumber, energyDeposit,
			       vertex_start,vertex_end,
			       npeFADC, npeATWD, npeBest,
			       nDOMsFADC, nDOMsATWD, nDOMsLaunch);
	// Set the reconsturcted geometry
	iceParticle.switchToReco();
	axisInIce3Reco.setR0(vertex_start); // Sets R0 as the start vertex
	double length = (vertex_end.getZ()-vertex_start.getZ())/axisInIce3Reco.getDirection().getZ();
	axisInIce3Reco.setAxisLength(length); // Sets the axis length so that it points to the end vertex posiion
	iceParticle.setParticleAxisInIceCubeCoordinate(axisInIce3Reco);
	iceParticle.transformParticleAxisToEarthCenterCoordinate();

	// set the (reco) energy
	iceParticle.putRecoEnergy(energyDeposit);

	// set the MC primary spectrum weight
	iceParticle.switchToMCTruth();

	return iceParticle;

    }	

  public static void main(String[] args) throws IOException{

	boolean isMCTruth = false;

	String fileName = null;
	I3Particle iceTrack;

        if(args.length!=2){
            System.out.println("Usage: I3ParticleStreamDemo file-name MCtruth(yes 1 no 0)");
	    System.exit(0);
        }else{
             fileName = args[0];
             if(Integer.valueOf(args[1]).intValue() == 1) isMCTruth = true;
        }

	/** I3ParticleBuilder */
	I3ParticleMillipedeBuilder builder =  new I3ParticleMillipedeBuilder(isMCTruth);

	/** Open Output data stream */
	FileOutputStream out = new FileOutputStream(fileName);
	DataInputStream input = new DataInputStream(System.in);

	/** Builds the I3Particle from the data via input.*/
	builder.process(input, out);
	out.close( );
	System.err.println("Building The I3Particle  done.");

	/** Check the serialized I3Particle objects stored in disk */

	FileInputStream in = new FileInputStream(fileName);
	I3Particle iceParticle = null;
	while((iceParticle = I3ParticleInputStream.inputI3Particle(in))
	      !=null){

	    double npeATWD = iceParticle.getIceCubeData().getNpeATWD();
	    double npeFADC = iceParticle.getIceCubeData().getNpeFADC();
	    double npeBest = iceParticle.getIceCubeData().getBestNpe();
	    if(npeBest>1.0e4){
		System.out.print("- EventNumber is " + 
				 iceParticle.getIceCubeData().getEventNumber());
		System.out.println(" The I3Particle Name is " + 
				   iceParticle.particleName(iceParticle.getFlavor(), 
							    iceParticle.getDoublet()));

		System.out.print("- The I3Particle Energy (" + 
				 iceParticle.getEnergy( ) + ") [GeV] ");
		System.out.println("- The I3Particle Reco Energy (" + 
				   iceParticle.getRecoEnergy( ) + ") [GeV] ");
		J3UnitVector n = iceParticle.getDirectionInIceCubeCoordinate();
		System.out.println("- The I3Particle direction in ice3(" + 
				   n.getX() + ", " + n.getY() + ", " + n.getZ() +
				   ")");
		System.out.print("- Npe by ATWD = " + 
				   iceParticle.getIceCubeData().getNpeATWD());
		System.out.print(" Npe by FADC = " + 
				   iceParticle.getIceCubeData().getNpeFADC());
		System.out.print(" Best Npe = " + 
				   iceParticle.getIceCubeData().getBestNpe());
		double fitQuality = iceParticle.getFirstGuessQuality();
		System.out.println("  FG velocity(" + fitQuality + ")");
	    }
	}

	in.close();

	System.err.println(" everything done --");
    }

}

