package iceCube.uhe.event;

import java.io.*;

public class RunJuliet {

    /** Define the type of run **/
    static int     numberOfEvent;
    static int     typeOfEvent;

    /** Main method. In order to run JULIET, it is need that RunManager object 
        and DataOutPutStream are generated. */
    public static void main(String[] args) throws IOException {


	// generate RunManager object
	JulietEventGenerator generator = new JulietEventGenerator(); 

	DataInputStream input = new DataInputStream(System.in);
	BufferedReader  d     = new BufferedReader(new InputStreamReader(input));
	String buffer;

	System.out.print("x [cm] in the ice3 coordinate->");
	buffer   = d.readLine();
	double x = Double.valueOf(buffer).doubleValue();
	System.out.print("y [cm] in the ice3 coordinate->");
	buffer   = d.readLine();
	double y = Double.valueOf(buffer).doubleValue();
	System.out.print("z [cm] in the ice3 coordinate->");
	buffer   = d.readLine();
	double z = Double.valueOf(buffer).doubleValue();
	System.out.print("nadir angle [deg] in the ice3 coordinate->");
	buffer   = d.readLine();
	double nadirAngleInDeg = Double.valueOf(buffer).doubleValue();
	System.out.print("Azimuth angle [deg] in the ice3 coordinate->");
	buffer   = d.readLine();
	double azimuthAngleInDeg = Double.valueOf(buffer).doubleValue();

	generator.definePropagationGeometry(x,y,z,nadirAngleInDeg,azimuthAngleInDeg);
	generator.configurePropagationGeometry();

	System.out.print("Set the neutrino interaction weight? yes(1)->");
	buffer   = d.readLine();
	if(Integer.valueOf(buffer).intValue()==1){
	    System.out.print("Enter the neutrino interaction weight->");
	    buffer   = d.readLine();
	    int weight = Integer.valueOf(buffer).intValue();
	    generator.setNeutrinoInteractionWeight(weight);
	    System.out.println("Has set the weight of " +
			       generator.getNeutrinoInteractionWeight());
	}


	String  fileName; 

	System.out.print("Name of File ->");
	buffer   = d.readLine();
	fileName = String.valueOf(buffer).trim();
	System.out.print("Number of Event ->");	
	buffer        = d.readLine();
	numberOfEvent = Integer.valueOf(buffer).intValue();

	System.out.print(
	"Type of Event (0:single energy (full data) 1:multi energy (Matrix ) ) ->");
	buffer      = d.readLine();
	typeOfEvent = Integer.valueOf(buffer).intValue();



	// make output stream
	DataOutputStream out = new DataOutputStream(new FileOutputStream(fileName));

	if(typeOfEvent==0){
	    for(int n=0; n<numberOfEvent; n++){
		generator.runSingleEvent();
		int numberOfCollisions = generator.getListedEvents(out);
		System.out.println("Number of collisions (" + 
				   numberOfCollisions + ")");
	    }
	}else if(typeOfEvent==1) {
	    generator.runEventOnMatrix(numberOfEvent, out);
	}
    	out.close();
    }

}
