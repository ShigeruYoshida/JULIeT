package iceCube.uhe.event;

import iceCube.uhe.particles.*;
import iceCube.uhe.points.*;
import geometry.*;

import java.awt.*;
import java.awt.geom.*;
import java.awt.event.*;
import javax.swing.*;
import java.util.*;

public class ParticleTravelGraphics extends JComponent {

    /** Common logarithm */
    static final double ln10 = Math.log(10.0);

    private JulietTask task;
    private int offset = 20;
    private double dMax = 1.0e5;
    private double eMax = 100.0;
    private double gridWidth;
    private double gridHeight;
    private double propagationDistance = 0.0;
    private double logEminForGraphics = 1.0; // Log(10 [GeV])
    final static Color fg = Color.black;

    /** Constructor */
    ParticleTravelGraphics(JulietTask task){
	this.task = task;
	setBackground(Color.WHITE);
	setForeground(fg);
	setOpaque(true);
    }


    protected void paintComponent(Graphics g) {
        Graphics2D g2 = (Graphics2D) g;
	Dimension d = getSize();

	//Paint background if we're opaque.
	if (isOpaque()) {
	    g2.setColor(getBackground());
	    g2.fillRect(0, 0, getWidth(), getHeight());
	}

	//Paint the Coordinate Axis.
	GradientPaint blueToWhite = new GradientPaint(0,0,Color.WHITE,
						       0, d.height,Color.BLUE);
        g2.setPaint(blueToWhite);
        g2.fill(new Rectangle2D.Double(0, 0, d.width, d.height));
	if(task.generator!= null){
	    if(task.generator.particleList!= null){
		g2.setColor(Color.BLACK);
		drawAxis(g2);
		g2.setColor(Color.CYAN);
		drawTrajectory(g2);
	    }
	}
    }


    private void drawAxis(Graphics2D g2){
	Dimension d = getSize();
        g2.drawString("Trajectory", 0, offset-10);
	g2.drawString("0.0 km", 0, offset);

	J3Vector rStart = task.generator.wherePrimaryParticleStartsInIceCubeCoordinate();
	J3Vector rEnd = task.generator.wherePrimaryParticleEndsInIceCubeCoordinate();
	J3Vector rStartToEnd = J3Vector.subtract(rEnd,rStart);
	propagationDistance = rStartToEnd.getLength();//[cm]
	dMax = propagationDistance*1.02;
        gridWidth = (double )(d.width - offset)/ dMax;
        g2.setPaint(fg);
        g2.draw(new Line2D.Double(offset, offset, dMax*gridWidth+offset, offset));
	String propagationDistanceString = Double.toString(propagationDistance/1.0e5);//[km]
	g2.drawString(propagationDistanceString, (int )(dMax*gridWidth)-30, offset);
	rStart = null; rEnd = null; rStartToEnd = null;

	double primaryEnergy = task.generator.primaryEnergy;
	double logPrimaryEnergy = Math.log(primaryEnergy)/ln10;
	eMax = (double )((int )(logPrimaryEnergy-logEminForGraphics));
        gridHeight = (double )(d.height - offset)/ eMax;
        g2.draw(new Line2D.Double(offset, offset, offset, eMax*gridHeight+offset));
	String s = Integer.toString((int )logPrimaryEnergy);
	String logEString = s.concat(" Log[GeV]");
	g2.drawString(logEString, 0, (int )(eMax*gridHeight)+offset);

	double logE = 1.0 + logEminForGraphics;
	while(logE < logPrimaryEnergy){ // draw scale axis
	    g2.draw(new Line2D.Double(offset, 
				      (logE-logEminForGraphics)*gridHeight+offset, 
				      dMax*gridWidth+offset, 
				      (logE-logEminForGraphics)*gridHeight+offset));
	    logEString = Double.toString(logE);
	    g2.drawString(logEString, 0, 
			  (int)((logE-logEminForGraphics)*gridHeight)+offset);
	    logE += 1.0;
	}

    }

    private void drawTrajectory(Graphics2D g2){

	ListIterator locationIterator = task.generator.getLocationIce3Iterator();
	ListIterator particleIterator = task.generator.getParticleIterator();

	Dimension d = getSize();

	if(locationIterator.hasNext()){
	    double lStart = 0.0;
	    double deltaL = propagationDistance*1.0e-2;
	    J3Vector rStart = task.generator.wherePrimaryParticleStartsInIceCubeCoordinate();
	    while(locationIterator.hasNext()){
		Particle cascadeParticle = (Particle )(particleIterator.next());
		double logEnergy = cascadeParticle.getLogEnergy();
		J3Vector cascadeLocation = (J3Vector )(locationIterator.next());
		double pathLength = J3Vector.subtract(cascadeLocation,rStart).getLength();
		g2.fill(new Rectangle2D.Double(offset+(int)(pathLength*gridWidth), offset, 
				       (int )(gridWidth*deltaL), 
					       (int)(gridHeight*(logEnergy-logEminForGraphics))));
	    }

	}

    }

 }
