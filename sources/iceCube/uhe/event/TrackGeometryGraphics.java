package iceCube.uhe.event;

import iceCube.uhe.particles.*;
import iceCube.uhe.points.*;
import geometry.*;

import java.awt.*;
import java.awt.geom.*;
import java.awt.event.*;
import javax.swing.*;
import java.util.*;

public class TrackGeometryGraphics extends JComponent {

    private JulietTask task;
    private double offset = 1.0e5;
    private double ice3Dimension = 1.0e5; // 1.0e5 [cm] = 1 [km]
    private double center = offset + 0.5*ice3Dimension;
    private double xWidth = 3.0e5; // 3.0e5 [cm] = 3 [km]
    private double yWidth = 3.0e5; // 3.0e5 [cm] = 3 [km]
    private double gridWidth;
    private double gridHeight;
    final static Color fg = Color.blue;
    final static BasicStroke stroke = new BasicStroke(2.0f);
    final static BasicStroke wideStroke = new BasicStroke(4.0f);

    /** Constructor */
    TrackGeometryGraphics(JulietTask task){
	this.task = task;
	setBackground(Color.WHITE);
	setForeground(fg);
	setOpaque(true);
    }


    protected void paintComponent(Graphics g) {
        Graphics2D g2 = (Graphics2D) g;
	Dimension d = getSize();
        gridWidth = 0.5*(double )(d.width)/xWidth;
        gridHeight = (double )(d.height)/yWidth;

	//Paint background if we're opaque.
	if (isOpaque()) {
	    g2.setColor(getBackground());
	    g2.fillRect(0, 0, getWidth(), getHeight());
	}

	//Paint the Coordinate Axis.
	g2.setColor(Color.BLACK);
	drawAxis(g2);
	g2.setColor(Color.CYAN);
	if(task.generator!= null){
	    if(task.generator.particleList!= null) drawTrack(g2);
	}

    }


    private void drawAxis(Graphics2D g2){
	Dimension d = getSize();
	GradientPaint blueToRed = new GradientPaint(0,0,Color.RED,
						       0, d.height,Color.BLUE);
        g2.setPaint(blueToRed);
        g2.fill(new Rectangle2D.Double(0, 0, d.width, d.height));

        g2.setPaint(Color.WHITE);
	double ice3x = -0.5*ice3Dimension;
	double ice3y = -0.5*ice3Dimension;
	// X-Z plain
        g2.fill(new Rectangle2D.Double((int)((ice3x+center)*gridWidth), 
				       (int)((ice3y+center)*gridHeight), 
				       (int)(ice3Dimension*gridWidth), 
				       (int)(ice3Dimension*gridHeight)));
	// Y-Z plain
        g2.fill(new Rectangle2D.Double((int)((xWidth+ice3x+center)*gridWidth), 
				       (int)((ice3y+center)*gridHeight), 
				       (int)(ice3Dimension*gridWidth), 
				       (int)(ice3Dimension*gridHeight)));
    }

    private void drawTrack(Graphics2D g2){
	Dimension d = getSize();
        g2.setPaint(Color.YELLOW);
	g2.setStroke(wideStroke);
	J3Vector rStart = task.generator.wherePrimaryParticleStartsInIceCubeCoordinate();
	J3Vector rEnd = task.generator.wherePrimaryParticleEndsInIceCubeCoordinate();
	// Projection to X-Z plain
        g2.draw(new Line2D.Double((int )((rStart.getX()+center)*gridWidth),
				  (int )((rStart.getZ()+center)*gridHeight),
				  (int )((rEnd.getX()+center)*gridWidth),
				  (int )((rEnd.getZ()+center)*gridHeight)));
	// Projection to Y-Z plain
        g2.draw(new Line2D.Double((int )((xWidth+rStart.getY()+center)*gridWidth),
				  (int )((rStart.getZ()+center)*gridHeight),
				  (int )((xWidth+rEnd.getY()+center)*gridWidth),
				  (int )((rEnd.getZ()+center)*gridHeight)));

    }

 }
