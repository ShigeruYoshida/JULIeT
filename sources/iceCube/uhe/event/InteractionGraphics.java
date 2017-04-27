package iceCube.uhe.event;

import iceCube.uhe.particles.*;
import iceCube.uhe.points.*;

import java.awt.*;
import java.awt.geom.*;
import java.awt.event.*;
import javax.swing.*;
import java.util.*;

public class InteractionGraphics extends JComponent {

    private JulietTask task;
    private int offset = 20;
    private double interactionMax = 20.0;
    private double gridWidth;
    private double gridHeight;
    final static Color fg = Color.black;

    /** Constructor */
    InteractionGraphics(JulietTask task){
	this.task = task;
	setBackground(Color.WHITE);
	setForeground(fg);
	setOpaque(true);
    }


    protected void paintComponent(Graphics g) {
        Graphics2D g2 = (Graphics2D) g;
	Dimension d = getSize();
        gridWidth = (double )(d.width - offset)/7.0;
        gridHeight = (double )((d.height - offset)/interactionMax);

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
	    if(task.generator.particleList!= null) drawHist(g2);
	}

    }


    private void drawAxis(Graphics2D g2){
	Dimension d = getSize();
	GradientPaint GreenToWhite = new GradientPaint(0,0,Color.WHITE,
						       0, d.height,Color.GREEN);
        g2.setPaint(GreenToWhite);
        g2.fill(new Rectangle2D.Double(0, 0, d.width, d.height));

        g2.setPaint(fg);
        g2.draw(new Line2D.Double(offset, offset, d.width+offset, offset));
        g2.drawString("Interactions", 0, offset-10);
        g2.drawString("Frequency", 0, d.height-5);
        g2.draw(new Line2D.Double(offset, offset, offset, d.height+offset));
        gridWidth = (double )(d.width - offset)/7.0;
        gridHeight = (double )((d.height - offset)/interactionMax);
        g2.drawString("Emg Cas", 1*(int )gridWidth, offset);
        g2.drawString("Hadron Cas", 4*(int )gridWidth, offset);
	//        g2.drawString("Photo-Nucl", 5*(int )gridWidth, offset);

    }

    private void drawHist(Graphics2D g2){
	Dimension d = getSize();
        gridWidth = (double )(d.width - offset)/7.0;
        gridHeight = (double )((d.height - offset)/interactionMax);
        g2.setPaint(Color.YELLOW);
	ListIterator particleIterator = task.generator.getParticleIterator();
	int hadCascade = 0;
	int eCascade = 0; 
	while(particleIterator.hasNext()){
	    Particle cascadeParticle = (Particle )(particleIterator.next());
	    if(cascadeParticle.getFlavor()==0) eCascade++;
	    else if(cascadeParticle.getFlavor()==3) hadCascade++;
	}
	int eHist= eCascade;
	int eHist2 = eCascade-(int)interactionMax;
	int eHist3 = eCascade-2*(int)interactionMax;
	int eHist4 = eCascade-3*(int)interactionMax;

        g2.fill(new Rectangle2D.Double(offset, offset, 
				       (int )(gridWidth), 
				       (int)(gridHeight*eHist)));
	if(eHist2>0)
        g2.fill(new Rectangle2D.Double(offset+(int )gridWidth, offset, 
				       (int )(gridWidth), 
				       (int)(gridHeight*eHist2)));
	if(eHist3>0)
        g2.fill(new Rectangle2D.Double(offset+2*(int )gridWidth, offset, 
				       (int )(gridWidth), 
				       (int)(gridHeight*eHist3)));
	if(eHist4>0)
        g2.fill(new Rectangle2D.Double(offset+3*(int )gridWidth, offset, 
				       (int )(gridWidth), 
				       (int)(gridHeight*eHist3)));
        //g2.setPaint(Color.ORANGE);
        g2.setPaint(Color.RED);
        g2.fill(new Rectangle2D.Double(offset+4*(int )gridWidth, offset, 
				       (int )(gridWidth), 
				       (int)(gridHeight*hadCascade)));
        //g2.setPaint(Color.RED);
        //g2.fill(new Rectangle2D.Double(offset+5*(int )gridWidth, offset, 
	//				       (int )(gridWidth), 
	//				       (int)(gridHeight*task.photoN)));

    }

 }
