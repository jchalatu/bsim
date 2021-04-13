package BSimCrossProtection;

import bsim.BSim;
import bsim.capsule.BSimCapsuleBacterium;
import bsim.draw.BSimDrawer;
import bsim.winter2021.P3DDrawer;
import processing.core.PConstants;
import processing.core.PGraphics3D;

import java.awt.*;
import java.util.ArrayList;

public class BasicDrawer extends P3DDrawer{

    // Bacterium should be a sub-class of BSimCapsuleBacterium
    final ArrayList<CrossProtectionBacterium> bacA;
    final ArrayList<CrossProtectionBacterium> bacB;
    final double simX;
    final double simY;
    final double window_height;
    final double window_width;
    
    public static ChemicalField antibioticA;
    public static ChemicalField antibioticB;
    public static double c; 
    
    /** Flag to draw fields on two separate simulation screens. */
    final boolean TWO_SCREENS;
	/** Steady state concentration level. */
	private double initial_conc = 310;
    
    // Two ways to show the antibiotic fields on a single screen
    int SINGLE_SCREEN;
    final int CHECKER_BOARD = 1;
    final int MIXED_CONC = 2;
    
    public BasicDrawer(BSim sim, double simX, double simY, int window_width, int window_height,
    		ArrayList<CrossProtectionBacterium> bac_to_drawA, ArrayList<CrossProtectionBacterium> bac_to_drawB,
    		ChemicalField antibioticA, ChemicalField antibioticB, double c, boolean TWO_SCREENS, int SINGLE_SCREEN) {
        super(sim, window_width, window_height);
        this.simX = simX;
        this.simY = simY;
        this.window_height = window_height;
        this.window_width = window_width;
        bacA = bac_to_drawA;
        bacB = bac_to_drawB;
        
        BasicDrawer.antibioticA = antibioticA;
        BasicDrawer.antibioticB = antibioticB;
        BasicDrawer.c = c;
        this.TWO_SCREENS = TWO_SCREENS;
        this.SINGLE_SCREEN = SINGLE_SCREEN;
    }

    /**
     * Draw the default cuboid boundary of the simulation as a partially transparent box
     * with a wireframe outline surrounding it.
     */
    @Override
    public void boundaries() {
        p3d.noFill();
        p3d.stroke(128, 128, 255);
        p3d.pushMatrix();
        p3d.translate((float)boundCentre.x,(float)boundCentre.y,(float)boundCentre.z);
        p3d.box((float)bound.x, (float)bound.y, (float)bound.z);
        p3d.popMatrix();
        p3d.noStroke();
    }
    
    /** Draws a colored box to screen. */
    public void legend( Color color, String title, int [] boxParams, int textX, int textY ) {
    	// Box
    	p3d.fill( color.getRed(), color.getGreen(), color.getBlue() );
    	p3d.stroke(128, 128, 255);
    	p3d.rect(boxParams[0], boxParams[1], boxParams[2], boxParams[3]);
    	
    	// Text
    	p3d.fill(50);
    	p3d.text(title, textX, textY);
    }

    @Override
    public void draw(Graphics2D g) {
        p3d.beginDraw();

        if(!cameraIsInitialised){
            // camera(eyeX, eyeY, eyeZ, centerX, centerY, centerZ, upX, upY, upZ)
            p3d.camera((float)bound.x*0.5f, (float)bound.y*0.5f,
                    // Set the Z offset to the largest of X/Y dimensions for a reasonable zoom-out distance:
                    simX > simY ? (float)simX : (float)simY,
                    (float)bound.x*0.5f, (float)bound.y*0.5f, 0,
                    0,1,0);
            cameraIsInitialised = true;
        }

        p3d.textFont(font);
        p3d.textMode(PConstants.SCREEN);

        p3d.sphereDetail(10);
        p3d.noStroke();
        p3d.background(255, 255,255);
        
        // Draw two simulation boxes if two screens are selected
        if ( TWO_SCREENS ) {
        	
        	// Apply light settings
            lights(p3d);
        	
            // Draw simulation with first field
            p3d.pushMatrix();
            p3d.translate((float) -45, 0, 0);
            scene(p3d, antibioticA, Color.MAGENTA);
            boundaries();
            time();
            p3d.popMatrix();
            
            // Draw simulation with second field
            p3d.pushMatrix();
            p3d.translate((float) 45, 0, 0);
            scene(p3d, antibioticB, Color.BLUE);
            boundaries();
            time();
            p3d.popMatrix();

            // Draw reference
            drawSingleFieldRef(-45, 60, 720, 3, Color.MAGENTA, Color.WHITE);
            drawLabel(65, (int) (window_height - 55), initial_conc, 0);
            drawSingleFieldRef(45, 60, 720, 3, Color.BLUE.brighter(), Color.WHITE);
            drawLabel(614, (int) (window_height - 55), initial_conc, 0);
            
        }
        // Draw only one simulation box
        else {
            scene(p3d);
            boundaries();
            if ( SINGLE_SCREEN == CHECKER_BOARD ) {
            	legend( Color.MAGENTA, "Toxin A", new int [] {-17, 0, 3, 3}, 57, 113 );
            	legend( Color.BLUE, "Toxin B", new int [] {-17, 6, 3, 3}, 57, 157 );
            }
            else if ( SINGLE_SCREEN == MIXED_CONC ) {
            	concDifferenceRef(-15, 4, 3, 525/*2.5f*/, 100, 50, "Field A", "Field B");
            }
            time();
        }

        p3d.endDraw();
        g.drawImage(p3d.image, 0,0, null);
    }

    /**
     * Draw the formatted simulation time to screen.
     */
    @Override
    public void time() {
        p3d.fill(0);
        //p3d.text(sim.getFormattedTimeHours(), 50, 50);
        p3d.text(sim.getFormattedTime(), 50, 50);
    }
    
    /** Applies light settings to scene. */
    public void lights( PGraphics3D p3d ) {
        p3d.ambientLight(128, 128, 128);
        p3d.directionalLight(128, 128, 128, 1, 1, -1);
        //p3d.lightFalloff(1, 1, 0);
    }
    
    /** Draws chemical fields in separate boxes. */
    public void scene(PGraphics3D p3d, ChemicalField field, Color field_color) {

        // Draw the antibiotic field
        draw2D(field, field_color, (float)(255/c));
        
        // Draw sub-population A
        for ( int i = 0; i < bacA.size(); i ++ ) {
        	draw(bacA.get(i), bacA.get(i).isAboveThreshold() ? Color.RED : Color.GREEN);
        }
        
        // Draw sub-population B
        for ( int i = 0; i < bacB.size(); i ++ ) {
        	draw(bacB.get(i), bacB.get(i).isAboveThreshold() ? Color.RED : Color.ORANGE);
        }

    }

    @Override
    public void scene(PGraphics3D p3d) {
        p3d.ambientLight(128, 128, 128);
        p3d.directionalLight(128, 128, 128, 1, 1, -1);
        
        // Draw the antibiotic field
        if ( SINGLE_SCREEN == CHECKER_BOARD ) {
        	draw2DGrid(sim, antibioticA, antibioticB, Color.MAGENTA, Color.BLUE, (float)(255/c));
        }
        else if ( SINGLE_SCREEN == MIXED_CONC ) {
        	int max_conc_difference = 100;
        	drawConcDifference2D(antibioticA, antibioticB, (float)(255/c), max_conc_difference);	
        }
		
        // Draw sub-population A
        for ( int i = 0; i < bacA.size(); i ++ ) {
        	draw(bacA.get(i), bacA.get(i).isAboveThreshold() ? Color.RED : Color.GREEN);
        }
        
        // Draw sub-population B 
        for ( int i = 0; i < bacB.size(); i ++ ) {
        	draw(bacB.get(i), bacB.get(i).isAboveThreshold() ? Color.RED : Color.ORANGE);
        }
    }
}
