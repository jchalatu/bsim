package BSimCrossProtection;

import bsim.BSim;
import bsim.BSimChemicalField;
import bsim.capsule.BSimCapsuleBacterium;
import bsim.draw.BSimDrawer;
import bsim.draw.BSimP3DDrawer;
import processing.core.PConstants;
import processing.core.PGraphics3D;

import java.awt.*;
import java.util.ArrayList;

public class BasicDrawer extends BSimP3DDrawer{

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
    
	/**
	 * Draws the concentration difference relative to field A based on its defined parameters, with custom transparency (alpha) parameters.
	 * @param field_A	The chemical field structure to be rendered.
	 * @param field_B	The chemical field structure to be rendered.
	 * @param c			Desired colour of the chemical field.
	 * @param alphaGrad	Alpha per unit concentration of the field.
	 * @param alphaMax	Maximum alpha value (enables better viewing).
	 */
	public void drawConcDifference2D(BSimChemicalField field_A, BSimChemicalField field_B, double alphaGrad, double alphaMax, int max_conc_difference) {
		int[] boxes;						// Number of boxes
		double[] boxSize;					// Size of each box
		double alpha = 0.0f;
		
		// Choose the minimum number of boxes
		if ( field_A.getBoxes().length < field_B.getBoxes().length ) {
			boxes = field_A.getBoxes();
			boxSize = field_A.getBox();
		}
		else {
			boxes = field_B.getBoxes();
			boxSize = field_B.getBox();
		}
		
		for(int i = 0; i < boxes[0]; i++) {
			for(int j = 0; j < boxes[1]; j++) {					
				p3d.pushMatrix();					
				p3d.translate((float)(boxSize[0]*i+boxSize[0]/2), (float)(boxSize[1]*j+boxSize[1]/2));
				
				alpha = alphaGrad*field_A.getConc(i, j, 0) + alphaGrad*field_B.getConc(i, j, 0);
				if (alpha > alphaMax) 
					alpha = alphaMax;
				
				//System.out.println("conc: " + (field_A.getConc(i, j, 0) - field_B.getConc(i, j, 0)));
				Color c = getColor( field_A.getConc(i, j, 0), field_B.getConc(i, j, 0), max_conc_difference);
				
				p3d.fill(c.getRed(), c.getGreen(), c.getBlue(), (float)alpha);
				p3d.box((float)boxSize[0], (float)boxSize[1], 0); 	// Only draw the x and y dimensions
				p3d.popMatrix();
			}
		}
	}
	
	/** Returns the associated color from the difference of two concentrations relative to concA. */
	public static Color getColor( double conc_A, double conc_B, double max_conc_difference ) {
		double conc = conc_A - conc_B;			// Concentration difference is relative to field A
		float value = 0.0f;						// Between 0 and 1
		
		// The floor of the hue is multiplied by 360 to get hue angle from HSB color model
		float max_hue = 0.1f;
		float min_hue = 0.9f;
		float mid_hue = (max_hue + min_hue) / 2;
		
		if ( conc >= max_conc_difference ) {
			value = 1/2f;// Blue 		
		}
		else if ( conc <= -max_conc_difference ) {
			value = 2/3f;// Green		
		}
		else {
			value = (float) ( (mid_hue)*(conc/max_conc_difference) );
		}
		
		float hue = value * max_hue + ( 1 - value ) * min_hue;
		return new Color( Color.HSBtoRGB(hue, 1f, 1f) );
	}
	
	/**
	 * Draw a chemical field structure based on its defined parameters (default alpha).
	 * @param field_A The chemical field to be drawn.
	 * @param field_B The chemical field to be drawn.
	 * @param c The desired colour.
	 * @param alphaGrad The alpha-per-unit-concentration.
	 */
	public void drawConcDifference2D(BSimChemicalField field_A, BSimChemicalField field_B, float alphaGrad, int max_conc_difference) {
		drawConcDifference2D(field_A, field_B, alphaGrad, 255, max_conc_difference);
	}
	
	/** Draws the reference for the mixed concentration color map. */
	public void concDifferenceRef(float x, float y, float w, float h, double max_conc_difference, int scale, String fieldA, String fieldB) {
		p3d.fill(50);
		int boxes = 21;
		int labelNum = 5;
		double step = (max_conc_difference / boxes);
		
		// Field labels
		p3d.text(fieldB, 35, 110);
		p3d.text(fieldA, 35, 520);
		
		// Draw text
		for ( int i = 0; i < labelNum; i ++ ) {
			p3d.text(Math.abs((int)max_conc_difference - scale * i), 70, 140 + (435/labelNum)*i );
		}
		
		// Draw colors
    	for ( int i = 0; i < h; i++ ) {
    		float value = 0.75f + i/(float)h;	
    		float max_hue = 0.1f;							
    		float min_hue = 0.9f;	
    		
    		float hue = value * max_hue + ( 1 - value ) * min_hue;
    		
        	Color c = new Color( Color.HSBtoRGB(hue, 1f, 1f) );
        	
			p3d.stroke(c.hashCode());
			p3d.line(x, y+i*0.1f, x+w, y+i*0.1f);
        	//p3d.fill(c.hashCode());
    		//p3d.rect(x, y + i * h, w, h);
    	}
    	p3d.noStroke();
	}
	
	/** Draws a reference for a single chemical field to screen. */
	public void drawSingleFieldRef( float x, float y, float w, float h, Color c_start, Color c_end ) {
		p3d.pushStyle();
		for ( int i = 0; i < w; i ++ ) {
			int c = p3d.lerpColor(c_start.hashCode(), c_end.hashCode(), i/(float)w);
			p3d.stroke(c);
			p3d.line(x+i*0.1f, y, x+i*0.1f, y+h);
		}
		p3d.popStyle();
	}
	/** Draws the label for the reference for a single chemical field to screen. */
	public void drawLabel( int x, int y, double conc_max, double conc_min ) {
		p3d.fill(50);
		p3d.text((int)conc_max, x, y);
		p3d.text((int)(conc_max + conc_min)*2/3, x + 171, y);
		p3d.text((int)(conc_max + conc_min)*1/3, x + 342, y);
		p3d.text((int)conc_min, x + 515, y);
	}
	
	/**
	 * Draws two chemical field structures in a checkerboard pattern based on its defined parameters, with custom transparency (alpha) parameters.
	 * @param field_A	The chemical field structure to be rendered.
	 * @param field_B	The chemical field structure to be rendered.
	 * @param c_A		Desired colour of chemical field A.
	 * @param c_A		Desired colour of chemical field B.
	 * @param alphaGrad	Alpha per unit concentration of the field.
	 * @param alphaMax	Maximum alpha value (enables better viewing).
	 */
	public void draw2DGrid(BSim sim, BSimChemicalField field_A, BSimChemicalField field_B, Color c_A, Color c_B, double alphaGrad, double alphaMax) {
		int[] boxes;
		double[] boxSize;
		double alpha = 0.0f;
		
		Color c = c_A;
		
		// Choose the minimum number of boxes
		if ( field_A.getBoxes().length < field_B.getBoxes().length ) {
			boxes = field_A.getBoxes();
			boxSize = field_A.getBox();
		}
		else {
			boxes = field_B.getBoxes();
			boxSize = field_B.getBox();
		}
		
		// Draw the two fields in a checkerboard pattern
		for ( int i = 0; i < boxes[0]; i ++ ) {
			for( int j = 0; j < boxes[1]; j++ ) {
				p3d.pushMatrix();					
				p3d.translate((float)(boxSize[0]*i+boxSize[0]/2), (float)(boxSize[1]*j+boxSize[1]/2));
				
				if ( (i + j) % 2 == 0 ) {
					alpha = alphaGrad*field_A.getConc(i, j, 0);
					c = c_A;
				}
				else {
					alpha = alphaGrad*field_B.getConc(i, j, 0);
					c = c_B;
				}
				
				if (alpha > alphaMax) 
					alpha = alphaMax;
				
				p3d.fill(c.getRed(), c.getGreen(), c.getBlue(), (float)alpha);
				p3d.box((float)boxSize[0], (float)boxSize[1], 0); 	// Only draw the x and y dimensions
				p3d.popMatrix();
			}
		}
	}
	
	/**
	 * Draws two chemical field structures in a checkerboard pattern based on its defined parameters (default alpha).
	 * @param field_A The chemical field to be drawn.
	 * @param field_B The chemical field to be drawn.
	 * @param c_A The desired colour for field A.
	 * @param c_B The desired colour for field B.
	 * @param alphaGrad The alpha-per-unit-concentration.
	 */
	public void draw2DGrid(BSim sim, BSimChemicalField field_A, BSimChemicalField field_B, Color c_A, Color c_B, float alphaGrad) {
		draw2DGrid(sim, field_A, field_B, c_A, c_B, alphaGrad, 255);
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
