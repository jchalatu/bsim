package bsim.winter2021;

import java.awt.Color; 

import bsim.BSim;
import bsim.BSimChemicalField;
import bsim.draw.BSimP3DDrawer;
import processing.core.PGraphics3D;

public class P3DDrawer extends BSimP3DDrawer {

	public P3DDrawer(BSim sim, int width, int height) {
		super(sim, width, height);
		// TODO Auto-generated constructor stub
	}

	@Override
	public void scene(PGraphics3D p3d) {
		// TODO Auto-generated method stub
		
	}
	
	/**
	 * Draws a 2D chemical field structure based on its defined parameters, with custom transparency (alpha) parameters.
	 * @param field		The chemical field structure to be rendered.
	 * @param c			Desired colour of the chemical field.
	 * @param alphaGrad	Alpha per unit concentration of the field.
	 * @param alphaMax	Maximum alpha value (enables better viewing).
	 */
	public void draw2D(BSimChemicalField field, Color c, double alphaGrad, double alphaMax) {
		int[] boxes = field.getBoxes();
		double[] boxSize = field.getBox();
		double alpha = 0.0f;
		
		for( int i = 0; i < boxes[0]; i++ ) {
			for( int j = 0; j < boxes[1]; j++ ) {						
				p3d.pushMatrix();					
				p3d.translate((float)(boxSize[0]*i+boxSize[0]/2), (float)(boxSize[1]*j+boxSize[1]/2));
				
				alpha = alphaGrad*field.getConc(i, j, 0);
				if (alpha > alphaMax) alpha = alphaMax;
				
				p3d.fill(c.getRed(), c.getGreen(), c.getBlue(),(float)alpha);
				
				p3d.box((float)boxSize[0], (float)boxSize[1], 0); 
				p3d.popMatrix();
			}
		}
	}
	
	/**
	 * Draw a chemical field structure based on its defined parameters (default alpha).
	 * @param field The chemical field to be drawn.
	 * @param c The desired colour.
	 * @param alphaGrad The alpha-per-unit-concentration.
	 */
	public void draw2D(BSimChemicalField field, Color c, float alphaGrad) {
		draw2D(field, c, alphaGrad, 255);
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
	

	
	
	
	
}




