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
	
	/*
	 * Updated for 2D simulation of phage field
	 */
	
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
	 * Draws two chemical field structures based on its defined parameters, with custom transparency (alpha) parameters.
	 * @param field_A	The chemical field structure to be rendered.
	 * @param field_B	The chemical field structure to be rendered.
	 * @param c			Desired colour of the chemical field.
	 * @param alphaGrad	Alpha per unit concentration of the field.
	 * @param alphaMax	Maximum alpha value (enables better viewing).
	 */
	public void drawMixedConc2D(BSimChemicalField field_A, BSimChemicalField field_B, Color c_A, Color c_B, double alphaGrad, double alphaMax) {
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
				
				Color c = getColor( field_A.getConc(i, j, 0), field_B.getConc(i, j, 0));
				
				p3d.fill(c.getRed(), c.getGreen(), c.getBlue(), (float)alpha);
				p3d.box((float)boxSize[0], (float)boxSize[1], 0); 	// Only draw the x and y dimensions
				p3d.popMatrix();
			}
		}
	}
	
	public static Color getColor( double conc_A, double conc_B ) {
		double conc = conc_A - conc_B;
		float value = 0.0f;						// Between 0 and 1
		
		double conc_threshold = 1000;
		
		if ( conc >= conc_threshold ) {
			value = 1.0f;
		}
		else if ( conc <= -conc_threshold ) {
			value = -1.0f;
		}
		else {
			value = (float) (conc/conc_threshold);
		}
		
		float max_hue = 0f;//2f;								
		//float min_hue = Color.RGBtoHSB(base_color.getRed(), base_color.getGreen(), base_color.getBlue(), null)[0];						// Red
		float min_hue = 1f;					// Red
		
		//System.out.println(value);
		
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
	public void drawMixedConc2D(BSimChemicalField field_A, BSimChemicalField field_B, Color c_A, Color c_B, float alphaGrad) {
		drawMixedConc2D(field_A, field_B, c_A, c_B, alphaGrad, 255);
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




