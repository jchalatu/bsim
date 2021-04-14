package BSimAsymGrowth;

import bsim.BSim;
import bsim.capsule.BSimCapsuleBacterium;
import bsim.draw.BSimDrawer;
import bsim.draw.BSimP3DDrawer;
import bsim.winter2021.Bacterium;
import processing.core.PConstants;
import processing.core.PGraphics3D;

import java.awt.*;
import java.util.ArrayList;

public class BasicDrawer extends BSimP3DDrawer{

    // Bacterium should be a sub-class of BSimCapsuleBacterium
    final ArrayList<Bacterium> bac;
    final double simX;
    final double simY;

    public BasicDrawer(BSim sim, int width_pixels, int height_pixels, double pixel_to_um_ratio, ArrayList<Bacterium> bac_to_draw) {
        super(sim, width_pixels, height_pixels);
        simX = width_pixels / pixel_to_um_ratio;
        simY = height_pixels / pixel_to_um_ratio;
        bac = bac_to_draw;
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

        scene(p3d);
        boundaries();
        time();

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

    @Override
    public void scene(PGraphics3D p3d) {
        p3d.ambientLight(128, 128, 128);
        p3d.directionalLight(128, 128, 128, 1, 1, -1);
		
		// Draw the infected bacteria in red and the non-infected bacteria in green
		for(Bacterium b : bac) {
			draw(b, Color.GREEN);
		}		

    }
}
