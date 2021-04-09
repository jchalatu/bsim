package BSimPhageLogger;

import java.io.IOException; 
import java.util.ArrayList;

import bsim.BSim;
import bsim.export.BSimLogger;

public class PhageFieldLogger extends BSimLogger {
	
	ArrayList<PhageFieldBacterium> bac;

	public PhageFieldLogger(BSim sim, String filename, ArrayList<PhageFieldBacterium> bac) {
		super(sim, filename);
		this.bac = bac;
	}
    
    @Override
    public void before() {
        super.before();
        String CellProfilerFields = "ImageNumber,ObjectNumber,"
        		+ "AreaShape_MajorAxisLength,AreaShape_MinorAxisLength,AreaShape_Center_X,AreaShape_Center_Y,"
        		+ "AreaShape_Orientation,AreaShape_Area\n";
        write(CellProfilerFields);
    }

	@Override
	public void during() {
        String buffer = "";
        for(PhageFieldBacterium b : bac) {
        	double area = (Math.PI * b.radius * b.radius + 2 * b.radius * b.L) * 13.89 * 13.89;
            buffer += ((int) (sim.getTimestep() * sim.getDt() / dt) + 1) + "," +
                    (b.id + 1) + "," + 
            		((b.L + 2 * b.radius) * BSimPhageField.pixel_to_um_ratio) + "," +
            		(2 * b.radius) * BSimPhageField.pixel_to_um_ratio + "," +
            		(b.position.x * BSimPhageField.pixel_to_um_ratio) + "," +
            		(b.position.y * BSimPhageField.pixel_to_um_ratio) + "," +
            		b.cell_profiler_angle() + "," +
            		area +
                    "\n";
        }
        write(buffer);
	}
	
    @Override
    public void write(String text) {
        try {
            bufferedWriter.write(text);
        } catch (IOException e) {
            e.printStackTrace();
        }
    }
	
}

