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
        String InfectionFields = "ImageNumber,ObjectNumber,SimulationTime,Population,InitialLength,"
        		+ "Length,TotalGrowth,InfStatus\n";
        write(InfectionFields);
    }

	@Override
	public void during() {
        String buffer = "";
        for(PhageFieldBacterium b : bac) {
            buffer += ((int) (sim.getTimestep() * sim.getDt() / dt) + 1) + "," +
                    (b.id + 1) + "," + 
            		(sim.getTimestep() * sim.getDt()) + "," +
            		bac.size() + "," +
            		(b.L_initial * BSimPhageField.pixel_to_um_ratio) + "," +
            		(b.L * BSimPhageField.pixel_to_um_ratio) + "," +
            		((b.L - b.L_initial) * BSimPhageField.pixel_to_um_ratio) + "," +
            		b.isInfected() + 
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