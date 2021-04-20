package bsim.winter2021;

import BasicSimulation2D.Parameters;
import bsim.BSim;
import bsim.export.BSimLogger;

import javax.vecmath.Vector3d;

public class MetaLogger extends BSimLogger {
    Parameters parameters;
    int initialPopulation;

    public MetaLogger(BSim sim, String filePath, Parameters parameters, int initialPopulation) {
        super(sim, filePath);
        this.parameters = parameters;
        this.initialPopulation = initialPopulation;
    }

    @Override
    public void before() {
        super.before();
        write("Simulation Metadata.");
        Vector3d sim_dimensions = sim.getBound();
        write("Simulation dimensions (in microns): (" + sim_dimensions.x + ", " + sim_dimensions.y + ", " + sim_dimensions.z + ")");
        write("Simulation Timestep: "+ sim.getDt());
        write("Initial Population Size: "+ initialPopulation);
        // can add: how the simulation was initialized (raw positions, CellProfiler data, or from a pop size),
        // simulation parameters (from the parameters field),
        // the start time and run time of the simulation
    }

    @Override
    public void during() {

    }
}
