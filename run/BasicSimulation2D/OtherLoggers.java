package BasicSimulation2D;

import bsim.BSim;
import bsim.capsule.BSimCapsuleBacterium;
import bsim.export.BSimLogger;

import java.text.DecimalFormat;
import java.text.DecimalFormatSymbols;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.Locale;

public class OtherLoggers {
    public class MetaLogger extends BSimLogger {

        Parameters parameters;

        public MetaLogger(BSim sim, String filePath, Parameters parameters) {
            super(sim, filePath);
            this.parameters = parameters;
        }

        @Override
        public void before() {
            super.before();
            write("Simulation metadata.");
            write("Catie Terrey Fall 2020."); //change name when new person :)
            List<Double> sim_dim = new ArrayList<>(Arrays.asList(new Double[]{198.0, 159.0, 1.})); // parameters.simDimensions;
            write("Simulation dimensions: (" + sim_dim.get(0) + ", " + sim_dim.get(1) + ", " + sim_dim.get(2) + ")");
            write("Initial population: "+ parameters.initialPopulation);
            write("Ratio " + parameters.populationRatio);



            if(parameters.fixedBounds){
                write("Boundaries: fixed");
            } else {
                write("Boundaries: leaky");
            }
        }

        @Override
        public void during() {

        }
    }

    public class PosLogger extends BSimLogger {

        ArrayList<BSimCapsuleBacterium> bacteriaAll;
        DecimalFormat formatter = new DecimalFormat("###.##", DecimalFormatSymbols.getInstance( Locale.ENGLISH ));

        public PosLogger(BSim sim, String filePath, ArrayList<BSimCapsuleBacterium> bacteriaAll) {
            super(sim, filePath);
            this.bacteriaAll = bacteriaAll;
        }

        @Override
        public void before() {
            super.before();
            write("per Act; per Rep; id, p1x, p1y, p1z, p2x, p2y, p2z, growth_rate,directions");
        }

        @Override
        public void during() {
            String buffer = new String();

            buffer += sim.getFormattedTime() + "\n";
            write(buffer);

            write("acts");

            buffer = "";
            for(BSimCapsuleBacterium b : bacteriaAll) {
                buffer += b.id + "," + formatter.format(b.x1.x)
                        + "," + formatter.format(b.x1.y)
                        + "," + formatter.format(b.x1.z)
                        + "," + formatter.format(b.x2.x)
                        + "," + formatter.format(b.x2.y)
                        + "," + formatter.format(b.x2.z)
                        + "," + formatter.format(b.getK_growth())
                        + "\n";
            }

            write(buffer);


        }
    }


    public class SumLogger extends BSimLogger {
        ArrayList<BSimCapsuleBacterium> bacteriaAll;

        public SumLogger(BSim sim, String filePath, ArrayList<BSimCapsuleBacterium> bacteriaAll) {
            super(sim, filePath);
            this.bacteriaAll = bacteriaAll;
        }


        @Override
        public void before() {
            super.before();
            write("time,id, p1x, p1y, p1z, p2x, p2y, p2z, px, py, pz, growth_rate, directions");
        }

        @Override
        public void during() {
            String buffer = new String();
            buffer = "";
            for(BSimCapsuleBacterium b : bacteriaAll) {
                buffer += sim.getFormattedTime()+","+b.id
                        + "," + b.x1.x
                        + "," + b.x1.y
                        + "," + b.x1.z
                        + "," + b.x2.x
                        + "," + b.x2.y
                        + "," + b.x2.z
                        + "," + b.position.x
                        + "," + b.position.y
                        + "," + b.position.z
                        + "," + b.getK_growth()
                        + "," + b.direction()
                        + "\n";
            }

            write(buffer);
        }
    }

}
