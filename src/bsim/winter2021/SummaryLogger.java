package bsim.winter2021;

import bsim.BSim;
import bsim.capsule.BSimCapsuleBacterium;
import bsim.export.BSimLogger;

import java.util.ArrayList;

public class SummaryLogger extends BSimLogger {
    ArrayList<? extends BSimCapsuleBacterium> bacteriaAll;

    public SummaryLogger(BSim sim, String filePath, ArrayList<? extends BSimCapsuleBacterium> bacteriaAll) {
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
