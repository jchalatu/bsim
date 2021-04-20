package bsim.winter2021;

import bsim.BSim;
import bsim.capsule.BSimCapsuleBacterium;
import bsim.export.BSimLogger;

import java.text.DecimalFormat;
import java.text.DecimalFormatSymbols;
import java.util.ArrayList;
import java.util.Locale;

public class PositionLogger extends BSimLogger {
    ArrayList<? extends BSimCapsuleBacterium> bacteriaAll;
    DecimalFormat formatter = new DecimalFormat("###.##", DecimalFormatSymbols.getInstance( Locale.ENGLISH ));

    public PositionLogger(BSim sim, String filePath, ArrayList<? extends BSimCapsuleBacterium> bacteriaAll) {
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
