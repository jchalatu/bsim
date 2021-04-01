package BSimCrossProtection;

import bsim.BSim;
import bsim.BSimChemicalField;

public class ChemicalField extends BSimChemicalField{

	public ChemicalField(BSim sim, int[] boxes, double diffusivity, double decayRate) {
		super(sim, boxes, diffusivity, decayRate);
	}
	
	/** Adds concentration to areas of low concentration. */
	public void fillArea( double conc, double initial_state ) {
		int[] index = {0,0,0};
		
        for(index[0] = 0; index[0]<boxes[0]; index[0]++)
			for(index[1] = 0; index[1]<boxes[1]; index[1]++)
				for(index[2] = 0; index[2]<boxes[2]; index[2]++) {
					double box_conc = getConc(index[0], index[1], index[2]);
					if ( box_conc < initial_state) {
						if ( box_conc + conc * sim.getDt() > initial_state ) {
							double concB = initial_state - box_conc;
							addQuantity(index[0], index[1], index[2], concB * sim.getDt());
						}
						else {
							addQuantity(index[0], index[1], index[2], conc * sim.getDt());
						}
					}
				}
	}
	
	/** Adds concentration to every point in the field in the direction specified by axis. */
	public void fill( int axis, double startConc, double endConc ) {
		assert ((axis >= 0) && (axis <= 2)) :
			"Chemical field gradient - axis selection out of range [0,2]\n" +
			"Check axis is one of: x-axis = '0', y = '1', or z = '2'";
		
		int[] index = {0,0,0};
		
		double grad = (endConc - startConc)/(double)boxes[axis];
        
        for(index[0] = 0; index[0]<boxes[0]; index[0]++)
			for(index[1] = 0; index[1]<boxes[1]; index[1]++)
				for(index[2] = 0; index[2]<boxes[2]; index[2]++)
					addQuantity(index[0], index[1], index[2], startConc + index[axis]*grad);
	}
	
	
	
	
	
	
	
}