package igtools.common.charts.palette;

import java.awt.Color;

import igtools.common.nucleotide.*;

public interface B3NucleotidePalette {
	
	public Color getColor(B3Nucleotide n);
	public Color getColor(int code);
	
	public int[] getColorComponents(int code);

}
