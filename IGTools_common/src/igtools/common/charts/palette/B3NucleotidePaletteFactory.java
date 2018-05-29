package igtools.common.charts.palette;

import igtools.common.nucleotide.B3Nucleotide;

import java.awt.Color;

public class B3NucleotidePaletteFactory {

	
	
	public static B3NucleotidePalette createPalette(
			final Color a_color,
			final Color c_color,
			final Color g_color,
			final Color t_color,
			final Color n_color
			){
		
		final Color[] colors = new Color[5];
		colors[B3Nucleotide.A.code()] = a_color;
		colors[B3Nucleotide.C.code()] = c_color;
		colors[B3Nucleotide.G.code()] = g_color;
		colors[B3Nucleotide.T.code()] = t_color;
		colors[B3Nucleotide.N.code()] = n_color;
		
		final int[][] colorsComponents = new int[5][4];
		for(int i=0; i<5; i++){
			colorsComponents[i][0] = colors[i].getRed();
			colorsComponents[i][1] = colors[i].getGreen();
			colorsComponents[i][2] = colors[i].getBlue();
			colorsComponents[i][3] = 255;
		}
		
		return new B3NucleotidePalette() {
			
			@Override
			public int[] getColorComponents(int code) {
				return colorsComponents[code];
			}
			
			@Override
			public Color getColor(int code) {
				return colors[code];
			}
			
			@Override
			public Color getColor(B3Nucleotide n) {
				return colors[n.code()];
			}
		};
	}
	
	
	
	
	
	
}
