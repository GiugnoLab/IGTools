package igtools.common.charts.palette;

import igtools.common.nucleotide.B3Nucleotide;

import java.awt.Color;

public class RedB3NucleotidePalette implements B3NucleotidePalette{
	//255, 128, 64, 32

	public static final Color A_COLOR = new Color(255, 0, 0);
	public static final Color C_COLOR = new Color(128, 0, 0);
	public static final Color G_COLOR = new Color(64, 0, 0);
	public static final Color T_COLOR = new Color(32, 0, 0);
	public static final Color BAD_COLOR = Color.GRAY;
	public static final Color NULL_COLOR = Color.PINK;
	
	private static Color[] colors;
	static{
		colors = new Color[5];
		colors[B3Nucleotide.A.code()] = A_COLOR;
		colors[B3Nucleotide.C.code()] = C_COLOR;
		colors[B3Nucleotide.G.code()] = G_COLOR;
		colors[B3Nucleotide.T.code()] = T_COLOR;
		colors[B3Nucleotide.N.code()] = BAD_COLOR;
	}
	
	private static int[][] colorsComponents;
	static{
		colorsComponents = new int[5][4];
		for(int i=0; i<5; i++){
			colorsComponents[i][0] = colors[i].getRed();
			colorsComponents[i][1] = colors[i].getGreen();
			colorsComponents[i][2] = colors[i].getBlue();
			colorsComponents[i][3] = 255;
		}
	}

	@Override
	public Color getColor(B3Nucleotide n) {
		return colors[n.code()];
	}

	@Override
	public Color getColor(int code) {
		return colors[code];
	}
	
	@Override
	public int[] getColorComponents(int code){
		return colorsComponents[code];
	}
}
