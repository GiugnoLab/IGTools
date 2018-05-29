package igtools.common.charts.palette;

import igtools.common.nucleotide.B3Nucleotide;

import java.awt.Color;

public class BandFilterB3NucleotidePalette implements B3NucleotidePalette{
	
	public static final Color A_COLOR = Color.BLUE;
	public static final Color C_COLOR = Color.RED;
	public static final Color G_COLOR = Color.GREEN;
	public static final Color T_COLOR = Color.ORANGE;
	public static final Color BAD_COLOR = Color.GRAY;
	public static final Color NULL_COLOR = Color.PINK;
	
	private  Color[] colors;
//	static{
//		colors = new Color[5];
//		colors[B3Nucleotide.A.code()] = A_COLOR;
//		colors[B3Nucleotide.C.code()] = C_COLOR;
//		colors[B3Nucleotide.G.code()] = G_COLOR;
//		colors[B3Nucleotide.T.code()] = T_COLOR;
//		colors[B3Nucleotide.BAD.code()] = BAD_COLOR;
//	}
	
	private int[][] colorsComponents;
//	static{
//		colorsComponents = new int[5][4];
//		for(int i=0; i<5; i++){
//			colorsComponents[i][0] = colors[i].getRed();
//			colorsComponents[i][1] = colors[i].getGreen();
//			colorsComponents[i][2] = colors[i].getBlue();
//			colorsComponents[i][3] = 255;
//		}
//	}
	
	
	
	public BandFilterB3NucleotidePalette(B3Nucleotide band){
		
		colors = new Color[5];
		colors[B3Nucleotide.A.code()] = Color.BLACK;
		colors[B3Nucleotide.C.code()] = Color.BLACK;
		colors[B3Nucleotide.G.code()] = Color.BLACK;
		colors[B3Nucleotide.T.code()] = Color.BLACK;
		colors[B3Nucleotide.N.code()] = BAD_COLOR;
		
		
		switch(band){
			case A:
				colors[B3Nucleotide.A.code()] = A_COLOR;
				break;
			case C:
				colors[B3Nucleotide.C.code()] = C_COLOR;
				break;
			case G:
				colors[B3Nucleotide.G.code()] = G_COLOR;
				break;
			case T:
				colors[B3Nucleotide.T.code()] = T_COLOR;
				break;
			case N:
				colors[B3Nucleotide.A.code()] = A_COLOR;
				colors[B3Nucleotide.C.code()] = C_COLOR;
				colors[B3Nucleotide.G.code()] = G_COLOR;
				colors[B3Nucleotide.T.code()] = T_COLOR;
				colors[B3Nucleotide.N.code()] = Color.WHITE;
				break;
		}
		
		
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
