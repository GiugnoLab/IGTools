package igtools.common.charts;

import java.awt.Color;
import java.awt.Font;
import java.awt.FontMetrics;
import java.awt.Graphics2D;
import java.awt.RenderingHints;
import java.awt.image.BufferedImage;
import java.lang.reflect.Field;
import java.text.DecimalFormat;
import java.util.Arrays;
import java.util.Map;

import igtools.common.charts.palette.DefaultB3NucleotidePalette;


/**
 * TO be used only for positive ([0 ... Double.MAX_VALUE]) values.
 * 
 * 
 * @author vbonnici
 *
 */
public class KHeatMap {
	
	public int margin = 20;
	public Color backgroundColor = Color.white;
	
	public int lineWidth = 4096;
	public int lineHeight = 100;
	public int lineSpace = 0;
	
	public boolean printTopRuler = true;
	public boolean printBottomRuler = true;
	public int rulerLineHeight = 20;
	public int rulerSpace = 10;
	
	public boolean printLegend  = true;
	public boolean printMultiLegend  = true;
	public int legendHeight = 50;
	public Font font = new Font(Font.MONOSPACED, Font.PLAIN, 50);
	public Color fontColor = Color.black;
	
	public int yAxisSpace = 2;
	
	public boolean normalizeByRow = true;
	public boolean fromZeroGradient = true;
	
	public double max[] = null;
	public double min[] = null;
	public double rmax[] = null;
	public double rmin[] = null;
	
	public double[][] values = null;
	public String[] labels = null;
	public int labelsSpace = 20;
	public boolean printLabels = true;
	
	public boolean printValues = true;
	public DecimalFormat iformat = new DecimalFormat("###,###,###,###");
	public DecimalFormat dformat = new DecimalFormat("0.###");
	
	public BufferedImage img = null;
	public Graphics2D g2d = null;
	
	private int colourValueDistance = 0;
	public Color highValueColor = Color.blue;
	public Color medianValueColor = Color.cyan;
	public Color lowValueColor = Color.white;
	public double hlColorLimit = 1.0;
	//public Color highValueColour = Color.cyan;
	//public Color lowValueColour = Color.black;
	//public Color medianValueColour = Color.blue;
	//public Color medianValueColour = null;
	public double colourScale = 1.0;
	/*
     * Sets the scale that is currently in use to map z-value to colour. A 
     * value of 1.0 will give a <strong>linear</strong> scale, which will 
     * spread the distribution of colours evenly amoungst the full range of 
     * represented z-values. A value of greater than 1.0 will give an 
     * <strong>exponential</strong> scale that will produce greater emphasis 
     * for the separation between higher values and a value between 0.0 and 1.0
     * will provide a <strong>logarithmic</strong> scale, with greater 
     * separation of low values. Values of 0.0 or less are illegal.
     * gh
     * <p>
     * Defaults to a linear scale value of 1.0.
	*/
	
	public KHeatMap(){
	}
	
	
	public void draw(){
		int height = 0;
		int width = 0;
		
		
		//get height
		height += margin * 2;
		if(printTopRuler){
			height += rulerLineHeight * 3;
			height += rulerSpace;
		}
		for(int i=0; i<values.length; i++){
			height += lineHeight;
			if(i != values.length)
				height += lineSpace;
		}
		if(printBottomRuler){
			height += rulerSpace;
			height += rulerLineHeight;
		}
		if(printLegend){
			height += margin;
			if(printMultiLegend){
				height += (values.length * legendHeight) + ((values.length -1 ) * lineSpace);
			}
			else{
				height += legendHeight;
			}
		}
		
		//get width
		width += margin * 2;
		width += lineWidth; 
		
		//get labels space
		img = new BufferedImage(1, 1, BufferedImage.TYPE_INT_RGB);
		g2d = (Graphics2D)img.getGraphics();
		g2d.setRenderingHint(RenderingHints.KEY_ANTIALIASING, RenderingHints.VALUE_ANTIALIAS_ON);
		
		
		int lstart = margin;
		
		
		g2d.setFont(font);
		FontMetrics metrics = g2d.getFontMetrics(font);
		int fontHeight = 0;
		int fw = 0;
		int maxLabel = 0;
		if(printLabels && labels != null){
			fontHeight = metrics.getHeight();
			if(fontHeight > lineHeight){
				fontHeight = lineHeight;
				font = new Font(font.getFamily(), font.getStyle(), lineHeight);
				g2d.setFont(font);
				metrics = g2d.getFontMetrics(font);
			}
					
			for(String s : labels){
				fw = metrics.stringWidth(s);
				if( fw > maxLabel)
					maxLabel = fw;
			}
			
			width += labelsSpace;
			width += maxLabel;
			lstart += labelsSpace;
			lstart += maxLabel;
		}
		
		
		//build image
		img = new BufferedImage(width, height, BufferedImage.TYPE_INT_RGB);
		g2d = (Graphics2D)img.getGraphics();
		g2d.setRenderingHint(RenderingHints.KEY_ANTIALIASING, RenderingHints.VALUE_ANTIALIAS_ON);
		g2d.setFont(font);
		metrics = g2d.getFontMetrics(font);
		g2d.setColor(backgroundColor);
		g2d.fillRect(0, 0, width, height);
		
		int y = 0;
		y += margin;
		int x = lstart;
		
		//print top ruler
		Color[] ncolors = new Color[4];
		ncolors[0] = DefaultB3NucleotidePalette.A_COLOR;
		ncolors[1] = DefaultB3NucleotidePalette.C_COLOR;
		ncolors[2] = DefaultB3NucleotidePalette.G_COLOR;
		ncolors[3] = DefaultB3NucleotidePalette.T_COLOR;
		
		if(printTopRuler){
			for(int k=1; k<2; k++){
				x = lstart;
				int nof = (int)Math.pow(4, k);
				int w = lineWidth / nof;
				
				int c = 0;
				for(int i=0; i<nof; i++){
					g2d.setColor(ncolors[c]);
					c++; if(c>3) c=0;
					
					g2d.fillRect(x, y, w, rulerLineHeight);
					x += w;
				}
				 y += rulerLineHeight;
			}
			y += rulerSpace;
		}
		
		
		//print k-lines
		updateColourDistance();
		
		int ly = y;
		
		boolean printValue = true;
		boolean printValueRow = true;
		String valueString = "";
		Color vcolor;
		Color compvcolor;
		
		max = new double[values.length];
		min = new double[values.length];
		rmax = new double[values.length];
		rmin = new double[values.length];
		
		if(normalizeByRow){
			for(int i=0; i<values.length; i++){
				max[i] = rmax[i] = values[i][0];
				min[i] = rmin[i] = values[i][0];
				for(int j=1; j<values[i].length; j++){
					if(max[i] < values[i][j]) max[i] = rmax[i] = values[i][j];
					if(min[i] > values[i][j]) min[i] = rmin[i] = values[i][j];
				}
			}
		}
		else{
			max[0] = values[0][0];
			min[0] = values[0][0];
			for(int i=0; i<values.length; i++){
				rmin[i] = values[i][0];
				rmax[i] = values[i][0];
				for(int j=0; j<values[i].length; j++){
					if(max[0] < values[i][j]) max[0] = values[i][j];
					if(min[0] > values[i][j]) min[0] = values[i][j];
					if(rmin[i] > values[i][j]) rmin[i] = values[i][j];
					if(rmax[i] < values[i][j]) rmax[i] = values[i][j];
				}
			}
			for(int i=1; i<values.length; i++){
				min[i] = min[0];
				max[i] = max[0];
			}
		}
		if(fromZeroGradient){
			for(int i=0; i<values.length; i++){
				min[i] = 0;
			}
		}
			
		
//		if(normalizeByRow){
//			max = new double[values.length];
//			min = new double[values.length];
//			for(int i=0; i<values.length; i++){
//				max[i] = values[i][0];
//				min[i] = values[i][0];
//				for(int j=1; j<values[i].length; j++){
//					if(max[i] < values[i][j]) max[i] = values[i][j];
//					if(min[i] > values[i][j]) min[i] = values[i][j];
//				}
//			}
			
			for(int l=0; l<values.length; l++){
				x = lstart;
				int w = (int)Math.ceil(lineWidth / (double)values[l].length);
				
				for(int j=0; j<values[l].length; j++){
						vcolor = linear1DCcolor(values[l][j], min[l], max[l]);
					g2d.setColor(vcolor);
					g2d.fillRect(x, y, w, lineHeight);
					
					if(printValue){
						if(values[l][j] < 0)
							valueString = dformat.format(values[l][j]);
						else
							valueString = iformat.format(values[l][j]);
						fw = metrics.stringWidth(valueString);
						if(fw < w){
							compvcolor = complementaryColor(vcolor);
							g2d.setColor(compvcolor);
							g2d.drawString(valueString, x + ((w - fw)/2 ), y + fontHeight + ((lineHeight - fontHeight)/2));
						}
						else{
							printValueRow  =false;
						}
					}
					
					x += w;
				}
				if(printValueRow == false)
					printValue = false;
				
				y += lineHeight;
				
				if(l < values.length-1)
					y += lineSpace;
			}
//		}
//		else{
//			max = new double[1];
//			min = new double[1];
//			max[0] = values[0][0];
//			min[0] = values[0][0];
//			for(int i=0; i<values.length; i++){
//				for(int j=0; j<values[i].length; j++){
//					if(max[0] < values[i][j]) max[0] = values[i][j];
//					if(min[0] > values[i][j]) min[0] = values[i][j];
//				}
//			}
//			
//			for(int l=0; l<values.length; l++){
//				x = lstart;
//				int w = (int)Math.ceil(lineWidth / (double)values[l].length);
//				
//				for(int j=0; j<values[l].length; j++){
//					if(fromZeroGradient)
//						vcolor = linear1DCcolor(values[l][j], 0, max[0]);
//					else
//						vcolor = linear1DCcolor(values[l][j], min[0], max[0]);
//					g2d.setColor(vcolor);
//					g2d.fillRect(x, y, w, lineHeight);
//					
//					if(printValue){
//						if(values[l][j] < 0)
//							valueString = dformat.format(values[l][j]);
//						else
//							valueString = iformat.format(values[l][j]);
//						fw = metrics.stringWidth(valueString);
//						if(fw < w){
//							compvcolor = complementaryColor(vcolor);
//							g2d.setColor(compvcolor);
//							g2d.drawString(valueString, x + ((w - fw)/2 ), y + fontHeight + ((lineHeight - fontHeight)/2));
//						}
//						else{
//							printValueRow  =false;
//						}
//					}
//					
//					
//					x += w;
//				}
//				
//				if(printValueRow == false)
//					printValue = false;
//				
//				y += lineHeight;
//				
//				if(l < values.length-1)
//					y += lineSpace;
//			}
//		}
		
//		for(int l=0; l<values.length; l++){
//			x = lstart;
//			g2d.setColor(Color.black);
//			g2d.fillRect(x, y, lineWidth, lineHeight);
//			
//			y += lineHeight;
//			
//			if(l < values.length-1)
//				y += lineSpace;
//		}
		
		//print labels
		if(printLabels && labels != null){
			g2d.setColor(fontColor);
			for(int i=0; i<(labels.length <= values.length ? labels.length : values.length); i++){
				fw = metrics.stringWidth(labels[i]);
				
				g2d.drawString(labels[i], margin + (maxLabel-fw), ly + (lineHeight - ((lineHeight - fontHeight)/2)));
				
				ly += lineHeight + lineSpace;
			}
		}
		
		//print bottom ruler
		if(printBottomRuler){
			y += rulerSpace;
			for(int k=3; k>0; k--){
				x = lstart;
				int nof = (int)Math.pow(4, k);
				int w = lineWidth / nof;
				
				int c = 0;
				for(int i=0; i<nof; i++){
					g2d.setColor(ncolors[c]);
					c++; if(c>3) c=0;
					
					g2d.fillRect(x, y, w, rulerLineHeight);
					x += w;
				}
				 y += rulerLineHeight;
			}
		}
		
		//gradient legend
		if(printLegend){
			//fontHeight = metrics.getHeight();
			//fontHeight =  (int)Math.ceil(g2d.getFont().getMaxCharBounds(g2d.getFontRenderContext()).getHeight());
			fontHeight = (int)Math.ceil(g2d.getFont().createGlyphVector(g2d.getFontRenderContext(),"0,123456789").getVisualBounds().getHeight());
			if(fontHeight > legendHeight){
				fontHeight = legendHeight;
				font = new Font(font.getFamily(), font.getStyle(), legendHeight);
				g2d.setFont(font);
				metrics = g2d.getFontMetrics(font);
			}
			
			x = lstart;
			y += margin;
			
			if(printMultiLegend){
				ly = y;
				int iv;
				for(double v = 0; v<values.length; v++){
					iv = (int)v;
					x = lstart;
					
//					double range = max - min;
//		            double position = data - min;
//		            double percentPosition = position / range;
					
//					System.out.println("0 rmin["+rmin+"] rmax["+rmax+"] x["+((rmin[iv]) + (0*(max[iv]-min[iv])/(double)lineWidth))+"]");
//					System.out.println("4096 rmin["+rmin+"] rmax["+rmax+"] x["+((rmin[iv]) + (4096*(max[iv]-min[iv])/(double)lineWidth))+"]");
//					
					for(double i=0; i<lineWidth; i++){
						//g2d.setColor(linear1DCcolor(i, 0.0, lineWidth));
						g2d.setColor(linear1DCcolor(rmin[iv] + (i*(max[iv]-min[iv])/(double)lineWidth), min[iv], max[iv]));
						g2d.fillRect(x, y, 1, legendHeight);
						x += 1;
					}
					
					compvcolor = complementaryColor(linear1DCcolor(rmin[iv] + (0*(max[iv]-min[iv])/(double)lineWidth), min[iv], max[iv]));
					g2d.setColor(compvcolor);
					valueString = formatValue(rmin[iv]);
					g2d.drawString(valueString, lstart, y + legendHeight - (legendHeight - fontHeight));
					valueString = formatValue(rmax[iv]);
					fw = metrics.stringWidth(valueString);
					compvcolor = complementaryColor(linear1DCcolor(rmin[iv] + (4096*(max[iv]-min[iv])/(double)lineWidth), min[iv], max[iv]));
					g2d.setColor(compvcolor);
					g2d.drawString(valueString, lstart + lineWidth - fw, y + legendHeight - (legendHeight - fontHeight));
					
					
					y += legendHeight;
					y += lineSpace;
				}
				
				if(printLabels && labels != null){
					g2d.setColor(fontColor);
					for(int i=0; i<(labels.length <= values.length ? labels.length : values.length); i++){
						fw = metrics.stringWidth(labels[i]);
						
						g2d.drawString(labels[i], lineSpace + (maxLabel-fw), ly + (legendHeight - ((legendHeight - fontHeight)/2)));
						
						ly += legendHeight + lineSpace;
					}
				}
			}
			else{
				for(double i=0; i<lineWidth; i++){
					//g2d.setColor(getCellColour(lineWidth - i, 0.0, lineWidth));
					g2d.setColor(linear1DCcolor(i, 0.0, lineWidth));
					g2d.fillRect(x, y, 1, legendHeight);
					x += 1;
				}
				
				
				if(normalizeByRow){
					compvcolor = complementaryColor(linear1DCcolor(0.0, 0.0, lineWidth));
					g2d.setColor(compvcolor);
					g2d.drawString("0", lstart, y + legendHeight - (legendHeight - fontHeight));
					fw = metrics.stringWidth("max");
					compvcolor = complementaryColor(linear1DCcolor(lineWidth -1, 0.0, lineWidth));
					g2d.setColor(compvcolor);
					g2d.drawString("max", lstart + lineWidth - fw, y + legendHeight - (legendHeight - fontHeight));
				}
				else{
					compvcolor = complementaryColor(linear1DCcolor(0.0, 0.0, lineWidth));
					g2d.setColor(compvcolor);
					valueString = formatValue(min[0]);
					g2d.drawString(valueString, lstart, y + legendHeight - (legendHeight - fontHeight));
					valueString = formatValue(max[0]);
					fw = metrics.stringWidth(valueString);
					compvcolor = complementaryColor(linear1DCcolor(lineWidth -1, 0.0, lineWidth));
					g2d.setColor(compvcolor);
					g2d.drawString(valueString, lstart + lineWidth - fw, y + legendHeight - (legendHeight - fontHeight));
				}
			}
		}
	}
	
	
	private String formatValue(double value){
		if(value < 0)
			return dformat.format(value);
		return iformat.format(value);
	}
	private Color complementaryColor(Color c){
		return new Color(255 - c.getRed(), 255 - c.getGreen(), 255 - c.getBlue());
	}
	
	/*
     * Determines what colour a heat map cell should be based upon the cell 
     * values.
     */
    private Color getCellColour(double data, double min, double max) {              
            double range = max - min;
            double position = data - min;

            // What proportion of the way through the possible values is that.
            double percentPosition = position / range;
            
            // Which colour group does that put us in.
            int colourPosition = getColourPosition(percentPosition);
            
            int r = highValueColor.getRed();
            int g = highValueColor.getGreen();
            int b = highValueColor.getBlue();
            
            // Make n shifts of the colour, where n is the colourPosition.
            for (int i=0; i<colourPosition; i++) {
                    int rDistance = r - lowValueColor.getRed();
                    int gDistance = g - lowValueColor.getGreen();
                    int bDistance = b - lowValueColor.getBlue();
                    
                    if ((Math.abs(rDistance) >= Math.abs(gDistance))
                                            && (Math.abs(rDistance) >= Math.abs(bDistance))) {
                            // Red must be the largest.
                            r = changeColourValue(r, rDistance);
                    } else if (Math.abs(gDistance) >= Math.abs(bDistance)) {
                            // Green must be the largest.
                            g = changeColourValue(g, gDistance);
                    } else {
                            // Blue must be the largest.
                            b = changeColourValue(b, bDistance);
                    }
            }
            
            return new Color(r, g, b);
    }
    
    /*
     * Returns how many colour shifts are required from the lowValueColour to 
     * get to the correct colour position. The result will be different 
     * depending on the colour scale used: LINEAR, LOGARITHMIC, EXPONENTIAL.
     */
    private int getColourPosition(double percentPosition) {
            return (int) Math.round(colourValueDistance * Math.pow(percentPosition, colourScale));
    }
    
    private int changeColourValue(int colourValue, int colourDistance) {
            if (colourDistance < 0) {
                    return colourValue+1;
            } else if (colourDistance > 0) {
                    return colourValue-1;
            } else {
                    // This shouldn't actually happen here.
                    return colourValue;
            }
    }
    
    private void updateColourDistance() {
        int r1 = highValueColor.getRed();
        int g1 = highValueColor.getGreen();
        int b1 = highValueColor.getBlue();
        int r2 = lowValueColor.getRed();
        int g2 = lowValueColor.getGreen();
        int b2 = lowValueColor.getBlue();
        
        colourValueDistance = Math.abs(r1 - r2);
        colourValueDistance += Math.abs(g1 - g2);
        colourValueDistance += Math.abs(b1 - b2);
    }
    
    private int getAColor(double imax, double i){
//		int G = (int)((imax - i)*255/imax);
//		return 255<<16 | G<<8; 
		
		
		//int R = (int)((i)*255/imax);
		//return R<<16;
		
//		if(i<(imax/2)){
//			int R = (int)((i)*255/(imax/2));
//			return R;
//		}
//		else{
//			int R = (int)((i-(imax/2))*255/(imax/2));
//			R = R<<8;
//			R = R | (255);
//			return R;
//		}
		
		//int R = 0x0, G = 0x0, B = 0x0;
		int C = 0x0;
		
		if(i < imax/3){
			C = (int)((i)*255/(imax/3));
		}
		else{
			if(i < imax*2/3){
				C = (int)((i - (imax/3) )*255/(imax/3));
				C = (C << 8) | 0xff;
			}
			else{
				C = (int)((i - (imax*2/3) )*255/(imax/3));
				C = (C << 16) | 0xffff;
			}
		}
		
		return C;
	}
    
    
    private Color linear1DCcolor(double value, double min, double max){
    	
    	int r,g,b;
		if(medianValueColor != null){
			double v = (value - min) / (max - min);
//			System.out.println("min,max: "+min+","+max);
//			System.out.println("vmean: "+((min + ((max -min) / 2)* hlColorLimit)));
			double mean = (min + ((max -min) / 2)* hlColorLimit);
//			System.out.println("mean: "+mean);
			if(value < mean){
				value = Math.pow((value - min) / (mean - min), colourScale);
				if(value > 1) value = 1;
		    	if(value < 0) value = 0;
				r = lowValueColor.getRed() + (int)(value*(medianValueColor.getRed() - lowValueColor.getRed()));
	    		g = lowValueColor.getGreen() + (int)(value*(medianValueColor.getGreen() - lowValueColor.getGreen()));
	    		b = lowValueColor.getBlue() + (int)(value*(medianValueColor.getBlue() - lowValueColor.getBlue()));
			}
			else{
				value = Math.pow((value - mean) / (max - mean), colourScale);
				if(value > 1) value = 1;
		    	if(value < 0) value = 0;
				r = medianValueColor.getRed() + (int)(value*(highValueColor.getRed() - medianValueColor.getRed()));
	    		g = medianValueColor.getGreen() + (int)(value*(highValueColor.getGreen() - medianValueColor.getGreen()));
	    		b = medianValueColor.getBlue() + (int)(value*(highValueColor.getBlue() - medianValueColor.getBlue()));
			}
	    		
		}
		else{
	    	value = Math.pow((value - min) / (max - min), colourScale);
	    	if(value > 1)
	    		value = 1;
	    	if(value < 0)
	    		value = 0;
	    	
			r = lowValueColor.getRed() + (int)(value*(highValueColor.getRed() - lowValueColor.getRed()));
			g = lowValueColor.getGreen() + (int)(value*(highValueColor.getGreen() - lowValueColor.getGreen()));
			b = lowValueColor.getBlue() + (int)(value*(highValueColor.getBlue() - lowValueColor.getBlue()));
		}
    	return new Color(r,g,b);
    }
    
    
    
    public void configByMap(Map<String,String> config) throws Exception{
    	String value;
    	
    	value = config.get("margin");
    	if(value != null)
    		margin = Integer.parseInt(value);
    	
    	value = config.get("backgroundColor");
    	if(value != null){		
    		Field field = Class.forName("java.awt.Color").getField(value);
    		backgroundColor = (Color)field.get(null);
    	}
    	
    	
    	value = config.get("lineWidth");
    	if(value != null)
    		lineWidth = Integer.parseInt(value);
    	
    	value = config.get("lineHeight");
    	if(value != null)
    		lineHeight = Integer.parseInt(value);
    	
    	value = config.get("lineSpace");
    	if(value != null)
    		lineSpace = Integer.parseInt(value);
    	
    	value = config.get("printTopRuler");
    	if(value != null)
    		printTopRuler = Boolean.parseBoolean(value);
    	
    	value = config.get("printBottomRuler");
    	if(value != null)
    		printBottomRuler = Boolean.parseBoolean(value);
    	
    	value = config.get("rulerLineHeight");
    	if(value != null)
    		rulerLineHeight = Integer.parseInt(value);
    	
    	value = config.get("rulerSpace");
    	if(value != null)
    		rulerSpace = Integer.parseInt(value);
    	
    	
    	value = config.get("printLegend");
    	if(value != null)
    		printLegend = Boolean.parseBoolean(value);
    	
    	value = config.get(printMultiLegend);
    	if(value != null)
    		printMultiLegend = Boolean.parseBoolean(value);
    	
    	value = config.get("legendHeight");
    	if(value != null)
    		legendHeight = Integer.parseInt(value);
    	
    	//font
    	value = config.get("fontColor");
    	if(value != null){		
    		Field field = Class.forName("java.awt.Color").getField(value);
    		fontColor = (Color)field.get(null);
    	}
    
    	value = config.get("yAxisSpace");
    	if(value != null)
    		yAxisSpace = Integer.parseInt(value);
    	
    	
    	value = config.get("normalizeByRow");
    	if(value != null)
    		normalizeByRow = Boolean.parseBoolean(value);
    	
    	value = config.get("fromZeroGradient");
    	if(value != null)
    		fromZeroGradient = Boolean.parseBoolean(value);
    	
    	
    	value = config.get("labelsSpace");
    	if(value != null)
    		labelsSpace = Integer.parseInt(value);
    	
    	value = config.get("printLabels");
    	if(value != null)
    		printLabels = Boolean.parseBoolean(value);
    	
    	value = config.get("iformat");
    	if(value != null)
    		iformat = new DecimalFormat(value);
    	
    	value = config.get("dformat");
    	if(value != null)
    		dformat = new DecimalFormat(value);
    	
    	
    	value = config.get("colourScale");
    	if(value != null)
    		colourScale = Double.parseDouble(value);
    	
    	value = config.get("hlColorLimit");
    	if(value != null)
    		hlColorLimit = Double.parseDouble(value);
    	
    	value = config.get("highValueColor");
    	if(value != null){		
    		Field field = Class.forName("java.awt.Color").getField(value);
    		highValueColor = (Color)field.get(null);
    	}
    	value = config.get("lowValueColor");
    	if(value != null){		
    		Field field = Class.forName("java.awt.Color").getField(value);
    		lowValueColor = (Color)field.get(null);
    	}
    	value = config.get("medianValueColor");
    	if(value != null){
    		if(value.compareTo("null") == 0){
    			medianValueColor = null;
    		}
    		else{
    			Field field = Class.forName("java.awt.Color").getField(value);
    			medianValueColor = (Color)field.get(null);
    		}
    	}
    }
}
