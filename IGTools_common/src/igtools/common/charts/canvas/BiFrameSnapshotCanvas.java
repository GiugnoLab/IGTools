package igtools.common.charts.canvas;

import igtools.common.charts.palette.B3NucleotidePalette;
import igtools.common.charts.palette.DefaultB3NucleotidePalette;
import igtools.common.nucleotide.B3Nucleotide;

import java.awt.Color;
import java.awt.Graphics2D;
import java.awt.image.BufferedImage;
import java.awt.image.WritableRaster;
import java.util.LinkedList;
import java.util.List;

public class BiFrameSnapshotCanvas {
	
	private Color backgroundColor = Color.BLACK;
	private Color diffColor = Color.WHITE;
	private int[] diffColorComponents;
	
	private B3NucleotidePalette palette = new DefaultB3NucleotidePalette();
	
	private int k;
	private int width;
	
	private BufferedImage img; 
	private Graphics2D g2d;
	private WritableRaster raster;
	
	private BufferedImage diff_img; 
	private Graphics2D diff_g2d;
	private WritableRaster diff_raster;
	
	private List<String> ordlist = null;
	private List<String> prev_ordlist = null;
	
	
	public BiFrameSnapshotCanvas(int k, int width){
		this.k = k;
		this.width = width;
		
		this.img = new BufferedImage(width, k,  BufferedImage.TYPE_INT_RGB);
		g2d = (Graphics2D)this.img.getGraphics();
		raster = img.getRaster();
		
		this.diff_img = new BufferedImage(width, k,  BufferedImage.TYPE_INT_RGB);
		diff_g2d = (Graphics2D)this.diff_img.getGraphics();
		diff_raster = diff_img.getRaster();
		
		diffColorComponents = new int[4];
		diffColorComponents[0] = diffColor.getRed();
		diffColorComponents[1] = diffColor.getGreen();
		diffColorComponents[2] = diffColor.getBlue();
		diffColorComponents[3] = 255;
		
		ordlist = new LinkedList<String>();
		
	}	
	
	public BufferedImage img(){
		return this.img;
	}
	public BufferedImage diff_img(){
		return this.diff_img;
	}
	
	public void setBackgroundColor(Color c){
		this.backgroundColor = c;
	}
	public void setPalette(B3NucleotidePalette p){
		this.palette = p;
	}
	
	public void push(String s){
		ordlist.add(s);
	}
	
	public void clear(){
		g2d.setColor(backgroundColor);
		g2d.fillRect(0, 0, width, k);
		
		diff_g2d.setColor(backgroundColor);
		diff_g2d.fillRect(0, 0, width, k);
		
		if(prev_ordlist != null)
			prev_ordlist.clear();
		prev_ordlist = ordlist;
		ordlist = new LinkedList<String>();
	}
	
	private int hamming(String a, String b){
		int l = (a.length() <= b.length()) ? a.length() : b.length();
		int m = (a.length() >= b.length()) ? a.length() : b.length();
		int dist = 0;
		for(int i=0; i<l; i++){
			if(a.charAt(i) != b.charAt(i)) dist++;
		}
		dist += m-l;
		return dist;
	}
	
	
	
	public void print(){
		if(prev_ordlist != null){
			
			List<String> temp_list = new LinkedList<String>();
			
			boolean[] choosen = new boolean[ordlist.size()]; 
			int score;
			int min_score = Integer.MAX_VALUE;
			String min_s = "";
			for(String sp : prev_ordlist){
				min_score = Integer.MAX_VALUE;
				min_s = "";
				int min_i = 0;
				
				int i = 0;
				for(String s : ordlist){
					score = hamming(sp, s);
					
					if(!choosen[i] && score < min_score){
						min_score = score;
						min_s = s;
						min_i = i;
					}
					i++;
				}
				
				
				temp_list.add(min_s);
				choosen[min_i] = true;
			}
			
			ordlist.clear();
			ordlist = temp_list;
			
			B3Nucleotide n;
			int i = 0;
			for(String s : ordlist){
				String ps = prev_ordlist.get(i);
				
//				System.out.println("-");
//				System.out.println(s);
//				System.out.println(ps);
				
				for(int j=0; j<s.length(); j++){
					
					n = B3Nucleotide.by((char)(s.charAt(j)));
					raster.setPixel(i, j, palette.getColorComponents(n.code()));
					
					if(s.charAt(j) == ps.charAt(j))
						diff_raster.setPixel(i, j, diffColorComponents);
					
				}
				i++;
			}
			
		}
		else{
			B3Nucleotide n;
			int i = 0;
			for(String s : ordlist){
				for(int j=0; j<s.length(); j++){
					
					n = B3Nucleotide.by((char)(s.charAt(j)));
					raster.setPixel(i, j, palette.getColorComponents(n.code()));
					
				}
				i++;
			}
		}
	}

}
