/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package igtools.gui2.sequence;

import igtools.common.charts.palette.B3NucleotidePalette;
import igtools.common.charts.palette.DefaultB3NucleotidePalette;
import igtools.common.nucleotide.B3Nucleotide;
import igtools.common.sequence.B3Sequence;
import java.awt.Color;
import java.awt.Graphics;
import javax.swing.JPanel;

/**
 *
 * @author vbonnici
 */
public class GenomePanel extends JPanel{
    private Color[] colors = {
        //Color.blue,         //0     A
        DefaultB3NucleotidePalette.A_COLOR,
        //Color.red,          //1     C
        DefaultB3NucleotidePalette.C_COLOR,
        //Color.green,        //2     G
        DefaultB3NucleotidePalette.G_COLOR,
        //Color.yellow,       //3     T
        DefaultB3NucleotidePalette.T_COLOR,
        //Color.lightGray,    //4     N
        DefaultB3NucleotidePalette.BAD_COLOR,
        
        Color.white,        //5     BACKGROUND
        Color.orange        //6     CONTINUE
    };
    
    
    private B3Sequence b3seq = null;
    private int b3seqLength = 0;
    private int sp = 0;
    private int to = 0;
    private int nW = 1, nH = 1;
    
    private GenomeFrame gframe = null;
    
    
    private int horPadding = 10;
    private int vertPadding = 5;
    
    
    private int upr = 1;
    private int rows  = 1;
    
    public GenomePanel(B3Sequence b3seq, GenomeFrame gframe){
        super();
        this.b3seq = b3seq;
        this.b3seqLength = b3seq.length();
        this.gframe = gframe;
    }
    
    
    
    public void setColorPalette(B3NucleotidePalette palette){
        this.colors[0] = palette.getColor(B3Nucleotide.A);
        this.colors[1] = palette.getColor(B3Nucleotide.C);
        this.colors[2] = palette.getColor(B3Nucleotide.G);
        this.colors[3] = palette.getColor(B3Nucleotide.T);
        this.colors[4] = palette.getColor(B3Nucleotide.N);
        
    }
    
    
    public int getUnitPerRow(){
        return upr;
    }
    public int rowCount(){
        return rows;
    }
    
    public void setOptions(int nW, int nH){
        this.nW = nW;
        this.nH = nH;
    }
    
    public int seqPointer(){
        return sp;
    }
    
    public void pageDown(){
        if(to != b3seqLength)
            sp = to;
        this.repaint();
    }
    
    public void pageUp(){
        if(sp > 0){
            int width = this.getWidth();
            int height = this.getHeight();
            int rows = (height- vertPadding*2) / nH;
            int unitPerRow = (width - horPadding*2) / nW;
            sp -= rows * unitPerRow;
            if(sp < 0){
                sp = 0;
            }
//            to = sp + (rows * unitPerRow);
//            if(to > b3seqLength)
//                to = b3seqLength;
            this.repaint();
        }
    }
    
    
    public void down(int step){
        if(to != b3seqLength){
            if(sp + step < b3seqLength){
                sp += step;
                this.repaint();
            }
        }
    }
    
    
    public void up(int step){
        if(sp > 0){
            sp -= step;
            if(sp < 0)
                sp = 0;
            this.repaint();
        }
    }
    
    
    @Override
    public void paintComponent(Graphics g) {
        super.paintComponent(g);
        
//        System.out.println("["+this.getWidth()+","+this.getHeight()+"]");

        g.setColor(colors[5]);
        g.fillRect(0, 0, this.getWidth(), this.getHeight());
        
        int width = this.getWidth();
        int height = this.getHeight();
        //B3Nucleotide[] ns;
        int nn;
        
        
        if(width>0 && height>0){
        
            int rows = (height- vertPadding*2) / nH;
            int unitPerRow = (width - horPadding*2) / nW;
            
            this.rows = rows;
            this.upr = unitPerRow;
            
            to = sp + (rows * unitPerRow);
            if(to > b3seqLength)
                to = b3seqLength;
            
//            System.out.println("nW: " + nW+"; nH: "+nH);
//            System.out.println("width: " + width+"; height: "+height);
//            System.out.println("rows: " + rows+"; unitPerRow: "+unitPerRow);
//            System.out.println("from: "+sp+"; to: "+to);
            
            int x=horPadding, y=vertPadding;
            
            int i = 0;
            for(int tsp = sp; tsp < to; ){
                
//                System.out.println(tsp+"\t"+i+"\t"+x+"\t"+y);
                
                g.setColor(colors[b3seq.getB3(tsp)]);
                g.fillRect(x, y, nW, nH);
                
                tsp++;
                i++;
                if(i == unitPerRow){
                    y += nH;
                    x = horPadding;
                    i=0;
                }
                else{
                    x += nW;
                }
            }
            
            
            gframe.setInfo(sp, to, upr, this.rows);
        }
    }
}
