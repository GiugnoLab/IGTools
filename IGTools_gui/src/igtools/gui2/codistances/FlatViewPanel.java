/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package igtools.gui2.codistances;

import igtools.common.nucleotide.B3Nucleotide;
import java.awt.Color;
import java.awt.Dimension;
import java.awt.Graphics;
import java.awt.Rectangle;
import javax.swing.JPanel;
import javax.swing.JScrollPane;
import javax.swing.Scrollable;

/**
 *
 * @author vbonnici
 */
public class FlatViewPanel extends JPanel{
    private final static Color[] colors = {
        Color.blue,         //0     A,
        Color.red,          //1     C
        Color.green,        //2     G
        Color.yellow,       //3     T
        Color.black        //4     N
    };
    private final static Color backgroundColor = Color.white;
    
    
    private String[] kmers = null;
    private int ikmer = 0;
    
    private int nW = 4, nH = 1;
    private int leftOffset = 10;
    private int prefW = 200, prefH = 200;
    private int dPrefH = 200;
    
    
    public FlatViewPanel(String[] kmers){
        setKmers(kmers);
    }
    
    public void setData(int nW, int nH){
        this.nW = nW;
        this.nH = nH;
    }
    public void setKmers(String[] kmers){
        this.kmers = kmers;
        prefSize();
    }
    
    public void prefSize(){
        prefW =leftOffset + (kmers[0].length() * nW);
        prefH = kmers.length * nH;
        
        this.setSize(new Dimension(prefW,dPrefH));
        this.setMinimumSize(new Dimension(prefW,dPrefH));
        this.setMaximumSize(new Dimension(prefW,prefH));
        this.setPreferredSize(new Dimension(prefW,dPrefH));
        
        this.revalidate();
    }
        
    
    @Override
    public void paintComponent(Graphics g) {
        super.paintComponent(g);
        
        

        g.setColor(backgroundColor);
        g.fillRect(0, 0, this.getWidth(), this.getHeight());
        
        int x = leftOffset, y = 0;
        
        char[] cs;
        for(int i=ikmer; i<kmers.length; i++){
            cs = kmers[i].toCharArray();
            
            x = leftOffset;
            for(int j=0; j<cs.length; j++){
                g.setColor(colors[B3Nucleotide.codeFor(cs[j])]);
                g.fillRect(x, y, nW, nH);
                x += nW;
            }
            
            y += nH;
            if(y>this.getHeight())
                break;
        }
        
    }
    
    
    
    public int rowCount(){
        return this.getHeight() / nH;
    }
    
    public void up(int step){
        if(ikmer > 0){
            ikmer -= step;
            if(ikmer < 0)
                ikmer = 0;
            this.repaint();
        }
    }
    
    public void down(int step){
        ikmer += step;
        if(ikmer >= kmers.length - (this.getHeight()/nH)){
            ikmer = (kmers.length - 1) - (this.getHeight()/nH);
            if(ikmer < 0)
                ikmer = 0;
        }
        this.repaint();
    }
     
}
