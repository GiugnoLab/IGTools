/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package igtools.gui2.dictionaries;

import igtools.common.charts.palette.B3NucleotidePalette;
import igtools.common.charts.palette.DefaultB3NucleotidePalette;
import igtools.common.nucleotide.B3Nucleotide;
import igtools.dictionaries.elsa.IELSAIterator;
import igtools.dictionaries.elsa.INELSA;
import java.awt.Color;
import java.awt.Dimension;
import java.awt.Graphics;
import javax.swing.JPanel;

/**
 *
 * @author vbonnici
 */
public class FlatDictionaryPanel extends JPanel{
    
    
    private INELSA nelsa = null;
    
    private IELSAIterator it = null;
    private IELSAIterator it_end = null;
    
    private int pointWith = 0;
    private int pointHeight = 0;
    
    private boolean byMulti = false;
    boolean extendToLCP = true;
    
    private int rows = 0;
    private IELSAIterator it_last = null;
    
    
    public boolean printMultsBar = true;
    public static final Color barColor = Color.BLACK;
    public static final Color leftBarColor = new Color(33, 66, 99);
    public static final Color tenBarColor = new Color(64, 128, 191);
    public static final Color hapaxBarColor = Color.RED;
    private int mbar_length = 40;
    private int ambar_length = 20;
    private int mbar_r_space = 10;
    public double log_max_multi = -1;
    public int pk = -1;
    
    
    private Color[] colors = {
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
    
    public FlatDictionaryPanel(INELSA nelsa){
        super();
        this.nelsa =  nelsa;
    }
    
    public void setData(
            IELSAIterator it,
            IELSAIterator it_end,
            int pointWidth,
            int pointHeight,
            boolean heightByMultiplicity,
            boolean extendToLCP
            ){
        this.it = it;
        this.it_end = it_end;
        this. pointWith = pointWidth;
        this.pointHeight = pointHeight;
        this.byMulti = heightByMultiplicity;
        this.extendToLCP = extendToLCP;
    }
    
    public void setColorPalette(B3NucleotidePalette palette){
        this.colors[0] = palette.getColor(B3Nucleotide.A);
        this.colors[1] = palette.getColor(B3Nucleotide.C);
        this.colors[2] = palette.getColor(B3Nucleotide.G);
        this.colors[3] = palette.getColor(B3Nucleotide.T);
        this.colors[4] = palette.getColor(B3Nucleotide.N);
        
    }
    
    public IELSAIterator lastIterator(){
        return this.it_last;
    }
    
    
    private int max(int a, int b){
        if(a>b)
            return a;
        return b;
    }
    
    private int min(int a, int b){
        if(a<b)
            return a;
        return b;
    }
    
    @Override
    public void paintComponent(Graphics g) {
        super.paintComponent(g);
        
        if(printMultsBar){
            if(this.it!=null){
                if(this.it.k()!=pk){
                    pk = this.it.k();
                    log_max_multi = 0;
                    IELSAIterator it = nelsa.begin(this.it.k());
                    while(it.next()){
                        if(log_max_multi < it.multiplicity())
                            log_max_multi = it.multiplicity();
                    }
                    System.out.println("max_multi: "+log_max_multi);
                    log_max_multi = Math.log(log_max_multi);
                    System.out.println("log_max_multi: "+log_max_multi);
                }
            }
        }
        
        
//        System.out.println("["+this.getWidth()+","+this.getHeight()+"]");

        g.setColor(colors[5]);
        g.fillRect(0, 0, this.getWidth(), this.getHeight());
            
        if(this.getHeight() > 0 && this.pointHeight>0){
//            g.setColor(colors[6]);
//            g.fillRect(0, 0, this.getWidth(), this.getHeight());
            
            int mers_width = this.getWidth();
            int mbar_width = 0;
            if(printMultsBar){
                mbar_width = mbar_length + ambar_length + mbar_r_space;
                mers_width -= mbar_width;
            }
            

            int rows = this.getHeight() / pointHeight;
            int cols = mers_width / pointWith;

            IELSAIterator it = this.it.clone();
            int row = 0;
            int  k = it.k();
            B3Nucleotide[] ns = new B3Nucleotide[it.k()];
            int lcp, i_s, i_e, lcp_s, lcp_e,lcp_i, lcp_l;
            B3Nucleotide[] lcp_ns;
            
//            System.out.println(it.good());
            Color abar_color;
            double abar_perc;
            int l;
            int n;
            int m, rm;
            while(row < rows && it.good()){
                it.kmer(ns);

                rm = it.multiplicity();
                m = 1;
                if(byMulti)
                    m = rm;
                for(n=0; n<k; n++){
                    g.setColor(colors[ns[n].code()]);
                    g.fillRect(mbar_width + n*pointWith, row*pointHeight, pointWith, pointHeight * m);
                }
                
                if(printMultsBar){
                    if(rm <= 10){
                            if(rm == 1){
                                    abar_perc = 1;
                                    abar_color = hapaxBarColor;
                            }
                            else{
                                    abar_perc = 1;
                                    abar_color = tenBarColor;
                            }	
                    }
                    else{
                            abar_perc =  1 - (Math.log(rm) / log_max_multi);
                            abar_color = leftBarColor;
                    }
                    
                    
                    /*if(mbar_length - (int)Math.floor(((double)mbar_length) * ( (Math.log(rm) / log_max_multi))) < 0){
                        System.out.println("pos <0 ");
                        System.out.println("rm "+rm);
                        System.out.println("log_max "+log_max_multi);
                        System.out.println("fract "+(Math.log(rm) / log_max_multi));
                        System.out.println("res "+((double)mbar_length) * ( (Math.log(rm) / log_max_multi)));
                    }
                    if(mbar_length - (int)Math.floor(((double)mbar_length) * ( (Math.log(rm) / log_max_multi))) > 200){
                        System.out.println("pos > length ");
                        System.out.println("rm "+rm);
                        System.out.println("log_max "+log_max_multi);
                        System.out.println("fract "+(Math.log(rm) / log_max_multi));
                        System.out.println("res "+((double)mbar_length) * ( (Math.log(rm) / log_max_multi)));
                    }
                    if((int)Math.floor(((double)mbar_length) * ( (Math.log(rm) / log_max_multi))) < 0){
                        System.out.println("w <0 ");
                        System.out.println("rm "+rm);
                        System.out.println("log_max "+log_max_multi);
                        System.out.println("fract "+(Math.log(rm) / log_max_multi));
                        System.out.println("res "+((double)mbar_length) * ( (Math.log(rm) / log_max_multi)));
                    }
                    if((int)Math.floor(((double)mbar_length) * ( (Math.log(rm) / log_max_multi))) > 200){
                        System.out.println("w > length ");
                        System.out.println("rm "+rm);
                        System.out.println("log_max "+log_max_multi);
                        System.out.println("fract "+(Math.log(rm) / log_max_multi));
                        System.out.println("res "+((double)mbar_length) * ( (Math.log(rm) / log_max_multi)));
                    }*/
                    
                    
                    g.setColor(barColor);
                    g.fillRect(
                            mbar_length - (int)Math.floor(((double)mbar_length) * ( (Math.log(rm) / log_max_multi))), 
                            row*pointHeight, 
                            (int)Math.floor(((double)mbar_length) * ( (Math.log(rm) / log_max_multi))), 
                            pointHeight * m);
                    g.setColor(abar_color);
                    g.fillRect(
                            mbar_length, 
                            row*pointHeight, 
                            (int)Math.floor(ambar_length * abar_perc), 
                            pointHeight * m);
                }
                
                
                if(extendToLCP){
                    i_s = it.istart();
                    i_e = it.iend();
                    
                    lcp = min(nelsa.lcp()[i_s], nelsa.ns()[i_s]);
                    lcp_i = nelsa.sa()[i_s];
                    for(l = i_s+1; l<=i_e; l++){
                        lcp_l = min(nelsa.lcp()[l], nelsa.ns()[l]);
                        if(lcp_l > lcp){
                            lcp = lcp_l;
                            lcp_i = nelsa.sa()[l];
                        }
                    }
                    
                    
//                    lcp_s = min(nelsa.lcp()[i_s], nelsa.ns()[i_s]);
//                    lcp_e = min(nelsa.lcp()[i_e], nelsa.ns()[i_e]);
//                    if(lcp_s >  lcp_e){
//                        lcp = lcp_s;
//                        lcp_i = i_s;
//                    }
//                    else{
//                        lcp = lcp_e;
//                        lcp_i = i_e;
//                    }
//                    
//                    System.out.println(lcp+"\t"+cols);
                    
                    if(lcp > cols)
                        lcp = cols;
                    
                    if(lcp > k){
                     lcp_ns = new B3Nucleotide[lcp];
                     nelsa.b3seq().getB3(lcp_i, lcp_ns);
                     for(n=0; n<lcp; n++){
                            g.setColor(colors[lcp_ns[n].code()]);
                            g.fillRect(mbar_width + n*pointWith, row*pointHeight, pointWith, pointHeight * m);
                        }
                    }
                }
                
                
                //row++;
                row += m;
                
                if(!it.next()){
                    it = null;
                    break;
                }
            }
            
            this.rows = row;
            this.it_last = it;
        }
    }
    
    public String getPickedKmer(int y){
        if(this.getHeight() > 0 && this.pointHeight>0){

            int pickedrow = y / pointHeight;
            //int rows = this.getHeight() / pointHeight;

            IELSAIterator it = this.it.clone();
            int row = 0;
            int  k = it.k();
            B3Nucleotide[] ns = new B3Nucleotide[it.k()];

//            System.out.println(it.good());
            
            int n;
            int m;
            while(row < pickedrow && it.good()){
//                it.kmer(ns);
//
//                for(n=0; n<k; n++){
//                    g.setColor(colors[ns[n].code()]);
//                    g.fillRect(n*pointWith, row*pointHeight, pointWith, pointHeight);
//                    
//                }
                
                //row++;
                m = 1;
                if(byMulti)
                    m = it.multiplicity();
                row += m;
                
                if(!it.next())
                    break;
            }
            
            if(it.good()){
                ns = it.kmer();
                return B3Nucleotide.toString(ns);
            }
        }
        return "";
    }
    
    
    
    
    public IELSAIterator  pageDownIterator(){
        if(this.getHeight() > 0 && this.pointHeight>0){
            int rows = this.getHeight() / pointHeight;
            IELSAIterator it = this.it.clone();
            int row = 0;
            int m;
            while(row < rows && it.good()){
                m = 1;
                if(byMulti)
                    m = it.multiplicity();
                row += m;
                
                if(!it.next()){
                    it = null;
                    break;
                }
            }
            return it;
        }
        return null;
    }
    
    public IELSAIterator  pageUpIterator(){
        if(this.getHeight() > 0 && this.pointHeight>0){
            int rows = this.getHeight() / pointHeight;
            IELSAIterator it = this.it.clone();
            int row = 0;
            int m;
            while(row < rows && it.good()){
                m = 1;
                if(byMulti)
                    m = it.multiplicity();
                row += m;
                
                if(!it.prev()){
                    it = null;
                    break;
                }
            }
            return it;
        }
        return null;
    }
    
    
    
}
