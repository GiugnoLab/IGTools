/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package igtools.gui2.codistances;

import java.awt.Color;
import java.awt.Component;
import java.awt.Dimension;
import java.awt.Font;
import java.awt.FontMetrics;
import java.awt.Graphics;
import javax.swing.JComponent;
import javax.swing.JTable;
import javax.swing.table.TableCellRenderer;

/**
 *
 * @author vbonnici
 */
public class KmerLabel extends JComponent  implements TableCellRenderer{

    private String text = "";
    private int row = 0;
    private boolean selected = false;
//    private FontMetrics fm = null;
    
    public KmerLabel(){
        super();
        Font f = new java.awt.Font("Andale Mono", 0, 15);
        this.setFont(f);
//        fm = f.get.getFontMetrics();
    }
    
    
    public void setText(String text){
        this.text = text;
        FontMetrics fm = this.getFontMetrics(this.getFont());
        int w = fm.charWidth('A');
        int h = this.getFont().getSize();
        this.setPreferredSize(new Dimension(w * text.length(),h));
        this.setMinimumSize(new Dimension(w * text.length(),h));

        this.row = row;
    }
    @Override
    public Component getTableCellRendererComponent(JTable table, Object value, boolean isSelected, boolean hasFocus, int row, int column) {
//            this.text = value.toString();
//
//            FontMetrics fm = this.getFontMetrics(this.getFont());
//            int w = fm.charWidth('A');
//            int h = this.getFont().getSize();
//            this.setPreferredSize(new Dimension(w * text.length(),h));
//            this.setMinimumSize(new Dimension(w * text.length(),h));
//            
//            this.row = row;
        this.setText(value.toString());
            this.selected = isSelected;
            return this;
    }
    
    @Override
    public void paintComponent(Graphics g) {
        super.paintComponent(g);       
        
        char[] cc = this.text.toCharArray();
        //FontMetrics fm = this.getFontMetrics(this.getFont());
        FontMetrics fm = g.getFontMetrics();
        
        g.setColor(this.row%2==0 ? Color.white : new Color(224,224,224));
        if(this.selected) g.setColor(new Color(33,66,99));
        //g.fillRect(1,1,this.getWidth()-1,this.getHeight()-1);
        g.fillRect(1,0,this.getWidth(),this.getHeight());
        
        int y = this.getHeight();
        int h = this.getFont().getSize();
        int x = 1;
        int w ;
        for(int i=0; i<cc.length; i++){
            switch (cc[i]){
                case 'A':
                    g.setColor(Color.blue);
                    break;
                case 'C':
                    g.setColor(Color.red);
                    break;
                case 'G':
                    g.setColor(Color.green);
                    break;
                case 'T':
                    g.setColor(Color.yellow);
                    break;
                case 'N':
                    g.setColor(Color.lightGray);
                    break;
                case '-':
                    g.setColor(Color.black);
                    break;  
                default:
                    g.setColor(Color.white);
            }
            w = fm.charWidth(cc[i]);
//            System.out.println("("+x+","+ y+","+ w+","+ h+")");
            g.fillRect(x, 1, w, h);
            x += w;
        }
       
        
        
        
        g.setColor(this.getForeground());
        g.drawString(this.text,1, this.getFont().getSize()-1);
        
    }
    
}
