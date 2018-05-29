/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package igtools.gui2;

import java.awt.Color;
import java.awt.Component;
import javax.swing.JLabel;
import javax.swing.JList;
import javax.swing.ListCellRenderer;

/**
 *
 * @author vbonnici
 */
public class WSSeqCellRenderer extends JLabel implements ListCellRenderer<WSSequence>{
  private static final Color HIGHLIGHT_COLOR = new Color(0, 0, 128);

  public WSSeqCellRenderer() {
    setOpaque(true);
    setIconTextGap(12);
  }

    @Override
    public Component getListCellRendererComponent(JList<? extends WSSequence> list, WSSequence value, int index, boolean isSelected, boolean cellHasFocus) {
        String text = value.getName();
        if(value.getB3seq()!=null || value.getNELSA()!=null){
            text =text + " (";
            if(value.getB3seq()!=null && value.getNELSA()!=null){
                text = text + "S,I";
            }
            else{
                if(value.getB3seq()!=null)
                    text = text + "S";
                if(value.getNELSA()!=null)
                    text = text + "S,I";
            }
            text = text + ")";
        }
        setText(text);
        
        if (isSelected) {
          setBackground(HIGHLIGHT_COLOR);
          setForeground(Color.white);
        } else {
          setBackground(Color.white);
          setForeground(Color.black);
        }
        return this;
    }
}
