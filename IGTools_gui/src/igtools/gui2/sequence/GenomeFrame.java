/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package igtools.gui2.sequence;

import igtools.common.charts.palette.B3NucleotidePaletteFactory;
import igtools.common.sequence.B3Sequence;
import igtools.gui2.WSSequence;
import java.awt.BorderLayout;
import java.awt.Color;
import java.awt.event.ActionEvent;
import java.awt.event.KeyEvent;
import java.awt.event.MouseWheelEvent;
import java.awt.event.MouseWheelListener;
import java.text.DecimalFormat;
import java.util.Set;
import java.util.TreeSet;
import javax.swing.AbstractAction;
import javax.swing.ActionMap;
import javax.swing.DefaultComboBoxModel;
import javax.swing.InputMap;
import javax.swing.JComponent;
import javax.swing.KeyStroke;

/**
 *
 * @author vbonnici
 */
public class GenomeFrame extends javax.swing.JFrame {
    private final static DecimalFormat df = new DecimalFormat("###,###,###,###");
    
    
    private B3Sequence b3seq = null;
    private int b3seqLength = 0;
    private GenomePanel gPanel = null;
    
    
    private class ColorItem{
        public String name;
        public Color color;
        public ColorItem(Color color, String name){
            this.color = color;
            this.name = name;
        }
        @Override
        public String toString(){
            return name;
        }
    }

    /**
     * Creates new form GenomeFrame
     */
    public GenomeFrame(WSSequence wsseq) {
       
        initComponents();
        
        this.setTitle("Sequence viewer: " + wsseq.getName());
        this.b3seq = wsseq.getB3seq();
        
        this.b3seqLength = b3seq.length();
        
        gPanel = new GenomePanel(b3seq, this);
        this.centralPanel.add(gPanel, BorderLayout.CENTER);
        setPanelOption();
        
        
        
        
        
        
        this.centralPanel.addMouseWheelListener(
                new MouseWheelListener() {

            @Override
            public void mouseWheelMoved(MouseWheelEvent e) {
                int notches = e.getWheelRotation();
                if (e.getScrollType() == MouseWheelEvent.WHEEL_UNIT_SCROLL) {
                    if(e.getUnitsToScroll()<0){
                        prev((int)(gPanel.rowCount() *0.25)*gPanel.getUnitPerRow());
                    }
                    else{
                        next((int)(gPanel.rowCount() *0.25)*gPanel.getUnitPerRow());
                    }
                }
             }
        });
        
        InputMap inputMap; 
        ActionMap actionMap;
        AbstractAction pageup_action;
        AbstractAction pagedown_action;

        inputMap  = this.centralPanel.getInputMap(JComponent.WHEN_IN_FOCUSED_WINDOW);
        actionMap = this.centralPanel.getActionMap();

        inputMap.put(KeyStroke.getKeyStroke(KeyEvent.VK_PAGE_UP, 0), "genomeViewPageUp");
        actionMap.put("genomeViewPageUp", new AbstractAction(){
           @Override
           public void actionPerformed(ActionEvent e){pageUp();}
        });
        inputMap.put(KeyStroke.getKeyStroke(KeyEvent.VK_PAGE_DOWN, 0), "genomeViewPageDown");
        actionMap.put("genomeViewPageDown", new AbstractAction(){
           @Override
           public void actionPerformed(ActionEvent e){pageDown();}
        });
        inputMap.put(KeyStroke.getKeyStroke(KeyEvent.VK_UP, 0), "genomeViewLineUp");
        actionMap.put("genomeViewLineUp", new AbstractAction(){
           @Override
           public void actionPerformed(ActionEvent e){lineUp();}
        });
        inputMap.put(KeyStroke.getKeyStroke(KeyEvent.VK_DOWN, 0), "genomeViewLineDown");
        actionMap.put("genomeViewLineDown", new AbstractAction(){
           @Override
           public void actionPerformed(ActionEvent e){lineDown();}
        });
        
        
        
    }
    
    
    
    private void setPanelOption(){
        int nW = Integer.parseInt(this.jTextField2.getText());
        int nH = Integer.parseInt(this.jTextField1.getText());
        this.gPanel.setOptions(nW, nH);
    }
    
    public void setInfo(int from, int to, int cols, int rows){
        this.fromField.setText(df.format(from));
        this.toField.setText(df.format(to));
        this.nofColumnsField.setText(df.format(cols));
        this.nofRowsField.setText(df.format(rows));
    }
    
    
    private void prev(int step){
        this.setEnabled(false);
        gPanel.up(step);
        this.setEnabled(true);
    }
    
    private void next(int step){
        this.setEnabled(false);
        gPanel.down(step);
        this.setEnabled(true);
    }
    
    private void pageUp(){
        this.setEnabled(false);
        gPanel.pageUp();
        this.setEnabled(true);
    }
    private void pageDown(){
        this.setEnabled(false);
        gPanel.pageDown();
        this.setEnabled(true);
    }
    private void lineUp(){
        this.prev(1);
    }
    private void lineDown(){
        this.next(1);
    }
    

    /**
     * This method is called from within the constructor to initialize the form.
     * WARNING: Do NOT modify this code. The content of this method is always
     * regenerated by the Form Editor.
     */
    @SuppressWarnings("unchecked")
    // <editor-fold defaultstate="collapsed" desc="Generated Code">//GEN-BEGIN:initComponents
    private void initComponents() {

        jPanel3 = new javax.swing.JPanel();
        jButton1 = new javax.swing.JButton();
        jTextField2 = new javax.swing.JTextField();
        jLabel2 = new javax.swing.JLabel();
        jLabel1 = new javax.swing.JLabel();
        jTextField1 = new javax.swing.JTextField();
        fromField = new javax.swing.JLabel();
        toField = new javax.swing.JLabel();
        jLabel5 = new javax.swing.JLabel();
        jLabel6 = new javax.swing.JLabel();
        jSeparator1 = new javax.swing.JSeparator();
        jLabel7 = new javax.swing.JLabel();
        nofColumnsField = new javax.swing.JLabel();
        nofRowsField = new javax.swing.JLabel();
        jLabel8 = new javax.swing.JLabel();
        jSeparator2 = new javax.swing.JSeparator();
        jLabel9 = new javax.swing.JLabel();
        aComboBox = new javax.swing.JComboBox();
        jLabel10 = new javax.swing.JLabel();
        cComboBox = new javax.swing.JComboBox();
        jLabel11 = new javax.swing.JLabel();
        gComboBox = new javax.swing.JComboBox();
        jLabel12 = new javax.swing.JLabel();
        jLabel13 = new javax.swing.JLabel();
        tComboBox = new javax.swing.JComboBox();
        jLabel14 = new javax.swing.JLabel();
        nComboBox = new javax.swing.JComboBox();
        jButton2 = new javax.swing.JButton();
        jSeparator3 = new javax.swing.JSeparator();
        centralPanel = new javax.swing.JPanel();

        setDefaultCloseOperation(javax.swing.WindowConstants.EXIT_ON_CLOSE);
        setMinimumSize(new java.awt.Dimension(200, 400));

        jButton1.setText("Refresh");
        jButton1.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                jButton1ActionPerformed(evt);
            }
        });

        jTextField2.setText("4");
        jTextField2.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                jTextField2ActionPerformed(evt);
            }
        });

        jLabel2.setText("w");

        jLabel1.setText("h");

        jTextField1.setText("1");

        fromField.setText("0");

        toField.setText("0");

        jLabel5.setText("From");

        jLabel6.setText("To");

        jLabel7.setText("Columns");

        nofColumnsField.setText("0");

        nofRowsField.setText("0");

        jLabel8.setText("Rows");

        jLabel9.setText("Color Palette");

        aComboBox.setModel(getColorListModel());
        aComboBox.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                aComboBoxActionPerformed(evt);
            }
        });

        jLabel10.setText("A");

        cComboBox.setModel(getColorListModel());
        cComboBox.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                cComboBoxActionPerformed(evt);
            }
        });

        jLabel11.setText("C");

        gComboBox.setModel(getColorListModel());
        gComboBox.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                gComboBoxActionPerformed(evt);
            }
        });

        jLabel12.setText("G");

        jLabel13.setText("T");

        tComboBox.setModel(getColorListModel());
        tComboBox.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                tComboBoxActionPerformed(evt);
            }
        });

        jLabel14.setText("N");

        nComboBox.setModel(getColorListModel());
        nComboBox.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                nComboBoxActionPerformed(evt);
            }
        });

        jButton2.setText("Apply");
        jButton2.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                jButton2ActionPerformed(evt);
            }
        });

        javax.swing.GroupLayout jPanel3Layout = new javax.swing.GroupLayout(jPanel3);
        jPanel3.setLayout(jPanel3Layout);
        jPanel3Layout.setHorizontalGroup(
            jPanel3Layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGroup(jPanel3Layout.createSequentialGroup()
                .addContainerGap()
                .addGroup(jPanel3Layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
                    .addComponent(jSeparator1)
                    .addGroup(jPanel3Layout.createSequentialGroup()
                        .addGroup(jPanel3Layout.createParallelGroup(javax.swing.GroupLayout.Alignment.TRAILING)
                            .addGroup(jPanel3Layout.createSequentialGroup()
                                .addComponent(jLabel2)
                                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                                .addComponent(jTextField2, javax.swing.GroupLayout.PREFERRED_SIZE, 50, javax.swing.GroupLayout.PREFERRED_SIZE))
                            .addComponent(jLabel5, javax.swing.GroupLayout.Alignment.LEADING)
                            .addComponent(jLabel6, javax.swing.GroupLayout.Alignment.LEADING))
                        .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED, javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE)
                        .addGroup(jPanel3Layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
                            .addGroup(jPanel3Layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
                                .addGroup(jPanel3Layout.createSequentialGroup()
                                    .addComponent(jLabel1)
                                    .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                                    .addComponent(jTextField1, javax.swing.GroupLayout.PREFERRED_SIZE, 50, javax.swing.GroupLayout.PREFERRED_SIZE))
                                .addComponent(fromField, javax.swing.GroupLayout.Alignment.TRAILING))
                            .addComponent(toField, javax.swing.GroupLayout.Alignment.TRAILING)))
                    .addGroup(javax.swing.GroupLayout.Alignment.TRAILING, jPanel3Layout.createSequentialGroup()
                        .addGroup(jPanel3Layout.createParallelGroup(javax.swing.GroupLayout.Alignment.TRAILING)
                            .addComponent(jLabel7, javax.swing.GroupLayout.Alignment.LEADING)
                            .addComponent(jLabel8, javax.swing.GroupLayout.Alignment.LEADING))
                        .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED, javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE)
                        .addGroup(jPanel3Layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
                            .addComponent(nofColumnsField)
                            .addComponent(nofRowsField, javax.swing.GroupLayout.Alignment.TRAILING)))
                    .addComponent(jSeparator2)
                    .addGroup(jPanel3Layout.createSequentialGroup()
                        .addComponent(jLabel10)
                        .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                        .addComponent(aComboBox, 0, javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE))
                    .addGroup(jPanel3Layout.createSequentialGroup()
                        .addComponent(jLabel11)
                        .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                        .addComponent(cComboBox, 0, javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE))
                    .addGroup(jPanel3Layout.createSequentialGroup()
                        .addComponent(jLabel12)
                        .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                        .addComponent(gComboBox, 0, javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE))
                    .addGroup(jPanel3Layout.createSequentialGroup()
                        .addComponent(jLabel13)
                        .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                        .addComponent(tComboBox, 0, javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE))
                    .addGroup(jPanel3Layout.createSequentialGroup()
                        .addComponent(jLabel14)
                        .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                        .addComponent(nComboBox, 0, javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE))
                    .addGroup(javax.swing.GroupLayout.Alignment.TRAILING, jPanel3Layout.createSequentialGroup()
                        .addGap(0, 0, Short.MAX_VALUE)
                        .addComponent(jButton2))
                    .addComponent(jSeparator3)
                    .addGroup(jPanel3Layout.createSequentialGroup()
                        .addGroup(jPanel3Layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
                            .addComponent(jButton1)
                            .addComponent(jLabel9))
                        .addGap(0, 0, Short.MAX_VALUE)))
                .addContainerGap())
        );
        jPanel3Layout.setVerticalGroup(
            jPanel3Layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGroup(jPanel3Layout.createSequentialGroup()
                .addContainerGap()
                .addComponent(jButton1)
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                .addGroup(jPanel3Layout.createParallelGroup(javax.swing.GroupLayout.Alignment.BASELINE)
                    .addComponent(jTextField1, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE)
                    .addComponent(jLabel1)
                    .addComponent(jTextField2, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE)
                    .addComponent(jLabel2))
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                .addComponent(jSeparator2, javax.swing.GroupLayout.PREFERRED_SIZE, 10, javax.swing.GroupLayout.PREFERRED_SIZE)
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                .addGroup(jPanel3Layout.createParallelGroup(javax.swing.GroupLayout.Alignment.BASELINE)
                    .addComponent(fromField)
                    .addComponent(jLabel5))
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                .addGroup(jPanel3Layout.createParallelGroup(javax.swing.GroupLayout.Alignment.BASELINE)
                    .addComponent(toField)
                    .addComponent(jLabel6))
                .addGap(18, 18, 18)
                .addGroup(jPanel3Layout.createParallelGroup(javax.swing.GroupLayout.Alignment.BASELINE)
                    .addComponent(nofColumnsField)
                    .addComponent(jLabel7))
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                .addGroup(jPanel3Layout.createParallelGroup(javax.swing.GroupLayout.Alignment.BASELINE)
                    .addComponent(nofRowsField)
                    .addComponent(jLabel8))
                .addGap(18, 18, 18)
                .addComponent(jSeparator1, javax.swing.GroupLayout.PREFERRED_SIZE, 10, javax.swing.GroupLayout.PREFERRED_SIZE)
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                .addComponent(jLabel9)
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                .addGroup(jPanel3Layout.createParallelGroup(javax.swing.GroupLayout.Alignment.BASELINE)
                    .addComponent(aComboBox, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE)
                    .addComponent(jLabel10))
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                .addGroup(jPanel3Layout.createParallelGroup(javax.swing.GroupLayout.Alignment.BASELINE)
                    .addComponent(cComboBox, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE)
                    .addComponent(jLabel11))
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                .addGroup(jPanel3Layout.createParallelGroup(javax.swing.GroupLayout.Alignment.BASELINE)
                    .addComponent(gComboBox, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE)
                    .addComponent(jLabel12))
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                .addGroup(jPanel3Layout.createParallelGroup(javax.swing.GroupLayout.Alignment.BASELINE)
                    .addComponent(tComboBox, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE)
                    .addComponent(jLabel13))
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                .addGroup(jPanel3Layout.createParallelGroup(javax.swing.GroupLayout.Alignment.BASELINE)
                    .addComponent(nComboBox, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE)
                    .addComponent(jLabel14))
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                .addComponent(jButton2)
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                .addComponent(jSeparator3, javax.swing.GroupLayout.PREFERRED_SIZE, 10, javax.swing.GroupLayout.PREFERRED_SIZE)
                .addContainerGap(108, Short.MAX_VALUE))
        );

        getContentPane().add(jPanel3, java.awt.BorderLayout.LINE_START);

        centralPanel.setMinimumSize(new java.awt.Dimension(20, 200));
        centralPanel.setPreferredSize(new java.awt.Dimension(600, 200));
        centralPanel.setLayout(new java.awt.BorderLayout());
        getContentPane().add(centralPanel, java.awt.BorderLayout.CENTER);

        pack();
    }// </editor-fold>//GEN-END:initComponents

    private void jTextField2ActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_jTextField2ActionPerformed
        // TODO add your handling code here:
    }//GEN-LAST:event_jTextField2ActionPerformed

    private void jButton1ActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_jButton1ActionPerformed
        setPanelOption();
        gPanel.repaint();
    }//GEN-LAST:event_jButton1ActionPerformed

    private void aComboBoxActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_aComboBoxActionPerformed
        // TODO add your handling code here:
    }//GEN-LAST:event_aComboBoxActionPerformed

    private void cComboBoxActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_cComboBoxActionPerformed
        // TODO add your handling code here:
    }//GEN-LAST:event_cComboBoxActionPerformed

    private void gComboBoxActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_gComboBoxActionPerformed
        // TODO add your handling code here:
    }//GEN-LAST:event_gComboBoxActionPerformed

    private void tComboBoxActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_tComboBoxActionPerformed
        // TODO add your handling code here:
    }//GEN-LAST:event_tComboBoxActionPerformed

    private void nComboBoxActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_nComboBoxActionPerformed
        // TODO add your handling code here:
    }//GEN-LAST:event_nComboBoxActionPerformed

    private void jButton2ActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_jButton2ActionPerformed
       this.gPanel.setColorPalette(
               B3NucleotidePaletteFactory.createPalette(
               ((ColorItem)aComboBox.getSelectedItem()).color,
               ((ColorItem)cComboBox.getSelectedItem()).color,
               ((ColorItem)gComboBox.getSelectedItem()).color,
               ((ColorItem)tComboBox.getSelectedItem()).color,
               ((ColorItem)nComboBox.getSelectedItem()).color
               )
               );
       gPanel.repaint();
    }//GEN-LAST:event_jButton2ActionPerformed

    
    
    
    private DefaultComboBoxModel<ColorItem> getColorListModel(){
        DefaultComboBoxModel<ColorItem> model = new DefaultComboBoxModel<>();
        Set<String> cset = new TreeSet<String>();
        java.lang.reflect.Field[] fields = java.awt.Color.class.getDeclaredFields();
        for (java.lang.reflect.Field f : fields) {
           try {
              if (! java.lang.reflect.Modifier.isPublic(f.getModifiers())) break;
              if (! java.lang.reflect.Modifier.isStatic(f.getModifiers())) break;
              if (! java.lang.reflect.Modifier.isFinal(f.getModifiers())) break;
              Object o = f.get(null);
              if (o instanceof java.awt.Color) { 
                 String name = f.getName();
                 name = name.substring(0, 1).toUpperCase() +  
                        name.substring(1).toLowerCase(); // convert to title case
                 name = name.replace("_",""); // de-duplicate Light_gray and Lightgray
                 if(!cset.contains(name)){
                    model.addElement(new ColorItem((java.awt.Color) o, name));
                    cset.add(name);
                 }
              }
           } catch (Exception e) {
              e.printStackTrace();
           }
        }
        return model;
    }
    
    
    
    
    // Variables declaration - do not modify//GEN-BEGIN:variables
    private javax.swing.JComboBox aComboBox;
    private javax.swing.JComboBox cComboBox;
    private javax.swing.JPanel centralPanel;
    private javax.swing.JLabel fromField;
    private javax.swing.JComboBox gComboBox;
    private javax.swing.JButton jButton1;
    private javax.swing.JButton jButton2;
    private javax.swing.JLabel jLabel1;
    private javax.swing.JLabel jLabel10;
    private javax.swing.JLabel jLabel11;
    private javax.swing.JLabel jLabel12;
    private javax.swing.JLabel jLabel13;
    private javax.swing.JLabel jLabel14;
    private javax.swing.JLabel jLabel2;
    private javax.swing.JLabel jLabel5;
    private javax.swing.JLabel jLabel6;
    private javax.swing.JLabel jLabel7;
    private javax.swing.JLabel jLabel8;
    private javax.swing.JLabel jLabel9;
    private javax.swing.JPanel jPanel3;
    private javax.swing.JSeparator jSeparator1;
    private javax.swing.JSeparator jSeparator2;
    private javax.swing.JSeparator jSeparator3;
    private javax.swing.JTextField jTextField1;
    private javax.swing.JTextField jTextField2;
    private javax.swing.JComboBox nComboBox;
    private javax.swing.JLabel nofColumnsField;
    private javax.swing.JLabel nofRowsField;
    private javax.swing.JComboBox tComboBox;
    private javax.swing.JLabel toField;
    // End of variables declaration//GEN-END:variables
}
