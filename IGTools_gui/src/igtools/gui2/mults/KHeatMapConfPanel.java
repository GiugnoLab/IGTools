/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package igtools.gui2.mults;

import igtools.common.charts.KHeatMap;
import java.awt.Font;
import javax.swing.JColorChooser;

/**
 *
 * @author vbonnici
 */
public class KHeatMapConfPanel extends javax.swing.JPanel {

    private KHeatMap heatmap = new KHeatMap();
    private String[] values = null;
    /**
     * Creates new form KHeatMapConfPanel
     */
    public KHeatMapConfPanel(KHeatMap heatmap) {
        initComponents();
        this.heatmap = heatmap;
        putValues();
    }
    
    private String toStringValue(int i){
        return (new Integer(i)).toString();
    }
    private String toStringValue(double i){
        return (new Double(i)).toString();
    }
    private Integer toInteger(String s){
        return Integer.parseInt(s);
    }
    private Double toDouble(String s){
        return Double.parseDouble(s);
    }
    private void putValues(){
        if(heatmap != null){
          this.marginsTextField.setText((new Integer(heatmap.margin)).toString());
          this.backgroundColorButton.setBackground(heatmap.backgroundColor);
          
          this.fontColorButton.setBackground(heatmap.fontColor);
          this.fontSizeTextField.setText(toStringValue(heatmap.font.getSize()));
          this.topRulerCheckbox.setSelected(heatmap.printTopRuler);
          this.bottomRulerCheckbox.setSelected(heatmap.printBottomRuler);
          
          this.printLegendCheckbox.setSelected(heatmap.printLegend);
          this.printMultiLegendCheckbox.setSelected(heatmap.printMultiLegend);
          this.legendHeightTextField.setText(toStringValue(heatmap.legendHeight));
          
          this.lineWidthTextField.setText(toStringValue(heatmap.lineWidth));
          this.lineHeightTextField.setText(toStringValue(heatmap.lineHeight));
          this.lineVspaceTextField.setText(toStringValue(heatmap.lineSpace));
          
          this.normalizeByRowCheckbox.setSelected(heatmap.normalizeByRow);
          this.fromZeroValuesCheckbox.setSelected(heatmap.fromZeroGradient);
          this.printLabelsCheckbox.setSelected(heatmap.printLabels);
          this.printValueCheckbox.setSelected(heatmap.printValues);
          
          this.lowColorButton.setBackground(heatmap.lowValueColor);
          this.meanColorCheckbox.setSelected(heatmap.medianValueColor != null);
          if(heatmap.medianValueColor != null)
            this.meanColorButton.setBackground(heatmap.medianValueColor);
          this.highColorButton.setBackground(heatmap.highValueColor);
          this.colorFactorTextField.setText(toStringValue(heatmap.colourScale));
          this.colorPositionTextField.setText(toStringValue(heatmap.hlColorLimit));
        }
    }
    public void getValues(){
        heatmap.margin = toInteger(this.marginsTextField.getText());
        heatmap.backgroundColor = this.backgroundColorButton.getBackground();
        
        heatmap.fontColor = this.fontColorButton.getBackground();
        heatmap.font = new Font(heatmap.font.getFamily(), Font.PLAIN, toInteger(this.fontSizeTextField.getText()));
        heatmap.printTopRuler = this.topRulerCheckbox.isSelected();
        heatmap.printBottomRuler = this.bottomRulerCheckbox.isSelected();
        
        heatmap.printLegend = this.printLegendCheckbox.isSelected();
        heatmap.printMultiLegend = this.printMultiLegendCheckbox.isSelected();
        heatmap.legendHeight = toInteger(this.legendHeightTextField.getText());
        
        heatmap.lineWidth = toInteger(this.lineWidthTextField.getText());
        heatmap.lineHeight = toInteger(this.lineHeightTextField.getText());
        heatmap.lineSpace =  toInteger(this.lineVspaceTextField.getText());
        
        heatmap.normalizeByRow = this.normalizeByRowCheckbox.isSelected();
        heatmap.fromZeroGradient = this.fromZeroValuesCheckbox.isSelected();
        this.heatmap.printLabels = this.printLabelsCheckbox.isSelected();
        heatmap.printValues = this.printValueCheckbox.isSelected();
        heatmap.lowValueColor = this.lowColorButton.getBackground();
        if(this.meanColorCheckbox.isSelected())
            heatmap.medianValueColor = this.meanColorButton.getBackground();
        else
            heatmap.medianValueColor = null;
        heatmap.highValueColor = this.highColorButton.getBackground();
        heatmap.colourScale = toDouble(this.colorFactorTextField.getText());
        heatmap.hlColorLimit = toDouble(this.colorPositionTextField.getText());
        
        
        
        
    }

    /**
     * This method is called from within the constructor to initialize the form.
     * WARNING: Do NOT modify this code. The content of this method is always
     * regenerated by the Form Editor.
     */
    @SuppressWarnings("unchecked")
    // <editor-fold defaultstate="collapsed" desc="Generated Code">//GEN-BEGIN:initComponents
    private void initComponents() {

        jTabbedPane1 = new javax.swing.JTabbedPane();
        jPanel1 = new javax.swing.JPanel();
        jLabel2 = new javax.swing.JLabel();
        marginsTextField = new javax.swing.JTextField();
        backgroundColorButton = new javax.swing.JButton();
        jLabel7 = new javax.swing.JLabel();
        fontColorButton = new javax.swing.JButton();
        jLabel11 = new javax.swing.JLabel();
        jSeparator4 = new javax.swing.JSeparator();
        topRulerCheckbox = new javax.swing.JCheckBox();
        bottomRulerCheckbox = new javax.swing.JCheckBox();
        fontSizeTextField = new javax.swing.JTextField();
        jSeparator5 = new javax.swing.JSeparator();
        jPanel4 = new javax.swing.JPanel();
        printLegendCheckbox = new javax.swing.JCheckBox();
        printMultiLegendCheckbox = new javax.swing.JCheckBox();
        legendHeightTextField = new javax.swing.JTextField();
        jLabel6 = new javax.swing.JLabel();
        jPanel2 = new javax.swing.JPanel();
        jLabel3 = new javax.swing.JLabel();
        jLabel4 = new javax.swing.JLabel();
        jLabel5 = new javax.swing.JLabel();
        lineWidthTextField = new javax.swing.JTextField();
        lineHeightTextField = new javax.swing.JTextField();
        lineVspaceTextField = new javax.swing.JTextField();
        jSeparator1 = new javax.swing.JSeparator();
        fromZeroValuesCheckbox = new javax.swing.JCheckBox();
        normalizeByRowCheckbox = new javax.swing.JCheckBox();
        jSeparator2 = new javax.swing.JSeparator();
        printLabelsCheckbox = new javax.swing.JCheckBox();
        printValueCheckbox = new javax.swing.JCheckBox();
        jSeparator3 = new javax.swing.JSeparator();
        lowColorButton = new javax.swing.JButton();
        jLabel8 = new javax.swing.JLabel();
        jLabel9 = new javax.swing.JLabel();
        meanColorButton = new javax.swing.JButton();
        jLabel10 = new javax.swing.JLabel();
        highColorButton = new javax.swing.JButton();
        meanColorCheckbox = new javax.swing.JCheckBox();
        colorFactorTextField = new javax.swing.JTextField();
        jLabel12 = new javax.swing.JLabel();
        colorPositionTextField = new javax.swing.JTextField();
        jLabel13 = new javax.swing.JLabel();

        setMaximumSize(new java.awt.Dimension(32767, 116));

        jTabbedPane1.setBorder(new javax.swing.border.SoftBevelBorder(javax.swing.border.BevelBorder.RAISED));
        jTabbedPane1.setMaximumSize(new java.awt.Dimension(32767, 116));

        jLabel2.setText("margins");

        marginsTextField.setText("jTextField1");

        backgroundColorButton.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                backgroundColorButtonActionPerformed(evt);
            }
        });

        jLabel7.setText("background");

        fontColorButton.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                fontColorButtonActionPerformed(evt);
            }
        });

        jLabel11.setText("font");

        jSeparator4.setOrientation(javax.swing.SwingConstants.VERTICAL);

        topRulerCheckbox.setText("top ruler");

        bottomRulerCheckbox.setText("bottom ruler");

        fontSizeTextField.setText("jTextField6");

        jSeparator5.setOrientation(javax.swing.SwingConstants.VERTICAL);

        javax.swing.GroupLayout jPanel1Layout = new javax.swing.GroupLayout(jPanel1);
        jPanel1.setLayout(jPanel1Layout);
        jPanel1Layout.setHorizontalGroup(
            jPanel1Layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGroup(jPanel1Layout.createSequentialGroup()
                .addContainerGap()
                .addGroup(jPanel1Layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
                    .addGroup(jPanel1Layout.createSequentialGroup()
                        .addComponent(marginsTextField, javax.swing.GroupLayout.PREFERRED_SIZE, 50, javax.swing.GroupLayout.PREFERRED_SIZE)
                        .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                        .addComponent(jLabel2))
                    .addComponent(jLabel7)
                    .addGroup(jPanel1Layout.createSequentialGroup()
                        .addGap(12, 12, 12)
                        .addComponent(backgroundColorButton, javax.swing.GroupLayout.PREFERRED_SIZE, 50, javax.swing.GroupLayout.PREFERRED_SIZE)))
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                .addComponent(jSeparator4, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE)
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                .addGroup(jPanel1Layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
                    .addComponent(jLabel11)
                    .addGroup(jPanel1Layout.createSequentialGroup()
                        .addComponent(fontColorButton, javax.swing.GroupLayout.PREFERRED_SIZE, 50, javax.swing.GroupLayout.PREFERRED_SIZE)
                        .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                        .addComponent(fontSizeTextField, javax.swing.GroupLayout.PREFERRED_SIZE, 50, javax.swing.GroupLayout.PREFERRED_SIZE)))
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                .addComponent(jSeparator5, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE)
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                .addGroup(jPanel1Layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
                    .addComponent(topRulerCheckbox)
                    .addComponent(bottomRulerCheckbox))
                .addContainerGap(385, Short.MAX_VALUE))
        );
        jPanel1Layout.setVerticalGroup(
            jPanel1Layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGroup(javax.swing.GroupLayout.Alignment.TRAILING, jPanel1Layout.createSequentialGroup()
                .addContainerGap()
                .addGroup(jPanel1Layout.createParallelGroup(javax.swing.GroupLayout.Alignment.TRAILING)
                    .addComponent(jSeparator5, javax.swing.GroupLayout.DEFAULT_SIZE, 145, Short.MAX_VALUE)
                    .addComponent(jSeparator4, javax.swing.GroupLayout.Alignment.LEADING)
                    .addGroup(javax.swing.GroupLayout.Alignment.LEADING, jPanel1Layout.createSequentialGroup()
                        .addGroup(jPanel1Layout.createParallelGroup(javax.swing.GroupLayout.Alignment.TRAILING)
                            .addGroup(javax.swing.GroupLayout.Alignment.LEADING, jPanel1Layout.createSequentialGroup()
                                .addGroup(jPanel1Layout.createParallelGroup(javax.swing.GroupLayout.Alignment.BASELINE)
                                    .addComponent(jLabel2)
                                    .addComponent(marginsTextField, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE))
                                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.UNRELATED)
                                .addComponent(jLabel7)
                                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                                .addComponent(backgroundColorButton, javax.swing.GroupLayout.PREFERRED_SIZE, 30, javax.swing.GroupLayout.PREFERRED_SIZE))
                            .addGroup(javax.swing.GroupLayout.Alignment.LEADING, jPanel1Layout.createSequentialGroup()
                                .addComponent(jLabel11)
                                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                                .addGroup(jPanel1Layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING, false)
                                    .addComponent(fontSizeTextField)
                                    .addComponent(fontColorButton, javax.swing.GroupLayout.DEFAULT_SIZE, 30, Short.MAX_VALUE)))
                            .addGroup(javax.swing.GroupLayout.Alignment.LEADING, jPanel1Layout.createSequentialGroup()
                                .addComponent(topRulerCheckbox)
                                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                                .addComponent(bottomRulerCheckbox)))
                        .addGap(0, 0, Short.MAX_VALUE)))
                .addContainerGap())
        );

        jTabbedPane1.addTab("General", jPanel1);

        printLegendCheckbox.setText("print legend");

        printMultiLegendCheckbox.setText("print multi-legend");

        legendHeightTextField.setText("jTextField5");

        jLabel6.setText("height");

        javax.swing.GroupLayout jPanel4Layout = new javax.swing.GroupLayout(jPanel4);
        jPanel4.setLayout(jPanel4Layout);
        jPanel4Layout.setHorizontalGroup(
            jPanel4Layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGroup(jPanel4Layout.createSequentialGroup()
                .addContainerGap()
                .addGroup(jPanel4Layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
                    .addComponent(printLegendCheckbox)
                    .addGroup(jPanel4Layout.createSequentialGroup()
                        .addComponent(legendHeightTextField, javax.swing.GroupLayout.PREFERRED_SIZE, 50, javax.swing.GroupLayout.PREFERRED_SIZE)
                        .addGap(2, 2, 2)
                        .addComponent(jLabel6)))
                .addGap(18, 18, 18)
                .addComponent(printMultiLegendCheckbox)
                .addContainerGap(490, Short.MAX_VALUE))
        );
        jPanel4Layout.setVerticalGroup(
            jPanel4Layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGroup(jPanel4Layout.createSequentialGroup()
                .addContainerGap()
                .addGroup(jPanel4Layout.createParallelGroup(javax.swing.GroupLayout.Alignment.BASELINE)
                    .addComponent(printLegendCheckbox)
                    .addComponent(printMultiLegendCheckbox))
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                .addGroup(jPanel4Layout.createParallelGroup(javax.swing.GroupLayout.Alignment.BASELINE)
                    .addComponent(legendHeightTextField, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE)
                    .addComponent(jLabel6))
                .addContainerGap(100, Short.MAX_VALUE))
        );

        jTabbedPane1.addTab("Legend", jPanel4);

        jLabel3.setText("width");

        jLabel4.setText("height");

        jLabel5.setText("vspace");

        lineWidthTextField.setText("jTextField2");

        lineHeightTextField.setText("jTextField2");

        lineVspaceTextField.setText("jTextField2");

        jSeparator1.setOrientation(javax.swing.SwingConstants.VERTICAL);

        fromZeroValuesCheckbox.setText("from zero value");

        normalizeByRowCheckbox.setText("normalize by row");

        printLabelsCheckbox.setText("print labels");

        printValueCheckbox.setText("print values");

        lowColorButton.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                lowColorButtonActionPerformed(evt);
            }
        });

        jLabel8.setText("low color");

        jLabel9.setText("mean color");

        meanColorButton.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                meanColorButtonActionPerformed(evt);
            }
        });

        jLabel10.setText("high color");

        highColorButton.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                highColorButtonActionPerformed(evt);
            }
        });

        meanColorCheckbox.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                meanColorCheckboxActionPerformed(evt);
            }
        });

        colorFactorTextField.setText("jTextField2");

        jLabel12.setText("scale");

        colorPositionTextField.setText("jTextField2");

        jLabel13.setText("position");

        javax.swing.GroupLayout jPanel2Layout = new javax.swing.GroupLayout(jPanel2);
        jPanel2.setLayout(jPanel2Layout);
        jPanel2Layout.setHorizontalGroup(
            jPanel2Layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGroup(jPanel2Layout.createSequentialGroup()
                .addContainerGap()
                .addGroup(jPanel2Layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
                    .addGroup(jPanel2Layout.createParallelGroup(javax.swing.GroupLayout.Alignment.TRAILING, false)
                        .addGroup(javax.swing.GroupLayout.Alignment.LEADING, jPanel2Layout.createSequentialGroup()
                            .addComponent(lineHeightTextField, javax.swing.GroupLayout.PREFERRED_SIZE, 50, javax.swing.GroupLayout.PREFERRED_SIZE)
                            .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED, javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE)
                            .addComponent(jLabel4))
                        .addGroup(javax.swing.GroupLayout.Alignment.LEADING, jPanel2Layout.createSequentialGroup()
                            .addComponent(lineWidthTextField, javax.swing.GroupLayout.PREFERRED_SIZE, 50, javax.swing.GroupLayout.PREFERRED_SIZE)
                            .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                            .addComponent(jLabel3)))
                    .addGroup(jPanel2Layout.createSequentialGroup()
                        .addComponent(lineVspaceTextField, javax.swing.GroupLayout.PREFERRED_SIZE, 50, javax.swing.GroupLayout.PREFERRED_SIZE)
                        .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                        .addComponent(jLabel5)))
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                .addComponent(jSeparator1, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE)
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                .addGroup(jPanel2Layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
                    .addComponent(jSeparator2)
                    .addComponent(jSeparator3, javax.swing.GroupLayout.DEFAULT_SIZE, 619, Short.MAX_VALUE)
                    .addGroup(jPanel2Layout.createSequentialGroup()
                        .addGroup(jPanel2Layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
                            .addGroup(jPanel2Layout.createSequentialGroup()
                                .addComponent(normalizeByRowCheckbox)
                                .addGap(18, 18, 18)
                                .addComponent(fromZeroValuesCheckbox))
                            .addGroup(jPanel2Layout.createSequentialGroup()
                                .addComponent(printLabelsCheckbox)
                                .addGap(58, 58, 58)
                                .addComponent(printValueCheckbox))
                            .addGroup(jPanel2Layout.createSequentialGroup()
                                .addGroup(jPanel2Layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
                                    .addComponent(jLabel8)
                                    .addComponent(lowColorButton, javax.swing.GroupLayout.PREFERRED_SIZE, 50, javax.swing.GroupLayout.PREFERRED_SIZE))
                                .addGap(18, 18, 18)
                                .addGroup(jPanel2Layout.createParallelGroup(javax.swing.GroupLayout.Alignment.TRAILING, false)
                                    .addComponent(jLabel9)
                                    .addGroup(jPanel2Layout.createSequentialGroup()
                                        .addComponent(meanColorCheckbox)
                                        .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED, javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE)
                                        .addComponent(meanColorButton, javax.swing.GroupLayout.PREFERRED_SIZE, 50, javax.swing.GroupLayout.PREFERRED_SIZE)))
                                .addGap(18, 18, 18)
                                .addGroup(jPanel2Layout.createParallelGroup(javax.swing.GroupLayout.Alignment.TRAILING)
                                    .addComponent(jLabel10)
                                    .addComponent(highColorButton, javax.swing.GroupLayout.PREFERRED_SIZE, 50, javax.swing.GroupLayout.PREFERRED_SIZE))
                                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                                .addComponent(colorFactorTextField, javax.swing.GroupLayout.PREFERRED_SIZE, 50, javax.swing.GroupLayout.PREFERRED_SIZE)
                                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                                .addComponent(jLabel12)
                                .addGap(18, 18, 18)
                                .addComponent(colorPositionTextField, javax.swing.GroupLayout.PREFERRED_SIZE, 50, javax.swing.GroupLayout.PREFERRED_SIZE)
                                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                                .addComponent(jLabel13)))
                        .addGap(0, 0, Short.MAX_VALUE)))
                .addContainerGap())
        );
        jPanel2Layout.setVerticalGroup(
            jPanel2Layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGroup(jPanel2Layout.createSequentialGroup()
                .addContainerGap()
                .addGroup(jPanel2Layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
                    .addGroup(jPanel2Layout.createSequentialGroup()
                        .addGroup(jPanel2Layout.createParallelGroup(javax.swing.GroupLayout.Alignment.BASELINE)
                            .addComponent(normalizeByRowCheckbox)
                            .addComponent(fromZeroValuesCheckbox))
                        .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                        .addComponent(jSeparator2, javax.swing.GroupLayout.PREFERRED_SIZE, 10, javax.swing.GroupLayout.PREFERRED_SIZE)
                        .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                        .addGroup(jPanel2Layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
                            .addComponent(printValueCheckbox)
                            .addComponent(printLabelsCheckbox, javax.swing.GroupLayout.Alignment.TRAILING))
                        .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                        .addComponent(jSeparator3, javax.swing.GroupLayout.PREFERRED_SIZE, 10, javax.swing.GroupLayout.PREFERRED_SIZE)
                        .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                        .addGroup(jPanel2Layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
                            .addGroup(jPanel2Layout.createSequentialGroup()
                                .addComponent(jLabel8)
                                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                                .addComponent(lowColorButton, javax.swing.GroupLayout.PREFERRED_SIZE, 30, javax.swing.GroupLayout.PREFERRED_SIZE))
                            .addGroup(jPanel2Layout.createSequentialGroup()
                                .addComponent(jLabel9)
                                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                                .addGroup(jPanel2Layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
                                    .addComponent(meanColorButton, javax.swing.GroupLayout.PREFERRED_SIZE, 30, javax.swing.GroupLayout.PREFERRED_SIZE)
                                    .addComponent(meanColorCheckbox)))
                            .addGroup(jPanel2Layout.createSequentialGroup()
                                .addComponent(jLabel10)
                                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                                .addGroup(jPanel2Layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
                                    .addComponent(highColorButton, javax.swing.GroupLayout.PREFERRED_SIZE, 30, javax.swing.GroupLayout.PREFERRED_SIZE)
                                    .addGroup(jPanel2Layout.createParallelGroup(javax.swing.GroupLayout.Alignment.BASELINE)
                                        .addComponent(colorFactorTextField, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE)
                                        .addComponent(jLabel12)
                                        .addGroup(jPanel2Layout.createParallelGroup(javax.swing.GroupLayout.Alignment.BASELINE)
                                            .addComponent(colorPositionTextField, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE)
                                            .addComponent(jLabel13))))))
                        .addGap(0, 12, Short.MAX_VALUE))
                    .addGroup(jPanel2Layout.createSequentialGroup()
                        .addGroup(jPanel2Layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
                            .addComponent(jSeparator1)
                            .addGroup(jPanel2Layout.createSequentialGroup()
                                .addGroup(jPanel2Layout.createParallelGroup(javax.swing.GroupLayout.Alignment.BASELINE)
                                    .addComponent(lineWidthTextField, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE)
                                    .addComponent(jLabel3))
                                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                                .addGroup(jPanel2Layout.createParallelGroup(javax.swing.GroupLayout.Alignment.BASELINE)
                                    .addComponent(lineHeightTextField, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE)
                                    .addComponent(jLabel4))
                                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                                .addGroup(jPanel2Layout.createParallelGroup(javax.swing.GroupLayout.Alignment.BASELINE)
                                    .addComponent(lineVspaceTextField, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE)
                                    .addComponent(jLabel5))
                                .addGap(0, 0, Short.MAX_VALUE)))
                        .addContainerGap())))
        );

        jTabbedPane1.addTab("Line", jPanel2);

        javax.swing.GroupLayout layout = new javax.swing.GroupLayout(this);
        this.setLayout(layout);
        layout.setHorizontalGroup(
            layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addComponent(jTabbedPane1, javax.swing.GroupLayout.DEFAULT_SIZE, 793, Short.MAX_VALUE)
        );
        layout.setVerticalGroup(
            layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGroup(layout.createSequentialGroup()
                .addComponent(jTabbedPane1, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE)
                .addGap(0, 0, 0))
        );
    }// </editor-fold>//GEN-END:initComponents

    private void meanColorCheckboxActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_meanColorCheckboxActionPerformed
        // TODO add your handling code here:
    }//GEN-LAST:event_meanColorCheckboxActionPerformed

    private void highColorButtonActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_highColorButtonActionPerformed
        highColorButton.setBackground(
            JColorChooser.showDialog(
                this,
                "Choose High Color",
                highColorButton.getBackground()));
    }//GEN-LAST:event_highColorButtonActionPerformed

    private void meanColorButtonActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_meanColorButtonActionPerformed
        meanColorButton.setBackground(
            JColorChooser.showDialog(
                this,
                "Choose Mean Color",
                meanColorButton.getBackground()));
    }//GEN-LAST:event_meanColorButtonActionPerformed

    private void lowColorButtonActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_lowColorButtonActionPerformed
        lowColorButton.setBackground(
            JColorChooser.showDialog(
                this,
                "Choose Low Color",
                lowColorButton.getBackground()));
    }//GEN-LAST:event_lowColorButtonActionPerformed

    private void fontColorButtonActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_fontColorButtonActionPerformed
        fontColorButton.setBackground(
            JColorChooser.showDialog(
                this,
                "Choose Font Color",
                fontColorButton.getBackground()));
    }//GEN-LAST:event_fontColorButtonActionPerformed

    private void backgroundColorButtonActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_backgroundColorButtonActionPerformed
        backgroundColorButton.setBackground(
            JColorChooser.showDialog(
                this,
                "Choose Background Color",
                backgroundColorButton.getBackground()));
    }//GEN-LAST:event_backgroundColorButtonActionPerformed

    // Variables declaration - do not modify//GEN-BEGIN:variables
    private javax.swing.JButton backgroundColorButton;
    private javax.swing.JCheckBox bottomRulerCheckbox;
    private javax.swing.JTextField colorFactorTextField;
    private javax.swing.JTextField colorPositionTextField;
    private javax.swing.JButton fontColorButton;
    private javax.swing.JTextField fontSizeTextField;
    private javax.swing.JCheckBox fromZeroValuesCheckbox;
    private javax.swing.JButton highColorButton;
    private javax.swing.JLabel jLabel10;
    private javax.swing.JLabel jLabel11;
    private javax.swing.JLabel jLabel12;
    private javax.swing.JLabel jLabel13;
    private javax.swing.JLabel jLabel2;
    private javax.swing.JLabel jLabel3;
    private javax.swing.JLabel jLabel4;
    private javax.swing.JLabel jLabel5;
    private javax.swing.JLabel jLabel6;
    private javax.swing.JLabel jLabel7;
    private javax.swing.JLabel jLabel8;
    private javax.swing.JLabel jLabel9;
    private javax.swing.JPanel jPanel1;
    private javax.swing.JPanel jPanel2;
    private javax.swing.JPanel jPanel4;
    private javax.swing.JSeparator jSeparator1;
    private javax.swing.JSeparator jSeparator2;
    private javax.swing.JSeparator jSeparator3;
    private javax.swing.JSeparator jSeparator4;
    private javax.swing.JSeparator jSeparator5;
    private javax.swing.JTabbedPane jTabbedPane1;
    private javax.swing.JTextField legendHeightTextField;
    private javax.swing.JTextField lineHeightTextField;
    private javax.swing.JTextField lineVspaceTextField;
    private javax.swing.JTextField lineWidthTextField;
    private javax.swing.JButton lowColorButton;
    private javax.swing.JTextField marginsTextField;
    private javax.swing.JButton meanColorButton;
    private javax.swing.JCheckBox meanColorCheckbox;
    private javax.swing.JCheckBox normalizeByRowCheckbox;
    private javax.swing.JCheckBox printLabelsCheckbox;
    private javax.swing.JCheckBox printLegendCheckbox;
    private javax.swing.JCheckBox printMultiLegendCheckbox;
    private javax.swing.JCheckBox printValueCheckbox;
    private javax.swing.JCheckBox topRulerCheckbox;
    // End of variables declaration//GEN-END:variables
}