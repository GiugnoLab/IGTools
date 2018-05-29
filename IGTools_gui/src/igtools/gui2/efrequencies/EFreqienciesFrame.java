/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package igtools.gui2.efrequencies;

import igtools.analyses.eprob.ExpectedFrequency;
import igtools.analyses.eprob.ExpectedFrequencyByMM;
import igtools.analyses.eprob.ExpectedFrequencyBySubFactors;
import igtools.common.nucleotide.B3Nucleotide;
import igtools.dictionaries.elsa.CompleteIterator;
import igtools.dictionaries.elsa.IELSAIterator;
import igtools.dictionaries.elsa.NELSA;
import igtools.gui2.WSSequence;
import java.awt.Color;
import java.awt.Dimension;
import javax.swing.BoxLayout;
import org.jfree.chart.ChartFactory;
import org.jfree.chart.ChartMouseEvent;
import org.jfree.chart.ChartMouseListener;
import org.jfree.chart.ChartPanel;
import org.jfree.chart.JFreeChart;
import org.jfree.chart.axis.NumberAxis;
import org.jfree.chart.entity.CategoryItemEntity;
import org.jfree.chart.plot.CategoryPlot;
import org.jfree.chart.plot.PlotOrientation;
import org.jfree.chart.plot.XYPlot;
import org.jfree.chart.renderer.category.BarRenderer;
import org.jfree.chart.renderer.category.StandardBarPainter;
import org.jfree.chart.renderer.xy.StandardXYBarPainter;
import org.jfree.chart.renderer.xy.XYBarRenderer;
import org.jfree.data.category.DefaultCategoryDataset;
import org.jfree.data.xy.XYSeries;
import org.jfree.data.xy.XYSeriesCollection;

/**
 *
 * @author vbonnici
 */
public class EFreqienciesFrame extends javax.swing.JFrame {

    private WSSequence wsseq;
    private NELSA nelsa = null;
    
    
    private enum ETYPE {FL,SUB};
    private ETYPE etype = ETYPE.FL;
    
    
    private final int nofChartRows = 4;
    private final int nofChartCols = 1;
    
    
     private ExpectedFrequency ef;
    private ExpectedFrequency efFL;
    private ExpectedFrequency efSub;
    
    /**
     * Creates new form EFreqienciesFrame
     */
    public EFreqienciesFrame(WSSequence wsseq) {
       this.wsseq = wsseq;
       this.nelsa = wsseq.getNELSA();
       
       this.ef = null;
       this.efFL = new ExpectedFrequencyByMM(this.nelsa);
       this.efSub = new ExpectedFrequencyBySubFactors(this.nelsa);
        
       initComponents();
       this.setTitle("Expected Frequencies: " + wsseq.getName());
       
       this.chartsPanel.setLayout(new BoxLayout(this.chartsPanel, BoxLayout.Y_AXIS));
    }
    
    
    
    
    private void getEType(){
        String ee = this.etypeComboBox.getSelectedItem().toString();
        if(ee.compareTo("Prefix-Suffix") == 0)
            this.etype = ETYPE.FL;
        else 
            this.etype = ETYPE.SUB;
    }
    
    private int getK(){
        int i = 1;
        try{
            i = Integer.parseInt(this.kTextField.getText());
        }catch(Exception e){
            i = 1;
            this.kTextField.setText("1");
        }
        
        if(i<1){
            i = 1;
            this.kTextField.setText("1");
        }
        return i;
    }
    
    private int getOrder(){
        return Integer.parseInt(this.orderComboBox.getSelectedItem().toString());
    }
    
    private void setOrders(){
        int k = getK();
        getEType();
        
        this.orderComboBox.removeAllItems();
        
        if(k == 1)
            this.orderComboBox.addItem("0");
        else{
            for(int i=0;i<k-1; i++)
                this.orderComboBox.addItem(""+i);
        }
    }
    
    
    private void getEF(){
        getEType();
        if(this.etype == ETYPE.FL)
            this.ef = this.efFL;
        else
            this.ef = this.efSub;
    }

    /**
     * This method is called from within the constructor to initialize the form.
     * WARNING: Do NOT modify this code. The content of this method is always
     * regenerated by the Form Editor.
     */
    @SuppressWarnings("unchecked")
    // <editor-fold defaultstate="collapsed" desc="Generated Code">//GEN-BEGIN:initComponents
    private void initComponents() {

        jPanel1 = new javax.swing.JPanel();
        etypeComboBox = new javax.swing.JComboBox();
        kPrevButton = new javax.swing.JButton();
        kTextField = new javax.swing.JTextField();
        kNextjButton = new javax.swing.JButton();
        multiplicityCheckBox = new javax.swing.JCheckBox();
        gapCheckBox = new javax.swing.JCheckBox();
        realCheckBox = new javax.swing.JCheckBox();
        expectedCheckBox = new javax.swing.JCheckBox();
        drawButton = new javax.swing.JButton();
        jLabel1 = new javax.swing.JLabel();
        jLabel2 = new javax.swing.JLabel();
        oderPrevButton = new javax.swing.JButton();
        oderNextButton = new javax.swing.JButton();
        jSeparator1 = new javax.swing.JSeparator();
        jSeparator2 = new javax.swing.JSeparator();
        orderComboBox = new javax.swing.JComboBox();
        completeIteratorCheckBox = new javax.swing.JCheckBox();
        chartsPanel = new javax.swing.JPanel();

        setDefaultCloseOperation(javax.swing.WindowConstants.EXIT_ON_CLOSE);

        etypeComboBox.setModel(new javax.swing.DefaultComboBoxModel(new String[] { "Prefix-Suffix", "Sub-factors" }));

        kPrevButton.setText("-");
        kPrevButton.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                kPrevButtonActionPerformed(evt);
            }
        });

        kTextField.setText("1");

        kNextjButton.setText("+");
        kNextjButton.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                kNextjButtonActionPerformed(evt);
            }
        });

        multiplicityCheckBox.setText("multiplicity");

        gapCheckBox.setSelected(true);
        gapCheckBox.setText("gap");

        realCheckBox.setText("real");

        expectedCheckBox.setSelected(true);
        expectedCheckBox.setText("expected");

        drawButton.setText("Draw");
        drawButton.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                drawButtonActionPerformed(evt);
            }
        });

        jLabel1.setText("k");

        jLabel2.setText("order");

        oderPrevButton.setText("-");
        oderPrevButton.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                oderPrevButtonActionPerformed(evt);
            }
        });

        oderNextButton.setText("+");
        oderNextButton.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                oderNextButtonActionPerformed(evt);
            }
        });

        jSeparator1.setOrientation(javax.swing.SwingConstants.VERTICAL);

        jSeparator2.setOrientation(javax.swing.SwingConstants.VERTICAL);

        orderComboBox.setModel(new javax.swing.DefaultComboBoxModel(new String[] { "0" }));

        completeIteratorCheckBox.setText("4^k");

        javax.swing.GroupLayout jPanel1Layout = new javax.swing.GroupLayout(jPanel1);
        jPanel1.setLayout(jPanel1Layout);
        jPanel1Layout.setHorizontalGroup(
            jPanel1Layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGroup(jPanel1Layout.createSequentialGroup()
                .addContainerGap()
                .addComponent(drawButton)
                .addGap(18, 18, 18)
                .addComponent(etypeComboBox, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE)
                .addGap(18, 18, 18)
                .addComponent(jLabel1)
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                .addComponent(kPrevButton)
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                .addComponent(kTextField, javax.swing.GroupLayout.PREFERRED_SIZE, 30, javax.swing.GroupLayout.PREFERRED_SIZE)
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                .addComponent(kNextjButton)
                .addGap(18, 18, 18)
                .addComponent(jLabel2)
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                .addComponent(oderPrevButton)
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                .addComponent(orderComboBox, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE)
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                .addComponent(oderNextButton)
                .addGap(18, 18, 18)
                .addComponent(jSeparator1, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE)
                .addGap(18, 18, 18)
                .addComponent(completeIteratorCheckBox)
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                .addComponent(multiplicityCheckBox)
                .addGap(18, 18, 18)
                .addComponent(jSeparator2, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE)
                .addGap(18, 18, 18)
                .addComponent(expectedCheckBox)
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                .addComponent(realCheckBox)
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                .addComponent(gapCheckBox)
                .addContainerGap(70, Short.MAX_VALUE))
        );
        jPanel1Layout.setVerticalGroup(
            jPanel1Layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGroup(jPanel1Layout.createSequentialGroup()
                .addContainerGap()
                .addGroup(jPanel1Layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
                    .addComponent(jSeparator2)
                    .addGroup(jPanel1Layout.createSequentialGroup()
                        .addGroup(jPanel1Layout.createParallelGroup(javax.swing.GroupLayout.Alignment.BASELINE)
                            .addComponent(etypeComboBox, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE)
                            .addComponent(kPrevButton)
                            .addComponent(kTextField, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE)
                            .addComponent(kNextjButton)
                            .addComponent(multiplicityCheckBox)
                            .addComponent(gapCheckBox)
                            .addComponent(realCheckBox)
                            .addComponent(expectedCheckBox)
                            .addComponent(drawButton)
                            .addComponent(jLabel1)
                            .addComponent(jLabel2)
                            .addComponent(oderPrevButton)
                            .addComponent(oderNextButton)
                            .addComponent(orderComboBox, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE)
                            .addComponent(completeIteratorCheckBox))
                        .addGap(0, 0, Short.MAX_VALUE))
                    .addComponent(jSeparator1, javax.swing.GroupLayout.Alignment.TRAILING))
                .addContainerGap())
        );

        getContentPane().add(jPanel1, java.awt.BorderLayout.PAGE_START);

        chartsPanel.setBackground(java.awt.Color.white);

        javax.swing.GroupLayout chartsPanelLayout = new javax.swing.GroupLayout(chartsPanel);
        chartsPanel.setLayout(chartsPanelLayout);
        chartsPanelLayout.setHorizontalGroup(
            chartsPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGap(0, 1007, Short.MAX_VALUE)
        );
        chartsPanelLayout.setVerticalGroup(
            chartsPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGap(0, 398, Short.MAX_VALUE)
        );

        getContentPane().add(chartsPanel, java.awt.BorderLayout.CENTER);

        pack();
    }// </editor-fold>//GEN-END:initComponents

    private void drawButtonActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_drawButtonActionPerformed
        drawCharts();
    }//GEN-LAST:event_drawButtonActionPerformed

    private void kPrevButtonActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_kPrevButtonActionPerformed
        int k = getK();
        if(k>1)
        this.kTextField.setText(""+(k-1));
        setOrders();
        drawCharts();
    }//GEN-LAST:event_kPrevButtonActionPerformed

    private void kNextjButtonActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_kNextjButtonActionPerformed
        int k = getK();
        this.kTextField.setText(""+(k+1));
        setOrders();
        drawCharts();
    }//GEN-LAST:event_kNextjButtonActionPerformed

    private void oderPrevButtonActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_oderPrevButtonActionPerformed
        int i = this.orderComboBox.getSelectedIndex();
        if(i>0)
            this.orderComboBox.setSelectedIndex(i-1);
        drawCharts();
    }//GEN-LAST:event_oderPrevButtonActionPerformed

    private void oderNextButtonActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_oderNextButtonActionPerformed
        int i = this.orderComboBox.getSelectedIndex();
        if(i < this.orderComboBox.getItemCount() - 1)
            this.orderComboBox.setSelectedIndex(i+1);
        drawCharts();
    }//GEN-LAST:event_oderNextButtonActionPerformed

    
    
    private void drawCharts(){
        getEType();
        getEF();
        int k = getK();
        
        IELSAIterator wit = null;
        if(this.completeIteratorCheckBox.isSelected()){
            wit = new CompleteIterator(k);
        }
        else{
            wit = nelsa.begin(k);
        }
        
        int order = this.orderComboBox.getSelectedIndex();
        
        
        clearChartsPanel();
        
        
        if(this.expectedCheckBox.isSelected()){
            System.out.println("echart");
            final ChartPanel chart = makeEChart(wit.clone(), order, this.multiplicityCheckBox.isSelected());
            this.chartsPanel.add(chart);
        }
        
        if(this.realCheckBox.isSelected()){
            final ChartPanel chart = makeRChart(wit.clone(), order, this.multiplicityCheckBox.isSelected());
            this.chartsPanel.add(chart);
        }
        
        if(this.gapCheckBox.isSelected()){
            final ChartPanel chart = makeGChart(wit.clone(), order, this.multiplicityCheckBox.isSelected());
            this.chartsPanel.add(chart);
        }
        
        forceChartsPanel();
    }
    
    
    private int mult(IELSAIterator wit){
        IELSAIterator  it = this.nelsa.find(wit.kmer());
        if(it == null)
                return 0;
        return it.multiplicity();
    }
    
    
    private ChartPanel makeEChart(IELSAIterator wit, int order, boolean mults){
        System.out.println("echart "+wit.k()+" "+order);
        DefaultCategoryDataset dataset = new DefaultCategoryDataset();
        
        if(mults){
            int k = wit.k();
            double evalue;
            while(wit.next()){
                evalue = ef.expFreq(wit, order);
                dataset.addValue(ef.toMultiplicity(evalue, k),"series1", B3Nucleotide.toString(wit.kmer()));
            }
        }
        else{
            double evalue;
            while(wit.next()){
                evalue = ef.expFreq(wit, order);
                dataset.addValue(evalue,"series1", B3Nucleotide.toString(wit.kmer()));
            }
        }
        
        String yLabel = "frequency";
        if(mults)
            yLabel = "multiplicity";
        
        ChartPanel pp = makeRawChart(dataset, "kmers",yLabel);
        pp.getChart().setTitle("Expected values");
        return pp;
    }
    private ChartPanel makeRChart(IELSAIterator wit, int order, boolean mults){
        DefaultCategoryDataset dataset = new DefaultCategoryDataset();
        
        if(mults){
            int k = wit.k();
            double evalue;
            while(wit.next()){
                dataset.addValue(mult(wit),"series1", B3Nucleotide.toString(wit.kmer()));
            }
        }
        else{
            double evalue;
            while(wit.next()){
                dataset.addValue(ef.toFrequency(mult(wit), wit.k()),"series1", B3Nucleotide.toString(wit.kmer()));
            }
        }  
        
        String yLabel = "frequency";
        if(mults)
            yLabel = "multiplicity";
        
        ChartPanel pp = makeRawChart(dataset, "kmers",yLabel);
        pp.getChart().setTitle("Real values");
        return pp;
    }
    private ChartPanel makeGChart(IELSAIterator wit, int order, boolean mults){
        DefaultCategoryDataset dataset = new DefaultCategoryDataset();
        
        if(mults){
            int k = wit.k();
            double evalue;
            while(wit.next()){
                evalue = ef.expFreq(wit, order);
                dataset.addValue(mult(wit) - ef.toMultiplicity(evalue, k),"series1", B3Nucleotide.toString(wit.kmer()));
            }
        }
        else{
            double evalue;
            while(wit.next()){
                evalue = ef.expFreq(wit, order);
                dataset.addValue(ef.toFrequency(mult(wit), wit.k()) - evalue,"series1", B3Nucleotide.toString(wit.kmer()));
            }
        }
        
        String yLabel = "frequency";
        if(mults)
            yLabel = "multiplicity";
        
        ChartPanel pp = makeRawChart(dataset, "kmers",yLabel);
        pp.getChart().setTitle("Gap (Real - Expected) values");
        return pp;
}
    
    
    private void clearChartsPanel(){
        this.chartsPanel.removeAll();
        this.chartsPanel.invalidate();
        this.chartsPanel.repaint();
    }
    
    private void forceChartsPanel(){
        this.chartsPanel.revalidate();
        this.chartsPanel.repaint();
    }
    
    
    private ChartPanel makeRawChart(
            final DefaultCategoryDataset dataset,
            final String xLabel, 
            final String yLabel)
    {
        
        System.out.println("chart "+xLabel+" "+yLabel);
        System.out.println((chartsPanel.getWidth() / nofChartCols) +" "+(chartsPanel.getHeight() / nofChartRows));
        
        final JFreeChart chart = ChartFactory.createBarChart(
                        null,         // chart title
                        xLabel,               // domain axis label
                        yLabel,                  // range axis label
                        dataset,                  // data
                        PlotOrientation.VERTICAL, // orientation
                        false,                     // include legend
                        true,                     // tooltips?
                        false                     // URLs?
                    );
        final ChartPanel chartPanel = new ChartPanel(chart){
            @Override
            public Dimension getPreferredSize() {
                return new Dimension(chartsPanel.getWidth() / nofChartCols,  chartsPanel.getHeight() / nofChartRows);
            }
        };
        chart.setBackgroundPaint(Color.white);
//                //chartPanel.setMouseZoomable(false);
        chartPanel.setMouseWheelEnabled(true);

        chartPanel.setMinimumDrawWidth( 0 );
        chartPanel.setMinimumDrawHeight( 0 );
        chartPanel.setMaximumDrawWidth( 1920 );
        chartPanel.setMaximumDrawHeight( 1200 );

        final CategoryPlot plot = chart.getCategoryPlot();
        plot.setBackgroundPaint(Color.white);
        ((BarRenderer)plot.getRenderer()).setBarPainter(new StandardBarPainter());

        // set the range axis to display integers only...
        final NumberAxis rangeAxis = (NumberAxis) plot.getRangeAxis();
        // disable bar outlines...
        final BarRenderer renderer = (BarRenderer) plot.getRenderer();
        plot.setRangePannable(true); 
        plot.setRangeGridlinesVisible(true);  
        plot.setRangeGridlinePaint(Color.gray);  
        return chartPanel;
    }
    
    
    
    
    
    
    // Variables declaration - do not modify//GEN-BEGIN:variables
    private javax.swing.JPanel chartsPanel;
    private javax.swing.JCheckBox completeIteratorCheckBox;
    private javax.swing.JButton drawButton;
    private javax.swing.JComboBox etypeComboBox;
    private javax.swing.JCheckBox expectedCheckBox;
    private javax.swing.JCheckBox gapCheckBox;
    private javax.swing.JLabel jLabel1;
    private javax.swing.JLabel jLabel2;
    private javax.swing.JPanel jPanel1;
    private javax.swing.JSeparator jSeparator1;
    private javax.swing.JSeparator jSeparator2;
    private javax.swing.JButton kNextjButton;
    private javax.swing.JButton kPrevButton;
    private javax.swing.JTextField kTextField;
    private javax.swing.JCheckBox multiplicityCheckBox;
    private javax.swing.JButton oderNextButton;
    private javax.swing.JButton oderPrevButton;
    private javax.swing.JComboBox orderComboBox;
    private javax.swing.JCheckBox realCheckBox;
    // End of variables declaration//GEN-END:variables
}
