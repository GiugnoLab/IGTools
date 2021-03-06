/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package igtools.gui2.mults;

import igtools.dictionaries.elsa.IELSAIterator;
import igtools.dictionaries.elsa.NELSA;
import igtools.gui2.WSSequence;
import java.awt.BorderLayout;
import java.awt.Color;
import java.text.DecimalFormat;
import java.util.Map;
import java.util.TreeMap;
import org.jfree.chart.ChartFactory;
import org.jfree.chart.ChartMouseEvent;
import org.jfree.chart.ChartMouseListener;
import org.jfree.chart.ChartPanel;
import org.jfree.chart.JFreeChart;
import org.jfree.chart.axis.LogAxis;
import org.jfree.chart.axis.NumberAxis;
import org.jfree.chart.entity.XYItemEntity;
import org.jfree.chart.plot.Plot;
import org.jfree.chart.plot.PlotOrientation;
import org.jfree.chart.plot.XYPlot;
import org.jfree.chart.renderer.xy.StandardXYBarPainter;
import org.jfree.chart.renderer.xy.XYBarRenderer;
import org.jfree.data.xy.XYDataItem;
import org.jfree.data.xy.XYSeries;
import org.jfree.data.xy.XYSeriesCollection;

/**
 *
 * @author vbonnici
 */
public class CoMultFrame extends javax.swing.JFrame {
     
    
    
    
    
    private DecimalFormat df = new DecimalFormat("###,###,###,###");
    
    private WSSequence wsseq;
    private NELSA nelsa;
    private ChartPanel chartPanel = null;

    /**
     * Creates new form ProperCoDistancesFrame
     */
    public CoMultFrame(WSSequence wsseq) {
        this.wsseq  = wsseq;
        this.nelsa = wsseq.getNELSA();
        
        initComponents();
        this.setTitle("Co-Multiplicity: "+wsseq.getName());
       
        center_panel.setLayout(new BorderLayout());
        this.pickedPairLabel.setText("");
    
    }
    
    
    public CoMultFrame(WSSequence wsseq, int k) {
        this(wsseq);
        kField.setText(""+k);
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
        kField = new javax.swing.JTextField();
        drawButton = new javax.swing.JButton();
        pickedPairLabel = new javax.swing.JLabel();
        decreaseButton = new javax.swing.JButton();
        increaseButton = new javax.swing.JButton();
        logxButton = new javax.swing.JToggleButton();
        logyButton = new javax.swing.JToggleButton();
        center_panel = new javax.swing.JPanel();

        setDefaultCloseOperation(javax.swing.WindowConstants.EXIT_ON_CLOSE);
        setMinimumSize(new java.awt.Dimension(600, 350));

        kField.setColumns(4);
        kField.setText("6");

        drawButton.setText("draw");
        drawButton.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                drawButtonActionPerformed(evt);
            }
        });

        pickedPairLabel.setText("jLabel1");

        decreaseButton.setText("-");
        decreaseButton.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                decreaseButtonActionPerformed(evt);
            }
        });

        increaseButton.setText("+");
        increaseButton.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                increaseButtonActionPerformed(evt);
            }
        });

        logxButton.setText("log X");
        logxButton.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                logxButtonActionPerformed(evt);
            }
        });

        logyButton.setText("log Y");
        logyButton.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                logyButtonActionPerformed(evt);
            }
        });

        javax.swing.GroupLayout jPanel1Layout = new javax.swing.GroupLayout(jPanel1);
        jPanel1.setLayout(jPanel1Layout);
        jPanel1Layout.setHorizontalGroup(
            jPanel1Layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGroup(jPanel1Layout.createSequentialGroup()
                .addContainerGap()
                .addComponent(kField, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE)
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                .addComponent(decreaseButton)
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                .addComponent(increaseButton)
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                .addComponent(drawButton)
                .addGap(64, 64, 64)
                .addComponent(logxButton)
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                .addComponent(logyButton)
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED, 197, Short.MAX_VALUE)
                .addComponent(pickedPairLabel)
                .addContainerGap())
        );
        jPanel1Layout.setVerticalGroup(
            jPanel1Layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGroup(jPanel1Layout.createSequentialGroup()
                .addContainerGap()
                .addGroup(jPanel1Layout.createParallelGroup(javax.swing.GroupLayout.Alignment.BASELINE)
                    .addComponent(kField, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE)
                    .addComponent(drawButton)
                    .addComponent(pickedPairLabel)
                    .addComponent(decreaseButton)
                    .addComponent(increaseButton)
                    .addComponent(logxButton)
                    .addComponent(logyButton))
                .addContainerGap(javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE))
        );

        getContentPane().add(jPanel1, java.awt.BorderLayout.PAGE_START);

        center_panel.setBackground(new java.awt.Color(254, 254, 254));
        center_panel.setPreferredSize(new java.awt.Dimension(600, 300));

        javax.swing.GroupLayout center_panelLayout = new javax.swing.GroupLayout(center_panel);
        center_panel.setLayout(center_panelLayout);
        center_panelLayout.setHorizontalGroup(
            center_panelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGap(0, 600, Short.MAX_VALUE)
        );
        center_panelLayout.setVerticalGroup(
            center_panelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGap(0, 297, Short.MAX_VALUE)
        );

        getContentPane().add(center_panel, java.awt.BorderLayout.CENTER);

        pack();
    }// </editor-fold>//GEN-END:initComponents

    private void drawButtonActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_drawButtonActionPerformed
        makeChart();
    }//GEN-LAST:event_drawButtonActionPerformed

    private void decreaseButtonActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_decreaseButtonActionPerformed
        try{
            int value = Integer.parseInt(kField.getText());
            if(value > 1)
                value--;
            kField.setText(""+value);
            makeChart();
            
        }catch(Exception e){  
        }
    }//GEN-LAST:event_decreaseButtonActionPerformed

    private void increaseButtonActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_increaseButtonActionPerformed
        try{
            int value = Integer.parseInt(kField.getText());
            value++;
            kField.setText(""+value);
            makeChart();
            
        }catch(Exception e){  
        }
    }//GEN-LAST:event_increaseButtonActionPerformed

    
    private void checkLogXScale(){
        if(this.chartPanel != null){
            if(this.logxButton.isSelected()){
                ((XYPlot)(this.chartPanel.getChart().getPlot())).setDomainAxis(new LogAxis("Multiplicity"));
            }
            else{
               ((XYPlot)(this.chartPanel.getChart().getPlot())).setDomainAxis(new NumberAxis("Multiplicity"));  
            }
        }
    }
    private void logxButtonActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_logxButtonActionPerformed
       checkLogXScale();
    }//GEN-LAST:event_logxButtonActionPerformed

    private void checkLogYScale(){
        if(this.chartPanel != null){
            if(this.logyButton.isSelected()){
                ((XYPlot)(this.chartPanel.getChart().getPlot())).setRangeAxis(new LogAxis("Co-Multiplicity"));
            }
            else{
               ((XYPlot)(this.chartPanel.getChart().getPlot())).setRangeAxis(new NumberAxis("Co-Multiplicity"));  
            }
        }
    }
    private void logyButtonActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_logyButtonActionPerformed
        checkLogYScale();
    }//GEN-LAST:event_logyButtonActionPerformed

   
    // Variables declaration - do not modify//GEN-BEGIN:variables
    private javax.swing.JPanel center_panel;
    private javax.swing.JButton decreaseButton;
    private javax.swing.JButton drawButton;
    private javax.swing.JButton increaseButton;
    private javax.swing.JPanel jPanel1;
    private javax.swing.JTextField kField;
    private javax.swing.JToggleButton logxButton;
    private javax.swing.JToggleButton logyButton;
    private javax.swing.JLabel pickedPairLabel;
    // End of variables declaration//GEN-END:variables

    
    
    
    
    public void makeChart(){
        if(this.nelsa != null){
            pickedPairLabel.setText("");
            
            this.chartPanel = null;
            
            kField.setEditable(false);
            drawButton.setEnabled(false);
            decreaseButton.setEnabled(false);
            increaseButton.setEnabled(false);
            
            center_panel.removeAll();
            center_panel.setBackground(Color.RED);
            center_panel.invalidate();
            center_panel.repaint();
            
           try{
        
                int k = Integer.parseInt(kField.getText());
                
                
                final Map<Integer,Integer> comults = new TreeMap<>();
                
                Integer mult;

                IELSAIterator it = nelsa.begin(k);
                while(it.next()){
                    mult = comults.get(it.multiplicity());
                    if(mult == null){
                        comults.put(it.multiplicity(), 1);
                    }
                    else{
                        comults.put(it.multiplicity(), mult + 1);
                    }
                }





                final XYSeries series = new XYSeries("");
                for(Map.Entry<Integer, Integer> entry : comults.entrySet()){
//                    System.out.println(entry.getKey() +"\t"+ entry.getValue());
                    series.add(entry.getKey(), entry.getValue());
                }


        //        System.out.println("["+it.istart()+","+it.iend()+"]");
                final XYSeriesCollection dataset = new XYSeriesCollection(series);


                final JFreeChart chart = ChartFactory.createXYBarChart(
                                                null, 
                                                "Multiplicity", 
                                                false, 
                                                "Co-Multiplicity", 
                                                dataset, 
                                                PlotOrientation.VERTICAL, 
                                                false, 
                                                true, 
                                                false);
                
//                final JFreeChart chart  = ChartFactory.createXYLineChart(
//                    "",
//                    "k", 
//                    "Dictionary size", 
//                    dataset,
//                    PlotOrientation.VERTICAL,
//                    false,
//                    true,
//                    false
//                );

                final ChartPanel chartPanel = new ChartPanel(chart);
                chartPanel.setPreferredSize(center_panel.getPreferredSize());
                chart.setBackgroundPaint(Color.white);
                //chartPanel.setMouseZoomable(false);
                chartPanel.setMouseWheelEnabled(true);

                chartPanel.setMinimumDrawWidth( 0 );
                chartPanel.setMinimumDrawHeight( 0 );
                chartPanel.setMaximumDrawWidth( 1920 );
                chartPanel.setMaximumDrawHeight( 1200 );

                XYPlot plot = (XYPlot) chart.getPlot();
                plot.setBackgroundPaint(Color.WHITE);
                ((XYBarRenderer) plot.getRenderer()).setBarPainter(new StandardXYBarPainter());
                plot.getRenderer().setSeriesPaint( 0, Color.BLUE);
                
                //((XYBarRenderer) plot.getRenderer()).setMargin(1);
                
                
                plot.setDomainPannable(true);
                plot.setRangePannable(true);
                
                plot.setDomainGridlinesVisible(true);  
                plot.setRangeGridlinesVisible(true);  
                plot.setRangeGridlinePaint(Color.gray);  
                plot.setDomainGridlinePaint(Color.gray);
                


                chartPanel.addChartMouseListener(new ChartMouseListener() {

                    @Override
                    public void chartMouseClicked(ChartMouseEvent cme) {

                        Plot p = cme.getChart().getPlot();
                        if(p instanceof XYPlot){                    
                            if(cme.getEntity() instanceof XYItemEntity){                        
                                int seriesIndex = ((XYItemEntity)cme.getEntity()).getSeriesIndex();
                                int item = ((XYItemEntity)cme.getEntity()).getItem();
                                XYSeries series = ((XYSeriesCollection)dataset).getSeries(seriesIndex);
                                XYDataItem xyItem = series.getDataItem(item);
//                                System.out.println(xyItem);
                                
                                pickedPairLabel.setText("("+df.format(xyItem.getX().intValue())+" ; "+df.format(xyItem.getY().intValue())+")");
                                
                                return;
                            }
                        }
                        
                        pickedPairLabel.setText("");
                    }

                    @Override
                    public void chartMouseMoved(ChartMouseEvent cme) {
                    }

                });
        
                this.chartPanel = chartPanel;
                
                checkLogXScale();
                checkLogYScale();
        
                center_panel.add(chartPanel, BorderLayout.CENTER);
                center_panel.revalidate();
                center_panel.repaint();
        
           }catch(Exception e){
                       
           }
           
           
            kField.setEditable(true);
            drawButton.setEnabled(true);
            decreaseButton.setEnabled(true);
            increaseButton.setEnabled(true);
        }
    }
    
}
