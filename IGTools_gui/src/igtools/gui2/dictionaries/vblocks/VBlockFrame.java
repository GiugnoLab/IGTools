/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package igtools.gui2.dictionaries.vblocks;

import igtools.common.alignment.B3NSmithWaterman;
import igtools.common.alignment.StringSmithWaterman;
import igtools.common.nucleotide.B3Nucleotide;
import igtools.common.sequence.B3LLSequence;
import igtools.dictionaries.elsa.IELSAIterator;
import igtools.gui2.codistances.KmerLabel;

import java.awt.FlowLayout;
import java.awt.Font;
import java.awt.FontMetrics;
import java.awt.font.FontRenderContext;
import java.awt.geom.AffineTransform;
import java.util.Set;
import java.util.Stack;
import javax.swing.JLabel;
import javax.swing.JTable;
import opal.OpalMSAAPI;

/**
 *
 * @author vbonnici
 */
public class VBlockFrame extends javax.swing.JFrame {
    
    public enum RecognitionType{
        SC,
        MSA
    }

    private VBlockTableModel tableModel = new VBlockTableModel();
    
    private double scoreTollerance;
    private double rowTollPerc;
    private double rowTollMin;
    private boolean remNs;
    private IELSAIterator pickedIt;
    private RecognitionType rtype;
    
    
    private final Stack<String> alignSequences = new Stack();
    private final Stack<Integer> alignScores = new Stack();
    
    /**
     * Creates new form VBlockFrame
     */
    public VBlockFrame(String title) {
        initComponents();
        this.setTitle(title);
        kmersTable.setDefaultRenderer(String.class, new KmerLabel());
        
//        kmerLabeljScrollPane.setBackground(kmerLabeljScrollPane.getParent().getBackground());
        
        kmersTable.setAutoResizeMode(JTable.AUTO_RESIZE_OFF);
    }
    

    
    public void recognize(
            double scoreTollerance,
            double rowTollPerc,
            double rowTollMin,
            boolean remNs,
            IELSAIterator pickedIt,
            RecognitionType rtype
            ){
        
        this.setEnabled(false);
        
        this.scoreTollerance = scoreTollerance;
        this.rowTollPerc = rowTollPerc;
        this.rowTollMin = rowTollMin;
        this.remNs = remNs;
        this.pickedIt = pickedIt.clone();
        this.rtype = rtype;
        
        
        if(pickedIt!=null && (pickedIt.hasNext() || pickedIt.hasPrev())){
            getKmersShape();
//            System.out.println("as:"+alignSequences.size());
        }
        
        if(this.alignSequences.size() == 0 || this.alignSequences.size()==1){
            this.setVisible(false);
            this.dispose();
        }
        
        else{
            if(rtype == RecognitionType.SC){
                fillTable();
                
                this.pack();
                this.revalidate();
                this.repaint();
            
                this.cs();
            }
            else if(rtype == RecognitionType.MSA){
                this.msa();
            }

//            System.out.println(alignSequences.size());
            this.pack();
            this.revalidate();
            this.repaint();
            this.setEnabled(true);
        }
    }
    
    
    
    private void getKmersShape(){
        B3Nucleotide[] pickedNs = pickedIt.kmer();
        String pickedString = B3Nucleotide.toString(pickedNs);

        B3NSmithWaterman falign = new B3NSmithWaterman(pickedNs, pickedNs);
        falign.process();
        double inScore = falign.score();
        inScore = inScore - (inScore * scoreTollerance);

        IELSAIterator cIt;
        B3Nucleotide[] cNs;
        String cString;


        String s, as1,as2;


        double cScore;
        int missRows = 0;

        int ii = 0;

        if(pickedIt != null && pickedIt.hasPrev()){
           cIt = pickedIt.clone();

           while(cIt.prev()){
               cNs = cIt.kmer();
               cString = B3Nucleotide.toString(cNs);

               B3NSmithWaterman align = new B3NSmithWaterman(pickedNs, cNs);
               align.process();
               cScore = align.score();

               if(cScore >= inScore){
//                   System.out.println(">seq_"+ii);
                   ii++;
//                   System.out.println(cString);
//                       alignSequences.push(cString);
//                       as2 = B3Nucleotide.toString(align.getAlignedS2());
//                       alignSequences.push(as2);
                   alignSequences.push(cString);
                   alignScores.push((int)cScore);


                   missRows = 0;
               }
               else{
                   missRows++;
               }
               if(missRows > rowTollMin){
                   break;
               }
           }
        }


        missRows = 0;

        if(pickedIt != null && pickedIt.hasNext()){
           cIt = pickedIt.clone();

           while(cIt.next()){
               cNs = cIt.kmer();
               cString = B3Nucleotide.toString(cNs);

               B3NSmithWaterman align = new B3NSmithWaterman(pickedNs, cNs);
               align.process();
               cScore = align.score();

               if(cScore >= inScore){
                   System.out.println(">seq_"+ii);
                   ii++;
                   System.out.println(cString);
//                       alignSequences.push(cString);
//                       as2 = B3Nucleotide.toString(align.getAlignedS2());
//                       alignSequences.push(as2);
                   alignSequences.push(cString);
                   alignScores.push((int)cScore);


                   missRows = 0;
               }
               else{
                   missRows++;
               }
               if(missRows > rowTollMin){
                   break;
               }
           }
        }
    }
    
    private void cs(){
        int[][] ascores = new int[this.alignSequences.size()][this.alignSequences.size()];
        int[] sums = new int[this.alignSequences.size()];
        
        
        int maxScore = Integer.MIN_VALUE;
        int maxI = -1;
        int maxJ = -1;
        
        
        for(int i=0; i<ascores.length; i++){
            for(int j=i+1; j<ascores[i].length; j++){
                StringSmithWaterman align = new StringSmithWaterman(alignSequences.get(i), alignSequences.get(j));
                ascores[i][j] = align.score();
                ascores[j][i] = ascores[i][j];
                
                sums[i] += ascores[i][j];
                sums[j] += ascores[j][i];
            }
        }
        
        
        
        for(int i=0; i<ascores.length; i++){
            if(sums[i]  > maxScore){
                maxScore = sums[i];
                maxI = i;
            }
        }
        maxScore = Integer.MIN_VALUE;
        for(int  j=0; j<ascores[maxI].length; j++){
            if(ascores[maxI][j] > maxScore){
                maxScore = ascores[maxI][j];
                maxJ = j;
            }
        }
        
        StringSmithWaterman align = new StringSmithWaterman(alignSequences.get(maxI), alignSequences.get(maxJ));
        align.process();
        
        String[] aligned = {align.getAlignedS1(), align.getAlignedS2()};
        
        
        
        getPattern(aligned);
        
        
    }
    
    private void msa(){
        String[] selected = new String[alignSequences.size()];
        alignSequences.toArray(selected);
        final String[] aligned = OpalMSAAPI.align(selected, new String[]{});

        if(remNs){
            double[] scores = new double[aligned[0].length()];
            int gapsCount = 0;
            int totCount = aligned[0].length();

            for(int i=0; i<aligned.length; i++){
                for(int j=0; j<aligned[i].length(); j++){
                    if(aligned[i].charAt(j) == '-'){
                        scores[j]++;
                    }
                }
            }

            for(int i=0; i<scores.length;i++){
                if(scores[i] / (double)aligned.length > 0.60){
                    scores[i] = 1;
                    gapsCount++;
                }
                else{
                    scores[i] = 0;
                }
            }

            if(gapsCount > 0){
                for(int i=0; i<aligned.length; i++){
                    char[] newS = new char[totCount - gapsCount];
                    for(int j=0, k=0; j<totCount; j++){
                        if(scores[j] == 0){
                            newS[k] = aligned[i].charAt(j);
                            k++;
                        }
                    }
                    aligned[i] = new String(newS);
                }
            }
        }
        
        
        Integer[] scores = new Integer[alignScores.size()];
        alignScores.toArray(scores);

        this.fillTable(aligned, scores);
        
        getPattern(aligned);
    }
    
    
    
    private int imax(int[] a){
        int max = Integer.MIN_VALUE;
        int imax = 0;
        for(int i=0; i<a.length; i++)
            if(a[i] > max){
                max = a[i];
                imax = i;
            }
        return imax;
    }
    private void getPattern(String[] seqs){
        int count = seqs.length;
        
        int[][] scores = new int[seqs[0].length()][5];
        
        char c;
        for(int i=0; i<seqs.length; i++){
            for(int j=0; j<seqs[i].length(); j++){
               c = seqs[i].charAt(j);
               switch(c){
                   case 'A':
                       scores[j][0]++;
                       break;
                   case 'C':
                       scores[j][1]++;
                       break;
                   case 'G':
                       scores[j][2]++;
                       break;
                   case 'T':
                       scores[j][3]++;
                       break;
                   case 'N':
                       scores[j][4]++;
                       break;
                   case '-':
                       scores[j][4]++;
                       break;
                   default:
                       scores[j][4]++;
                       break;
               }
            }
        }
        
        char[] pattern;
         
        if(remNs){
            int countSkip = 0;
            int max;
            for(int i=0; i<scores.length; i++){
                max = imax(scores[i]);
                if(max == 4)
                    countSkip++;
            }
            
            pattern = new char[scores.length - countSkip];
            
            for(int i=0, ip=0; i<pattern.length; i++){
                max = imax(scores[i]);
                if(max != 4){
                    switch(max){
                   case 0:
                       pattern[ip] = 'A';
                       break;
                   case 1:
                       pattern[ip] = 'C';
                       break;
                   case 2:
                       pattern[ip] = 'G';
                       break;
                   case 3:
                       pattern[ip] = 'T';
                       break;
                    }
                    
                    ip++;
                }
            }
        }
        else{
            pattern = new char[scores.length];
            
            int max;
            for(int i=0; i<pattern.length; i++){
                max = imax(scores[i]);
                switch(max){
                   case 0:
                       pattern[i] = 'A';
                       break;
                   case 1:
                       pattern[i] = 'C';
                       break;
                   case 2:
                       pattern[i] = 'G';
                       break;
                   case 3:
                       pattern[i] = 'T';
                       break;
                   case 4:
                       pattern[i] = 'N';
                       break;
               }
            }
        }
        
        System.out.println("pattern: "+ new String(pattern));
        this.patternField.setText(new String(pattern));
        this.kmerLabel.setText(new String(pattern));
    }
    
    private void fillTable(){
        String[] aseqs = new String[alignSequences.size()];
        aseqs = alignSequences.toArray(aseqs);
        
        Integer[] asc = new Integer[alignScores.size()];
        asc = alignScores.toArray(asc);
        
        
        if(aseqs.length > 0 )
        fillTable(aseqs, asc);
    }
    
    public  void fillTable(String[] sequences, Integer[] scores){
        tableModel.setData(sequences, scores);
        
        Font f = new java.awt.Font("Andale Mono", 0, 15);
        FontMetrics fm = this.getFontMetrics(this.getFont());
        int w = fm.charWidth('A');
//        int h = this.getFont().getSize();
//        this.setPreferredSize(new Dimension(w * text.length(),h));
//        kmersTable.getColumnModel().getColumn(0).setMinWidth(w * sequences[0].length());
//        this.setMinimumSize(new Dimension((w * sequences[0].length()) +50, 600));
        AffineTransform affinetransform = new AffineTransform();     
        FontRenderContext frc = new FontRenderContext(affinetransform,true,true);     
        kmersTable.getColumnModel().getColumn(0).setMinWidth((int)f.getStringBounds(sequences[0], frc).getWidth());
//        kmersTable.setMinimumSize(new java.awt.Dimension((int)f.getStringBounds(sequences[0], frc).getWidth()+75, this.getHeight()));
//        kmersTable.setPreferredSize(new java.awt.Dimension((int)f.getStringBounds(sequences[0], frc).getWidth()+75, this.getHeight()));
//        
//        kmersTable.setPreferredScrollableViewportSize(kmersTable.getMinimumSize());
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
        patternField = new javax.swing.JTextField();
        jLabel1 = new javax.swing.JLabel();
        kmerLabeljScrollPane = new javax.swing.JScrollPane();
        kmerLabel = new igtools.gui2.codistances.KmerLabel();
        jPanel2 = new javax.swing.JPanel();
        jScrollPane1 = new javax.swing.JScrollPane();
        kmersTable = new javax.swing.JTable();

        setDefaultCloseOperation(javax.swing.WindowConstants.EXIT_ON_CLOSE);

        patternField.setEditable(false);
        patternField.setFont(new java.awt.Font("Andale Mono", 0, 15)); // NOI18N

        jLabel1.setText("Pattern");

        kmerLabeljScrollPane.setOpaque(false);
        kmerLabeljScrollPane.setViewportView(kmerLabel);

        javax.swing.GroupLayout jPanel1Layout = new javax.swing.GroupLayout(jPanel1);
        jPanel1.setLayout(jPanel1Layout);
        jPanel1Layout.setHorizontalGroup(
            jPanel1Layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGroup(jPanel1Layout.createSequentialGroup()
                .addContainerGap()
                .addGroup(jPanel1Layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
                    .addComponent(patternField)
                    .addGroup(jPanel1Layout.createSequentialGroup()
                        .addComponent(jLabel1)
                        .addGap(920, 920, 920))
                    .addComponent(kmerLabeljScrollPane))
                .addContainerGap())
        );
        jPanel1Layout.setVerticalGroup(
            jPanel1Layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGroup(jPanel1Layout.createSequentialGroup()
                .addContainerGap()
                .addComponent(jLabel1)
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                .addComponent(patternField, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE)
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                .addComponent(kmerLabeljScrollPane, javax.swing.GroupLayout.PREFERRED_SIZE, 40, javax.swing.GroupLayout.PREFERRED_SIZE)
                .addContainerGap())
        );

        getContentPane().add(jPanel1, java.awt.BorderLayout.PAGE_END);

        jPanel2.setLayout(new java.awt.BorderLayout());

        jScrollPane1.setFont(new java.awt.Font("Andale Mono", 0, 15)); // NOI18N

        kmersTable.setModel(tableModel);
        jScrollPane1.setViewportView(kmersTable);

        jPanel2.add(jScrollPane1, java.awt.BorderLayout.CENTER);

        getContentPane().add(jPanel2, java.awt.BorderLayout.CENTER);

        pack();
    }// </editor-fold>//GEN-END:initComponents

    /**
     * @param args the command line arguments
     */
    public static void main(String args[]) {
        /* Set the Nimbus look and feel */
        //<editor-fold defaultstate="collapsed" desc=" Look and feel setting code (optional) ">
        /* If Nimbus (introduced in Java SE 6) is not available, stay with the default look and feel.
         * For details see http://download.oracle.com/javase/tutorial/uiswing/lookandfeel/plaf.html 
         */
        try {
            for (javax.swing.UIManager.LookAndFeelInfo info : javax.swing.UIManager.getInstalledLookAndFeels()) {
                if ("Nimbus".equals(info.getName())) {
                    javax.swing.UIManager.setLookAndFeel(info.getClassName());
                    break;
                }
            }
        } catch (ClassNotFoundException ex) {
            java.util.logging.Logger.getLogger(VBlockFrame.class.getName()).log(java.util.logging.Level.SEVERE, null, ex);
        } catch (InstantiationException ex) {
            java.util.logging.Logger.getLogger(VBlockFrame.class.getName()).log(java.util.logging.Level.SEVERE, null, ex);
        } catch (IllegalAccessException ex) {
            java.util.logging.Logger.getLogger(VBlockFrame.class.getName()).log(java.util.logging.Level.SEVERE, null, ex);
        } catch (javax.swing.UnsupportedLookAndFeelException ex) {
            java.util.logging.Logger.getLogger(VBlockFrame.class.getName()).log(java.util.logging.Level.SEVERE, null, ex);
        }
        //</editor-fold>

        /* Create and display the form */
        java.awt.EventQueue.invokeLater(new Runnable() {
            public void run() {
                new VBlockFrame("").setVisible(true);
            }
        });
    }
    // Variables declaration - do not modify//GEN-BEGIN:variables
    private javax.swing.JLabel jLabel1;
    private javax.swing.JPanel jPanel1;
    private javax.swing.JPanel jPanel2;
    private javax.swing.JScrollPane jScrollPane1;
    private igtools.gui2.codistances.KmerLabel kmerLabel;
    private javax.swing.JScrollPane kmerLabeljScrollPane;
    private javax.swing.JTable kmersTable;
    private javax.swing.JTextField patternField;
    // End of variables declaration//GEN-END:variables
}
