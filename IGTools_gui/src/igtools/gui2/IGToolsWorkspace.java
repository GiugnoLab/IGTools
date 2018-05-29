/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package igtools.gui2;

import igtools.gui2.codistances.DistalRecurrencesFrame;
import igtools.gui2.codistances.DistalRecurrencesFrame2;
import igtools.gui2.codistances.distributions.BiEDistalDistributionsFrame;
import igtools.gui2.codistances.distributions.DistalDistributionsFrame;
//import igtools.gui2.codistances.distributions.DistrCoefficientsFrame;
//import igtools.gui2.codistances.distributions.DistrCoefficientsFrame2;
import igtools.gui2.codistances.distributions.ProperEDistalDistributionsFrame;
import igtools.gui2.codistances.distributions.WordDistance;
import igtools.gui2.dictionaries.FlatDictionaryFrame;
import igtools.gui2.dictionaries.GeneralTrendsFrame;
import igtools.gui2.dictionaries.KDictsSizesFrame;
import igtools.gui2.efrequencies.EFreqienciesFrame;
import igtools.gui2.mults.CoMultFrame;
import igtools.gui2.mults.K6HeatMapFrame;
import igtools.gui2.mults.LabeledMultiplicitiesFrame;
import igtools.gui2.mults.MultiplicitiesFrame;
//import igtools.gui2.parikh.ParikhDistributionsFrame;
//import igtools.gui2.parikh.ParikhSizesFrame;
import igtools.gui2.sequence.GenomeFrame;
import igtools.gui2.words.FileWordProperties;

import java.awt.event.MouseAdapter;
import java.awt.event.MouseEvent;
import java.awt.event.WindowAdapter;
import java.awt.event.WindowEvent;
import java.io.File;
import java.nio.file.Files;
import java.text.DecimalFormat;
import java.util.Arrays;
import java.util.Enumeration;
import java.util.List;
import javax.swing.DefaultListModel;
import javax.swing.JCheckBox;
import javax.swing.JList;
import javax.swing.JOptionPane;
import static javax.swing.WindowConstants.DISPOSE_ON_CLOSE;

/**
 *
 * @author vbonnici
 */
public class IGToolsWorkspace extends javax.swing.JFrame {
    
    private DecimalFormat iformat = new DecimalFormat("###,###,###,###.####");

     private WSConfiguration config = null;
     private ResourcesManager resourcesManager = null;
     private List<WSSequence> wssequences = null;
     private WSSeqListModel listModel = null;
     private String lineSep = "--------------------------------------------------";
     
    /**
     * Creates new form IGToolsWorkspace
     */
    public IGToolsWorkspace(WSConfiguration config) {
        this.config = config;
        
        wssequences = config.getWSSequences();
        listModel = new WSSeqListModel(wssequences);
        resourcesManager = new ResourcesManager(wssequences);
        
        initComponents();
       
        
        jList1.addMouseListener(new MouseAdapter() {
            public void mouseClicked(MouseEvent evt) {
                JList list = (JList)evt.getSource();
                if (evt.getClickCount() == 2) {
                    int index = list.locationToIndex(evt.getPoint());
                    //if(index == 0){
                    //    addNewSequence();
                    //}
                    //else{
                        editSequence(index);
                    //}
                }
            }
        });
        
        
        addWindowListener(new WindowAdapter() {
            public void windowClosing(WindowEvent e) {
              exit();
            }
          });
    }
    
    
    private void addNewSequence(){
        WSSequence seq = new WSSequence();
        WSSequenceEditDialog dialog = new WSSequenceEditDialog(this,true,seq);
        dialog.setVisible(true);
        if(dialog.getAnswer()){
            listModel.add(seq);
            //listModel.addElement(seq);
            //jList1.invalidate();
        }
    }
    
    private void editSequence(int index){
        WSSequenceEditDialog dialog = new WSSequenceEditDialog(this,true,(WSSequence)listModel.getElementAt(index));
        dialog.setVisible(true);
    }
    
    
    private void removeSelectedSequences(){
        int[] indexes = this.jList1.getSelectedIndices();
        if(indexes!=null && indexes.length !=0){
            //if(indexes.length == 1 && indexes[0] == 0){
            //}
            //else{
            
                /*JOptionPane.showInputDialog(
                            this, 
                            "Are you sure to delete selected sequences?", 
                            "Delete sequences and data", 
                            JOptionPane.CANCEL_OPTION
                        );*/
                
                JCheckBox checkbox = new JCheckBox("Delete also data.");  
                String message = "Are you sure you want to delete the selected sequences?";  
                Object[] params = {message, checkbox};  
                int n = JOptionPane.showConfirmDialog(this, params, "Disconnect Products", JOptionPane.YES_NO_OPTION);  
                boolean deleteData = checkbox.isSelected();
                if(n == JOptionPane.OK_OPTION){
                    if(deleteData)
                        System.out.println("deleting sequences and data...");
                    else
                        System.out.println("deleting sequences...");
                    
                    
                    if(deleteData){
                        for(int i=0; i<indexes.length; i++){
                           // if(indexes[i] != 0){
                                WSSequence seq = (WSSequence)listModel.getElementAt(indexes[i]);
                                if(seq.getB3seqFile() != null){
                                    try{
                                        Files.deleteIfExists(new File(seq.getB3seqFile()).toPath());
                                        seq.setB3seq(null);
                                    }catch(Exception ex){
                                        ex.printStackTrace();
                                    }
                                }
                                if(seq.getNELSAFile() != null){
                                    try{
                                        Files.deleteIfExists(new File(seq.getNELSAFile()).toPath());
                                        seq.setNELSA(null);
                                    }catch(Exception ex){
                                        ex.printStackTrace();
                                    }
                                }
                            //}
                        }
                    }
                    
                    System.out.println("done.");
                }
                
                
                //Arrays.sort(indexes);
                //int minus = 0;
                for(int i=0; i<indexes.length; i++){
                    listModel.remove(indexes[i]);
                    //if(indexes[i] != 0){
                        //this.listModel.removeElementAt(indexes[i] - minus);
                        //this.wssequences.remove(indexes[i] - minus -1);
                        //minus++;
                    //}
                }
                
            //}
        }
    }
    
    
    
    private void exit(){
        int confirmed = JOptionPane.showConfirmDialog(null, 
                  "Are you sure you want to exit the program?", "Exit Program Message Box",
                  JOptionPane.YES_NO_OPTION);

        
        System.out.println(confirmed+" "+JOptionPane.YES_OPTION+" " +JOptionPane.NO_OPTION);
        if (confirmed == JOptionPane.YES_OPTION) {
          try{
              if(wssequences != null){
                  config.save(wssequences);
              }
          }catch(Exception ex){
              ex.printStackTrace();
          }
          
          System.out.println("exit");
          setDefaultCloseOperation(EXIT_ON_CLOSE);
          dispose();
        }
        else{
            setDefaultCloseOperation(DO_NOTHING_ON_CLOSE);
        }
    }

    /**
     * This method is called from within the constructor to initialize the form.
     * WARNING: Do NOT modify this code. The content of this method is always
     * regenerated by the Form Editor.
     */
    @SuppressWarnings("unchecked")
    // <editor-fold defaultstate="collapsed" desc="Generated Code">//GEN-BEGIN:initComponents
    private void initComponents() {

        jScrollPane1 = new javax.swing.JScrollPane();
        jList1 = new javax.swing.JList();
        jMenuBar1 = new javax.swing.JMenuBar();
        jMenu1 = new javax.swing.JMenu();
        jMenuItem4 = new javax.swing.JMenuItem();
        toolsMenu = new javax.swing.JMenu();
        jMenuItem7 = new javax.swing.JMenuItem();
        jSeparator2 = new javax.swing.JPopupMenu.Separator();
        jMenuItem8 = new javax.swing.JMenuItem();
        jSeparator4 = new javax.swing.JPopupMenu.Separator();
        jMenuItem9 = new javax.swing.JMenuItem();
        jMenuItem15 = new javax.swing.JMenuItem();
        jMenuItem22 = new javax.swing.JMenuItem();
        jSeparator5 = new javax.swing.JPopupMenu.Separator();
        jMenuItem21 = new javax.swing.JMenuItem();
        GeneralTrendsMenuItem = new javax.swing.JMenuItem();
        jSeparator10 = new javax.swing.JPopupMenu.Separator();
        jMenuItem10 = new javax.swing.JMenuItem();
        jMenuItem18 = new javax.swing.JMenuItem();
        jMenuItem12 = new javax.swing.JMenuItem();
        jMenuItem13 = new javax.swing.JMenuItem();
        jSeparator6 = new javax.swing.JPopupMenu.Separator();
        jMenuItem11 = new javax.swing.JMenuItem();
        jSeparator7 = new javax.swing.JPopupMenu.Separator();
        mmeprobsMenu = new javax.swing.JMenuItem();
        jSeparator8 = new javax.swing.JPopupMenu.Separator();
        jSeparator9 = new javax.swing.JPopupMenu.Separator();
        jMenuItem19 = new javax.swing.JMenuItem();
        jMenuItem23 = new javax.swing.JMenuItem();
        jMenuItem24 = new javax.swing.JMenuItem();
        jMenuItem26 = new javax.swing.JMenuItem();
        jMenuItem27 = new javax.swing.JMenuItem();
        distrCoefficientsMenuItem = new javax.swing.JMenuItem();
        jSeparator11 = new javax.swing.JPopupMenu.Separator();
        jMenuItem25 = new javax.swing.JMenuItem();
        jSeparator12 = new javax.swing.JPopupMenu.Separator();
        jMenuItem6 = new javax.swing.JMenuItem();
        jMenuItem5 = new javax.swing.JMenuItem();
        jSeparator13 = new javax.swing.JPopupMenu.Separator();
        jMenuItem28 = new javax.swing.JMenuItem();
        jMenu4 = new javax.swing.JMenu();
        jMenuItem20 = new javax.swing.JMenuItem();
        jMenuItem1 = new javax.swing.JMenuItem();
        jSeparator3 = new javax.swing.JPopupMenu.Separator();
        jSeparator1 = new javax.swing.JPopupMenu.Separator();
        jMenuItem3 = new javax.swing.JMenuItem();
        jMenuItem2 = new javax.swing.JMenuItem();

        setDefaultCloseOperation(javax.swing.WindowConstants.EXIT_ON_CLOSE);
        setTitle("InfoGenomics Tools - Workspace");

        jList1.setModel(listModel.getModel());
        jList1.setCellRenderer(new WSSeqCellRenderer());
        jScrollPane1.setViewportView(jList1);

        jMenu1.setText("File");

        jMenuItem4.setText("Exit");
        jMenuItem4.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                jMenuItem4ActionPerformed(evt);
            }
        });
        jMenu1.add(jMenuItem4);

        jMenuBar1.add(jMenu1);

        toolsMenu.setText("Tools");

        jMenuItem7.setText("Sequence info");
        jMenuItem7.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                jMenuItem7ActionPerformed(evt);
            }
        });
        toolsMenu.add(jMenuItem7);
        toolsMenu.add(jSeparator2);

//        jMenuItem8.setText("Informational analysis");
//        toolsMenu.add(jMenuItem8);
//        toolsMenu.add(jSeparator4);

        jMenuItem9.setText("Sequence viewer");
        jMenuItem9.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                jMenuItem9ActionPerformed(evt);
            }
        });
        toolsMenu.add(jMenuItem9);

        jMenuItem15.setText("Dictionary viewer");
        jMenuItem15.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                jMenuItem15ActionPerformed(evt);
            }
        });
        toolsMenu.add(jMenuItem15);

//        jMenuItem22.setText("LCP viewer");
//        toolsMenu.add(jMenuItem22);
//        toolsMenu.add(jSeparator5);

        jMenuItem21.setText("Size trends");
        jMenuItem21.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                jMenuItem21ActionPerformed(evt);
            }
        });
        toolsMenu.add(jMenuItem21);

        GeneralTrendsMenuItem.setText("General trends");
        GeneralTrendsMenuItem.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                GeneralTrendsMenuItemActionPerformed(evt);
            }
        });
        toolsMenu.add(GeneralTrendsMenuItem);
        toolsMenu.add(jSeparator10);

        jMenuItem10.setText("Multiplicity");
        jMenuItem10.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                jMenuItem10ActionPerformed(evt);
            }
        });
        toolsMenu.add(jMenuItem10);

        jMenuItem18.setText("Labeled Multiplicity");
        jMenuItem18.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                jMenuItem18ActionPerformed(evt);
            }
        });
        toolsMenu.add(jMenuItem18);

//        jMenuItem12.setText("Multiplicity heatmap");
//        toolsMenu.add(jMenuItem12);

        jMenuItem13.setText("[1,6] tree heatmap");
        jMenuItem13.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                jMenuItem13ActionPerformed(evt);
            }
        });
        toolsMenu.add(jMenuItem13);
        toolsMenu.add(jSeparator6);

        jMenuItem11.setText("Co-Multiplicity");
        jMenuItem11.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                jMenuItem11ActionPerformed(evt);
            }
        });
        toolsMenu.add(jMenuItem11);
        toolsMenu.add(jSeparator7);

        mmeprobsMenu.setText("Expected frequencies");
        mmeprobsMenu.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                mmeprobsMenuActionPerformed(evt);
            }
        });
        toolsMenu.add(mmeprobsMenu);
        toolsMenu.add(jSeparator8);
        toolsMenu.add(jSeparator9);

        jMenuItem19.setText("Pair-RDD");
        jMenuItem19.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                jMenuItem19ActionPerformed(evt);
            }
        });
        toolsMenu.add(jMenuItem19);

//        jMenuItem23.setText("Distal rec2");
//        jMenuItem23.setEnabled(false);
//        jMenuItem23.addActionListener(new java.awt.event.ActionListener() {
//            public void actionPerformed(java.awt.event.ActionEvent evt) {
//                jMenuItem23ActionPerformed(evt);
//            }
//        });
//        toolsMenu.add(jMenuItem23);
//
//        jMenuItem24.setText("Distal distribution");
//        jMenuItem24.setEnabled(false);
//        jMenuItem24.addActionListener(new java.awt.event.ActionListener() {
//            public void actionPerformed(java.awt.event.ActionEvent evt) {
//                jMenuItem24ActionPerformed(evt);
//            }
//        });
//        toolsMenu.add(jMenuItem24);

        jMenuItem26.setText("Self-RDD");
        jMenuItem26.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                jMenuItem26ActionPerformed(evt);
            }
        });
        toolsMenu.add(jMenuItem26);

//        jMenuItem27.setText("Bi est distr");
//        jMenuItem27.setEnabled(false);
//        jMenuItem27.addActionListener(new java.awt.event.ActionListener() {
//            public void actionPerformed(java.awt.event.ActionEvent evt) {
//                jMenuItem27ActionPerformed(evt);
//            }
//        });
//        toolsMenu.add(jMenuItem27);

//        distrCoefficientsMenuItem.setText("Distr. Coefficients");
//        distrCoefficientsMenuItem.addActionListener(new java.awt.event.ActionListener() {
//            public void actionPerformed(java.awt.event.ActionEvent evt) {
//                distrCoefficientsMenuItemActionPerformed(evt);
//            }
//        });
//        toolsMenu.add(distrCoefficientsMenuItem);
//        toolsMenu.add(jSeparator11);

//        jMenuItem25.setText("WordList Prop.");
//        jMenuItem25.setEnabled(false);
//        jMenuItem25.addActionListener(new java.awt.event.ActionListener() {
//            public void actionPerformed(java.awt.event.ActionEvent evt) {
//                jMenuItem25ActionPerformed(evt);
//            }
//        });
//        toolsMenu.add(jMenuItem25);
//        toolsMenu.add(jSeparator12);

//        jMenuItem6.setText("Parikh class sizes");
//        jMenuItem6.addActionListener(new java.awt.event.ActionListener() {
//            public void actionPerformed(java.awt.event.ActionEvent evt) {
//                jMenuItem6ActionPerformed(evt);
//            }
//        });
//        toolsMenu.add(jMenuItem6);
//
//        jMenuItem5.setText("Pairkh distributions");
//        jMenuItem5.addActionListener(new java.awt.event.ActionListener() {
//            public void actionPerformed(java.awt.event.ActionEvent evt) {
//                jMenuItem5ActionPerformed(evt);
//            }
//        });
//        toolsMenu.add(jMenuItem5);
//        toolsMenu.add(jSeparator13);

//        jMenuItem28.setText("Word distance");
//        jMenuItem28.setEnabled(false);
//        jMenuItem28.addActionListener(new java.awt.event.ActionListener() {
//            public void actionPerformed(java.awt.event.ActionEvent evt) {
//                jMenuItem28ActionPerformed(evt);
//            }
//        });
//        toolsMenu.add(jMenuItem28);

        jMenuBar1.add(toolsMenu);

        jMenu4.setText("Resources");

        jMenuItem20.setText("Add sequence");
        jMenuItem20.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                jMenuItem20ActionPerformed(evt);
            }
        });
        jMenu4.add(jMenuItem20);

        jMenuItem1.setText("Delete sequences");
        jMenuItem1.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                jMenuItem1ActionPerformed(evt);
            }
        });
        jMenu4.add(jMenuItem1);
        jMenu4.add(jSeparator3);
        jMenu4.add(jSeparator1);

        jMenuItem3.setText("Release all data");
        jMenuItem3.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                jMenuItem3ActionPerformed(evt);
            }
        });
        jMenu4.add(jMenuItem3);

        jMenuItem2.setText("Release selected data");
        jMenuItem2.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                jMenuItem2ActionPerformed(evt);
            }
        });
        jMenu4.add(jMenuItem2);

        jMenuBar1.add(jMenu4);

        setJMenuBar(jMenuBar1);

        javax.swing.GroupLayout layout = new javax.swing.GroupLayout(getContentPane());
        getContentPane().setLayout(layout);
        layout.setHorizontalGroup(
            layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addComponent(jScrollPane1, javax.swing.GroupLayout.DEFAULT_SIZE, 400, Short.MAX_VALUE)
        );
        layout.setVerticalGroup(
            layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addComponent(jScrollPane1, javax.swing.GroupLayout.DEFAULT_SIZE, 296, Short.MAX_VALUE)
        );

        pack();
    }// </editor-fold>//GEN-END:initComponents

    private void jMenuItem1ActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_jMenuItem1ActionPerformed
       removeSelectedSequences();
    }//GEN-LAST:event_jMenuItem1ActionPerformed

    private void jMenuItem3ActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_jMenuItem3ActionPerformed
        Enumeration<WSSequence> seqs = listModel.getModel().elements();
        while(seqs.hasMoreElements()){
            WSSequence seq = seqs.nextElement();
            if(seq.getB3seq() != null){
                System.out.println("Releasing sequence "+ seq.getName());
                seq.setB3seq(null);
            }
            if(seq.getNELSA() != null){
                System.out.println("Releasing NELSA "+ seq.getName());
                seq.setNELSA(null);
            }
        }
        System.out.println("Calling GC...");
        System.gc();
        System.gc();
        System.gc();
        System.out.println("done.");
    }//GEN-LAST:event_jMenuItem3ActionPerformed

    private void jMenuItem4ActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_jMenuItem4ActionPerformed
        exit();
    }//GEN-LAST:event_jMenuItem4ActionPerformed

    private void jMenuItem7ActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_jMenuItem7ActionPerformed
        seqsInfo();
    }//GEN-LAST:event_jMenuItem7ActionPerformed

    private void jMenuItem20ActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_jMenuItem20ActionPerformed
        addNewSequence();
    }//GEN-LAST:event_jMenuItem20ActionPerformed

    private void jMenuItem9ActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_jMenuItem9ActionPerformed
        int[] indexes = this.jList1.getSelectedIndices();
        for(int i=0; i<indexes.length; i++){
            final WSSequence seq = (WSSequence)listModel.getElementAt(indexes[i]);
            resourcesManager.loadB3seq(indexes[i], false);
        
            java.awt.EventQueue.invokeLater(new Runnable() {
                public void run() {
                    GenomeFrame frame = new GenomeFrame(seq);
                    frame.setDefaultCloseOperation(DISPOSE_ON_CLOSE);
                    frame.setVisible(true);
                }
            });
        
        }
    }//GEN-LAST:event_jMenuItem9ActionPerformed

    private void jMenuItem15ActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_jMenuItem15ActionPerformed
       int[] indexes = this.jList1.getSelectedIndices();
        for(int i=0; i<indexes.length; i++){
            final WSSequence seq = (WSSequence)listModel.getElementAt(indexes[i]);
            resourcesManager.loadB3seq(indexes[i], false);
            resourcesManager.loadNELSA(indexes[i], false);
            
            java.awt.EventQueue.invokeLater(new Runnable() {
                public void run() {
                    FlatDictionaryFrame frame = new FlatDictionaryFrame(seq);
                    frame.setDefaultCloseOperation(DISPOSE_ON_CLOSE);
                    frame.setVisible(true);
                }
            });
        
        }
    }//GEN-LAST:event_jMenuItem15ActionPerformed

    private void jMenuItem21ActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_jMenuItem21ActionPerformed
        int[] indexes = this.jList1.getSelectedIndices();
        for(int i=0; i<indexes.length; i++){
            final WSSequence seq = (WSSequence)listModel.getElementAt(indexes[i]);
            resourcesManager.loadB3seq(indexes[i], false);
            resourcesManager.loadNELSA(indexes[i], false);
            
            java.awt.EventQueue.invokeLater(new Runnable() {
                public void run() {
                    KDictsSizesFrame frame = new KDictsSizesFrame(seq);
                    frame.setDefaultCloseOperation(DISPOSE_ON_CLOSE);
                    frame.setVisible(true);
                }
            });
        
        }
    }//GEN-LAST:event_jMenuItem21ActionPerformed

    private void jMenuItem2ActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_jMenuItem2ActionPerformed
        int[] indexes = this.jList1.getSelectedIndices();
        for(int i=0; i<indexes.length; i++){
            final WSSequence seq = (WSSequence)listModel.getElementAt(indexes[i]);
             if(seq.getB3seq() != null){
                System.out.println("Releasing sequence "+ seq.getName());
                seq.setB3seq(null);
            }
            if(seq.getNELSA() != null){
                System.out.println("Releasing NELSA "+ seq.getName());
                seq.setNELSA(null);
            }
        }
        System.out.println("Calling GC...");
        System.gc();
        System.gc();
        System.gc();
        System.out.println("done.");
    }//GEN-LAST:event_jMenuItem2ActionPerformed

    private void jMenuItem10ActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_jMenuItem10ActionPerformed
        int[] indexes = this.jList1.getSelectedIndices();
        for(int i=0; i<indexes.length; i++){
            final WSSequence seq = (WSSequence)listModel.getElementAt(indexes[i]);
            resourcesManager.loadB3seq(indexes[i], false);
            resourcesManager.loadNELSA(indexes[i], false);
            
            java.awt.EventQueue.invokeLater(new Runnable() {
                public void run() {
                    MultiplicitiesFrame frame = new MultiplicitiesFrame(seq);
                    frame.setDefaultCloseOperation(DISPOSE_ON_CLOSE);
                    frame.setVisible(true);
                }
            });
        }
    }//GEN-LAST:event_jMenuItem10ActionPerformed

    private void jMenuItem13ActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_jMenuItem13ActionPerformed
       int[] indexes = this.jList1.getSelectedIndices();
        for(int i=0; i<indexes.length; i++){
            final WSSequence seq = (WSSequence)listModel.getElementAt(indexes[i]);
            resourcesManager.loadB3seq(indexes[i], false);
            
            java.awt.EventQueue.invokeLater(new Runnable() {
                public void run() {
                    K6HeatMapFrame frame = new K6HeatMapFrame(seq);
                    frame.setDefaultCloseOperation(DISPOSE_ON_CLOSE);
                    frame.setVisible(true);
                }
            });
        }
    }//GEN-LAST:event_jMenuItem13ActionPerformed

    private void jMenuItem11ActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_jMenuItem11ActionPerformed
        int[] indexes = this.jList1.getSelectedIndices();
        for(int i=0; i<indexes.length; i++){
            final WSSequence seq = (WSSequence)listModel.getElementAt(indexes[i]);
            resourcesManager.loadB3seq(indexes[i], false);
            resourcesManager.loadNELSA(indexes[i], false);
            
            java.awt.EventQueue.invokeLater(new Runnable() {
                public void run() {
                    CoMultFrame frame = new CoMultFrame(seq);
                    frame.setDefaultCloseOperation(DISPOSE_ON_CLOSE);
                    frame.setVisible(true);
                }
            });
        }
    }//GEN-LAST:event_jMenuItem11ActionPerformed

    private void jMenuItem18ActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_jMenuItem18ActionPerformed
        int[] indexes = this.jList1.getSelectedIndices();
        for(int i=0; i<indexes.length; i++){
            final WSSequence seq = (WSSequence)listModel.getElementAt(indexes[i]);
            resourcesManager.loadB3seq(indexes[i], false);
            resourcesManager.loadNELSA(indexes[i], false);
            
            java.awt.EventQueue.invokeLater(new Runnable() {
                public void run() {
                    LabeledMultiplicitiesFrame frame = new LabeledMultiplicitiesFrame(seq);
                    frame.setDefaultCloseOperation(DISPOSE_ON_CLOSE);
                    frame.setVisible(true);
                }
            });
        }
    }//GEN-LAST:event_jMenuItem18ActionPerformed

    private void jMenuItem19ActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_jMenuItem19ActionPerformed
        int[] indexes = this.jList1.getSelectedIndices();
        for(int i=0; i<indexes.length; i++){
            final WSSequence seq = (WSSequence)listModel.getElementAt(indexes[i]);
            resourcesManager.loadB3seq(indexes[i], false);
            resourcesManager.loadNELSA(indexes[i], false);
            
            java.awt.EventQueue.invokeLater(new Runnable() {
                public void run() {
                    DistalRecurrencesFrame frame = new DistalRecurrencesFrame(seq);
                    //BiEDistalDistributionsFrame frame = new BiEDistalDistributionsFrame(seq);
                    frame.setDefaultCloseOperation(DISPOSE_ON_CLOSE);
                    frame.setVisible(true);
                }
            });
        }
    }//GEN-LAST:event_jMenuItem19ActionPerformed

    private void jMenuItem23ActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_jMenuItem23ActionPerformed
       int[] indexes = this.jList1.getSelectedIndices();
        for(int i=0; i<indexes.length; i++){
            final WSSequence seq = (WSSequence)listModel.getElementAt(indexes[i]);
            resourcesManager.loadB3seq(indexes[i], false);
            resourcesManager.loadNELSA(indexes[i], false);
            
            java.awt.EventQueue.invokeLater(new Runnable() {
                public void run() {
                    DistalRecurrencesFrame2 frame = new DistalRecurrencesFrame2(seq);
                    frame.setDefaultCloseOperation(DISPOSE_ON_CLOSE);
                    frame.setVisible(true);
                }
            });
        }
    }//GEN-LAST:event_jMenuItem23ActionPerformed

    private void GeneralTrendsMenuItemActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_GeneralTrendsMenuItemActionPerformed
        int[] indexes = this.jList1.getSelectedIndices();
        for(int i=0; i<indexes.length; i++){
            final WSSequence seq = (WSSequence)listModel.getElementAt(indexes[i]);
            resourcesManager.loadB3seq(indexes[i], false);
            resourcesManager.loadNELSA(indexes[i], false);
            
            java.awt.EventQueue.invokeLater(new Runnable() {
                public void run() {
                    GeneralTrendsFrame  frame = new GeneralTrendsFrame(seq);
                    frame.setDefaultCloseOperation(DISPOSE_ON_CLOSE);
                    frame.setVisible(true);
                }
            });
        
        }
    }//GEN-LAST:event_GeneralTrendsMenuItemActionPerformed

    private void jMenuItem24ActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_jMenuItem24ActionPerformed
        int[] indexes = this.jList1.getSelectedIndices();
        for(int i=0; i<indexes.length; i++){
            final WSSequence seq = (WSSequence)listModel.getElementAt(indexes[i]);
            resourcesManager.loadB3seq(indexes[i], false);
            resourcesManager.loadNELSA(indexes[i], false);
            
            java.awt.EventQueue.invokeLater(new Runnable() {
                public void run() {
                    DistalDistributionsFrame frame = new DistalDistributionsFrame(seq);
                    frame.setDefaultCloseOperation(DISPOSE_ON_CLOSE);
                    frame.setVisible(true);
                }
            });
        }
    }//GEN-LAST:event_jMenuItem24ActionPerformed

    private void jMenuItem25ActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_jMenuItem25ActionPerformed
       int[] indexes = this.jList1.getSelectedIndices();
        for(int i=0; i<indexes.length; i++){
            final WSSequence seq = (WSSequence)listModel.getElementAt(indexes[i]);
            resourcesManager.loadB3seq(indexes[i], false);
            resourcesManager.loadNELSA(indexes[i], false);
            
            java.awt.EventQueue.invokeLater(new Runnable() {
                public void run() {
                    FileWordProperties frame = new FileWordProperties(seq);
                    frame.setDefaultCloseOperation(DISPOSE_ON_CLOSE);
                    frame.setVisible(true);
                }
            });
        }
    }//GEN-LAST:event_jMenuItem25ActionPerformed

    private void distrCoefficientsMenuItemActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_distrCoefficientsMenuItemActionPerformed
         int[] indexes = this.jList1.getSelectedIndices();
        for(int i=0; i<indexes.length; i++){
            final WSSequence seq = (WSSequence)listModel.getElementAt(indexes[i]);
            resourcesManager.loadB3seq(indexes[i], false);
            resourcesManager.loadNELSA(indexes[i], false);
            
            java.awt.EventQueue.invokeLater(new Runnable() {
                public void run() {
                    //DistrCoefficientsFrame2 frame = new DistrCoefficientsFrame2(seq);
                    //frame.setDefaultCloseOperation(DISPOSE_ON_CLOSE);
                    //frame.setVisible(true);
                }
            });
        }
    }//GEN-LAST:event_distrCoefficientsMenuItemActionPerformed

    private void jMenuItem26ActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_jMenuItem26ActionPerformed
        int[] indexes = this.jList1.getSelectedIndices();
        for(int i=0; i<indexes.length; i++){
            final WSSequence seq = (WSSequence)listModel.getElementAt(indexes[i]);
            resourcesManager.loadB3seq(indexes[i], false);
            resourcesManager.loadNELSA(indexes[i], false);
            
            java.awt.EventQueue.invokeLater(new Runnable() {
                public void run() {
                    ProperEDistalDistributionsFrame frame = new ProperEDistalDistributionsFrame(seq);
                    frame.setDefaultCloseOperation(DISPOSE_ON_CLOSE);
                    frame.setVisible(true);
                }
            });
        }
    }//GEN-LAST:event_jMenuItem26ActionPerformed

    private void jMenuItem27ActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_jMenuItem27ActionPerformed
        int[] indexes = this.jList1.getSelectedIndices();
        for(int i=0; i<indexes.length; i++){
            final WSSequence seq = (WSSequence)listModel.getElementAt(indexes[i]);
            resourcesManager.loadB3seq(indexes[i], false);
            resourcesManager.loadNELSA(indexes[i], false);
            
            java.awt.EventQueue.invokeLater(new Runnable() {
                public void run() {
                    BiEDistalDistributionsFrame frame = new BiEDistalDistributionsFrame(seq);
                    frame.setDefaultCloseOperation(DISPOSE_ON_CLOSE);
                    frame.setVisible(true);
                }
            });
        }
    }//GEN-LAST:event_jMenuItem27ActionPerformed

    private void mmeprobsMenuActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_mmeprobsMenuActionPerformed
        int[] indexes = this.jList1.getSelectedIndices();
        for(int i=0; i<indexes.length; i++){
            final WSSequence seq = (WSSequence)listModel.getElementAt(indexes[i]);
            resourcesManager.loadB3seq(indexes[i], false);
            resourcesManager.loadNELSA(indexes[i], false);
            
            java.awt.EventQueue.invokeLater(new Runnable() {
                public void run() {
                    EFreqienciesFrame frame = new EFreqienciesFrame(seq);
                    frame.setDefaultCloseOperation(DISPOSE_ON_CLOSE);
                    frame.setVisible(true);
                }
            });
        }
    }//GEN-LAST:event_mmeprobsMenuActionPerformed

    private void jMenuItem5ActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_jMenuItem5ActionPerformed
        int[] indexes = this.jList1.getSelectedIndices();
        for(int i=0; i<indexes.length; i++){
            final WSSequence seq = (WSSequence)listModel.getElementAt(indexes[i]);
            resourcesManager.loadB3seq(indexes[i], false);
            resourcesManager.loadNELSA(indexes[i], false);
            
            java.awt.EventQueue.invokeLater(new Runnable() {
                public void run() {
                   // ParikhDistributionsFrame frame = new ParikhDistributionsFrame(seq);
                    //frame.setDefaultCloseOperation(DISPOSE_ON_CLOSE);
                    //frame.setVisible(true);
                }
            });
        }
    }//GEN-LAST:event_jMenuItem5ActionPerformed

    private void jMenuItem6ActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_jMenuItem6ActionPerformed
        int[] indexes = this.jList1.getSelectedIndices();
        for(int i=0; i<indexes.length; i++){
            final WSSequence seq = (WSSequence)listModel.getElementAt(indexes[i]);
            resourcesManager.loadB3seq(indexes[i], false);
            resourcesManager.loadNELSA(indexes[i], false);
            
            java.awt.EventQueue.invokeLater(new Runnable() {
                public void run() {
                    //ParikhSizesFrame frame = new ParikhSizesFrame(seq);
                    //frame.setDefaultCloseOperation(DISPOSE_ON_CLOSE);
                    //frame.setVisible(true);
                }
            });
        }
    }//GEN-LAST:event_jMenuItem6ActionPerformed

    private void jMenuItem28ActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_jMenuItem28ActionPerformed
       int[] indexes = this.jList1.getSelectedIndices();
        for(int i=0; i<indexes.length; i++){
            final WSSequence seq = (WSSequence)listModel.getElementAt(indexes[i]);
            resourcesManager.loadB3seq(indexes[i], false);
            resourcesManager.loadNELSA(indexes[i], false);
            
            java.awt.EventQueue.invokeLater(new Runnable() {
                public void run() {
                    WordDistance frame = new WordDistance(seq);
                    frame.setDefaultCloseOperation(DISPOSE_ON_CLOSE);
                    frame.setVisible(true);
                }
            });
        }
    }//GEN-LAST:event_jMenuItem28ActionPerformed

    // Variables declaration - do not modify//GEN-BEGIN:variables
    private javax.swing.JMenuItem GeneralTrendsMenuItem;
    private javax.swing.JMenuItem distrCoefficientsMenuItem;
    private javax.swing.JList jList1;
    private javax.swing.JMenu jMenu1;
    private javax.swing.JMenu jMenu4;
    private javax.swing.JMenuBar jMenuBar1;
    private javax.swing.JMenuItem jMenuItem1;
    private javax.swing.JMenuItem jMenuItem10;
    private javax.swing.JMenuItem jMenuItem11;
    private javax.swing.JMenuItem jMenuItem12;
    private javax.swing.JMenuItem jMenuItem13;
    private javax.swing.JMenuItem jMenuItem15;
    private javax.swing.JMenuItem jMenuItem18;
    private javax.swing.JMenuItem jMenuItem19;
    private javax.swing.JMenuItem jMenuItem2;
    private javax.swing.JMenuItem jMenuItem20;
    private javax.swing.JMenuItem jMenuItem21;
    private javax.swing.JMenuItem jMenuItem22;
    private javax.swing.JMenuItem jMenuItem23;
    private javax.swing.JMenuItem jMenuItem24;
    private javax.swing.JMenuItem jMenuItem25;
    private javax.swing.JMenuItem jMenuItem26;
    private javax.swing.JMenuItem jMenuItem27;
    private javax.swing.JMenuItem jMenuItem28;
    private javax.swing.JMenuItem jMenuItem3;
    private javax.swing.JMenuItem jMenuItem4;
    private javax.swing.JMenuItem jMenuItem5;
    private javax.swing.JMenuItem jMenuItem6;
    private javax.swing.JMenuItem jMenuItem7;
    private javax.swing.JMenuItem jMenuItem8;
    private javax.swing.JMenuItem jMenuItem9;
    private javax.swing.JScrollPane jScrollPane1;
    private javax.swing.JPopupMenu.Separator jSeparator1;
    private javax.swing.JPopupMenu.Separator jSeparator10;
    private javax.swing.JPopupMenu.Separator jSeparator11;
    private javax.swing.JPopupMenu.Separator jSeparator12;
    private javax.swing.JPopupMenu.Separator jSeparator13;
    private javax.swing.JPopupMenu.Separator jSeparator2;
    private javax.swing.JPopupMenu.Separator jSeparator3;
    private javax.swing.JPopupMenu.Separator jSeparator4;
    private javax.swing.JPopupMenu.Separator jSeparator5;
    private javax.swing.JPopupMenu.Separator jSeparator6;
    private javax.swing.JPopupMenu.Separator jSeparator7;
    private javax.swing.JPopupMenu.Separator jSeparator8;
    private javax.swing.JPopupMenu.Separator jSeparator9;
    private javax.swing.JMenuItem mmeprobsMenu;
    private javax.swing.JMenu toolsMenu;
    // End of variables declaration//GEN-END:variables

    
    
    
    
    private void seqsInfo(){
        int[] indexes = this.jList1.getSelectedIndices();
        WSSequence seq;
        for(int i=0; i<indexes.length; i++){
            seq = (WSSequence)listModel.getElementAt(indexes[i]);
            System.out.println(seq.getName());
            resourcesManager.loadB3seq(indexes[i], false);
            //System.out.println(seq.getB3seq());
            int length = seq.getB3seq().length();
            System.out.println("Total length: "+iformat.format(length));
            int bads = seq.getB3seq().countBads();
            System.out.println("Number of N: "+iformat.format(bads));
            System.out.println("Effective length: "+iformat.format((length - bads)));
            System.out.println("N ratio: "+((double)bads / (double)length));
            System.out.println("Effective ratio: "+((double)(length-bads) / (double)length));
            System.out.println(lineSep);
        }
    }


}
