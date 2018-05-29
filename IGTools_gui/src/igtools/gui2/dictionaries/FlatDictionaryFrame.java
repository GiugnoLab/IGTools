/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package igtools.gui2.dictionaries;

import igtools.common.alignment.B3NSmithWaterman;
import igtools.common.alignment.StringSmithWaterman;
import igtools.common.charts.palette.B3NucleotidePaletteFactory;
import igtools.common.nucleotide.B3Nucleotide;
import igtools.common.sequence.B3LLSequence;
import igtools.dictionaries.elsa.IELSAIterator;
import igtools.dictionaries.elsa.NELSA;
import igtools.gui2.WSSequence;
//import igtools.gui.codistances.BiCoDistancesFrame;
import igtools.gui2.dictionaries.vblocks.VBlockFrame;
import igtools.gui2.sequence.GenomeFrame;
import java.awt.BorderLayout;
import java.awt.Color;
import java.awt.event.ActionEvent;
import java.awt.event.KeyEvent;
import java.awt.event.MouseAdapter;
import java.awt.event.MouseEvent;
import java.awt.event.MouseWheelEvent;
import java.awt.event.MouseWheelListener;
import java.util.Set;
import java.util.Stack;
import java.util.TreeSet;
import java.util.Vector;
import javax.swing.AbstractAction;
import javax.swing.ActionMap;
import javax.swing.DefaultComboBoxModel;
import javax.swing.InputMap;
import javax.swing.JComponent;
import javax.swing.KeyStroke;
import static javax.swing.WindowConstants.DISPOSE_ON_CLOSE;
import opal.OpalMSAAPI;

/**
 *
 * @author vbonnici
 */
public class FlatDictionaryFrame extends javax.swing.JFrame {

    private WSSequence wsseq = null;
    private String title;
     private NELSA nelsa;
     private IELSAIterator first_it = null;
     private FlatDictionaryPanel dictPanel = null;
     
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
     * Creates new form FlatDictionaryView
     */
    public FlatDictionaryFrame(WSSequence wsseq) {
        this.wsseq = wsseq;
        this.title = "Dictionary viewer: "+ wsseq.getName();
        this.setTitle(title);
        this.nelsa = wsseq.getNELSA();
        
        initComponents();
        
        this.dictPanel = new FlatDictionaryPanel(nelsa);
        this.centerPanel.add(BorderLayout.CENTER, dictPanel);
        
        
        this.centerPanel.addMouseWheelListener(
                new MouseWheelListener() {

            @Override
            public void mouseWheelMoved(MouseWheelEvent e) {
//                String message;
                int notches = e.getWheelRotation();
//                System.out.println("wheeltoration: "+ notches);
                if (notches < 0) {
//                    message = "Mouse wheel moved UP "
//                                 + -notches + " notch(es)" + newline;
                } else {
//                    message = "Mouse wheel moved DOWN "
//                                 + notches + " notch(es)" + newline;
                }
                if (e.getScrollType() == MouseWheelEvent.WHEEL_UNIT_SCROLL) {
                    if(e.getUnitsToScroll()<0){
                        prev(100);
                    }
                    else{
                        next(100);
                    }
//                    System.out.println("unitscroll: "+e.getUnitsToScroll());
//                    System.out.println("scrollammount: "+e.getScrollAmount());
                    
                    
//                    message += "    Scroll type: WHEEL_UNIT_SCROLL" + newline;
//                    message += "    Scroll amount: " + e.getScrollAmount()
//                            + " unit increments per notch" + newline;
//                    message += "    Units to scroll: " + e.getUnitsToScroll()
//                            + " unit increments" + newline;
//                    message += "    Vertical unit increment: "
//                        + scrollPane.getVerticalScrollBar().getUnitIncrement(1)
//                        + " pixels" + newline;
                } else { //scroll type == MouseWheelEvent.WHEEL_BLOCK_SCROLL
//                    message += "    Scroll type: WHEEL_BLOCK_SCROLL" + newline;
//                    message += "    Vertical block increment: "
//                        + scrollPane.getVerticalScrollBar().getBlockIncrement(1)
//                        + " pixels" + newline;
                }
             }
        });
        
        
        this.centerPanel.addMouseListener(
                new MouseAdapter() {

                    @Override
                    public void mouseClicked(MouseEvent e) {
//                        System.out.println(e.getX()+","+e.getY());
                        pickedKmerField.setText(dictPanel.getPickedKmer(e.getY()));
                        pickedKmerField.setCaretPosition(0);
                        requestFocus();
                    }

                });
        
        
//        this.dictPanel.addKeyListener(new KeyAdapter() {
//            
//            @Override
//            public void keyPressed(KeyEvent e){
//                System.out.println(e.getKeyCode());
//            }
//            
//             @Override
//            public void keyTyped(KeyEvent e){
//                System.out.println(e.getKeyCode());
//            }
//        });
//        this.dictPanel.setFocusable(true);
        
        InputMap inputMap; 
        ActionMap actionMap;
        AbstractAction pageup_action;
        AbstractAction pagedown_action;

        inputMap  = this.centerPanel.getInputMap(JComponent.WHEN_IN_FOCUSED_WINDOW);
        actionMap = this.centerPanel.getActionMap();

        inputMap.put(KeyStroke.getKeyStroke(KeyEvent.VK_PAGE_UP, 0), "dictViewPageUp");
        actionMap.put("dictViewPageUp", new AbstractAction(){
           @Override
           public void actionPerformed(ActionEvent e){pageUp();}
        });
        inputMap.put(KeyStroke.getKeyStroke(KeyEvent.VK_PAGE_DOWN, 0), "dictViewPageDown");
        actionMap.put("dictViewPageDown", new AbstractAction(){
           @Override
           public void actionPerformed(ActionEvent e){pageDown();}
        });
        inputMap.put(KeyStroke.getKeyStroke(KeyEvent.VK_UP, 0), "dictViewLineUp");
        actionMap.put("dictViewLineUp", new AbstractAction(){
           @Override
           public void actionPerformed(ActionEvent e){lineUp();}
        });
        inputMap.put(KeyStroke.getKeyStroke(KeyEvent.VK_DOWN, 0), "dictViewLineDown");
        actionMap.put("dictViewLineDown", new AbstractAction(){
           @Override
           public void actionPerformed(ActionEvent e){lineDown();}
        });
        
    }

    /**
     * This method is called from within the constructor to initialize the form.
     * WARNING: Do NOT modify this code. The content of this method is always
     * regenerated by the Form Editor.
     */
    @SuppressWarnings("unchecked")
    // <editor-fold defaultstate="collapsed" desc="Generated Code">//GEN-BEGIN:initComponents
    private void initComponents() {

        jScrollPane3 = new javax.swing.JScrollPane();
        jPanel1 = new javax.swing.JPanel();
        jLabel2 = new javax.swing.JLabel();
        kField = new javax.swing.JTextField();
        startButton = new javax.swing.JButton();
        jSeparator1 = new javax.swing.JSeparator();
        prevButton = new javax.swing.JButton();
        nextButton = new javax.swing.JButton();
        stepField = new javax.swing.JTextField();
        jSeparator2 = new javax.swing.JSeparator();
        jButton4 = new javax.swing.JButton();
        jButton5 = new javax.swing.JButton();
        jSeparator3 = new javax.swing.JSeparator();
        jLabel3 = new javax.swing.JLabel();
        pointWidthField = new javax.swing.JTextField();
        jLabel4 = new javax.swing.JLabel();
        pointHeightField = new javax.swing.JTextField();
        jButton6 = new javax.swing.JButton();
        jTextField5 = new javax.swing.JTextField();
        heightByMultiCheck = new javax.swing.JCheckBox();
        jSeparator4 = new javax.swing.JSeparator();
        jButton3 = new javax.swing.JButton();
        scoreTolleranceField = new javax.swing.JTextField();
        jLabel1 = new javax.swing.JLabel();
        rowTollPercField = new javax.swing.JTextField();
        jLabel5 = new javax.swing.JLabel();
        rowTollMinField = new javax.swing.JTextField();
        remNsField = new javax.swing.JCheckBox();
        jButton7 = new javax.swing.JButton();
        lcpCheck = new javax.swing.JCheckBox();
        multsBarCheckBox = new javax.swing.JCheckBox();
        jLabel9 = new javax.swing.JLabel();
        aComboBox = new javax.swing.JComboBox();
        jLabel10 = new javax.swing.JLabel();
        cComboBox = new javax.swing.JComboBox();
        jLabel11 = new javax.swing.JLabel();
        jLabel12 = new javax.swing.JLabel();
        gComboBox = new javax.swing.JComboBox();
        tComboBox = new javax.swing.JComboBox();
        jLabel13 = new javax.swing.JLabel();
        jButton8 = new javax.swing.JButton();
        jSeparator5 = new javax.swing.JSeparator();
        jPanel2 = new javax.swing.JPanel();
        jPanel3 = new javax.swing.JPanel();
        pickedKmerField = new javax.swing.JTextField();
        goButton = new javax.swing.JButton();
        jButton1 = new javax.swing.JButton();
        jButton2 = new javax.swing.JButton();
        jScrollPane1 = new javax.swing.JScrollPane();
        centerPanel = new javax.swing.JPanel();

        setDefaultCloseOperation(javax.swing.WindowConstants.EXIT_ON_CLOSE);
        setMinimumSize(new java.awt.Dimension(200, 100));

        jScrollPane3.setHorizontalScrollBarPolicy(javax.swing.ScrollPaneConstants.HORIZONTAL_SCROLLBAR_NEVER);
        jScrollPane3.setPreferredSize(new java.awt.Dimension(220, 602));

        jPanel1.setMinimumSize(new java.awt.Dimension(200, 100));
        jPanel1.setPreferredSize(new java.awt.Dimension(200, 800));

        jLabel2.setText("k");

        kField.setText("100");

        startButton.setText("Start");
        startButton.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                startButtonActionPerformed(evt);
            }
        });

        prevButton.setText("<");
        prevButton.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                prevButtonActionPerformed(evt);
            }
        });

        nextButton.setText(">");
        nextButton.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                nextButtonActionPerformed(evt);
            }
        });

        stepField.setText("10000");

        jButton4.setText("||");
        jButton4.setEnabled(false);

        jButton5.setText(">");
        jButton5.setEnabled(false);

        jLabel3.setText("w");

        pointWidthField.setText("4");

        jLabel4.setText("h");

        pointHeightField.setText("1");

        jButton6.setText("go");
        jButton6.setEnabled(false);

        jTextField5.setEnabled(false);

        heightByMultiCheck.setText("multiplicity");

        jButton3.setText("MSA");
        jButton3.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                jButton3ActionPerformed(evt);
            }
        });

        scoreTolleranceField.setText("0.30");

        jLabel1.setText("Score tollerance");

        rowTollPercField.setEditable(false);
        rowTollPercField.setText("0.05");
        rowTollPercField.setEnabled(false);

        jLabel5.setText("Row tollerance");

        rowTollMinField.setText("10");

        remNsField.setText("Remove Gaps");
        remNsField.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                remNsFieldActionPerformed(evt);
            }
        });

        jButton7.setText("SC");
        jButton7.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                jButton7ActionPerformed(evt);
            }
        });

        lcpCheck.setText("LCP");

        multsBarCheckBox.setText("multiplicity bars");

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

        jLabel12.setText("G");

        gComboBox.setModel(getColorListModel());
        gComboBox.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                gComboBoxActionPerformed(evt);
            }
        });

        tComboBox.setModel(getColorListModel());
        tComboBox.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                tComboBoxActionPerformed(evt);
            }
        });

        jLabel13.setText("T");

        jButton8.setText("Apply");
        jButton8.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                jButton8ActionPerformed(evt);
            }
        });

        javax.swing.GroupLayout jPanel1Layout = new javax.swing.GroupLayout(jPanel1);
        jPanel1.setLayout(jPanel1Layout);
        jPanel1Layout.setHorizontalGroup(
            jPanel1Layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGroup(jPanel1Layout.createSequentialGroup()
                .addContainerGap()
                .addGroup(jPanel1Layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
                    .addComponent(jSeparator1, javax.swing.GroupLayout.Alignment.TRAILING)
                    .addGroup(javax.swing.GroupLayout.Alignment.TRAILING, jPanel1Layout.createSequentialGroup()
                        .addComponent(jLabel2)
                        .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                        .addComponent(kField, javax.swing.GroupLayout.PREFERRED_SIZE, 60, javax.swing.GroupLayout.PREFERRED_SIZE)
                        .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED, javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE)
                        .addComponent(startButton))
                    .addGroup(jPanel1Layout.createSequentialGroup()
                        .addComponent(prevButton)
                        .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                        .addComponent(stepField)
                        .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                        .addComponent(nextButton))
                    .addComponent(jSeparator2)
                    .addComponent(jSeparator3)
                    .addGroup(jPanel1Layout.createSequentialGroup()
                        .addComponent(jLabel3)
                        .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                        .addComponent(pointWidthField, javax.swing.GroupLayout.PREFERRED_SIZE, 50, javax.swing.GroupLayout.PREFERRED_SIZE)
                        .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED, 32, Short.MAX_VALUE)
                        .addComponent(jLabel4)
                        .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                        .addComponent(pointHeightField, javax.swing.GroupLayout.PREFERRED_SIZE, 50, javax.swing.GroupLayout.PREFERRED_SIZE))
                    .addGroup(jPanel1Layout.createSequentialGroup()
                        .addComponent(jButton6)
                        .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                        .addComponent(jTextField5))
                    .addComponent(jSeparator4)
                    .addGroup(jPanel1Layout.createSequentialGroup()
                        .addComponent(jLabel10)
                        .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                        .addComponent(aComboBox, 0, javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE))
                    .addGroup(jPanel1Layout.createSequentialGroup()
                        .addComponent(jLabel11)
                        .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                        .addComponent(cComboBox, 0, javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE))
                    .addGroup(jPanel1Layout.createSequentialGroup()
                        .addComponent(jLabel12)
                        .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                        .addComponent(gComboBox, 0, javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE))
                    .addGroup(jPanel1Layout.createSequentialGroup()
                        .addComponent(jLabel13)
                        .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                        .addComponent(tComboBox, 0, javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE))
                    .addGroup(javax.swing.GroupLayout.Alignment.TRAILING, jPanel1Layout.createSequentialGroup()
                        .addGap(0, 0, Short.MAX_VALUE)
                        .addComponent(jButton8))
                    .addComponent(jSeparator5, javax.swing.GroupLayout.Alignment.TRAILING, javax.swing.GroupLayout.DEFAULT_SIZE, 176, Short.MAX_VALUE)
                    .addGroup(jPanel1Layout.createSequentialGroup()
                        .addGroup(jPanel1Layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
                            .addGroup(jPanel1Layout.createSequentialGroup()
                                .addComponent(jButton4)
                                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                                .addComponent(jButton5))
                            .addComponent(multsBarCheckBox)
                            .addComponent(heightByMultiCheck)
                            .addComponent(lcpCheck)
                            .addComponent(jLabel9)
                            .addComponent(jLabel1)
                            .addComponent(scoreTolleranceField, javax.swing.GroupLayout.PREFERRED_SIZE, 45, javax.swing.GroupLayout.PREFERRED_SIZE)
                            .addGroup(jPanel1Layout.createParallelGroup(javax.swing.GroupLayout.Alignment.TRAILING, false)
                                .addGroup(javax.swing.GroupLayout.Alignment.LEADING, jPanel1Layout.createSequentialGroup()
                                    .addComponent(rowTollPercField, javax.swing.GroupLayout.PREFERRED_SIZE, 50, javax.swing.GroupLayout.PREFERRED_SIZE)
                                    .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                                    .addComponent(rowTollMinField))
                                .addComponent(jLabel5, javax.swing.GroupLayout.Alignment.LEADING))
                            .addGroup(jPanel1Layout.createSequentialGroup()
                                .addComponent(jButton7)
                                .addGap(18, 18, 18)
                                .addComponent(jButton3))
                            .addComponent(remNsField))
                        .addGap(0, 0, Short.MAX_VALUE)))
                .addContainerGap())
        );
        jPanel1Layout.setVerticalGroup(
            jPanel1Layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGroup(jPanel1Layout.createSequentialGroup()
                .addContainerGap()
                .addGroup(jPanel1Layout.createParallelGroup(javax.swing.GroupLayout.Alignment.BASELINE)
                    .addComponent(jLabel2)
                    .addComponent(kField, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE)
                    .addComponent(startButton))
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                .addComponent(jSeparator1, javax.swing.GroupLayout.PREFERRED_SIZE, 10, javax.swing.GroupLayout.PREFERRED_SIZE)
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                .addGroup(jPanel1Layout.createParallelGroup(javax.swing.GroupLayout.Alignment.BASELINE)
                    .addComponent(prevButton)
                    .addComponent(nextButton)
                    .addComponent(stepField, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE))
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                .addGroup(jPanel1Layout.createParallelGroup(javax.swing.GroupLayout.Alignment.BASELINE)
                    .addComponent(jButton6)
                    .addComponent(jTextField5, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE))
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                .addComponent(jSeparator2, javax.swing.GroupLayout.PREFERRED_SIZE, 10, javax.swing.GroupLayout.PREFERRED_SIZE)
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                .addGroup(jPanel1Layout.createParallelGroup(javax.swing.GroupLayout.Alignment.BASELINE)
                    .addComponent(jButton4)
                    .addComponent(jButton5))
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                .addComponent(jSeparator3, javax.swing.GroupLayout.PREFERRED_SIZE, 10, javax.swing.GroupLayout.PREFERRED_SIZE)
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                .addGroup(jPanel1Layout.createParallelGroup(javax.swing.GroupLayout.Alignment.BASELINE)
                    .addComponent(jLabel3)
                    .addComponent(pointWidthField, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE)
                    .addComponent(pointHeightField, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE)
                    .addComponent(jLabel4))
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                .addComponent(multsBarCheckBox)
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                .addComponent(heightByMultiCheck)
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                .addComponent(lcpCheck)
                .addGap(18, 18, 18)
                .addComponent(jSeparator4, javax.swing.GroupLayout.PREFERRED_SIZE, 10, javax.swing.GroupLayout.PREFERRED_SIZE)
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                .addComponent(jLabel9)
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                .addGroup(jPanel1Layout.createParallelGroup(javax.swing.GroupLayout.Alignment.BASELINE)
                    .addComponent(aComboBox, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE)
                    .addComponent(jLabel10))
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                .addGroup(jPanel1Layout.createParallelGroup(javax.swing.GroupLayout.Alignment.BASELINE)
                    .addComponent(cComboBox, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE)
                    .addComponent(jLabel11))
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                .addGroup(jPanel1Layout.createParallelGroup(javax.swing.GroupLayout.Alignment.BASELINE)
                    .addComponent(gComboBox, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE)
                    .addComponent(jLabel12))
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                .addGroup(jPanel1Layout.createParallelGroup(javax.swing.GroupLayout.Alignment.BASELINE)
                    .addComponent(tComboBox, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE)
                    .addComponent(jLabel13))
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                .addComponent(jButton8)
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                .addComponent(jSeparator5, javax.swing.GroupLayout.PREFERRED_SIZE, 10, javax.swing.GroupLayout.PREFERRED_SIZE)
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                .addComponent(jLabel1)
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                .addComponent(scoreTolleranceField, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE)
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.UNRELATED)
                .addComponent(jLabel5)
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                .addGroup(jPanel1Layout.createParallelGroup(javax.swing.GroupLayout.Alignment.BASELINE)
                    .addComponent(rowTollPercField, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE)
                    .addComponent(rowTollMinField, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE))
                .addGap(18, 18, 18)
                .addGroup(jPanel1Layout.createParallelGroup(javax.swing.GroupLayout.Alignment.BASELINE)
                    .addComponent(jButton7)
                    .addComponent(jButton3, javax.swing.GroupLayout.PREFERRED_SIZE, 29, javax.swing.GroupLayout.PREFERRED_SIZE))
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                .addComponent(remNsField)
                .addContainerGap(54, Short.MAX_VALUE))
        );

        jScrollPane3.setViewportView(jPanel1);

        getContentPane().add(jScrollPane3, java.awt.BorderLayout.WEST);

        jPanel2.setLayout(new java.awt.BorderLayout());

        goButton.setText("go");
        goButton.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                goButtonActionPerformed(evt);
            }
        });

        jButton1.setText("<");
        jButton1.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                jButton1ActionPerformed(evt);
            }
        });

        jButton2.setText("<+");
        jButton2.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                jButton2ActionPerformed(evt);
            }
        });

        javax.swing.GroupLayout jPanel3Layout = new javax.swing.GroupLayout(jPanel3);
        jPanel3.setLayout(jPanel3Layout);
        jPanel3Layout.setHorizontalGroup(
            jPanel3Layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGroup(javax.swing.GroupLayout.Alignment.TRAILING, jPanel3Layout.createSequentialGroup()
                .addComponent(goButton)
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                .addComponent(pickedKmerField, javax.swing.GroupLayout.DEFAULT_SIZE, 488, Short.MAX_VALUE)
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                .addComponent(jButton1)
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                .addComponent(jButton2))
        );
        jPanel3Layout.setVerticalGroup(
            jPanel3Layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGroup(javax.swing.GroupLayout.Alignment.TRAILING, jPanel3Layout.createSequentialGroup()
                .addGap(0, 0, Short.MAX_VALUE)
                .addGroup(jPanel3Layout.createParallelGroup(javax.swing.GroupLayout.Alignment.BASELINE)
                    .addComponent(pickedKmerField, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE)
                    .addComponent(goButton)
                    .addComponent(jButton1)
                    .addComponent(jButton2)))
        );

        jPanel2.add(jPanel3, java.awt.BorderLayout.PAGE_END);

        centerPanel.setBackground(java.awt.Color.white);
        centerPanel.setLayout(new java.awt.BorderLayout());
        jScrollPane1.setViewportView(centerPanel);

        jPanel2.add(jScrollPane1, java.awt.BorderLayout.CENTER);

        getContentPane().add(jPanel2, java.awt.BorderLayout.CENTER);

        pack();
    }// </editor-fold>//GEN-END:initComponents

    private void startButtonActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_startButtonActionPerformed
        start();
    }//GEN-LAST:event_startButtonActionPerformed

    private void prevButtonActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_prevButtonActionPerformed
        prev();
    }//GEN-LAST:event_prevButtonActionPerformed

    private void nextButtonActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_nextButtonActionPerformed
        next();
    }//GEN-LAST:event_nextButtonActionPerformed

    private void goButtonActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_goButtonActionPerformed
        try{
            
            this.setEnabled(false);
            IELSAIterator it = nelsa.find(new B3LLSequence(this.pickedKmerField.getText().trim()));
            if(it != null){
                for(int i=0; i<100 && it.hasPrev(); i++)
                    it.prev();
                this.first_it = it;
                this.kField.setText("" + it.k());
                this.pickedKmerField.setBackground(Color.white);
                int pointWidth = Integer.parseInt(this.pointWidthField.getText());
                int pointHeight = Integer.parseInt(this.pointHeightField.getText());
                this.dictPanel.setData(first_it.clone(), null, pointWidth, pointHeight, this.heightByMultiCheck.isSelected(), this.lcpCheck.isSelected());
                this.dictPanel.repaint();
            }
            else{
                 this.pickedKmerField.setBackground(Color.red);
            }
            
        }catch(Exception e){
        }
        this.setEnabled(true);
    }//GEN-LAST:event_goButtonActionPerformed

    private void jButton1ActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_jButton1ActionPerformed
        if(this.pickedKmerField.getText().trim().length() > 1){
            try{
            this.setEnabled(false);
            String sk = this.pickedKmerField.getText().trim();
            IELSAIterator it = nelsa.find(new B3LLSequence(sk.substring(1, sk.length())));
            if(it != null){
                for(int i=0; i<100 && it.hasPrev(); i++)
                    it.prev();
                this.first_it = it;
                this.kField.setText("" + it.k());
                this.pickedKmerField.setBackground(Color.white);
                int pointWidth = Integer.parseInt(this.pointWidthField.getText());
                int pointHeight = Integer.parseInt(this.pointHeightField.getText());
                this.dictPanel.setData(first_it.clone(), null, pointWidth, pointHeight, this.heightByMultiCheck.isSelected(), this.lcpCheck.isSelected());
                this.dictPanel.repaint();
                this.pickedKmerField.setText(sk.substring(1, sk.length()));
            }
            else{
                 this.pickedKmerField.setBackground(Color.red);
            }
            }catch(Exception e){
            }
            this.setEnabled(true);
        }
    }//GEN-LAST:event_jButton1ActionPerformed

    private void jButton2ActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_jButton2ActionPerformed
        // TODO add your handling code here:
    }//GEN-LAST:event_jButton2ActionPerformed

    private void jButton3ActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_jButton3ActionPerformed
        //recognize();
        rec(VBlockFrame.RecognitionType.MSA);
    }//GEN-LAST:event_jButton3ActionPerformed

    private void remNsFieldActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_remNsFieldActionPerformed
        // TODO add your handling code here:
    }//GEN-LAST:event_remNsFieldActionPerformed

    private void jButton7ActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_jButton7ActionPerformed
        rec(VBlockFrame.RecognitionType.SC);
    }//GEN-LAST:event_jButton7ActionPerformed

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

    private void jButton8ActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_jButton8ActionPerformed
        this.dictPanel.setColorPalette(
            B3NucleotidePaletteFactory.createPalette(
                ((ColorItem)aComboBox.getSelectedItem()).color,
                ((ColorItem)cComboBox.getSelectedItem()).color,
                ((ColorItem)gComboBox.getSelectedItem()).color,
                ((ColorItem)tComboBox.getSelectedItem()).color,
                Color.black
            )
        );
        dictPanel.repaint();
    }//GEN-LAST:event_jButton8ActionPerformed

    // Variables declaration - do not modify//GEN-BEGIN:variables
    private javax.swing.JComboBox aComboBox;
    private javax.swing.JComboBox cComboBox;
    private javax.swing.JPanel centerPanel;
    private javax.swing.JComboBox gComboBox;
    private javax.swing.JButton goButton;
    private javax.swing.JCheckBox heightByMultiCheck;
    private javax.swing.JButton jButton1;
    private javax.swing.JButton jButton2;
    private javax.swing.JButton jButton3;
    private javax.swing.JButton jButton4;
    private javax.swing.JButton jButton5;
    private javax.swing.JButton jButton6;
    private javax.swing.JButton jButton7;
    private javax.swing.JButton jButton8;
    private javax.swing.JLabel jLabel1;
    private javax.swing.JLabel jLabel10;
    private javax.swing.JLabel jLabel11;
    private javax.swing.JLabel jLabel12;
    private javax.swing.JLabel jLabel13;
    private javax.swing.JLabel jLabel2;
    private javax.swing.JLabel jLabel3;
    private javax.swing.JLabel jLabel4;
    private javax.swing.JLabel jLabel5;
    private javax.swing.JLabel jLabel9;
    private javax.swing.JPanel jPanel1;
    private javax.swing.JPanel jPanel2;
    private javax.swing.JPanel jPanel3;
    private javax.swing.JScrollPane jScrollPane1;
    private javax.swing.JScrollPane jScrollPane3;
    private javax.swing.JSeparator jSeparator1;
    private javax.swing.JSeparator jSeparator2;
    private javax.swing.JSeparator jSeparator3;
    private javax.swing.JSeparator jSeparator4;
    private javax.swing.JSeparator jSeparator5;
    private javax.swing.JTextField jTextField5;
    private javax.swing.JTextField kField;
    private javax.swing.JCheckBox lcpCheck;
    private javax.swing.JCheckBox multsBarCheckBox;
    private javax.swing.JButton nextButton;
    private javax.swing.JTextField pickedKmerField;
    private javax.swing.JTextField pointHeightField;
    private javax.swing.JTextField pointWidthField;
    private javax.swing.JButton prevButton;
    private javax.swing.JCheckBox remNsField;
    private javax.swing.JTextField rowTollMinField;
    private javax.swing.JTextField rowTollPercField;
    private javax.swing.JTextField scoreTolleranceField;
    private javax.swing.JButton startButton;
    private javax.swing.JTextField stepField;
    private javax.swing.JComboBox tComboBox;
    // End of variables declaration//GEN-END:variables

    private void changeSettings(){
        try{            
        }catch(Exception e){
            
        }
    }
    
    private void start(){
        this.setEnabled(false);
        try{
            
            int k = Integer.parseInt(this.kField.getText());
            this.first_it = nelsa.begin(k);
            int pointWidth = Integer.parseInt(this.pointWidthField.getText());
            int pointHeight = Integer.parseInt(this.pointHeightField.getText());
            
            if(this.first_it.next()){
                this.dictPanel.setData(first_it.clone(), null, pointWidth, pointHeight, this.heightByMultiCheck.isSelected(), this.lcpCheck.isSelected());
                this.dictPanel.printMultsBar = this.multsBarCheckBox.isSelected();
                this.dictPanel.repaint();
            }
                    
                    
        }catch(Exception e){
            e.printStackTrace();
        }
        this.setEnabled(true);
    }
    
    
    private void next(int step){
        this.setEnabled(false);
        try{
            int pointWidth = Integer.parseInt(this.pointWidthField.getText());
            int pointHeight = Integer.parseInt(this.pointHeightField.getText());
            boolean byMulti = this.heightByMultiCheck.isSelected();
            
            if(step > 0 && this.first_it.good()){
                for(int s = 0; s<step && this.first_it.hasNext(); s++){
                    if(byMulti)
                        s += this.first_it.multiplicity() - 1;
                    if(!this.first_it.next())
                        break;
                }
                
                this.dictPanel.setData(first_it.clone(), null, pointWidth, pointHeight, this.heightByMultiCheck.isSelected(), this.lcpCheck.isSelected());
                this.dictPanel.printMultsBar = this.multsBarCheckBox.isSelected();
                this.dictPanel.repaint();
            }
            
        }catch(Exception e){
            e.printStackTrace();
        }
        this.setEnabled(true);
    }
    
    private void next(){
        this.setEnabled(false);
        try{
            int step = Integer.parseInt(this.stepField.getText());
            next(step);
        }catch(Exception e){
            e.printStackTrace();
        }
        this.setEnabled(true);
    }
    
    
    private void prev(int step){
        this.setEnabled(false);
        try{
            int pointWidth = Integer.parseInt(this.pointWidthField.getText());
            int pointHeight = Integer.parseInt(this.pointHeightField.getText());
            boolean byMulti = this.heightByMultiCheck.isSelected();
            
            if(step > 0 && this.first_it.good()){
                for(int s = 0; s<step && this.first_it.hasPrev(); s++){
                    if(byMulti)
                        s += this.first_it.multiplicity() - 1;
                    if(!this.first_it.prev())
                        break;
                }
                
                this.dictPanel.setData(first_it.clone(), null, pointWidth, pointHeight, this.heightByMultiCheck.isSelected(), this.lcpCheck.isSelected());
                this.dictPanel.printMultsBar = this.multsBarCheckBox.isSelected();
                this.dictPanel.repaint();
            }
        }catch(Exception e){
            e.printStackTrace();
        }
        this.setEnabled(true);
    }
    
    private void prev(){
        this.setEnabled(false);
        try{
            int step = Integer.parseInt(this.stepField.getText());
            prev(step);
        }catch(Exception e){
            e.printStackTrace();
        }
        this.setEnabled(true);
    }
    
    
    private void pageUp(){
        if(dictPanel != null){
            try{
              IELSAIterator it = dictPanel.pageUpIterator();
              if(it != null){
                  dictPanel.setData(
                          it,
                          null, 
                          Integer.parseInt(pointWidthField.getText()), 
                          Integer.parseInt(pointHeightField.getText()), 
                          heightByMultiCheck.isSelected()
                          , this.lcpCheck.isSelected());
                  this.dictPanel.printMultsBar = this.multsBarCheckBox.isSelected();
                  dictPanel.repaint();
                  first_it = it;
              }
            }catch(Exception ex){    
          }
        }
    }
    
    private void pageDown(){
        if(dictPanel != null){
            try{
              IELSAIterator it = dictPanel.pageDownIterator();
              if(it != null){
                  dictPanel.setData(
                          it,
                          null, 
                          Integer.parseInt(pointWidthField.getText()), 
                          Integer.parseInt(pointHeightField.getText()), 
                          heightByMultiCheck.isSelected()
                          , this.lcpCheck.isSelected());
                  this.dictPanel.printMultsBar = this.multsBarCheckBox.isSelected();
                  dictPanel.repaint();
                  first_it = it;
              }
            }catch(Exception ex){    
          }
        }
    }
    
    private void lineUp(){
        if(dictPanel != null){
            try{
              IELSAIterator it = first_it.clone();
              if(it.prev()){
                  dictPanel.setData(
                          it,
                          null, 
                          Integer.parseInt(pointWidthField.getText()), 
                          Integer.parseInt(pointHeightField.getText()), 
                          heightByMultiCheck.isSelected()
                          , this.lcpCheck.isSelected());
                  this.dictPanel.printMultsBar = this.multsBarCheckBox.isSelected();
                  dictPanel.repaint();
                  first_it = it;
              }
            }catch(Exception ex){    
          }
        }
    }
    
    private void lineDown(){
        if(dictPanel != null){
            try{
              IELSAIterator it = first_it.clone();
              if(it.next()){
                  dictPanel.setData(
                          it,
                          null, 
                          Integer.parseInt(pointWidthField.getText()), 
                          Integer.parseInt(pointHeightField.getText()), 
                          heightByMultiCheck.isSelected()
                          , this.lcpCheck.isSelected());
                  this.dictPanel.printMultsBar = this.multsBarCheckBox.isSelected();
                  dictPanel.repaint();
                  first_it = it;
              }
            }catch(Exception ex){    
          }
        }
    }
    
    
    private int max(int a, int b){
        if(a>=b)
            return a;
        return b;
    }
    
    
    private void rec(final VBlockFrame.RecognitionType rtype){
        try{
            final double scoreTollerance = Double.parseDouble(this.scoreTolleranceField.getText());
            final double rowTollPerc = Double.parseDouble(this.rowTollPercField.getText());
            final double rowTollMin = Double.parseDouble(this.rowTollMinField.getText());
            
            final boolean remNs = this.remNsField.isSelected();
            
            final Stack<String> alignSequences = new Stack();
            final Stack<Integer> alignScores = new Stack();
            
            final IELSAIterator pickedIt = nelsa.find(new B3LLSequence(this.pickedKmerField.getText().trim()));
            
            if(pickedIt != null){
                java.awt.EventQueue.invokeLater(new Runnable() {
                    public void run() {
                        VBlockFrame frame = new VBlockFrame(title);
                        frame.setDefaultCloseOperation(DISPOSE_ON_CLOSE);
                        frame.setVisible(true);

                        frame.recognize(
                                scoreTollerance, 
                                rowTollPerc, 
                                rowTollMin, 
                                remNs, 
                                pickedIt, 
                                rtype);

                    }
                });
            }
            
        }catch(Exception e){
            
        }
    }
    
    private void recognize(){
        this.setEnabled(false);
        try{
            
            final double scoreTollerance = Double.parseDouble(this.scoreTolleranceField.getText());
            final double rowTollPerc = Double.parseDouble(this.rowTollPercField.getText());
            final double rowTollMin = Double.parseDouble(this.rowTollMinField.getText());
            
            final boolean remNs = this.remNsField.isSelected();
            
            final Stack<String> alignSequences = new Stack();
            final Stack<Integer> alignScores = new Stack();
            
            IELSAIterator pickedIt = nelsa.find(new B3LLSequence(this.pickedKmerField.getText().trim()));
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
            
            
//            int maxLength = Integer.MIN_VALUE;
//            int maxLengthIndex = -1;
//            String maxS;
//            
//            
//            for(int i=0; i< alignSequences.size(); i++){
//                maxS = alignSequences.get(i);
//                if(maxS.length() > maxLength){
//                    maxLength =  maxS.length();
//                    maxLengthIndex = i;
//                }
//            }
            
//            if(maxLengthIndex >= 0){
//                maxS = alignSequences.get(maxLengthIndex);
//                for(int i=0; i<alignSequences.size(); i++){
//                    if(i!= maxLengthIndex){
//                        s = alignSequences.get(i);
//                        StringSmithWaterman align = new StringSmithWaterman(maxS, s);
//                        align.process();
//                        
//                        as1 = align.getAlignedS1();
//                        as2 = align.getAlignedS2();
//                        
//                        if(as1.length() != maxLength){
//                            System.out.println("WARNING");
//                            System.out.println(maxS);
//                            System.out.println(as1);
//                        }
//                        
//                        alignSequences.set(i, as2);
//                    }
//                }
//            }
            
            
            java.awt.EventQueue.invokeLater(new Runnable() {
                public void run() {
                    VBlockFrame frame = new VBlockFrame(title);
//                    String[] seqs = new String[alignSequences.size()];
//                    alignSequences.toArray(seqs);
                    Integer[] scores = new Integer[alignScores.size()];
                    alignScores.toArray(scores);

                    frame.fillTable(aligned, scores);
                    frame.setDefaultCloseOperation(DISPOSE_ON_CLOSE);
                    frame.setVisible(true);
                }
            });
            
        }catch(Exception e){
        }
        this.setEnabled(true);
    }
    
    
    
    
    
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
}
