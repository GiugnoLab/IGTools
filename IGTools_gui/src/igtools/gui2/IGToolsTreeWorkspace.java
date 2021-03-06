/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package igtools.gui2;

import java.awt.event.MouseAdapter;
import java.awt.event.MouseEvent;
import java.awt.event.MouseListener;
import java.io.File;
import javax.swing.JButton;
import javax.swing.tree.DefaultMutableTreeNode;
import javax.swing.tree.DefaultTreeModel;
import javax.swing.tree.TreeNode;
import javax.swing.tree.TreePath;

/**
 *
 * @author vbonnici
 */
public class IGToolsTreeWorkspace extends javax.swing.JFrame {
    
    private static String ACTION_ADD_SEQUENCE = "add_sequence"; 

    private WSConfiguration config = null;
    private DefaultMutableTreeNode wsRootNode = new DefaultMutableTreeNode("Sequences");
    private DefaultTreeModel wsTreeModel = new DefaultTreeModel(wsRootNode); 
    
    private DefaultMutableTreeNode wsAddSequenceTreeNode = new DefaultMutableTreeNode("Add...");
    
    /**
     * Creates new form IGToolsWorkspace
     */
    public IGToolsTreeWorkspace(WSConfiguration config) {
        this.config = config;
        wsTreeModel.insertNodeInto(new DefaultMutableTreeNode(wsAddSequenceTreeNode), wsRootNode, 0);
        initComponents();
        jTree2.setModel(wsTreeModel);
        
        
        MouseListener ml = new MouseAdapter() {
            public void mousePressed(MouseEvent e) {
                int selRow = jTree2.getRowForLocation(e.getX(), e.getY());
                TreePath selPath = jTree2.getPathForLocation(e.getX(), e.getY());
                if(selPath != null){
                     System.out.println(selPath.getLastPathComponent());
                }
                if(selRow != -1) {
                    if(e.getClickCount() == 1) {
                       
                    }
                    else if(e.getClickCount() == 2) {
                    }
                }
            }
        };
        jTree2.addMouseListener(ml);
    }
    
    
    
    private void initSeqsTree(){
        
    }
    
    
    

    /**
     * This method is called from within the constructor to initialize the form.
     * WARNING: Do NOT modify this code. The content of this method is always
     * regenerated by the Form Editor.
     */
    @SuppressWarnings("unchecked")
    // <editor-fold defaultstate="collapsed" desc="Generated Code">//GEN-BEGIN:initComponents
    private void initComponents() {

        jScrollPane2 = new javax.swing.JScrollPane();
        jTree2 = new javax.swing.JTree();
        jMenuBar1 = new javax.swing.JMenuBar();
        jMenu1 = new javax.swing.JMenu();
        jMenu2 = new javax.swing.JMenu();
        jMenu3 = new javax.swing.JMenu();

        setDefaultCloseOperation(javax.swing.WindowConstants.EXIT_ON_CLOSE);
        setTitle("InfoGenomic Tools - Workspace");

        javax.swing.tree.DefaultMutableTreeNode treeNode1 = new javax.swing.tree.DefaultMutableTreeNode("Loading...");
        jTree2.setModel(new javax.swing.tree.DefaultTreeModel(treeNode1));
        jScrollPane2.setViewportView(jTree2);

        jMenu1.setText("File");
        jMenuBar1.add(jMenu1);

        jMenu2.setText("Edit");
        jMenuBar1.add(jMenu2);

        jMenu3.setText("Tools");
        jMenuBar1.add(jMenu3);

        setJMenuBar(jMenuBar1);

        javax.swing.GroupLayout layout = new javax.swing.GroupLayout(getContentPane());
        getContentPane().setLayout(layout);
        layout.setHorizontalGroup(
            layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addComponent(jScrollPane2, javax.swing.GroupLayout.DEFAULT_SIZE, 671, Short.MAX_VALUE)
        );
        layout.setVerticalGroup(
            layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addComponent(jScrollPane2, javax.swing.GroupLayout.DEFAULT_SIZE, 361, Short.MAX_VALUE)
        );

        pack();
    }// </editor-fold>//GEN-END:initComponents


    // Variables declaration - do not modify//GEN-BEGIN:variables
    private javax.swing.JMenu jMenu1;
    private javax.swing.JMenu jMenu2;
    private javax.swing.JMenu jMenu3;
    private javax.swing.JMenuBar jMenuBar1;
    private javax.swing.JScrollPane jScrollPane2;
    private javax.swing.JTree jTree2;
    // End of variables declaration//GEN-END:variables
}
