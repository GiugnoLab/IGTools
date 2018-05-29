/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package igtools.gui2.dictionaries.vblocks;

import javax.swing.table.AbstractTableModel;

/**
 *
 * @author vbonnici
 */
public class VBlockTableModel extends AbstractTableModel{
    
    private String[] columnNames = {"k-mer","score"};
    private String[] kmers = new String[0];
    private Integer[] scores = new Integer[0];
    
    
    public VBlockTableModel(){
    }
    
    public void setData(String[] kmers, Integer[] scores){
        this.kmers = kmers;
        this.scores = scores;
    }
    
    

    @Override
    public int getRowCount() {
        return kmers.length;
    }

    @Override
    public int getColumnCount() {
        return columnNames.length;
    }

    @Override
    public Object getValueAt(int rowIndex, int columnIndex) {
        if(rowIndex < kmers.length){
            if(columnIndex == 0){
                return kmers[rowIndex];
            }
            else if(columnIndex == 1){
                return scores[rowIndex];
            }
            else{
                return null;
            }
        }
        return null;
    }
    
    public String getColumnName(int col) {
        return columnNames[col];
    }
    
    public Class getColumnClass(int c){
        if(c == 0){
            return String.class;
        }
        else if(c == 1){
            return Integer.class;
        }
        else{
            return null;
        }
    }
}
