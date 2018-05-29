/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package igtools.gui2.codistances;

import javax.swing.table.AbstractTableModel;

/**
 *
 * @author vbonnici
 */
public class DistalWordListTableModel extends AbstractTableModel{
    
    private String[] columnNames = {"k-mer", "count"};
    private String[] kmers = new String[0];
    private int[] counts = new int[0];
    
    
    public DistalWordListTableModel(){
    }
    
    public void setData(String[] kmers, int[] counts){
        this.kmers = kmers;
        this.counts = counts;
    }
    
    public String[] kmers(){
        return this.kmers;
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
                return new Integer(counts[rowIndex]);
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
