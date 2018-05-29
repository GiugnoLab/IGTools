/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package igtools.gui2;

import java.util.List;
import javax.swing.AbstractListModel;
import javax.swing.DefaultListModel;

/**
 *
 * @author vbonnici
 */
public class WSSeqListModel{

    private List<WSSequence> wssequences = null;
    private DefaultListModel<WSSequence> realModel = null;
    
    public WSSeqListModel(List<WSSequence> wssequences){
        this.wssequences = wssequences;
        this.realModel = new DefaultListModel<WSSequence>();
        for(WSSequence wseq : wssequences)
            realModel.addElement(wseq);
        
    }
    
    
    public DefaultListModel<WSSequence> getModel(){
        return realModel;
    }
    
    
    public void add(WSSequence wsseq){
        wssequences.add(wsseq);
        realModel.addElement(wsseq);
    }
    
    public void remove(int index){
        if(index < wssequences.size()){
            wssequences.remove(index);
            realModel.removeElementAt(index);
        }
    }
    
    public void remove(WSSequence wsseq){
      wssequences.remove(wsseq);
      realModel.removeElement(wsseq);
    }
    
    public WSSequence getElementAt(int index){
        if(index<wssequences.size()){
            return wssequences.get(index);
        }
        return null;
    }
}
