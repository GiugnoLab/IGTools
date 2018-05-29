/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package igtools.gui2;

import igtools.common.sequence.B3LLSequence;
import java.util.List;

import igtools.common.util.Timer;
import igtools.dictionaries.elsa.NELSA;

/**
 *
 * @author vbonnici
 */
public class ResourcesManager {
    
    private Timer timer = new Timer();
    
    private List<WSSequence> sequences = null;
    
    
    public ResourcesManager(List<WSSequence> sequences){
        this.sequences = sequences;
    }
    
    
    public void loadB3seq(int index, boolean forceReload){
        try{
            if(index < sequences.size()){
                WSSequence wsseq = sequences.get(index);
                if(wsseq.getB3seq()==null || forceReload){
                    timer.reset();
                    System.out.println("Loading sequence...");
                    B3LLSequence b3seq = B3LLSequence.load(wsseq.getB3seqFile());
                    wsseq.setB3seq(b3seq);
                    System.out.println("done "+timer.getElapsedSecs() +" sec.\n");
                    //System.out.println(b3seq+"\t"+wsseq.getB3seq());
                }
            }
        }catch(Exception ex){
            ex.printStackTrace();
        }
    }
    
    
    public void loadNELSA(int index, boolean forceReload){
        try{
            if(index < sequences.size()){
                WSSequence wsseq = sequences.get(index);
                if(wsseq.getNELSA()==null || forceReload){
                    timer.reset();
                    System.out.println("Loading NELSA...");
                    NELSA nelsa = new NELSA();
                    nelsa.load(wsseq.getNELSAFile());
                    nelsa.setSequence(wsseq.getB3seq());
                    wsseq.setNELSA(nelsa);
                    System.out.println("done "+timer.getElapsedSecs() +" sec.\n");
                    //System.out.println(b3seq+"\t"+wsseq.getB3seq());
                }
            }
        }catch(Exception ex){
            ex.printStackTrace();
        }
    }
}
