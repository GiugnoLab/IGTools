/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package igtools.gui2;

import igtools.common.sequence.B3LLSequence;
import igtools.dictionaries.elsa.INELSA;
import igtools.dictionaries.elsa.NELSA;

/**
 *
 * @author vbonnici
 */
public class WSSequence {
    
    private String name = "";
    
    private String fastaFile = null;
    private String b3seqFile = null;
    private String nelsaFile = null;
    
    private B3LLSequence b3seq = null;
    private NELSA nelsa = null;
    
    
    public WSSequence(){
    }
    
    public WSSequence(
            String name,
            String fastaFile,
            String b3seqFile,
            String nelsaFile,
            B3LLSequence b3seq,
            NELSA nelsa
            ){
        this.name = name;
        this.fastaFile = fastaFile;
        this.b3seqFile = b3seqFile;
        this.nelsaFile = nelsaFile;
        this.b3seq = b3seq;
        this.nelsa = nelsa;
    }
    
    public WSSequence(WSSequence seq){
        this.name = seq.name;
        this.fastaFile = seq.fastaFile;
        this.b3seqFile = seq.b3seqFile;
        this.nelsaFile = seq.nelsaFile;
        this.b3seq = seq.b3seq;
        this.nelsa = seq.nelsa;
    }
    
    public String getName(){
        return name;
    }
    public void setName(String name){
        this.name = name;
    }
    
    public String getFASTAFile(){
        return fastaFile;
    }
    public void setFASTAFile(String fastaFile){
        this.fastaFile = fastaFile;
    }
    
    public String getB3seqFile(){
        return b3seqFile;
    }
    public void setB3seqFile(String b3seqFile){
        this.b3seqFile = b3seqFile;
    }
    
    public String getNELSAFile(){
        return nelsaFile;
    }
    public void setNELSAFile(String nelsaFile){
        this.nelsaFile = nelsaFile;
    }
    
    public B3LLSequence getB3seq(){
        return b3seq;
    }
    public void setB3seq(B3LLSequence b3seq){
        this.b3seq = b3seq;
    }
    
    public NELSA getNELSA(){
        return nelsa;
    }
    public void setNELSA(NELSA nelsa){
        this.nelsa = nelsa;
    }
    
    
    @Override
    public String toString(){
        return name;
    }
}
