
package igtools.dictionaries.elsa;

import igtools.common.kmer.b2.unit.B2UnitRLKmer;
import igtools.common.nucleotide.B3Nucleotide;

/**
 *
 * @author vbonnici
 */
public class CompleteUnitIterator implements IELSAIterator{
    
    private int k = 0;
    private int code = 0x0;
    private B2UnitRLKmer kmer = new B2UnitRLKmer(code);
    private B3Nucleotide[] ns = null;
    
    public CompleteUnitIterator(int k){
        this.k = k;
        ns = B3Nucleotide.toB3(kmer.nucleotides(k));
    }
    
    protected CompleteUnitIterator(
            int k,
            int code
            ){
        this.k = k;
        this.code = code;
        kmer.setCode(code);
    }
    
    public int maxK(){
        return B2UnitRLKmer.MAX_K;
    }

    public int code(){
        return code;
    }
    
    @Override
    public int istart() {
        throw new UnsupportedOperationException("Not supported yet."); //To change body of generated methods, choose Tools | Templates.
    }

    @Override
    public int iend() {
        throw new UnsupportedOperationException("Not supported yet."); //To change body of generated methods, choose Tools | Templates.
    }

    @Override
    public int k() {
        return k;
    }

    @Override
    public B3Nucleotide[] kmer() {
        //return B3Nucleotide.toB3(kmer.nucleotides(k));
        return ns;
    }
    @Override
    public void kmer(B3Nucleotide[] ns) {
            throw new UnsupportedOperationException("Not supported yet."); //To change body of generated methods, choose Tools | Templates.
    }

    @Override
    public int multiplicity() {
        throw new UnsupportedOperationException("Not supported yet."); //To change body of generated methods, choose Tools | Templates.
    }

    @Override
    public int[] positions() {
        throw new UnsupportedOperationException("Not supported yet."); //To change body of generated methods, choose Tools | Templates.
    }

    @Override
    public int[] sortedPositions() {
        throw new UnsupportedOperationException("Not supported yet."); //To change body of generated methods, choose Tools | Templates.
    }

    @Override
    public boolean isMinimalHapax() {
        throw new UnsupportedOperationException("Not supported yet."); //To change body of generated methods, choose Tools | Templates.
    }

    @Override
    public boolean isGlobalMaximalRepeat() {
        throw new UnsupportedOperationException("Not supported yet."); //To change body of generated methods, choose Tools | Templates.
    }

    @Override
    public boolean next() {
        if(code < B2UnitRLKmer.MAX_CODE[k]){
            code++;
            kmer.setCode(code);
            ns = B3Nucleotide.toB3(kmer.nucleotides(k));
            
            return true;
        }
        return false;
    }

    @Override
    public boolean prev() {
        if(code > 0){
            code--;
            kmer.setCode(code);
            ns = B3Nucleotide.toB3(kmer.nucleotides(k));
            
            return true;
        }
        return false;
    }

    @Override
    public boolean hasNext() {
        return code < B2UnitRLKmer.MAX_CODE[k];
    }

    @Override
    public boolean hasPrev() {
        return code > 0x0;
    }

    @Override
    public IELSAIterator clone() {
        return new CompleteUnitIterator(k, code);
    }

    @Override
    public int compare(IELSAIterator it) {
        B3Nucleotide[] itns = it.kmer();
        for(int i=0; i<(k <= it.k() ? k : it.k()); i++){
            if(ns[i] != itns[i]){
                return ((int)ns[i].code()) - (int)itns[i].code();
            }
        }
        if(k<=it.k())
            return 0;
        return -1;
    }

    @Override
    public ELSA elsa() {
        return null;
    }

    @Override
    public boolean good() {
        throw new UnsupportedOperationException("Not supported yet."); //To change body of generated methods, choose Tools | Templates.
    }
    
}
