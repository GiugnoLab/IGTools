/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package igtools.common.alignment;

import igtools.common.nucleotide.B3Nucleotide;
import igtools.common.sequence.B3LLSequence;
import igtools.common.sequence.B3Sequence;

/**
 *
 * @author vbonnici
 */
public class SmithWaterman {
    
    private final static int gapD = -1;
    
    private final static int[][] simMatrix = 
    {  
        {2,-1,-1,-1,-1},
        {-1,2,-1,-1,-1},
        {-1,-1,2,-1,-1},
        {-1,-1,-1,2,-1},
        {-1,-1,-1,-1,2},
    };
    
    
    private B3Sequence s1;
    private B3Sequence s2;
    
    private int[][] mF = null;
    
    private B3Sequence as1;
    private B3Sequence as2;
    
    
    public SmithWaterman(B3Sequence s1, B3Sequence s2){
        this.s1 = s1;
        this.s2 = s2;
    }
    
    public void process(){
        fillMf();
        getAligedSequences();
    }
    
    
    public int[][] getScoreMatrix(){
        return mF;
    }
    public B3Sequence getAlignedS1(){
        return as1;
    }
    public B3Sequence getAlignedS2(){
        return as2;
    }
    
    public int score(){
        if(mF == null)
            fillMf();
        
        
        int i = 1;
        int j = 1;
        int max = mF[i][j];

        for (int k = 1; k <= s1.length(); k++) {
                for (int l = 1; l <= s2.length(); l++) {
                        if (mF[k][l] > max) {
                                i = k;
                                j = l;
                                max = mF[k][l];
                        }
                }
        }

        return mF[i][j];
    }
    
    
    public void fillMf(){
        mF = new int[s1.length() +1][s2.length() +1];
        
        for(int i=0; i<mF.length; i++)
            mF[i][0] = 0;
        for(int j=0; j<mF[0].length; j++)
            mF[0][j] = 0;
        
        int match, delete, insert;
        
        for(int i=1; i<mF.length; i++){
            for(int j=1; j<mF[i].length; j++){
                
                match = mF[i-1][j-1] + simMatrix[s1.getB3(i-1)] [s2.getB3(j-1)];
                delete = mF[i-1][j] + gapD;
                insert = mF[i][j-1] + gapD;
                
                mF[i][j] = Math.max(match, Math.max(delete, insert));
                
            }
        }
    }
    
    
    public void getAligedSequences(){
        
        String mAlignmentSeqB = "";
        String mAlignmentSeqA = "";
        
//        int a_i;
//        int a_j;
        
        int i = 1;
        int j = 1;
        int max = mF[i][j];

        for (int k = 1; k <= s1.length(); k++) {
                for (int l = 1; l <= s2.length(); l++) {
                        if (mF[k][l] > max) {
                                i = k;
                                j = l;
                                max = mF[k][l];
                        }
                }
        }

//        mScore = mD[i][j];

        int k = s1.length();
        int l = s2.length();

        while (k > i) {
                mAlignmentSeqB += "N";
                mAlignmentSeqA += B3Nucleotide.by(s1.getB3(k - 1));
                k--;
        }
        while (l > j) {
                mAlignmentSeqA += "N";
                mAlignmentSeqB += B3Nucleotide.by(s2.getB3(l - 1));
                l--;
        }
        
        
       while (mF[i][j] != 0) {                 
                if (mF[i][j] == mF[i-1][j-1] +simMatrix[s1.getB3(i-1)][s2.getB3(j-1)]) {                          
                        mAlignmentSeqA += B3Nucleotide.by(s1.getB3(i-1));
                        mAlignmentSeqB += B3Nucleotide.by(s2.getB3(j-1));
                        i--;
                        j--;                            
                        continue;
                } else if (mF[i][j] == mF[i][j-1] - 1) {
                        mAlignmentSeqA += "N";
                        mAlignmentSeqB += B3Nucleotide.by(s2.getB3(j-1));
                        j--;
                        continue;
                } else {
                        mAlignmentSeqA += B3Nucleotide.by(s1.getB3(i-1));
                        mAlignmentSeqB += "N";
                        i--;
                        continue;
                }
        }

        while (i > 0) {
                mAlignmentSeqB += "N";
                mAlignmentSeqA += B3Nucleotide.by(s1.getB3(i - 1));
                i--;
        }
        while (j > 0) {
                mAlignmentSeqA += "N";
                mAlignmentSeqB += B3Nucleotide.by(s2.getB3(j - 1));
                j--;
        }
        
        mAlignmentSeqA = new StringBuffer(mAlignmentSeqA).reverse().toString();
        mAlignmentSeqB = new StringBuffer(mAlignmentSeqB).reverse().toString();
        
        as1 = new B3LLSequence(mAlignmentSeqA);
        as2 = new B3LLSequence(mAlignmentSeqB);
        
    }
    
    
    
    
    
    
    public static int score(B3Sequence s1, B3Sequence s2){
        int[][] mF = new int[s1.length() +1][s2.length() +1];
        
        for(int i=0; i<mF.length; i++)
            mF[i][0] = 0;
        for(int j=0; j<mF[0].length; j++)
            mF[0][j] = 0;
        
        int match, delete, insert;
        
        for(int i=1; i<mF.length; i++){
            for(int j=1; j<mF[i].length; j++){
                
                match = mF[i-1][j-1] + simMatrix[s1.getB3(i-1)] [s2.getB3(j-1)];
                delete = mF[i-1][j] + gapD;
                insert = mF[i][j-1] + gapD;
                
                mF[i][j] = Math.max(match, Math.max(delete, insert));
                
            }
        }
        
        int i = 1;
        int j = 1;
        int max = mF[i][j];

        for (int k = 1; k <= s1.length(); k++) {
                for (int l = 1; l <= s2.length(); l++) {
                        if (mF[k][l] > max) {
                                i = k;
                                j = l;
                                max = mF[k][l];
                        }
                }
        }

        return mF[i][j];
    }
    
    
    
    
    
    public static void main(String[] args){
        B3LLSequence s1 = new B3LLSequence("TAATC");
        B3LLSequence s2 = new B3LLSequence("TACGAGTATGGAGC");
        NeedlemanWunsch alignment = new NeedlemanWunsch(s1, s2);
        //NeedlemanWunsch alignment = new NeedlemanWunsch(new ReverseB3Sequence(s1), new ReverseB3Sequence(s2));
        
        System.out.println(alignment.getAlignedS1());
        System.out.println(alignment.getAlignedS2());
        
        int[][] mat = alignment.getScoreMatrix();
        for(int i=0; i<mat.length; i++){
            for(int j=0; j<mat[i].length; j++){
                System.out.print(mat[i][j]+"\t");
            }
            System.out.println();
        }
    }
}
