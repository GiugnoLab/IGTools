/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package igtools.common.alignment;

import igtools.common.nucleotide.B3Nucleotide;

/**
 *
 * @author vbonnici
 */
public class B3NSmithWaterman {
    
    private final static int gapD = -1;
    
//    private final static int[][] simMatrix = 
//    {  
//        {2,-1,-1,-1,-1},
//        {-1,2,-1,-1,-1},
//        {-1,-1,2,-1,-1},
//        {-1,-1,-1,2,-1},
//        {-1,-1,-1,-1,2},
//    };
    
    
    private B3Nucleotide[] s1;
    private B3Nucleotide[] s2;
    
    private int[][] mF = null;
    
    private B3Nucleotide[] as1;
    private B3Nucleotide[] as2;
    
    
    public B3NSmithWaterman(B3Nucleotide[] s1, B3Nucleotide[] s2){
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
    public B3Nucleotide[] getAlignedS1(){
        return as1;
    }
    public B3Nucleotide[] getAlignedS2(){
        return as2;
    }
    
    public int score(){
        if(mF == null)
            fillMf();
        
        
        int i = 1;
        int j = 1;
        int max = mF[i][j];

        for (int k = 1; k <= s1.length; k++) {
                for (int l = 1; l <= s2.length; l++) {
                        if (mF[k][l] > max) {
                                i = k;
                                j = l;
                                max = mF[k][l];
                        }
                }
        }

        return mF[i][j];
    }
    
    private int weight(int i, int j){
        if(s1[i] == s2[j])
            return 2;
        return -1;
    }
    
    public void fillMf(){
        mF = new int[s1.length +1][s2.length +1];
        
        for(int i=0; i<mF.length; i++)
            mF[i][0] = 0;
        for(int j=0; j<mF[0].length; j++)
            mF[0][j] = 0;
        
        int match, delete, insert;
        
        for(int i=1; i<mF.length; i++){
            for(int j=1; j<mF[i].length; j++){
                
                match = mF[i-1][j-1] + weight(i-1, j-1);
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

        for (int k = 1; k <= s1.length; k++) {
                for (int l = 1; l <= s2.length; l++) {
                        if (mF[k][l] > max) {
                                i = k;
                                j = l;
                                max = mF[k][l];
                        }
                }
        }

//        mScore = mD[i][j];

        int k = s1.length;
        int l = s2.length;

        while (k > i) {
                mAlignmentSeqB += "N";
                mAlignmentSeqA += s1[k-1];
                k--;
        }
        while (l > j) {
                mAlignmentSeqA += "N";
                mAlignmentSeqB += s2[l-1];
                l--;
        }
        
        
       while (mF[i][j] != 0) {                 
                if (mF[i][j] == mF[i-1][j-1] + weight(i-1,j-1)) {                          
                        mAlignmentSeqA += s1[i-1];
                        mAlignmentSeqB += s2[j-1];
                        i--;
                        j--;                            
                        continue;
                } else if (mF[i][j] == mF[i][j-1] - 1) {
                        mAlignmentSeqA += "N";
                        mAlignmentSeqB += s2[j-1];
                        j--;
                        continue;
                } else {
                        mAlignmentSeqA += s1[i-1];
                        mAlignmentSeqB += "N";
                        i--;
                        continue;
                }
        }

        while (i > 0) {
                mAlignmentSeqB += "N";
                mAlignmentSeqA += s1[i-1];
                i--;
        }
        while (j > 0) {
                mAlignmentSeqA += "N";
                mAlignmentSeqB += s2[j-1];
                j--;
        }
        
//        mAlignmentSeqA = new StringBuffer(mAlignmentSeqA).reverse().toString();
//        mAlignmentSeqB = new StringBuffer(mAlignmentSeqB).reverse().toString();
        
        int length = mAlignmentSeqA.length() - 1;
        as1 = new B3Nucleotide[mAlignmentSeqA.length()];
        for(i=0; i<as1.length; i++){
            as1[i] = B3Nucleotide.by(mAlignmentSeqA.charAt(length - i));
        }
        
        length = mAlignmentSeqB.length() - 1;
        as2 = new B3Nucleotide[mAlignmentSeqB.length()];
        for(i=0; i<as2.length; i++){
            as2[i] = B3Nucleotide.by(mAlignmentSeqB.charAt(length - i));
        }
        
    }
    
    
    
    
    
    private static int weight(B3Nucleotide[] s1, B3Nucleotide[] s2, int i, int j){
        if(s1[i] == s2[j])
            return 2;
        return -1;
    }
    public static int score(B3Nucleotide[] s1, B3Nucleotide[] s2){
        int[][] mF = new int[s1.length +1][s2.length +1];
        
        for(int i=0; i<mF.length; i++)
            mF[i][0] = 0;
        for(int j=0; j<mF[0].length; j++)
            mF[0][j] = 0;
        
        int match, delete, insert;
        
        for(int i=1; i<mF.length; i++){
            for(int j=1; j<mF[i].length; j++){
                
                match = mF[i-1][j-1] + weight(s1,s2,i-1,j-1);
                delete = mF[i-1][j] + gapD;
                insert = mF[i][j-1] + gapD;
                
                mF[i][j] = Math.max(match, Math.max(delete, insert));
                
            }
        }
        
        int i = 1;
        int j = 1;
        int max = mF[i][j];

        for (int k = 1; k <= s1.length; k++) {
                for (int l = 1; l <= s2.length; l++) {
                        if (mF[k][l] > max) {
                                i = k;
                                j = l;
                                max = mF[k][l];
                        }
                }
        }

        return mF[i][j];
    }
    
}
