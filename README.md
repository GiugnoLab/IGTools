# IGTools
>Infogenomics Tools: A Computational Suite for Informational Analyses of Genomes.


<hr />

### Description
Infogenomics is a project aimed at developing informational analyses of genomes adding a new perspective to the more common biological and biochemical investigation on genomes. InfoGenomics Tools, shortly IGTools, is a computational framework, which consists of a collection of interactive tools, designed to support typical analyses required in the context of Infogenomics. Modularity, interactivenes, data visualization, and low-cost computational requirements are the main goal of IGTools, where suffix arrays are a key point in the data structures representing genomes, and their algorithmic power allows us to speed-up genomic computations that would be prohibitive by using naive methods.
IGTools is completely developed in Java 8, and it is published under GNU General Public License 3.

IGTools works on top of data structures optimized for the genomic alphabet {A,C,G,T,N}. Public available genomic sequences are often distributed in FASTA format.
Before to run an IGTools analysis, you need to convert the FASTA file into a specific binary sequence format named 3bit. After the conversion, an indexing data structure, called NELSA, must be built. Please, see the CLI documentation page in order to understand how to convert from FASTA to 3bit and how to built a NELSA index and save it on file.
Moreover, if you are using the GUI application, please make sure that those two files have been made before to run analyses.

By default, the Java Virtual Environment uses a limited amount of memory (RAM) that is available on your computer. In some cases, you may need to increase the total amount used by the execution by setting some JVM directives.
 ```
java -server -d64 -Xmn2560M -Xms6144M -Xmx6144M -cp igtools.jar <command> <parameters>
 ```
 
The example shows optimal settings for a 64bit machine with 8Gb of RAM.

<hr />

### CLI - Command Line Interface
```
java -cp igtools.jar command parameters
```

Running commands without parameters will output the help.


##### Basic commands
Convert from FASTA to 3bit format.
```
java -cp igtools.jar igtools.cli.util.FASTATo3bit <isequence.fa> <osequence.3bit>
```
The tool takes in input a file in FASTA format storing a genomic sequence and convert it into the 3bit IGTools format.

Generate the NELSA structure and save it on file
```
java -cp igtools.jar igtools.cli.dictionaries.BuildNELSA <isequence.3bit> <onelsa.nelsa>
```

##### Informational indexes
Obtain some features of a genomic sequence G, such as sequence length and number of characters N , plus

+ MHL(G): minimal hapax length
+ MRL(G): maximal repeat length
+ MFL(G): minimal forbidden length
+ empirical entropy at MFL(G) and MFL(G) - 1

```
java -cp igtools.jar igtools.cli.GenomeStats <isequence.3bit> <onelsa.nelsa>
```
Obtain cardinality of the following information sets, with respect to a genomic sequence G, for a given range of word length k

+ |D_k| the set of k-mers in G
+ |H_k| the set of hapax k-mer in G
+ |R_k| the set of repeat k-mers in G
+ |T_k| the overall k-mer multiplicity
+ |E_k| the empirical entropy at word length k

```
java -cp igtools.jar igtools.cli.GenomeKStats <from_k> <to_k> <isequence.3bit> <onelsa.nelsa>
```

Get the coverage of H_k(G) and R_k(G) for several values of k.
```
java -cp igtools.jar igtools.cli.coverage.HRCoverage <isequence.3bit> <onelsa.nelsa> <k...>
```
Enumerate k-mers, for a given k, and their multiplicity.
```
java -cp igtools.jar igtools.cli.dictionaries.EnumerateKmers <isequence.3bit> <onelsa.nelsa> <k>
```

Show k-mer frequencies.
```
java -cp igtools.jar igtools.cli.dictionaries.EnumerateKmers <isequence.3bit> <onelsa.nelsa> <k>
```

##### Genomic distributions
Multiplicity/Co-Multiplicity distribution.
```
java -cp igtools.jar igtools.cli.distributions.MultiplicityDistribution <k> <isequence.3bit> <onelsa.nelsa>
```

##### Dictionaries in FASTA format
Words are listed in FASTA files. Every line is a word.
Add reverse complement (RC) of words.
```
java -cp igtools.jar igtools.cli.wordset.AddRC <seqs.fa>
```
Calculate the sequence coverage of a set of sequences (a dictionary) in a given genomic sequence (3bit + NELSA)
```
java -cp igtools.jar igtools.cli.wordset.Coverage <seqs.fa> <isequence.3bit> <inelsa.nelsa>
```
Calculate the following statistics for a set of input words over a genomic sequence (3bit + NEKLA):

+ total number of input words
+ um of their length
+ sequence coverage
+ coverage ratio
+ number of contiguous covered regions
+ length distribution
+ multiplicity distribution
+ number of hapaxes
+ number of minimal hapaxes
+ number of repeats
+ number of maximal repeats

```
java -cp igtools.jar igtools.cli.wordset.WordsProperties <seqs.fa> <isequence.3bit> <inelsa.nelsa>
```
Calculate only the multiplicity distribution
```
java -cp igtools.jar igtools.cli.wordset.MultiplicityDistribution <seqs.fa> <isequence.3bit> <inelsa.nelsa>
```


##### Dictionaries in FASTA format : set-theoretic operations
Words are listed in FASTA files. Every line is a word.
Select words of a given length.
```
java -cp igtools.jar igtools.cli.wordset.SelectByLength [eq|geq|leq] <k> <seqs.fa>
```

Print statistics about two distinct sets of words s1 and s2:

+ number of words in s1 and s2
+ number of words in s1 (s2) that are the reverse complement of other words in s1 (s2)
+ number of words in s2 that are in s1
+ cardinality of the union of the two input sets

```
java -cp igtools.jar igtools.cli.wordset.SharingStats <s1.fa> <s2.fa>
```
Set union for two sets of words s1 and s2:
```
java -cp igtools.jar igtools.cli.wordset.Difference <s1.fa> <s2.fa>
```
Set union for a list of input sets.
```
java -cp igtools.jar igtools.cli.wordset.Difference {files...}
```

Set intersection between two sets of words s1 and s2:

+ exact : list words that are shared between the two sets
+ prefix : list words of s1 if there exists at least a word in s2 such that the word in s1 is a prefix of the word in s2
+ suffix : list words of s1 if there exists at least a word in s2 such that the word in s1 is a suffix of the word in s2
+ sub : list words of s1 if there exists at least a word in s2 such that the word in s1 is a substring of the word in s2

```
java -cp igtools.jar igtools.cli.wordset.Intersection2 [exact|prefix|suffix|sub] <s1.fa> <s2.fa>
```

Set difference between two sets of words s1 and s2:
```
java -cp igtools.jar igtools.cli.wordset.Difference <s1.fa> <s2.fa>
```
Given a set of words, remove duplicates and sub-included words
```
java -cp igtools.jar igtools.cli.wordset.MinimalList <s1.fa>
```

<hr />

### For developers
Read a sequence from a 3bit format file.
```
 import igtools.common.sequence.B3LLSequence;
 String sequence_file_path = ...
 B3LLSequence b3seq = B3LLSequence.load(sequence_file_path);
```
 Load a NELSA structure from file.
```
 import igtools.common.sequence.B3LLSequence;
 import igtools.dictionaries.elsa.NELSA;
 ...

 String a_inelsa = ...
 NELSA nelsa = new NELSA();
 nelsa.load(a_inelsa);
 nelsa.setSequence(b3seq);
```

Enumerate k-mers and print their multiplcitiy.
```
 import igtools.common.sequence.B3LLSequence;
 import igtools.dictionaries.elsa.NELSA;
 import igtools.dictionaries.elsa.IELSAIterator;
 import igtools.common.nucleotide.B3Nucleotide;
 ...

 B3Nucleotide[] kmer;
 int multiplicity;
 IELSAIterator it = nelsa.begin(k);
 while(it.next()){
  kmer = it.kmer();
  multiplicity = it.multiplcity()
  print(B3Nucleotide.toString(kmer) +"\t"+ multiplicity);
}
```
...a more efficient solution
```
 B3Nucleotide[] kmer = new B3Nucleotide[k];
 int multiplicity;
 IELSAIterator it = nelsa.begin(k);
 while(it.next()){
  it.kmer(kmer);
  multiplicity = it.multiplcity()
  print(B3Nucleotide.toString(kmer) +"\t"+ multiplicity);
}
```
Retrieve k-mer positions.
```
B3Nucleotide[] kmer = new B3Nucleotide[k];
int positions[];
IELSAIterator it = nelsa.begin(k);
while(it.next()){
 positions = it.positions();
 //positions = it.sortedPositions();
}
```
Search a k-mer over a genomic sequence.
```
B3Nucleotide kmer[] = B3Nucleotide.toB3("ACGGTGGCT");
IELSAIterator it = nelsa.find(kmer);
if(it != null){
 ...print k-mer
}
```
**Informational indexes*** and others.
```
import igtools.analyses.*;
...

String sequence_file_path = ...
B3LLSequence b3seq = B3LLSequence.load(sequence_file_path);
String a_inelsa = ...
NELSA nelsa = new NELSA();
nelsa.load(a_inelsa);
nelsa.setSequence(b3seq);

System.out.println("MRL " + Inform.mrl(nelsa));
System.out.println("MHL " + Inform.mhl(nelsa));
System.out.println(Inform.kCompleteness(nelsa) + "is the greatest k for which the sequence is k-complete");
System.out.println("Empirical entropy at k=3, "+Inform.entropy(nelsa, 3));
```
Get the multiplicity distribution:
```
import igtools.analyses.*;
...

String sequence_file_path = ...
B3LLSequence b3seq = B3LLSequence.load(sequence_file_path);
String a_inelsa = ...
NELSA nelsa = new NELSA();
nelsa.load(a_inelsa);
nelsa.setSequence(b3seq);

int k = 3;

Map<Integer,Integer> distr = Infrom.multiplicityDistribution(nelsa, 3);
for(Map.Entry<Integer,Integer> ee : distr.entrySet()){
 System.out.println(ee.getKey() +"\t"+ ee.getValue());
}
```
Extract the proper recurrence distance distribution (RDD) of k-mers.
```
import igtools.common.sequence.B3LLSequence;
import igtools.dictionaries.elsa.NELSA;
import igtools.dictionaries.elsa.IELSAIterator;
import igtools.common.nucleotide.B3Nucleotide;
import igtools.analyses.recurrences.distances.*;
...

B3LLSequence seq = B3LLSequence.load(a_seq);
NELSA nelsa = new NELSA();
nelsa.load(a_nelsa);
nelsa.setSequence(seq);
ProperRecurrenceDistancesExtractor rddext = ProperMinimalRecurrenceDistancesExtractor.factory(true, true);
IELSAIterator it = nelsa.begin(k);
while(it.next()){
 System.out.println(seq.toString(it.kmer()));
 Map<Double,Double> it_rdd = rddext.recurrenceDistanceDistributionMap(it);
}
```
Extract the recurrence distance distribution of two distinct k-mers.
```
import igtools.common.sequence.B3LLSequence;
import igtools.dictionaries.elsa.NELSA;
import igtools.dictionaries.elsa.IELSAIterator;
import igtools.common.nucleotide.B3Nucleotide;
import igtools.analyses.recurrences.distances.*;
...

B3LLSequence seq = B3LLSequence.load(a_seq);
NELSA nelsa = new NELSA();
nelsa.load(a_nelsa);
nelsa.setSequence(seq);


BiMinialRecurrenceDistancesExtractor rddext = BiMinimalRecurrenceDistancesExtractor.factory(true, true);
IELSAIterator it_a = nelsa.find( B3Nucleotide.toB3("ACT") );
IELSAIterator it_b = nelsa.find( B3Nucleotide.toB3("TGC") );

if(it_a != null && it_b != null){
 Map<Double,Double> rdd = rddext.recurrenceDistanceDistributionMap(it_a, it_b);
}
```
 
 <hr />
 
 ### GUI - Graphical User Interface
How to run the GUI:
```
java -cp igtools.jar igtools.gui2.IGToolsGUI2
```

 <hr />
 
### People
Prof. Vincenzo Manca

Dr. Vincenzo Bonnici

Department of Computer Science, University of Verona, Italy.

<hr />

### Publications
Bonnici, V. & Manca, V. Infogenomics Tools: A Computational Suite for Informational Analyses of Genomes. (2015) J. Bioinfo. Proteomics Rev 1(1): 7- 14.
