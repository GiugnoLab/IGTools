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
 
 
