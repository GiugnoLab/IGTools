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
