#!/bin/bash


#java -server -d64 -Xmn2560M -Xms6144M -Xmx6144M -cp igtools.jar <command> <parameters>
#The example shows optimal settings for a 64bit machine with 8Gb of RAM.
#see also https://bitbucket.org/infogenomics/igtools/wiki/Documentation

java -cp igtools.jar igtools.cli.util.FASTATo3bit nanoa.fa nanoa.3bit
java -cp igtools.jar igtools.cli.dictionaries.BuildNELSA nanoa.3bit nanoa.nelsa
