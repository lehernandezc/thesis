#!/usr/bin/env python

from os import listdir
from os.path import isfile, join
from subprocess import call
from os import system


if __name__ == "__main__":
    mypath = "."
    fastqs = [ f for f in listdir(mypath) if isfile(join(mypath,f)) and f.endswith(".fastq") ]

    for f in fastqs:
        print "Processing file: {}".format(f)
        splittedName = f.split("_")
        baseName = "_".join(splittedName[:2])
        c1 = "fastx_trimmer -f 1 -Q 33 -i {} -o {}.fastq".format(f, baseName)
        system(c1)
        c2 = "cat {}.fastq | perl -e '$i=0;while(<>){{if(/^\@/&&$i==0){{s/^\@/\>/;print;}}elsif($i==1){{print;$i=-3}}$i++;}}' > {}.fasta".format(baseName, baseName)
        system(c2)

        c3 = "perl 2b_Extract.pl {}.fasta  .{{10}}CCAC.{{5}}TTG.{{10}} {}_R.fasta".format(baseName,baseName)
        system(c3)

        c4 = "perl 2b_Extract.pl {}.fasta  .{{10}}CAA.{{5}}GTGG.{{10}} {}_F.fasta".format(baseName, baseName)
        system(c4)

        c5 = "perl ReverseComplement.pl {}_R.fasta {}_R_reverse.fasta".format(baseName, baseName)
        system(c5)

        c6 = "cat {}_F.fasta {}_R_reverse.fasta > {}_merged.fasta".format(baseName, baseName, baseName)
        system(c6)

