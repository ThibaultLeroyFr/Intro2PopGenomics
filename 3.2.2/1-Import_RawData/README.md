Our objective was to be as simple as possible at the very beginning of this manual. As a consequence, we do not use loops. <br/>
It is however possible to download sequentially all files using a simple bash script with a loop such as:<br/>
<pre>
cd ./AfricanRice/Oryza_barthii
for i in ERR2008855 ERR2008856 ERR2008857 ... ; do #all SRA/EBI accession ID to download
  mkdir $i # create a directory for each run
  cd $i
  start=$(echo -n $i| head -c 6) # first 6 characters of the ID
  end=$(echo -n $i| tail -c 1) # last character of the ID
  path1=$(echo "ftp://ftp.sra.ebi.ac.uk/vol1/fastq/""$start""/00""$end""/""$i""/""$i""_1.fastq.gz")
  path2=$(echo "ftp://ftp.sra.ebi.ac.uk/vol1/fastq/""$start""/00""$end""/""$i""/""$i""_1.fastq.gz")
  wget $path1
  wget $path2
  cd ..
done
</pre>
