Several scripts was made available to parallelize the download, but you can also download sequentially all files using a simple script using a loop such as:<br/>
cd AfricanRice/Oryza_barthii<br/>
for i in {ERR2008855 ERR2008856 ERR2008857 .. }; #all SRA/EBI accession ID to download<br/>
  mkdir $i # create a directory for each run<br/>
  cd $i<br/>
  do start(echo -n $i| head -c 6) # first 6 characters of the ID<br/>
  end=$(echo -n $i| tail -c 1) # last character of the ID<br/>
  path1=$(echo "ftp://ftp.sra.ebi.ac.uk/vol1/fastq/""$start""/00""$end""/""$i""/""$i""_1.fastq.gz"<br/>
  path2=$(echo "ftp://ftp.sra.ebi.ac.uk/vol1/fastq/""$start""/00""$end""/""$i""/""$i""_1.fastq.gz"<br/>
  wget $path1<br/>
  wget $path2<br/>
  cd ..<br/>
done<br/>
