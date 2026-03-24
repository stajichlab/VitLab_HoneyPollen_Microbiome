pushd input
echo "sample,read_1,read_2" > ../samples.csv
ls -1 *_R1_* |  perl -p -e 'chomp; $a=$_; $b=$_; $b =~ s/^JS_//; $a =~ s/_R1_/_R2_/; $b=~s/(\S+)_S\d+_L\d+_R\d_001\.fastq\.gz/$1/; $_="$b,$_,$a\n"' >> ../samples.csv

