#!/bin/bash
input=$1
output=$2
name=$3
dir="$(cd "$(dirname "$0")" && pwd)"
mkdir -p $PWD/$name/files
f=$(file $input)
zip7Type="7-zip archive"
tarType="tar archive"
bzip2Type="bzip2 compressed"
gzipType="gzip compressed"
zipType="Zip archive"
rarType="RAR archive"
zxType="XZ compressed data"

if [[ "$f" == *"$zip7Type"* ]]; then
	echo "7-zip"
	echo "Trying: 7za e $input -o$PWD/files/"
	7za e $input -o$PWD/$name/files
fi

if [[ "$f" == *"$tarType"* ]]
then
	echo "tar archive"
	echo "Trying: tar xvf $input -C $PWD/files/"
	tar -xvf $input -C $PWD/$name/files
fi

if [[ "$f" == *"$bzip2Type"* ]]
then
	echo "bzip2 compressed data"
	echo "Trying: tar jxf $input -C $PWD/files/"
	tar -jxf $input -C $PWD/$name/files
fi

if [[ "$f" == *"$gzipType"* ]]
then
	echo "gzip compressed data"
	echo "Trying: tar xvzf $input -C $PWD/files/"
	tar -xvzf $input -C $PWD/$name/files
fi

if [[ "$f" == *"$zipType"* ]]
then
	echo "Zip archive"
	echo "Trying: unzip $input -d $PWD/files/"
	unzip $input -d $PWD/$name/files > $PWD/unziplog.log
fi

if [[ "$f" == *"$rarType"* ]]
then
	echo "RAR archive"
	echo "Trying: unrar e $input $PWD/files/"
	unrar e $input $PWD/$name/files
fi

if [[ "$f" == *"$zxType"* ]]
then
	echo "xz compressed data"
	echo "Trying: tar -xJf $input -C $PWD/files/"
	tar xJf $input -C $PWD/$name/files
fi
find $PWD/$name/files -iname "1_*" -exec cat {} + > $PWD/$name/summ.txt
find $PWD/$name/files -iname "3_*" -exec cat {} + > $PWD/$name/sequences.txt
find $PWD/$name/files -iname "4_*" -exec cat {} + > $PWD/$name/gapped_aa.txt
find $PWD/$name/files -iname "5_*" -exec cat {} + > $PWD/$name/aa.txt
find $PWD/$name/files -iname "6_*" -exec cat {} + > $PWD/$name/junction.txt

echo "summ.txt `cat $PWD/$name/summ.txt | wc -l`"
echo "aa.txt `cat $PWD/$name/aa.txt | wc -l`"

#python $dir/imgt_loader.py --summ $PWD/$name/summ.txt --aa $PWD/$name/aa.txt --junction $PWD/$name/junction.txt --output $output

Rscript --verbose $dir/imgt_loader.r $PWD/$name/summ.txt $PWD/$name/sequences.txt $PWD/$name/aa.txt $PWD/$name/junction.txt $PWD/$name/gapped_aa.txt $output 2>&1
