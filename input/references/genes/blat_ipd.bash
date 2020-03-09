#!/bin/bash

echo
cwd=$(echo "$(pwd)")

# docker run -it -v $(pwd):/inout/ --rm hiroko/blat-hg19:latest blat /db/hg19/hg19.2bit /inout/KIR2DL1S1S2_nuc.fasta /inout/blat/KIR2DL1S1S2_nuc.psi
for file in *_nuc.fasta
do
	# remove the file extension
	name=$(echo $file | cut -f 1 -d '.')
	echo "docker run -it -v "$cwd":/inout/ --rm hiroko/blat-hg19:latest blat /db/hg19/hg19.2bit /inout/"$file" /inout/blat/"$name".psi"
#	`docker run -it -v "$cwd":/inout/ --rm hiroko/blat-hg19:latest blat /db/hg19/hg19.2bit /inout/"$file" /inout/blat/"$name".psi`


done

# docker run -it --rm -v $(pwd):/inout cgwyx/augustus:latest blat2hints.pl --in=/inout/KIR2DL1S1S2_nuc.psi --out=/inout/blat/KIR2DL1S1S2_nuc.gff
for file in blat/*_nuc.psi
do
	# remove the file extension
	name=$(echo $file | cut -f 1 -d '.')
    echo "docker run -it --rm -v $(pwd):/inout chrishah/premaker-plus:18 blat2hints.pl  --in=/inout/"$file" --out=/inout/"$name".gff"
done


exit 0
