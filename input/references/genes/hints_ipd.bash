#!/bin/bash
for file in blat/*.psi
do
        # remove the file extension
        name=$(echo $file | cut -f 1 -d '.')
       echo "docker run -it --rm -v $(pwd):/inout chrishah/premaker-plus:18 blat2hints.pl  --in=/inout/"$file" --out=/inout/"$name".gff"
#       `docker run -it --rm -v $(pwd):/inout chrishah/premaker-plus:18 blat2hints.pl  --in=/inout/"$file" --out=/inout/"$name".gff`
done
