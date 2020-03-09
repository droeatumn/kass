#!/bin/bash
# for 2DP1, had to edit msa2prfl.pl and change die to warn on line 555
for file in *_prot_msa.fasta
do
        # remove the file extension
    name=$(echo $file | cut -f 1 -d '.')
    echo "#!/bin/bash"
    echo "$HOME/git/Augustus/scripts/msa2prfl.pl --qij $HOME/git/Augustus/config/profile/default.qij "$file" > "$name".prfl"
done
