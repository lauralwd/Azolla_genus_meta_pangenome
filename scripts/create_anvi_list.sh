#!/bin/bash

output=scripts/test.tab
MT=data/bin_master_table.tab

echo -e "name\tbin_id\tcollection_id\tprofile_db_path\tcontigs_db_path" > "$output"

while read line
do  bin_id=$(echo $line | cut -f 2 -d ' ')
    spec=$(echo $line | cut -f 1 -d ' ')
    name="$spec"_"$bin_id"
    col='refined'
    pdb=../data/MAG_anvi_dbs/"$spec"_PROFILE.db
    cdb=../data/MAG_anvi_dbs/"$spec"_contigs.db
    echo -e "$name\t$bin_id\t$col\t$pdb\t$cdb" >> "$output"
done < <(tail -n +2 "$MT")
