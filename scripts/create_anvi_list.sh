#!/bin/bash
# it a subset is set, call it subset
subset="$1"

# bin Master table
MT=data/bin_master_table.tab

# all is no subset
if   [ "$subset" == 'all' ]
then unset subset
fi

# use subset var in output name if a subset is set al all
if   [ "$subset" ]
then output=scripts/anvi_genomes_internal_"$subset".tab
else output=scripts/anvi_genomes_internal_all.tab
fi

# create temporary output file for subsetting purposes
MTt="$MT".temp
cut -f 1,2,27- -d ' ' "$MT" | grep "$subset" > "$MTt"

# first make a header in the output file
echo -e "name\tbin_id\tcollection_id\tprofile_db_path\tcontigs_db_path" > "$output"

# build up the output file line by line
while read line
do  bin_id=$(echo $line | cut -f 2 -d ' ')
    spec=$(echo $line | cut -f 1 -d ' ')
    name="$spec"_"$bin_id"
    col='refined'
    pdb=../data/MAG_anvi_dbs/"$spec"_PROFILE.db
    cdb=../data/MAG_anvi_dbs/"$spec"_contigs.db
    echo -e "$name\t$bin_id\t$col\t$pdb\t$cdb" >> "$output"
done < <(tail -n +2 "$MTt" )

rm "$MTt"
