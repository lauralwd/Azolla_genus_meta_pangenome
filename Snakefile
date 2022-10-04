MAG_HOSTS=['Azfil_lab','Azfil_wild','Azfil_minus_cyano','Azmex','Azmic','Aznil','Azrub','Azcar1','Azcar2']
ORDERS=['Burkholderiales', 'Caulobacterales', 'Nevskiales', 'Rhizobiales', 'Sphingomonadales']

############################### stage 1: collect MAGs and genomes from the Azolla genus metagenome project ###############################
rule gather_anvi_MAGs:
  output:
    expand("data/MAG_anvi_dbs/{mag_host}_contigs.db" ,mag_host=MAG_HOSTS),
    expand("data/MAG_anvi_dbs/{mag_host}_PROFILE.db" ,mag_host=MAG_HOSTS)
  shell:
    "bash ./scripts/collect_mag_dbs.sh"

# you may need to remove old hmm, kegg, cogg results, or choose to keep your current ones and ignore the next three rules
rule anvi_contigdb_runhmms:
  input:
    "{db}_contigs.db"
  output:
    touch("{db}_contigs.db.hmms")
  log:
    stderr="logs/anvi_annotate/{db}_hmms.stderr",
    stdout="logs/anvi_annotate/{db}_hmms.stdout"
  threads: 6
  shell:
    """
    anvi-run-hmms -c {input}   \
                  -T {threads} \
                  --just-do-it \
    > {log.stdout} 2> {log.stderr}
    """

rule anvi_contigdb_kegg:
  input:
    "{db}_contigs.db"
  output:
    touch("{db}_contigs.db.kegg")
  log:
    stderr="logs/anvi_annotate/{db}_kegg.stderr",
    stdout="logs/anvi_annotate/{db}_kegg.stdout"
  threads: 6
  shell:
    """
    anvi-run-kegg-kofams -c {input}   \
                         -T {threads} \
                         --just-do-it \
    > {log.stdout} 2> {log.stderr}
    """

rule anvi_contigdb_cogs:
  input:
    "{db}_contigs.db"
  output:
    touch("{db}_contigs.db.cogs")
  log:
    stderr="logs/anvi_annotate/{db}_cogs.stderr",
    stdout="logs/anvi_annotate/{db}_cogs.stdout"
  threads: 6
  shell:
    """
    anvi-run-ncbi-cogs -c {input}   \
                       -T {threads} \
                       --sensitive  \
    > {log.stdout} 2> {log.stderr}
    """

rule anvi_MAGs_to_fasta:
  input:
    "data/MAG_anvi_dbs/{mag_host}_contigs.db"
  output:
    "data/MAG_anvi_dbs/{mag_host}_contigs.fasta"
  log:
    stderr="logs/MAG_anvi_dbs/export_contigs_{mag_host}.stderr",
    stdout="logs/MAG_anvi_dbs/export_contigs_{mag_host}.stdout"
  shell:
    """
    anvi-export-contigs -c {input}  \
                        -o {output} \
    > {log.stdout} 2> {log.stderr}
    """

rule anvi_internal_list_all:
  input:
    "data/bin_master_table.tab"
  output:
    "scripts/anvi_genomes_internal_all.tab"
  shell:
    """
    scripts/create_anvi_list.sh
    """

rule anvi_internal_list_subset:
  input:
    "data/bin_master_table.tab"
  output:
    "scripts/anvi_genomes_internal_{subset}.tab"
  shell:
    """
    scripts/create_anvi_list.sh {wildcards.subset}
    """

rule create_pangenome_storage_internal:
  input:
    dbs=expand("data/MAG_anvi_dbs/{mag_host}_contigs.db"      ,mag_host=MAG_HOSTS                           ),
    ann=expand("data/MAG_anvi_dbs/{mag_host}_contigs.db.{ext}",mag_host=MAG_HOSTS,ext=['hmms','kegg','cogs']),
    txt="scripts/anvi_genomes_internal_{subset}.tab"
  output:
    "data/anvio_genomes_storage/{subset}_GENOMES.db"
  log:
    stdout="logs/pangenomics/anvi_create_pangenome_storage_{subset}.stdout",
    stderr="logs/pangenomics/anvi_create_pangenome_storage_{subset}.stderr"
  shell:
    """
    anvi-gen-genomes-storage \
      -i {input.txt} \
      -o {output}    \
      > {log.stdout} 2> {log.stderr}
    """

rule all_genome_storages:
  input:
    expand("data/anvio_genomes_storage/{subset}_GENOMES.db",subset=ORDERS),
    "data/anvio_genomes_storage/all_GENOMES.db"

############################### stage 2 create Azolla meta-pangenomes ###############################
ruleorder: create_pangenome_analysis_all > create_pangenome_analysis_subset

rule create_pangenome_analysis_all:
  input:
    "data/anvio_genomes_storage/all_GENOMES.db"
  output:
    "data/anvio_pangenomes/all/all_mcl{mcl}-PAN.db"
  log:
    stdout="logs/pangenomics/anvi_create_pangenome_all_mcl{mcl}.stdout",
    stderr="logs/pangenomics/anvi_create_pangenome_all_mcl{mcl}.stderr"
  threads: 12
  params:
    dir=lambda w: expand ("data/anvio_pangenomes/{subset}/",subset='all')
  shell:
    """
    anvi-pan-genome -g {input}                           \
                    --project-name all_mcl{wildcards.mcl}\
                    --output-dir {params.dir}            \
                    --num-threads  {threads}             \
                    --minbit 0.5                         \
                    --min-occurrence 3                   \
                    --mcl-inflation {wildcards.mcl}      \
                    --exclude-partial-gene-calls         \
                    --enforce-hierarchical-clustering    \
                    --sensitive                          \
    > {log.stdout} 2> {log.stderr}
    """

rule create_pangenome_analysis_subset:
  input:
    storage="data/anvio_genomes_storage/{subset}_GENOMES.db",
    txt="scripts/anvi_genomes_internal_{subset}.tab"
  output:
    "data/anvio_pangenomes/{subset}/{subset}_mcl{mcl}-PAN.db"
  log:
    stdout="logs/pangenomics/anvi_create_pangenome_{subset}_mcl{mcl}.stdout",
    stderr="logs/pangenomics/anvi_create_pangenome_{subset}_mcl{mcl}.stderr"
  threads: 12
  params:
    dir=lambda w: expand ("data/anvio_pangenomes/{subset}/",subset=w.subset),
    name=lambda w: expand ("{subset}_mcl{mcl}",
                            subset = w.subset,
                            mcl    = w.mcl)
  shell:
    """
    genomes=$(cut -f 1 {input.txt} | tail -n +2 | tr '\n' ',' | sed "s/,$//g" )
    anvi-pan-genome -g {input.storage}                   \
                    --project-name {params.name}         \
                    --genome-names $genomes              \
                    --output-dir {params.dir}            \
                    --num-threads  {threads}             \
                    --minbit 0.5                         \
                    --exclude-partial-gene-calls         \
                    --enforce-hierarchical-clustering    \
                    --min-occurrence 2                   \
                    --mcl-inflation {wildcards.mcl}      \
                    --sensitive                          \
    > {log.stdout} 2> {log.stderr}
    """

rule create_pangenome_ANI:
  input:
    pangenome="data/anvio_pangenomes/{subset}/{subset}_mcl{mcl}-PAN.db",
    txt="scripts/anvi_genomes_internal_{subset}.tab"
  output:
    directory("data/anvio_pangenomes/{subset}/{subset}_ANI_mcl{mcl}")
  log:
    stdout="logs/pangenomics/anvi_create_pangenome_ANI_{subset}_mcl{mcl}.stdout",
    stderr="logs/pangenomics/anvi_create_pangenome_ANI_{subset}_mcl{mcl}.stderr"
  threads: 12
  shell:
    """
    anvi-compute-genome-similarity --internal-genomes {input.txt} \
                                   --program pyANI                     \
                                   --output-dir {output}               \
                                   --num-threads {threads}             \
                                   --pan-db {input.pangenome}          \
     > {log.stdout} 2> {log.stderr}
    """

rule collect_PAN_ANI_for_mcl:
  input:
    expand("data/anvio_pangenomes/{subset}/{subset}_ANI_mcl{{mcl}}",subset=ORDERS),
    expand("data/anvio_pangenomes/{subset}/{subset}_ANI_mcl{{mcl}}",subset='all')
  output:
    touch("data/anvio_pangenomes/all_pangenomes_MCL{mcl}.touch")


rule extract_phylogenomic_fasta:
  input:
    pangenome="data/anvio_pangenomes/{subset}/{subset}_mcl{mcl}-PAN.db",
    genomestorage="data/anvio_genomes_storage/{subset}_GENOMES.db"
  output:
    fasta="data/anvio_pangenomes/{subset}_mcl{mcl}_phylogenomic_core.fasta",
    partition="data/anvio_pangenomes/{subset}_mcl{mcl}_phylogenomic_core.partitions"
  log:
    stdout="logs/phylogenomics/anvi_pangenome_fasta_{subset}_mcl{mcl}.stdout",
    stderr="logs/phylogenomics/anvi_pangenome_fasta_{subset}_mcl{mcl}.stderr"
  shell:
    """
    anvi-get-sequences-for-gene-clusters -p {input.pangenome}        \
                                         -g {input.genomestorage}    \
                                         -C default                  \
                                         -b phylogenomic_core        \
                                         --concatenate-gene-clusters \
                                         --report-DNA-sequences --just-do-it  \
                                         --align-with muscle         \
                                         --partition-file {output.partition} \
                                         -o {output.fasta}           \
    > {log.stdout} 2> {log.stderr}
    """

rule phylogenomic_tree_nonpar:
  input:
    fasta="data/anvio_pangenomes/{subset}_mcl{mcl}_phylogenomic_core.fasta",
    partition="data/anvio_pangenomes/{subset}_mcl{mcl}_phylogenomic_core.partitions"
  output:
    tree="data/anvio_pangenomes/{subset}_mcl{mcl}_phylogenomics/{subset}_nonparametric.treefile"
  params:
    pre=lambda w: expand ("data/anvio_pangenomes/{subset}_mcl{mcl}_phylogenomics/{subset}_nonparametric",
                          subset=w.subset,
                          mcl=w.mcl)
  threads: 6
  log:
    stdout="logs/phylogenomics/anvi_phylogenomic_{subset}_mcl{mcl}.stdout",
    stderr="logs/phylogenomics/anvi_phylogenomic_{subset}_mcl{mcl}.stderr"
  shell:
    """
    iqtree -s {input.fasta}     \
           -p {input.partition} \
           -m MFP+MERGE         \
           -b 100               \
           -nt {threads}        \
           -ntmax {threads}     \
           -pre {params.pre}    \
    > {log.stdout} 2> {log.stderr}
    """

rule phylogenomic_tree_big:
  input:
    fasta="data/anvio_pangenomes/{subset}_mcl{mcl}_phylogenomic_core.fasta",
    partition="data/anvio_pangenomes/{subset}_mcl{mcl}_phylogenomic_core.partitions"
  output:
    tree="data/anvio_pangenomes/{subset}_mcl{mcl}_phylogenomics/{subset}_UFbootstrap.treefile"
  params:
    pre=lambda w: expand ("data/anvio_pangenomes/{subset}_mcl{mcl}_phylogenomics/{subset}_UFbootstrap",
                          mcl=w.mcl,
                          subset=w.subset)
  threads: 6
  log:
    stdout="logs/phylogenomics/anvi_phylogenomic_{subset}_mcl{mcl}.stdout",
    stderr="logs/phylogenomics/anvi_phylogenomic_{subset}_mcl{mcl}.stderr"
  shell:
    """
    iqtree -s {input.fasta}     \
           -p {input.partition} \
           -m MFP+MERGE         \
           -bb 2000             \
           -alrt 2000           \
           -nt {threads}        \
           -ntmax {threads}     \
           -pre {params.pre}    \
    > {log.stdout} 2> {log.stderr}
    """
