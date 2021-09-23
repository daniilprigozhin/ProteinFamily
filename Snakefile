import pandas as pd
import pprint as pp
from snakemake.utils import validate, min_version
##### set minimum snakemake version #####
min_version("5.1.2")
configfile: "config.yaml"
samples = pd.read_table(config["samples"]).set_index("sample", drop=False)
models = pd.read_table(config["models"]).set_index("model", drop=False)
nbarc_clades = pd.read_table(config["nbarc_clades"]).set_index("clade", drop=False)
r1_nbarc_clades = pd.read_table(config["r1_nbarc_clades"]).set_index("clade", drop=False)
r2_nbarc_clades = pd.read_table(config["r2_nbarc_clades"]).set_index("clade", drop=False)
r3_nbarc_clades = pd.read_table(config["r3_nbarc_clades"]).set_index("clade", drop=False)

rule all:
    input:
        ## for name_correction
        #expand("Proteomes/corr.{line}.protein.fasta", line=samples["sample"])
        ## for hmm_search
        # expand("HMM_search_{model}/hmmsearch.{model}.{sample}.out", sample=samples["sample"], model = "pbNB-ARC")
        # expand("HMM_search_{model}/hmmsearch.{model}.{sample}.out", sample=samples["sample"], model = "RPW8")
        ## for merge_sto
        #expand("HMM_search_{model}/all_samples.{model}.sto", model = models["model"])
        ## for sto_to_afa_with_trim
        #expand("HMM_search_{model}/all_samples.{model}.afa", model = "pbNB-ARC")
        #expand("HMM_search_{model}/all_samples.{model}.afa", model = "RPW8")
        ## for afa_to_fasta
        #expand("HMM_search_{model}/{sample}_{model}.fasta",sample=samples["sample"], model = models["model"])
        ## for combine_fasta_hits
        #expand("HMM_search_{model}/all_samples.{model}.fasta", model = models["model"])
        ## for filter_isoforms:
        #expand("HMM_search_{model}/all_samples.{model}.selected.afa", model = models["model"])
        ##for build_raxml_tree
        #expand("RAxML_tree_{model}/RAxML_bipartitionsBranchLabels.{model}.Raxml.out", model = "pbNB-ARC")
        ## for hmmscan_subset
        #expand("HMM_search_{model}/all_samples.{model}.Pfam_scan.out", model = models["model"])
        ## for reduce_pfam
        #expand("HMM_search_{model}/all_samples.{model}.Pfam_scan_reduced.out", model = models["model"])
        ## for annotate_itol_domains
        #expand("Annotation/{model}.Bigtree.Domains.iTOL.txt", model = "pbNB-ARC")
        ## for split_big_tree
        #expand("RAxML_tree_{model}/Initial_Clades.{model}.tsv", model = "pbNB-ARC")
        ## for print_init_clades
        #expand("RAxML_tree_{model}/Initial_Clades.{model}.list", model = "pbNB-ARC")
        ## for clade fasta files
        #expand("RAxML_tree_{model}/{stage}/{clade}.fasta", stage = "Initial_Clades", model = "pbNB-ARC", clade = nbarc_clades["clade"])
        ## for align_mafft
        #expand("RAxML_tree_{model}/{stage}/{clade}.afa",stage = "Initial_Clades",model = "pbNB-ARC", clade = nbarc_clades["clade"])
        ## for build_raxml_tree
        # expand("RAxML_tree_{model}/{stage}/RAxML_bipartitionsBranchLabels.{clade}.Raxml.out",stage = "Initial_Clades", model = "pbNB-ARC", clade = nbarc_clades["clade"])
        ## for refinement_1
        #expand("RAxML_tree_{model}/Refinement_1.{model}.BigTable.RData", model = "pbNB-ARC")
        ## for print_refined_clades
        #expand("RAxML_tree_{model}/Refinement_1.{model}.Clades.list", model = "pbNB-ARC")
        ## for clade fasta files
        #expand("RAxML_tree_{model}/{stage}/{clade}.fasta", stage = "Clade_Refinement_1", model = "pbNB-ARC", clade = r1_nbarc_clades["clade"])
        ## for align_mafft
        #expand("RAxML_tree_{model}/{stage}/{clade}.afa",stage = "Clade_Refinement_1",model = "pbNB-ARC", clade = r1_nbarc_clades["clade"])
        ## for build_raxml_tree
        #expand("RAxML_tree_{model}/{stage}/RAxML_bipartitionsBranchLabels.{clade}.Raxml.out",stage = "Clade_Refinement_1", model = "pbNB-ARC", clade = r1_nbarc_clades["clade"])
        ## for refinement_2
        #expand("RAxML_tree_{model}/Refinement_2.{model}.BigTable.RData", model = "pbNB-ARC")
        ## for print_refined_clades_2
        #expand("RAxML_tree_{model}/Refinement_2.{model}.Clades.list", model = "pbNB-ARC")
        ## for align_mafft
        #expand("RAxML_tree_{model}/{stage}/{clade}.afa",stage = "Clade_Refinement_2",model = "pbNB-ARC", clade = r2_nbarc_clades["clade"])
        ## for build_raxml_tree
        #expand("RAxML_tree_{model}/{stage}/RAxML_bipartitionsBranchLabels.{clade}.Raxml.out",stage = "Clade_Refinement_2", model = "pbNB-ARC", clade = r2_nbarc_clades["clade"])
        ## for refinement_3
        #expand("RAxML_tree_{model}/Refinement_3.{model}.BigTable.RData", model = "pbNB-ARC")
        ## for print_refined_clades_3
        #expand("RAxML_tree_{model}/Refinement_3.{model}.Clades.list", model = "pbNB-ARC")
        ## for build_raxml_tree
        #expand("RAxML_tree_{model}/{stage}/RAxML_bipartitionsBranchLabels.{clade}.Raxml.out",stage = "Clade_Refinement_3", model = "pbNB-ARC", clade = r3_nbarc_clades["clade"])
        ## for refinement_4
        expand("RAxML_tree_{model}/Refinement_4.{model}.BigTable.RData", model = "pbNB-ARC")

#### It works much better to have the rules in reverse chronological order.
#### This way the new rules gets written directly under "rule all" used to set tasks.

## Refinement converged after 3 rounds. No good cut nodes were proposed in refinement 4.
rule refinement_4:
    threads:
        8
    input:
        afa = expand("RAxML_tree_{{model}}/Clade_Refinement_3/{clade}.afa", clade = r3_nbarc_clades["clade"]),
        tree = expand("RAxML_tree_{{model}}/Clade_Refinement_3/RAxML_bipartitionsBranchLabels.{clade}.Raxml.out", clade = r3_nbarc_clades["clade"]),
        annotation = "NLR_Map.csv",
    output:
        tsv = "RAxML_tree_{model}/Refinement_4.{model}.tsv",
        BigTable = "RAxML_tree_{model}/Refinement_4.{model}.BigTable.RData",
        cladestar = "RAxML_tree_{model}/Refinement_4.{model}.iTOL.cladestar.txt",
    params:
        MinGapFraction = 1,
        MinGapBlockWidth = 1,
        hvSiteEntCutoff = 1.5,
    script:
        "scripts/Zm_NLRome_Refinement_sm.R"

rule print_refined_clades_3:
    input:
        tsv = "RAxML_tree_{model}/Refinement_3.{model}.edited.tsv",
        BigTable = "RAxML_tree_{model}/Refinement_3.{model}.BigTable.RData",
    output:
        list = "RAxML_tree_{model}/Refinement_3.{model}.Clades.list",
        dir = directory("RAxML_tree_{model}/Clade_Refinement_3/")
    script:
        "scripts/Zm_NLRome_PrintRefinedClades.R"

rule refinement_3:
    threads:
        8
    input:
        afa = expand("RAxML_tree_{{model}}/Clade_Refinement_2/{clade}.afa", clade = r2_nbarc_clades["clade"]),
        tree = expand("RAxML_tree_{{model}}/Clade_Refinement_2/RAxML_bipartitionsBranchLabels.{clade}.Raxml.out", clade = r2_nbarc_clades["clade"]),
        annotation = "NLR_Map.csv",
    output:
        tsv = "RAxML_tree_{model}/Refinement_3.{model}.tsv",
        BigTable = "RAxML_tree_{model}/Refinement_3.{model}.BigTable.RData",
        cladestar = "RAxML_tree_{model}/Refinement_3.{model}.iTOL.cladestar.txt",
    params:
        MinGapFraction = 1,
        MinGapBlockWidth = 1,
        hvSiteEntCutoff = 1.5,
    script:
        "scripts/Zm_NLRome_Refinement_sm.R"

rule print_refined_clades_2:
    input:
        tsv = "RAxML_tree_{model}/Refinement_2.{model}.edited.tsv",
        BigTable = "RAxML_tree_{model}/Refinement_2.{model}.BigTable.RData",
    output:
        list = "RAxML_tree_{model}/Refinement_2.{model}.Clades.list",
        dir = directory("RAxML_tree_{model}/Clade_Refinement_2/")
    script:
        "scripts/Zm_NLRome_PrintRefinedClades.R"

rule refinement_2:
    threads:
        8
    input:
        afa = expand("RAxML_tree_{{model}}/Clade_Refinement_1/{clade}.afa", clade = r1_nbarc_clades["clade"]),
        tree = expand("RAxML_tree_{{model}}/Clade_Refinement_1/RAxML_bipartitionsBranchLabels.{clade}.Raxml.out", clade = r1_nbarc_clades["clade"]),
        annotation = "NLR_Map.csv",
    output:
        tsv = "RAxML_tree_{model}/Refinement_2.{model}.tsv",
        BigTable = "RAxML_tree_{model}/Refinement_2.{model}.BigTable.RData",
        cladestar = "RAxML_tree_{model}/Refinement_2.{model}.iTOL.cladestar.txt",
    params:
        MinGapFraction = 1,
        MinGapBlockWidth = 1,
        hvSiteEntCutoff = 1.5,
    script:
        "scripts/Zm_NLRome_Refinement_sm.R"

rule print_refined_clades:
    input:
        tsv = "RAxML_tree_{model}/Refinement_1.{model}.edited.tsv",
        BigTable = "RAxML_tree_{model}/Refinement_1.{model}.BigTable.RData",
    output:
        list = "RAxML_tree_{model}/Refinement_1.{model}.Clades.list",
        dir = directory("RAxML_tree_{model}/Clade_Refinement_1/")
    script:
        "scripts/Zm_NLRome_PrintRefinedClades.R"

rule refinement_1:
    threads:
        8
    input:
        afa = expand("RAxML_tree_{{model}}/Initial_Clades/{clade}.afa", clade = nbarc_clades["clade"]),
        tree = expand("RAxML_tree_{{model}}/Initial_Clades/RAxML_bipartitionsBranchLabels.{clade}.Raxml.out", clade = nbarc_clades["clade"]),
        annotation = "NLR_Map.csv",
    output:
        tsv = "RAxML_tree_{model}/Refinement_1.{model}.tsv",
        BigTable = "RAxML_tree_{model}/Refinement_1.{model}.BigTable.RData",
        cladestar = "RAxML_tree_{model}/Refinement_1.{model}.iTOL.cladestar.txt",
    params:
        MinGapFraction = 0.9,
        MinGapBlockWidth = 1,
        hvSiteEntCutoff = 1.5,
    script:
        "scripts/Zm_NLRome_Refinement_sm.R"

rule build_raxml_clade_trees:
    input:
        afa = "RAxML_tree_{model}/{stage}/{clade}.afa",
    output:
         bpbl = "RAxML_tree_{model}/{stage}/RAxML_bipartitionsBranchLabels.{clade}.Raxml.out",
    threads: 4
    ### Sequential RAxML (single thread) can be called at ~/Installers/standard-RAxML/raxmlHPC-AVX
    ### On iMac use ~/Installers/standard-RAxML/raxmlHPC-PTHREADS-AVX; \
    ### On Air use /Users/prigozhin/Installers/Git/standard-RAxML/raxmlHPC-PTHREADS-SSE3 \

    shell:
        "cd RAxML_tree_{wildcards.model}/{wildcards.stage};\
         /Users/prigozhin/Installers/Git/standard-RAxML/raxmlHPC-PTHREADS-SSE3 \
         -T {threads} \
         -n {wildcards.clade}.Raxml.out \
         -f a -x 12345 -p 12345 -# 100 -m PROTCATJTT \
         -s ../../{input.afa}"

rule align_mafft:
    input:
        fasta = "RAxML_tree_{model}/{stage}/{clade}.fasta"
    output:
        afa = "RAxML_tree_{model}/{stage}/{clade}.afa"
    shell:
        "mafft --genafpair --maxiterate 100 {input.fasta} > {output.afa}"

rule get_fasta:
    input:
        list = "RAxML_tree_{model}/{stage}/{clade}.txt",
        fasta = "HMM_search_{model}/all_samples.{model}.fasta"
    output:
        fasta = "RAxML_tree_{model}/{stage}/{clade}.fasta",
    run:
        shell("~/Dropbox/DP_KVK/Scripts/K-get_fasta_from_ids.pl \
                -i {input.list} \
                -f {input.fasta} > \
                {output.fasta}")

rule print_init_clades:
    input:
        tree = "RAxML_tree_{model}/RAxML_bipartitionsBranchLabels.{model}.Raxml.out",
        tsv = "RAxML_tree_{model}/Initial_Clades.{model}.tsv",
    output:
        list = "RAxML_tree_{model}/Initial_Clades.{model}.list",
        dir = directory("RAxML_tree_{model}/Initial_Clades/")
    script:
        "scripts/Zm_NLRome_PrintInitClades.R"

rule split_big_tree:
    threads:
        8
    input:
        tree = "RAxML_tree_{model}/RAxML_bipartitionsBranchLabels.{model}.Raxml.out"
    output:
        tsv = "RAxML_tree_{model}/Initial_Clades.{model}.tsv",
        cladestrip = "RAxML_tree_{model}/Initial_Clades.{model}.iTOL.cladestrip.txt",
        cladestar = "RAxML_tree_{model}/Initial_Clades.{model}.iTOL.cladestar.txt",
    params:
        min_seq = 25,
        max_seq = 500,
    script:
        "scripts/Zm_NLRome_InitialAssignment_sm.R"

rule annotate_itol_domains:
    input:
        pfam = "HMM_search_{model}/all_samples.{model}.Pfam_scan_reduced.out",
        lrrpred = "Annotation/all_samples.LRRpred.tsv",
        fasta = expand("Proteomes/corr.{sample}.protein.fasta",sample = samples["sample"]),
        afa = "HMM_search_{model}/all_samples.{model}.afa"
    output:
        domains = "Annotation/all_samples.{model}.iTOLdomains.tsv",
        bigtree = "Annotation/{model}.Bigtree.Domains.iTOL.txt",
        smalltree = "Annotation/{model}.Smalltree.Domains.iTOL.tsv",
    script:
        "scripts/DomainDiagrams_sm.R"

rule extract_lrrpred:
    input:
        lrrpred = "Annotation/LRRpredictor",
    output:
        lrrpred = "Annotation/all_samples.LRRpred.tsv",
    script:
        "scripts/Extract_LRRpred_sm.R"

rule reduce_pfam:
    params:
        e_val_max = 1e-3,
        hmm_frac_min = 0.3,
        max_overlap = 10
    input:
        pfam = "HMM_search_{model}/all_samples.{model}.Pfam_scan.out"
    output:
        pfam = "HMM_search_{model}/all_samples.{model}.Pfam_scan_reduced.out"
    script:
        "scripts/reduce_pfam.R"

rule hmmsearch_subset:
    input:
        fasta = "HMM_search_{model}/all_samples.{model}.fasta",
    output:
        tblout = "HMM_search_{model}/all_samples.{model}.Pfam_scan.out"
    threads:
        4
    shell:
        "hmmsearch \
        --cpu {threads} \
        --domtblout {output.tblout} \
        HMM_models/Pfam-A_Plus.hmm \
        {input.fasta}"

# rule build_raxml_tree:
#     input:
#         afa = "HMM_search_{model}/all_samples.{model}.selected.afa",
#     output:
#          bpbl = protected("RAxML_tree_{model}/RAxML_bipartitionsBranchLabels.{model}.Raxml.out"),
#     threads: 8
#     shell:
#         "cd RAxML_tree_{wildcards.model};\
#          ~/Installers/standard-RAxML/raxmlHPC-PTHREADS-AVX \
#          -T {threads} \
#          -n {wildcards.model}.Raxml.out \
#          -f a -x 12345 -p 12345 -# 100 -m PROTCATJTT \
#          -s ../{input.afa}"
#
# rule filter_isoforms:
#     input:
#         fasta = "HMM_search_{model}/all_samples.{model}.afa",
#     output:
#         all = "HMM_search_{model}/all_samples.{model}.allgene.list",
#         select = "HMM_search_{model}/all_samples.{model}.sel_gene.list",
#         fasta = "HMM_search_{model}/all_samples.{model}.selected.afa"
#     run:
#         shell("grep '>' {input.fasta}|tr -d '>' >{output.all}")
#         shell("cat {output.all} |grep -e '-R1/' >{output.select}")
#         shell("cat {output.all} |grep -e '_1/' >>{output.select}")
#         shell("~/Dropbox/DP_KVK/Scripts/K-get_fasta_from_ids.pl \
#                 -i {output.select} \
#                 -f {input.fasta} > \
#                 {output.fasta}")
#
# rule combine_fasta_hits:
#     input:
#         fasta = expand("HMM_search_{{model}}/{lines}_{{model}}.fasta", lines=samples["sample"]),
#     output:
#         fasta = "HMM_search_{model}/all_samples.{model}.fasta",
#     run:
#         shell("cat {input.fasta} > {output.fasta}")
#
#
# rule afa_to_fasta:
#     input:
#         fasta = "Proteomes/corr.{sample}.protein.fasta",
#         afa = "HMM_search_{model}/all_samples.{model}.afa",
#     output:
#         fasta = "HMM_search_{model}/{sample}_{model}.fasta",
#         list = "HMM_search_{model}/{sample}_{model}.list",
#     run:
#         shell("echo " " > {output.list}")
#         shell("grep '>' {input.afa} | tr -d '>' | cut -d '/' -f 1 > {output.list} || true")
#         shell("~/Dropbox/DP_KVK/Scripts/K-get_fasta_from_ids.pl \
#                 -i {output.list} \
#                 -f {input.fasta} > \
#                 {output.fasta}")
#
# rule sto_to_afa_with_trim:
#     input:
#         sto = "HMM_search_{model}/all_samples.{model}.sto",
#     output:
#         afa = "HMM_search_{model}/all_samples.{model}.afa",
#     params:
#         length = 237
#         #length = 100
#         #length = 30
#     shell:
#         """
#         tr a-z - < {input.sto} |\
#         esl-reformat --mingap afa -|\
#         esl-alimanip --lmin {params.length} - |\
#         esl-reformat afa -|\
#         cut -d ' ' -f 1 |tr -d ' '>{output.afa}
#         """
#
# rule merge_sto:
#     input:
#         sto = expand("HMM_search_{{model}}/hmmsearch.{{model}}.{lines}.sto",lines = samples["sample"]),
#     output:
#         list = temp("{model}.sto.list"),
#         sto = "HMM_search_{model}/all_samples.{model}.sto",
#     run:
#         shell("ls {input.sto} > {output.list}")
#         shell("esl-alimerge --list {output.list} > {output.sto}")
#
# rule hmmsearch:
#     params:
#         e=1e-5
#         #e=10
#     input:
#         fasta = "Proteomes/corr.{sample}.protein.fasta",
#         hmm = ancient("HMM_models/{model}.hmm")
#     output:
#         out = "HMM_search_{model}/hmmsearch.{model}.{sample}.out",
#         sto = "HMM_search_{model}/hmmsearch.{model}.{sample}.sto"
#     shell:
#         "hmmsearch --cpu 1 -E {params.e} -A {output.sto} --domtblout {output.out} {input.hmm} {input.fasta}"
#
# rule hmmfetch:
#     output:
#         hmm = "HMM_models/{model}.hmm"
#     shell:
#         "hmmfetch ~/BioInf/Pfam/Pfam-A.hmm {wildcards.model} > {output.hmm}"
#
# rule name_correction:
#     input:
#         fasta = "Proteomes/{sample}.aa.fa"
#
#     output:
#         fasta = "Proteomes/corr.{sample}.protein.fasta",
#
#     shell:
#         "cut -d ' ' -f 1 {input.fasta}|\
#         tr a-z A-Z|tr -d ' '|\
#         tr '.' '_' > {output.fasta}"
