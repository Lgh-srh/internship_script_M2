# Récupération de la profondeur de couverture &  des sequences consensus de nos 14 promoteurs
Contact: sarah.loughani@univ-tlse3.fr

## Table of Contents
- **[I - Récupération de la profondeur de couverture](#)**
  - Importation des données sur son repertoire bettik
  - Installation du package SAMTOOLS - Environnement Nix
  - On vérifie que les séquences à analyser soient prêtes à être utilisées (indexation)
    - Pour un cas 
    - L'appliquer à toutes les séquences
  - Récupération de la profondeur de couverture (Depth)
    - Pour un cas
    - Mettre au format csv & ajouter une annotation
    - Concatener tous les csv en un seul document
    - Importer les fichiers excel dans mon ordi en gardant une copie dans le repertoire
    - Créer un script pour l'appliquer à tous les genes - scipt final
- **[II- Creation des sequences consensus ](#)**
  - Pour un cas
  - Creation des sequences consensus Reunion

## Samples used
For each mosquitoe line, a pool of 100 mosquitoe were sequenced with the DNBseq technology. Raw sequences (fastq.gz files) are available at XXXX.
The genome version used is this paper is AALBF3.v2.

Mosquitoe lines:
- SRun
- NS
- Sel
- Saint-Paul
- Saint-Louis
- Sainte-Suzanne

# Tools used
samtools-1.18
#
chemin du dossier de travail : /bettik/loughans/promReads/

Avant-propos : 
Dahu est un cluster de calcul pour le calcul haute performance (HPC) et l’analyse de données (DA) (cluster de calcul) 
Dahu donne accès à : 
- Bettik: est un gestionnaire de donnée partagé (scratch) où sont réalisés les calculs (cluster de stockage)
- Mantis: est un nuage de stockage (cluster de stockage) -> ici en iROD
oarsub -I --project mosquitomics permet d'utiliser les ressources de calculs dans un noeud dédié au projet
!/bin/bash indique que le script suivant doit être exécuté en utilisant l'interpréteur de commandes Bash (Bourne Again Shell)

## I - Récupération de la profondeur de couverture
### 1) Importation des données sur son repertoire bettik 

    ssh dahu.access #connexion automatique grâce à la clé access 
    cp /bettik/bonnevje/PromoterAnalysis_albo/promReads/* /bettik/loughans/promReads/

### 2) Installation du package SAMTOOLS - Environnement Nix

Attention: Tous les packages doivent être installés depuis l'environnement Nix et appelés par sa commande nix-env !



    oarsub -I --project mosquitomics
    source /applis/site/nix.sh
    nix-env -i -A samtools # télécharge le package samtools
    nix-env -q # permet de voir les packages téléchargés    
    
> Output:  
hello-2.12.1  
samtools-1.18

Documentation SAMTOOLS: http://www.htslib.org/doc/samtools.html

### À faire à chaque session :
    oarsub -I --project mosquitomics
    source /applis/site/nix.sh

### 3) On vérifie que les séquences à analyser soient prêtes à être utilisées (indexation)
#### a) Pour un cas
Exemple de: LOC109412993_NS


Attention à bien indexer les données avant d'appliquer le consensus ou le vérifier.

On indexe tous les fichiers on applique les "étiquettes" pour les informations associées à nos sequences

"étiquette" : créer un index qui associe des positions spécifiques dans le génome aux emplacements physiques correspondants dans le fichier .BAM

    cd /bettik/loughans/promReads/
    samtools view LOC109412993_NS.bam
    samtools view LOC109412993_NS.bam chr1_pchr:131605860 # on vérifie --> pas indéxé

Dans ce cas on indexe:

    samtools index LOC109412993_NS.bam 
    samtools view LOC109412993_NS.bam 
    # on vérifie:
    chr1_pchr:131605860 # ok et on a bien un fichier .Bai 

si problème dans l'indxation, utiliser un fichier .csi:  

    samtools index -c #fichir en .csi

#### b) L'appliquer à toutes les séquences

    samtools index *.bam

### 4) Récupération de la profondeur de couverture (Depth)

#### a) Pour un cas 

Exemple de: LOC109412993_NS

    samtools depth /bettik/loughans/promReads/LOC109412993_NS.bam

> Example of output:    
chr1_pchr	131606005	1  
chr1_pchr	131606006	1  
chr1_pchr	131606007	1   
chr1_pchr	131606008	1  
chr1_pchr	131606009	1

avec : chromosome  position couverture

Atteention: on n'a pas exactement 2,5kb //
On precise : samtools depth -r chromosome:debut-fin gene[i].bam

    samtools depth -r chr1_pchr:131603368-131605868 /bettik/loughans/promReads/LOC109412993_NS.bam
    samtools depth -r chr1_pchr:131603368-131605868 /bettik/loughans/promReads/LOC109412993_NS.bam |wc -l #ok on a bien 2500
NB: les informations sur les positions pour nos genes se trouvent dans RegionList_candidates_M2_Sarah.xlxs

#### b) Mettre au format csv & ajouter une annotation
    samtools depth input.bam > output.txt
    awk 'BEGIN {FS="\t"; OFS=","} {print $1, $2, $3, $4}' output.txt > output.csv
    
    samtools depth -r chr1_pchr:131603368-131605868 /bettik/loughans/promReads/LOC109412993_NS.bam > /bettik/loughans/work_directory/coverage/LOC109412993_NS.txt
    
    awk 'BEGIN {FS="\t"; OFS=","} {$1="LOC109412993_NS"; print}' /bettik/loughans/work_directory/coverage/LOC109412993_NS.txt > /bettik/loughans/work_directory/coverage/LOC109412993_NS_with_gene.csv

#### c) Concatener tous les csv en un seul document
    awk 'FNR==1 && NR!=1 {next;}{print}' "${directory}"*.csv > "$output_file"
    
    paste -d',' directory > output

à rédiger sous cette forme : 
> paste -d',' /bettik/loughans/work_directory/coverage/*.csv > /bettik/loughans/work_directory/coverage/all_cov_assmbly_cb.csv

#### d) Importer les fichiers excel dans mon ordi en gardant une copie dans le repertoire
    monlogin@monpc:~$ scp f-dahu.ciment:my_file my_file_destination (copier depuis le reprtoire dahu vers mon ordinateur) // ouvrir un autre terminal
    ou 
    monlogin@monpc:~$ scp my_file f-dahu.ciment:my_file_destination (copier vers)
    
    scp f-dahu.access:/bettik/loughans/work_directory/coverage/LOC109412993_NS_with_gene.csv /Users/princess/Downloads/coverage_reader 
    
    scp f-dahu.access:/bettik/loughans/work_directory/coverage/*.csv /Users/princess/Downloads/coverage_reader #apppliquer a tous les fichier .csv

#### e) Créer un script pour l'appliquer à tous les genes - scipt final

#### Indexation de nos fichiers bam
    samtools index -c -M /bettik/loughans/promReads/*.bam # mettre -M pour pouvoir utiliser "*"
    
    nano /bettik/loughans/work_directory/nano_scripts/Depth_writter
#
Script nano :  

    ### Depth_writter ####
    #!/bin/bash
    #OAR -n Depth_worker_mosquitom_coverage
    #OAR --project mosquitomics
    ##OAR -t besteffort
    #OAR -l core=10,walltime=01:00:00
    #OAR -O OAR.%jobid%mapping.stdout
    #OAR -E OAR.%jobid%mapping.stderr

    source /applis/site/nix.sh
    #samtools index -c -M /bettik/loughans/promReads/*.bam #enlever le"#" si besoin

    genes=("LOC115259486" "LOC109412993" "LOC115257289" "LOC109413318" "LOC109408672" "LOC109410283" "LOC109410284" "LOC109426259" "LOC115259235"  "LOC109421559" "LOC109409164" "LOC109407008" "LOC115262073" "LOC109429691")
    lineages=("NS" "Sel" "Srun" "StL" "StP" "StS")

    for gene in "${genes[@]}"; do
        for lineage in "${lineages[@]}"; do

            gene_with_lineage="${gene}_${lineage}"
            header="dowstream_gene,Chromosome,start,Total_Reads_${lineage}"

            # Commande samtools depth avec redirection vers un fichier texte
            samtools depth "/bettik/loughans/promReads/${gene_with_lineage}.bam" > "/bettik/loughans/work_directory/coverage/${gene_with_lineage}.txt"

            # Commande awk pour ajouter le nom du gene comme premiere colonne et redirection vers un fichier CSV
            awk -v gene="$gene" -v header="$header" 'BEGIN {OFS=","; print header} {gsub("\t", ",", $0); print gene,$0}' "/bettik/loughans/work_directory/coverage/${gene_with_lineage}.txt" > "/bettik/loughans/work_directory/coverage/${gene_with_lineage}.csv"
        done
    done
#

    > bash /bettik/loughans/work_directory/nano_scripts/Depth_writter
    > ls |wc -l # pour vérifier: 168 

#### Rassembler les csv en un seul (par la ligne ou par la colonne) (bonus)

    awk 'FNR==1 && NR!=1 {next;}{print}' /bettik/loughans/work_directory/coverage/*.csv > /bettik/loughans/work_directory/coverage/all_cov_assmbly_rb.csv #concatène en un seul fichier csv
    
    paste -d',' /bettik/loughans/work_directory/coverage/*.csv > /bettik/loughans/work_directory/coverage/all_cov_assmbly_cb.csv

#### importation des donnees dans mon dossier personnel
    scp f-dahu.access:/bettik/loughans/work_directory/coverage/*.csv /Users/princess/Downloads/coverage_reader #apppliquer à tous les fichier .csv
    
    scp f-dahu.access:/bettik/loughans/work_directory/coverage/all_cov_assmbly_cb.csv /Users/princess/Downloads/coverage_reader

=> resultat de l'utilisation de ces données de couvertures brutes dans Rstudio pour le plot_genomique: R_studio_script_dataWork

## II- Creation des sequences consensus 
#### 1) Pour un cas

    samtools consensus /bettik/loughans/promReads/LOC115259486_NS.bam -o /bettik/loughans/work_directory/consensus_seq/LOC115259486_NS.fq

Documentation
http://www.htslib.org/doc/samtools-consensus.html : 
> samtools consensus -f fastq in.bam -o cons.fq  
> samtools consensus in.bam -o cons.fq  
> consensus samtools [-saAMq] [-r région] [-f format] [-l line-len] [-d min-depth] [-C cutoff] [-c call-fract] [-H het-fract] in.bam

#### 2) Creation des sequences consensus Reunion

Celle-ci s fait sur toutes les lignees & tous les gènes ensemble (ceux qui serviront pour le résultat final)

    nano /bettik/loughans/work_directory/nano_scripts/Consensus_Reu_writer
#
Script nano :

    #### Consensus_Reu_writer ####
    #!/bin/bash
    #OAR -n Consensus_Reu_worker_mosquitom
    #OAR --project mosquitomics
    ##OAR -t besteffort
    #OAR -l core=10,walltime=01:00:00
    #OAR -O OAR.%jobid%mapping.stdout
    #OAR -E OAR.%jobid%mapping.stderr


    ## Liste de nos 14 genes
    genes=("LOC115259486" "LOC109412993" "LOC115257289" "LOC109413318" "LOC109408672" "LOC109410283" "LOC109410284" "LOC109426259" "LOC115259235"  "LOC109421559" "LOC109409164" "LOC109407008" "LOC115262073" "LOC109429691")
    ## Coordonnees des genes : chr1_pchr:debut-fin
    ###Liste de nos coordonnees
    coordinate_start=("65321351" "131603368" "333648863" "357956239" "236223628" "300350485" "300355022" "300355290" "473844199" "525917704" "526016124" "547743964" "547811438" "554432583")
    coordinate_end=("65323851" "131605868" "333651363" "357958739" "236226128" "300352985" "300357522" "300357790" "473846699" "525920204" "526018624" "547746464" "547813938" "554435083")
    ###Liste du chromosome concerne
    chromosomes=("chr1_pchr" "chr1_pchr" "chr1_pchr" "chr1_pchr" "chr2_pchr" "chr2_pchr" "chr2_pchr" "chr2_pchr" "chr2_pchr" "chr2_pchr" "chr2_pchr" "chr2_pchr" "chr2_pchr" "chr2_pchr")
    echo " "
    # load du package sur le cluster:samtools
    source /applis/site/nix.sh

    # On cree un vecteur vide
    declare -a gene_coordinates

    #On realise la boucle
    for (( i=0; i<${#genes[@]}; i++ )); do 

        gene="${genes[i]}"

        #lignee"NS" "Sel" "Srun" "StL" "StP" "StS"
        
        samtools merge -o /bettik/loughans/work_directory/merged_seq/${gene}_ReuMerged.bam /bettik/loughans/promReads/${gene}_NS.bam /bettik/loughans/promReads/${gene}_Sel.bam /bettik/loughans/promReads/${gene}_Srun.bam /bettik/loughans/promReads/${gene}_StL.bam /bettik/loughans/promReads/${gene}_StP.bam /bettik/loughans/promReads/${gene}_StS.bam

        echo " "
        echo "_echo#! OK: creation de la sequence merge de ${gene}"
        echo " "

        
        coordinate_start="${coordinate_start[i]}"
        coordinate_end="${coordinate_end[i]}"
        chromosomes="${chromosomes[i]}"

        gene_coordinates="${chromosomes}:${coordinate_start}-${coordinate_end}"

        samtools index -c -M /bettik/loughans/work_directory/merged_seq/${gene}_ReuMerged.bam

        echo " "
        echo "_echo#! OK: creation de l'index du gene  ${gene} .csi"
        echo " "

        #modifier les seuils par defaut et la profondeur - consensus
        samtools consensus -A --mode simple --no-use-qual --call-fract 0.25 --min-depth 60 -r "$gene_coordinates" /bettik/loughans/work_directory/merged_seq/${gene}_ReuMerged.bam -o /bettik/loughans/work_directory/consensus_seq/${gene}_ReuMerged_consensus.fq

        echo " "
        echo "_echo#! OK: creation du consensus ${gene}"
        echo " "

    done

    echo "       OOOOooo work is done \_(O,O)_/  oooOOOO            "
#
    bash /bettik/loughans/work_directory/nano_scripts/Consensus_Reu_writer

    cd /bettik/loughans/work_directory/merged_seq
    cd /bettik/loughans/work_directory/consensus_seq

On verifie la taille d'une sequence:

    samtools view /bettik/loughans/work_directory/consensus_seq |wc -m
    
    samtools view /bettik/loughans/work_directory/consensus_seq |grep 'N'|wc -m

    scp f-dahu.access:/bettik/loughans/work_directory/consensus_seq/* /Users/princess/Downloads/consensus_c25_d60

=> resultat de l'utilisation de ces données dee consensus dans Rstudio pour l'édition des séquences: R_studio_script_editing

Documentation:

http://www.htslib.org/doc/samtools-merge.html 
http://www.htslib.org/doc/samtools-consensus.html
http://www.htslib.org/doc/samtools-index.html 



