

# Frises genomiques avec couvertures
Contact: sarah.loughani@univ-tlse3.fr

## Table of Contents
- **[I - plot génomique ](#)**
  - Vérification du tableau final : plot again
  - Création du plot génomique pour tous les gènes
  - Création du plot génomique pour les gènes chevauchants

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
Chemin du dossier de travail bash (cf. script: Bash_script_depth&consensus.md): /bettik/loughans/promReads/

#
Nettoyage de l'espace de travail:

    # clean up
    rm(list=ls())
    graphics.off()

Importation des tableaux: 

    setwd("~/Desktop/Documents-stage/Data_analyses_TIGERISK2/")#_Sarah

    library(readxl)
    #region_list<-read_excel("./INPUT_data/RegionList_candidates_M2_Sarah.xlsx") #JeanM
    region_list<-read_excel("~/Desktop/Documents-stage/Data_analyses_TIGERISK2/Data_StrandNGS/RegionList_candidates_M2_Sarah.xlsx")#_Sarah

    region_list <-  as.data.frame(region_list)
    head(region_list, n=7)


    setwd("~/Desktop/Documents-stage/Data_analyses_TIGERISK2/promoteursSeq_fasta")#_Sarah

    #SNP info YES/NA for the 2500 pb of each promoter
    #prbbyb <- read.table('~/Desktop/Documents-stage/Data_analyses_TIGERISK2/promoters_base_by_base.tab',
    prbbyb <- read.table('./promoters_base_by_base.tab',
                        sep='\t', h=T) # produced with SNPfr_vs_gene_exp_V6g.R
    prbbyb$position= NULL

    colnames(prbbyb)
    prbbyb <-  prbbyb[,-grep('Total_Reads_', colnames(prbbyb) )]
    head(prbbyb)


# I - plot génomique
###### 1) Vérification du tableau final : plot again 
    ## plots genomiques avec frise min/MAJ

    couls <-  cbind(gene_class=unique(prbbyb$gene_class), couleur=c('red', 'darkgreen', 'blue', 'purple'))
    rownames(couls) <- couls[,1]
    couls

    # graphics.off()
    # quartz(h=10, w=15)
    #par(mfrow=c(3,5))


    for (k in 1:length(unique(prbbyb$dowstream_gene))) {
    
    prbk <- prbbyb[prbbyb$dowstream_gene==unique(prbbyb$dowstream_gene)[k],]
    
    head(prbk, n=3)
    
    #seqDFk <- read.table(file = paste('seqDFk', unique(SNPp$dowstream_gene)[k], sep='_'))
    
    plot(x=prbk$Distance_to_prom,
        y=prbk$diffSel_NS,
        xlim=c(-2500,0),
        ylim=c(0,100),
        type = "h",
        col=prbk$couleur,
        main =unique(prbk$dowstream_gene),
        ylab = "diffSel_NS", xlab = "Distance to promoter",
        # yaxs = "i",xaxs = "i",
    )
    
    coul <-  couls[unique(prbk$gene_class), 'couleur']
    mtext(text = unique(prbk$gene_class), col = coul, cex=0.7, line=0.6)
    
    # # Ajout des boules/triangles
    points(x=prbk$Distance_to_prom,
            y=prbk$diffSel_NS, xlim=c(-2500,0),
            col=prbk$col_g, pch=prbk$pch_g)
    
    # # frise des minuscules/majuscules
    # points(x=as.numeric(prbk$Distance_to_prom), y=80-as.numeric(as.logical(prbk$masked))*5, type='l', col='grey')
    
    # if (unique(prbk$gene_strand)=="-") {
    #   #points(x=-as.numeric(seqDFk$position), y=80-as.numeric(as.logical(seqDFk$masked))*5, type='l', col='grey')
    #   #points(x=-as.numeric(prbk$position), y=80-as.numeric(as.logical(prbk$masked))*5, type='l', col='grey')
    #   points(x=as.numeric(prbk$Distance_to_prom), y=80-as.numeric(as.logical(prbk$masked))*5, type='l', col='grey')
    # }else{
    #   #prbk$position <- rev(as.numeric(prbk$position))
    #  positions<- as.numeric(prbk$Distance_to_prom)
    # 
    #   #points(x= -as.numeric(prbk$position), y=80-as.numeric(as.logical(prbk$masked))*5, type='l', col='grey')
    #   points(x= positions, y=80-as.numeric(as.logical(prbk$masked))*5, type='l', col='grey')
    # }
    }

    #dev.copy2pdf(file='./Graphics/genomicPlots_frises_v2.pdf ')

    #" verification conforme :)



###### 2) Création du plot génomique pour tous les gènes

    # Import coverage data ########

    rownames(region_list) <- region_list$`gene ID`

    region_list <-  region_list[unique(prbbyb$dowstream_gene),]

    if(! all(unique(prbbyb$dowstream_gene)==region_list$`gene ID`))
    stop('genes in ≠ order in prbbyb & region_list')

    graphics.off()
    quartz(h=10, w=15)
    par(mfrow=c(3,5))

    setwd("~/Desktop/Documents-stage/Data_analyses_TIGERISK2")#_Sarah

    # list.files('./coverage_reader/')


    #changement de la légende boule triangle : sur prbbyb

    #  calcul des couleurs et formes des boules

    g_pch_col_2 <- function(i, SNPpk){
    
    pch_g_2 <-  NA; col_g_2 <-  NA
    ## pas de boule (pch=NA) si adj_pval>0.05 ou si diffSel_NS≤15%
    ## boule pleine (pch=19) si  0.01 > adj_pval  &  diffSel_NS > 15%, rouge ou noire selon slope
    ## boule grise vide (pch=1) si  0.01 <adj_pval<0.05 & diffSel_NS > 15%
    ## triangle vide (pch=2) si  0.01 <adj_pval<0.05   &   25% <diffSel_NS < 15%, gris
    ## triangle plein (pch=17) si  adj_pval<0.01   &   25% <diffSel_NS < 15%, gris
    
    ## changement: que des boules pleine rouge ou noir si sup à 0.25 & 0.05
    
    ## filtrer sur : SNP$min_reads>20 & SNP$polymorphic & SNP$polyallelic==F
    
    if(!is.na(SNPpk$R2[i])){
        
        ## pas de boule (pch=NA) si adj_pval>0.05 ou si diffSel_NS≤15%
        if(SNPpk$adj_p_val[i] > 0.05 | SNPpk$diffSel_NS[i] < 25){
        pch_g_2 <-  NA; col_g_2 <-  NA
        }
        
        # ## boule pleine (pch=19) si  0.01 > adj_pval  &  diffSel_NS > 15%, rouge ou noire selon slope
        # #if(SNPpk$adj_p_val[i] < 0.01 | SNPpk$diffSel_NS[i] > 15){
        # if(SNPpk$adj_p_val[i] < 0.01 & SNPpk$diffSel_NS[i] > 15){
        #   pch_g <-  19;  if(SNPpk$frSNPp[i]>0) col_g <- "black" else col_g='red' #col_g <-  SNPpk$couleur[i]
        # }
        # 
        # boule pleine (pch=19) si  0.05 > adj_pval  &  diffSel_NS > 25%, rouge ou noire selon slope
        #if(SNPpk$adj_p_val[i] < 0.01 | SNPpk$diffSel_NS[i] > 15){
        if(SNPpk$adj_p_val[i] < 0.05 & SNPpk$diffSel_NS[i] > 25){
        pch_g_2 <-  19;  if(SNPpk$frSNPp[i]>0) col_g_2 <- "black" else col_g_2='red' #col_g <-  SNPpk$couleur[i]
        }
        
        
        # ## boule grise vide (pch=1) si  0.01 <adj_pval<0.05 & diffSel_NS > 15%
        # if(SNPpk$adj_p_val[i] < 0.05 & SNPpk$adj_p_val[i] > 0.01 &
        #    SNPpk$diffSel_NS[i] > 15){
        #   pch_g <-  1; col_g <-  'grey'
        # }
        
        # ## triangle vide (pch=2) si  0.01 <adj_pval<0.05   &   25% <diffSel_NS < 15%, gris
        # if( SNPpk$adj_p_val[i] < 0.05 & SNPpk$adj_p_val[i] > 0.01 &
        #     SNPpk$diffSel_NS[i] > 15 & SNPpk$diffSel_NS[i] < 25){
        #   
        #   pch_g <-  2; col_g <- 'grey'
        # }
        # 
        # 
        # ## triangle plein (pch=17) si  adj_pval<0.01   &   25% <diffSel_NS < 15%, gris
        # if(SNPpk$adj_p_val[i] < 0.01 &
        #    SNPpk$diffSel_NS[i] > 15 & SNPpk$diffSel_NS[i] < 25){
        #   pch_g <-  17; col_g <- 'grey'
        # }
        
    }
    
    ans <-  c(pch_g_2=pch_g_2, col_g_2=col_g_2)
    
    return(ans)
    
    }

    SNPp_Test4 <-  rbind()
    for (k in 1:length(unique(prbbyb$dowstream_gene))) {
    SNPpk <- prbbyb[prbbyb$dowstream_gene==unique(prbbyb$dowstream_gene)[k],]
    
    hans <- t(mapply(g_pch_col_2, i=1:nrow(SNPpk), MoreArgs = list (SNPpk=SNPpk)))
    
    SNPpk <- cbind(SNPpk, hans)
    SNPpk$pch_g_2 <-  as.numeric(SNPpk$pch_g_2)
    
    SNPp_Test4 <-  rbind(SNPp_Test4, SNPpk)
    }
    prbbyb <- SNPp_Test4









    # mise en couleur des SNPs
    prbbyb_Test <- rbind()

    couvPlot <- function(k, xlim=c(-2500, 0), PATH='~/Desktop/Documents-stage/Data_analyses_TIGERISK2/coverage_reader/', PLOT=TRUE){
    
    #genek <-  unique(prbbyb$dowstream_gene)[k]
    genek <-  region_list$`gene ID`[k]
    
    
    genek # "LOC115259486", chr1
    
    orientation <-  region_list$strand[k]
    prom_start <- region_list$`Prom start`[k]
    prom_end <- region_list$`Prom end`[k]
    
    #covk <- cbind()
    covNS <- read.csv(paste0(PATH, genek,'_NS', '.csv'))
    colnames(covNS)[3] <- 'position_AALBF3.v2'
    covk <-  covNS
    
    covSel <- read.csv(paste0(PATH, genek,'_Sel', '.csv'))
    colnames(covSel)[3] <- 'position_AALBF3.v2'
    head(covSel)
    covk <-  merge(covk, covSel[,3:4], by='position_AALBF3.v2', all=TRUE)
    
    covSrun <- read.csv(paste0(PATH, genek,'_Srun', '.csv'))
    colnames(covSrun)[3] <- 'position_AALBF3.v2'
    head(covSrun)
    covk <-  merge(covk, covSrun[,3:4], by='position_AALBF3.v2', all=TRUE)
    
    
    covStL <- read.csv(paste0(PATH, genek,'_StL', '.csv'))
    colnames(covStL)[3] <- 'position_AALBF3.v2'
    head(covStL)
    covk <-  merge(covk, covStL[,3:4], by='position_AALBF3.v2', all=TRUE)
    
    covStP <- read.csv(paste0(PATH, genek,'_StP', '.csv'))
    colnames(covStP)[3] <- 'position_AALBF3.v2'
    head(covStP)
    covk <-  merge(covk, covStP[,3:4], by='position_AALBF3.v2', all=TRUE)
    
    covStSz <- read.csv(paste0(PATH, genek,'_StS', '.csv'))
    colnames(covStSz)[3] <- 'position_AALBF3.v2'
    head(covStSz)
    covk <-  merge(covk, covStSz[,3:4], by='position_AALBF3.v2', all=TRUE)
    
    colors <- ifelse(prbbyb$diffSel_NS < 25 
                    ,'grey'#,'pink')
                    ,(prbbyb$frSNPp < 0)+1)
    
    
    prbbyb$couleur_2 <- colors
    
    prbk <- prbbyb[prbbyb$dowstream_gene==genek,]
    
    if(orientation=='-'){
        
        prbk$position_AALBF3.v2 <- (prom_start+2500 -1):prom_start
        
    } else{ 
        # gene in positive orientation
        
        prbk$position_AALBF3.v2 <- prom_start:(prom_end -1)
    }
    
    for(j in 2:9)
        for(i in 1:nrow(covk))
        if(is.na(covk[i,j]))
            covk[i,j] <-  0
    
    #NOTE:
    range(covk$position_AALBF3.v2) # is larger than
    range(prbk$position_AALBF3.v2) 
    #'  note larger range because reads are 150 nt long (on each strand),
    #'  but are retained by samtools view if as little as 1 pb is in the required
    #'  interval
    #'  
    
    prbk <- merge(prbk, covk[,c("position_AALBF3.v2", 
                                "Total_Reads_NS" , 
                                "Total_Reads_Sel" ,  
                                "Total_Reads_Srun"  ,
                                "Total_Reads_StL" , 
                                "Total_Reads_StP"  ,
                                "Total_Reads_StS" )], by='position_AALBF3.v2', all.x = TRUE)
    if(PLOT)  {
        
        plot(x=prbk$Distance_to_prom,
            y=prbk$diffSel_NS,
            xlim=xlim,
            ylim=c(0,100),
            type = "h",
            col=prbk$couleur_2,
            main =unique(prbk$dowstream_gene),
            ylab = "diffSel_NS", xlab = "Distance to promoter",
            # yaxs = "i",xaxs = "i",
        )
        
        coul <-  couls[unique(prbk$gene_class), 'couleur']
        mtext(text = unique(prbk$gene_class), col = coul, cex=0.7, line=0.6)
        
        # # Ajout des boules/triangles
        # points(x=prbk$Distance_to_prom,
        #        y=prbk$diffSel_NS, xlim=xlim,
        #        col=prbk$col_g, pch=prbk$pch_g)
        points(x=prbk$Distance_to_prom,
            y=prbk$diffSel_NS, xlim=xlim,
            col=prbk$col_g_2, pch=prbk$pch_g_2)
        
        # # frise des minuscules/majuscules
        # points(x=prbk$Distance_to_prom,
        #       y=50-as.numeric(as.logical(prbk$masked))*5, type='l', col='grey')

        # frise des couvertures
        
        range(prbk$Total_Reads_NS) # 103 280
        which(is.na(prbk$Total_Reads_NS))
        range(50+prbk$Total_Reads_NS/10)
        range(60+prbk$Total_Reads_Sel/10)
        range(50-as.numeric(as.logical(prbk$masked))*5)
        
        # points(x= prbk$Distance_to_prom, 
        #        y=50+prbk$Total_Reads_NS/10,  type='l', col='red')
        points(x= prbk$Distance_to_prom, 
            y=50+prbk$Total_Reads_NS/10,  type='l', col='red')
        
        points(x= prbk$Distance_to_prom, 
            y=60+prbk$Total_Reads_Sel/10,  type='l', col='darkgreen')
    
        
        # mtext(text = 'cov_NS',side = 4,at=70, col='red', cex=0.7, las=2)
        # mtext(text = 'cov_Sel',side = 4,at=85, col='darkgreen', cex=0.7, las=2)
        
        axis(4)
        
        abline(h=60, lty=3, col='grey')
        abline(h=50, lty=3, col='grey')
    }
    
    
    return(prbk)
    }

    #nos plot genomique

    graph_table <- data.frame(
    dowstream_gene = character(),
    trait_noir = numeric(),
    trait_rouge = numeric(),
    trait_gris = numeric(),
    points_noir = numeric(),
    points_rouge = numeric(),
    stringsAsFactors = FALSE
    )

    graphics.off()
    quartz(h=10, w=15)
    par(mfrow=c(3,5))


    for(k in 1:nrow(region_list)) {
    
    hans <- couvPlot(k)
    prbbyb_Test <- rbind(prbbyb_Test, hans)
    
    # #on stock les resultats
    # new_row <- data.frame(
    #   dowstream_gene = unique(prbbyb_Test$dowstream_gene),
    #   trait_noir = as.numeric(table(prbbyb_Test$couleur_2)['1']),
    #   trait_rouge = as.numeric(table(prbbyb_Test$couleur_2)['2']),
    #   trait_gris = as.numeric(table(prbbyb_Test$couleur_2)['grey']),
    #   points_noir = as.numeric(table(prbbyb_Test$col_g_2)['black']),
    #   points_rouge = as.numeric(table(prbbyb_Test$col_g_2)['red'])
    # )
    # 
    # graph_table <- rbind(graph_table, new_row)
    }

    dev.copy2pdf(file='./FrisesG_Couv_modify.pdf')

    prbbyb <- prbbyb_Test





    #####

    #' 
    #' k=5
    #' 
    #' zoom <- couvPlot(k=5, xlim=c(-2200, -2000))
    #' which(zoom$Distance_to_prom==-2100) #401
    #' head(zoom[,c(1,2, grep("Total_Reads", colnames(zoom)))])
    #' 
    #' zoom[405:415,c(1,2, 32,grep("Total_Reads", colnames(zoom)))]
    #' #' for this UDPGT, a 2 pb peak in Sel coverage, chr2:236224036
    #' 
    #' 
    #' z2 <- couvPlot(2,  xlim=c(-1800, -1400))
    #' which(z2$Distance_to_prom==-1500)
    #' head(z2[,c(1,2, grep("Total_Reads", colnames(z2)))])
    #' z2[1490:1510,c(1,2, 32,grep("Total_Reads", colnames(z2)))]
    #' # for this GST, a 2 pb peak (chr1; 131604865) in Sel coverage,
    #' #' but also seen in SRun



###### 3) Création du plot génomique pour les gènes chevauchants
 
Pour les gènes chevauchants : LOC...259 et LOC...284 

    SNPp_259 <- prbbyb[prbbyb$dowstream_gene==unique(prbbyb$dowstream_gene)[7],]

    SNPp_284 <- prbbyb[prbbyb$dowstream_gene==unique(prbbyb$dowstream_gene)[8],]
    SNPp_284$diffSel_NS <- -SNPp_284$diffSel_NS

    SNPp_loc <- rbind(SNPp_259, SNPp_284)

    (SNPp_loc$Start)-300350000
    head((SNPp_loc$Start), n=4)
    head((SNPp_loc$Start)-300350000, n=4)
    head((SNPp_loc$Start)/300000000, n=4)

    SNPp_loc$Start<-(SNPp_loc$Start-300350000)

    ## Plot
    graphics.off()
    # quartz()

    plot(x=SNPp_loc$Start, 
        y=SNPp_loc$diffSel_NS,
        #xlim=c(-2500,0),
        ylim=c(-70,70),
        type = "h",
        col=SNPp_loc$couleur_2,
        ylab = "diffSel_NS", xlab = "SNP position (bp)",
        main =paste(unique(prbbyb$dowstream_gene)[7],unique(prbbyb$dowstream_gene)[8]),
        # xaxt="n"
    )

    ## Ajout des boules/triangles
    points(x=SNPp_loc$Start,
        y=SNPp_loc$diffSel_NS, 
        #xlim=c(-2500,0),
        col=SNPp_loc$col_g_2, pch=SNPp_loc$pch_g_2,
        # xaxt="n"
        )

    # calcul de la distance eentre les deux gènes
    max(SNPp_259$position_AALBF3.v2)-min(SNPp_259$position_AALBF3.v2)+1

    max(SNPp_284$position_AALBF3.v2)-min(SNPp_284$position_AALBF3.v2)+1

    min(SNPp_259$position_AALBF3.v2)-min(SNPp_284$position_AALBF3.v2)+1
    max(SNPp_259$position_AALBF3.v2)-max(SNPp_284$position_AALBF3.v2)+1