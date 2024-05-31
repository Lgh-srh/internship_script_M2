
# Filtrage des données: la fréquence d'un SNP permet-il de prédire le niveau d'expression d'un gène ? (modèle linéaire)
Contact: sarah.loughani@univ-tlse3.fr

## Table of Contents
- **[Introduction ](#)**
  - Adapter les colonnes à R
  - Dedoublement LOC109410284 / LOC109426259
  - Distance SNP to transcription start 
  - Donner une couleur pour la catégorie de couleur
- **[I - Filtrage des données](#)**
  - Filtrage sur les couvertures
  - Filtrage des SNPs non polymorphes
  - Filtrage pval
  - Ajout diffSel_NS
- **[II - plot à partir du filtrage des données](#)**
  - biplots adj-pval ~ diffSel_NS avec les coloriages necessaires, x14 genes
  - plots  genomiques frises 1/2
  - plots  genomiques frises 2/2
  - Filtrage polyallelique :Reduction a une seule ligne des SNPs polyalleliques
  - Tableau Génomique


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
    #nb de reads RNA normalisés par gene
    #normReads <- read_xlsx('./INPUT_data/Data_Assembly_filtered.xlsx', sheet='Gene_RNA_norm_counts_filtered')#_JeanM
    normReads <- read_excel("Clean_ALL_Data_Assembly_filtered.xlsx",sheet =3)#_Sarah

    # head(normReads)
    normReads <-  as.data.frame(normReads)
    # head(normReads)

    # Les SNPs des promoteurs
    # SNPs0 <- read_xlsx('./INPUT_data/Data_Assembly_filtered.xlsx', sheet='SNP_freqs_filtered') #_JeanM
    SNPs0 <- read_excel("Clean_ALL_Data_Assembly_filtered.xlsx",sheet =2) #_Sarah
    #SNPs0 <- read_excel("Data_Assembly_filtered.xlsx",sheet =2) #_JeanM
    # head(SNPs0)
    # dim(SNPs0)

    SNPs0 <- as.data.frame(SNPs0)
    head(SNPs0)
    colnames(SNPs0)
    # ne garder que les colonnes d'interet
    SNP <- SNPs0[,c(1:9,
                    grep("Supporting", colnames(SNPs0)),
                    grep("Total Reads", colnames(SNPs0))
    )]

    head(SNP)

# Introduction
Paramètrer les tableaux:

###### 1) Adapter les noms des colonnes a R
    colnames(SNP) <- gsub(' ', "_", colnames(SNP))
    colnames(SNP) <- gsub('%', "PCT", colnames(SNP))

    colnames(SNP)

###### 2) Dedoublement LOC109410284 / LOC109426259 
    SNP0 <- SNP
    # operation dedoublement des lignes
    udpgt0 <- SNP[SNP$dowstream_gene=='LOC109410284 / LOC109426259',]

    head(udpgt0)

    udpgt1 <-  udpgt0
    udpgt2 <-  udpgt0
    udpgt1$dowstream_gene <- 'LOC109410284' # -
    udpgt1$gene_strand<-"-"
    udpgt2$dowstream_gene <- 'LOC109426259' # +
    udpgt2$gene_strand<-"+"

    SNP <- SNP[-which(SNP$dowstream_gene=='LOC109410284 / LOC109426259'),]
    SNP <-  rbind(SNP, udpgt1, udpgt2)

###### 3) Distance SNP to transcription start 

    setwd("~/Desktop/Documents-stage/Data_analyses_TIGERISK2/Data_StrandNGS/") #_Sarah

    #region_list<-read_excel("./INPUT_data/RegionList_candidates_M2_Sarah.xlsx")#_JM

    region_list<-read_excel("RegionList_candidates_M2_Sarah.xlsx")#_Sarah

    # principe: si c'est - alors on fait GeneStart= colEnd sinon GeneStart= colStart
    # et donc on aura GeneStart-position_SNP

    # dprom <-  function(i){
    #   SNP_start<-SNP[i,"Start"]
    #   
    #   if (SNP[i,'gene_strand']=='-'){
    #     # selectionne la position Gen end
    #     Gen_start<-as.numeric(region_list[region_list$`gene ID`==SNP[i,'dowstream_gene'],"Gene end" ])
    #     ans <-  Gen_start-SNP_start
    #     
    #   }  
    #   else if(SNP[i,'gene_strand']=='+'){
    #     # selectionne la position Gen end
    #     Gen_start<-as.numeric(region_list[region_list$`gene ID`==SNP[i,'dowstream_gene'],"Gene start" ])
    #     ans <-  SNP_start-Gen_start
    #   }
    #   else
    #     warning("Le marqueur génétique n'est ni '-' ni '+' pour la ligne ", i)
    #   
    #   return(ans)
    # }


    dprom2 <-  function(i){
    SNP_start<-SNP[i,"Start"]
    
    if (SNP[i,'gene_strand']=='-'){
        # selectionne la position Gen end
        Gen_start<-as.numeric(region_list[region_list$`gene ID`==SNP[i,'dowstream_gene'],"Gene end" ])
        #ans <-  Gen_start-SNP_start
        ans <-  Gen_start-SNP_start +1
        
    }  
    else if(SNP[i,'gene_strand']=='+'){
        # selectionne la position Gen end
        Gen_start<-as.numeric(region_list[region_list$`gene ID`==SNP[i,'dowstream_gene'],"Gene start" ])
        ans <-  SNP_start-Gen_start
        #ans <-  SNP_start-Gen_start -1
    }
    else
        warning("Le marqueur génétique n'est ni '-' ni '+' pour la ligne ", i)
    
    
    return(ans)
    }


    # dprom(4)
    #SNP$Distance_to_prom <- mapply(dprom, i=1:nrow(SNP))
    SNP$Distance_to_prom <- mapply(dprom2, i=1:nrow(SNP))

    hist(SNP$Distance_to_prom)  

###### 4) Give each gene a class  and a color 

    unique(SNP$dowstream_gene)
    #topMm14$dowstream_gene <-  factor(topMm14$dowstream_gene, levels=unique(topMm14$dowstream_gene))
    SNP$dowstream_gene <-  factor(SNP$dowstream_gene, levels=unique(SNP$dowstream_gene)) # ne change pas l'ordre des lignes

    unique(SNP$gene_description)

    uu <-  data.frame(gene_description=unique(SNP$gene_description), gene_class=NA)
    uu$gene_class[grep('P450', uu$gene_description)] <-  'P450'
    uu$gene_class[grep('GST', uu$gene_description)] <-  'GST'
    uu$gene_class[grep('UDPGT', uu$gene_description)] <-  'UDPGT'
    uu$gene_class[grep('ABC', uu$gene_description)] <-  'ABC transporter'

    # uu
    vv <-  data.frame(gene_class= unique(uu$gene_class), colG = NA)

    vv$colG <- rainbow(4)
    #
    # uu
    # vv

    uu <- merge(uu, vv, by='gene_class')

    # head(SNP)
    # head(uu)

    SNP <-  merge(SNP, uu, by='gene_description')
    SNP <-  SNP[order(SNP$Chromosome, SNP$Start),]





# I - Filtrage des données

En parallèle on crée un tableau de stockage qui gardera en mémoire les SNPs perdus lors du filtrage : 

    myTableFil <- rbind()
    table(SNP$dowstream_gene)

    myTableFil <-  rbind(myTableFil, all_SNPs=table(SNP$dowstream_gene))
    myTableFil

###### 1) Filtrage sur les couvertures 

Refiltrer le tableau sur la couverture afin d'avoir une estimation correcte des fréquences pour chaque échantillon et limiter les biais.couverture mini colonne total reads) pour chaque échantillon >20.

Visualisation avant filtrage:
    cov_minMax<-cbind(Distance_to_prom=SNP$Distance_to_prom,
                    Cov_min=apply(SNP[,16:21], 1, min),
                    Cov_max=apply(SNP[,16:21], 1, max))
    cov_minMax<-as.data.frame(cov_minMax)

    plot(cov_minMax$Distance_to_prom,
        cov_minMax$Cov_max,
        col='pink', cex=0.1, type='h',yaxs = "i",xaxs = "i",xlab = "Distance_to_gene_start", ylab = "Coverage")

    points(cov_minMax$Distance_to_prom,
        cov_minMax$Cov_min, 
        col='skyblue', cex=0.1, type='h',yaxs = "i",xaxs = "i")

    title(main = "couverture minMmax -all SNPs")

CAR ON VEUT AU MOINS 20 READS POUR CHAQUE LIGNEE ET CHAQUE SNP


    #SNP$Mean_reads<-apply(SNP[,grep("Total", colnames(SNP))],1,mean)
    SNP$min_reads<-apply(SNP[,grep("Total", colnames(SNP))],1,min)

    hist(SNP$min_reads, breaks=20)
    abline(v=20, col='red')
    length(SNP$min_reads<20)
    sum(SNP$min_reads<20) # on va perdre ~400 SNPs

    ## afficher sur le tableau de départ oui/non
    Cov_min=function(i){
    
    #g=SNP$Mean_reads[i]>20  NON !
    g=SNP$min_reads[i]>20 
    if (is.na(g)) Cov_min=NA
    else { 
        if (g) Cov_min=TRUE # sup
        
        else Cov_min=FALSE # inf
    }
    return(Cov_min)
    }

    SNP$covMin20reads<- mapply(Cov_min,i=c(1:nrow(SNP)))
    head(SNP)
    table(SNP$covMin20reads)


    myTableFil <- rbind(myTableFil, covMin20reads=table(SNP$dowstream_gene[(SNP$min_reads>20)]))
    myTableFil

    # ## on verifie
    rowSums(myTableFil)
    # all_SNPs covMin20reads 
    # 3943          3538 
    ## on a perdu 400 reads en filtrant sur 20 reads mini

    # # SNP_filtered<-SNP[SNP$Mean_reads>20,]
    # # length(SNP_filtered$gene_description) # = 3922
    # 
    # # unique(SNP$Cov_min)
    # # summary(SNP$Cov_min=="oui")# 3922 ok on est bon -> 21 en-dessous du seuil


###### 2) Filtrage des SNPs non polymorphes 

On refuse les SNPs pour lesquels tous les % de variants sont > 90% ou tous < 10%

    colnames(SNP)
    colnames(SNP)[grep('Supporting_Reads_PCT', colnames(SNP))]
    head(SNP[, grep('Supporting_Reads_PCT', colnames(SNP))])


    # i=1
    # SNP[i, grep('Supporting_Reads_PCT', colnames(SNP))]
    # mini <- min(SNP[i, grep('Supporting_Reads_PCT', colnames(SNP))])
    # maxi <- max(SNP[i, grep('Supporting_Reads_PCT', colnames(SNP))])
    # 
    # polymorphic <- mini<90 & maxi >5
    # ans <- c(mini=mini, maxi=maxi, polymorphic=polymorphic)
    # ans

    polym <- function(i, minT=90, maxT=5){
    mini <- min(SNP[i, grep('Supporting_Reads_PCT', colnames(SNP))])
    maxi <- max(SNP[i, grep('Supporting_Reads_PCT', colnames(SNP))])
    
    ans <- c(miniPCT=mini, maxiPCT=maxi
            #, polymorphic=polymorphic
    )
    return(ans)
    }

    polym(8)

    hans=t(mapply(polym, i=1:nrow(SNP)))

    head(hans)


    SNP <-  cbind(SNP, hans)

    sum(SNP$miniPCT>90) # loss= 62 SNPs
    sum(SNP$miniPCT>95) # loss= 41 SNPs

    sum(SNP$maxiPCT<0) # loss=0 SNP
    sum(SNP$maxiPCT<5) # loss =434 SNPs
    sum(SNP$maxiPCT<10) # loss 839 SNPs
    SNP$polymorphic <- SNP$miniPCT<90 & SNP$maxiPCT>10

    table(SNP$polymorphic)
    # FALSE  TRUE 
    # 912  3031         # loss of about 900 non polymorphic SNPs

    # on verifie
    head(SNP[!SNP$polymorphic, grep('Supporting_Reads_PCT', colnames(SNP))])
    tail(SNP[!SNP$polymorphic, grep('Supporting_Reads_PCT', colnames(SNP))], n=10)
    # semble OK

    table(SNP$dowstream_gene[SNP$polymorphic])
    table(SNP$dowstream_gene[SNP$min_reads>20 & SNP$polymorphic])


    #table(SNP$dowstream_gene[(SNP$min_reads>20)] )

    #myTableFil <- rbind(myTableFil, covMin20reads=table(SNP$dowstream_gene[(SNP$min_reads>20)]))
    myTableFil <- rbind(myTableFil, polymorphic=table(SNP$dowstream_gene[SNP$min_reads>20 & 
                                                                        SNP$polymorphic]))
    myTableFil

    rowSums(myTableFil)
    # all_SNPs covMin20reads   polymorphic 
    # 3943          3538          2668 

    # Elimnination des indesirables
    SNP_all <-  SNP # as a back up copy

    #SNP <-  SNP[SNP$min_reads>20 & SNP$polymorphic,]
    SNPp <-  SNP[SNP$min_reads>20 & SNP$polymorphic,]




###### 3) Filtrage pval 

Mettre en regard les frequqnce d'un SNPp dans les 6 lignees et les valeurs d'expression du gene: 4 replicats de mesure d'expression > 1 tableau  6 lignes x(4+1) colonnes

    corresp <- cbind(colnames(SNPp)[10:15], t(matrix(colnames(normReads)[2:25], nrow=4))) #tableau à remplir

    colnames(corresp) <- c('frSNPp', 'norm_expr1','norm_expr2','norm_expr3','norm_expr4')

    corresp <-  as.data.frame(corresp)
    corresp$line <- c('SRun', 'NS', 'Sel', 'SPau', 'StLou', 'Ssuz')


    ### Modeles lineaires expr=f(frSNPp) ########

    i=1
    R2__maker <-  function(i,SNPp=SNPp, PLOT=FALSE, LEGEND=FALSE){ 
    SNPp[i,]
    gene=SNPp[i,'dowstream_gene'] 
    # #gene
    # 
    # 
    expr <- t(matrix(normReads[normReads$GeneID==SNPp[i,'dowstream_gene'],2:25], nrow=4))
    #expr
    
    rownames(expr) <- corresp$line
    colnames(expr) <- colnames(corresp)[2:5]
    
    #expr # ressemble a un tableau
    #class(expr)
    
    # exporter le tableau de R pour le ré-importer
    write.table(expr, './temp')
    expr <- read.table('./temp', h=T)
    
    datai <- cbind( expr, frSNPp=unlist(as.vector(SNPp[i, 10:15])) )
    
    title <-  paste(SNPp[i, 'gene_description'], SNPp$dowstream_gene[i], 
                    #SNPp$Start[i],
                    SNPp$Distance_to_prom[i],
                    SNPp$Reference[i], SNPp$Variant_Allele[i], sep ='_')
    
    # Un tableau 24 lignes x 2 colonnes, necessaire pour un modele lineaire
    expr2 <- data.frame(frSNPp=rep(datai$frSNPp, each=4), 
                        norm_expr=as.vector(t(datai[,1:4])))
    rownames(expr2) <- paste(rep(rownames(datai), each=4), rep(1:4, 6), sep='_')
    
    
    modele1 <- lm(norm_expr~frSNPp, data = expr2)
    
    if(is.na(modele1$coefficients[2])){
        p_val = NA
        R2 <-  NA
    } else{
        
        p_val<-summary(modele1)$coefficients['frSNPp', 'Pr(>|t|)'] # p value for H0
        R2<-summary(modele1)$r.squared # R2
        
    }
    
    if(! is.na(modele1$coefficients[2]))
        if(PLOT){
        plot(norm_expr1~frSNPp, data = datai, xlab='variant SNPp (%)', 
            ylim=c(0, max(expr)), 
            xlim=c(0,100), #
            ylab='norm. expression', main=title, col=1:6, cex.axis=0.8, 
            cex.main=0.8)
        
        points(norm_expr2~frSNPp, data = datai, col=1:6)
        points(norm_expr3~frSNPp, data = datai, col=1:6)
        points(norm_expr4~frSNPp, data = datai, col=1:6)
        
        
        abline(a = modele1$coefficients[1], b = modele1$coefficients[2], lty=2)
        
        # limod <- paste("R2=", round(summary(modele1)$r.squared, 2), ", P-value=", 
        #                round(summary(modele1)$coefficients['frSNPp', 'Pr(>|t|)'],2))
        limod <- paste("R2=", round(summary(modele1)$r.squared, 2), ", Slope=", 
                        round(modele1$coefficients[2], 3))
        
        mtext(limod, cex=0.6, line=0.8)
        
        if(LEGEND){
            if( modele1$coefficients[2] >0)
            place='topleft' else place='topright'
            
            #legend('topright', legend = rownames(datai),pch=1, col=1:6, cex=0.8, ncol=6)
            legend(place, legend = rownames(datai),pch=1, col=1:6, cex=0.8, ncol=2)
        }
        
        
        }
    
    #result<-rbind(result,c(gene=gene,R2=R2,p_val=p_val))
    ans <-c(gene=gene,modele1$coefficients,R2=R2,p_val=p_val)
    
    return(ans)
    #result<-rbind(result,c(gene=gene,R2=R2,p_val=p_val))
    }

    hans2 <- mapply(R2__maker, i=1:nrow(SNPp), MoreArgs = list(PLOT=F, SNPp=SNPp))
    hans2 <- t(hans2)

    #SNPp <- cbind(SNPp, hans2[,-1])
    SNPp <- cbind(SNPp, hans2[,-1])

On plot all 14 R2X p-val plots avec ces données filtrés:

    SNPp$R2 <- as.numeric(SNPp$R2)
    SNPp$p_val <- as.numeric(SNPp$p_val)

    graphics.off()
    quartz(h=9,w=12)
    par(mfrow=c(3,5))

    SNPp_Test2 <- rbind() # JM

    for(k in 1: length(unique(SNPp$dowstream_gene))){
    SNPpk <- SNPp[SNPp$dowstream_gene==unique(SNPp$dowstream_gene)[k],]
    SNPpk$adj_p_val <- p.adjust(SNPpk$p_val, method = 'fdr')
    
    lowP <- nrow(SNPpk[SNPpk$adj_p_val<0.05,])
    
    titre <- paste(SNPpk$dowstream_gene[1],":", SNPpk$gene_description[1])
    
    plot(R2~adj_p_val, data=SNPpk, main=titre, cex=0.3, cex.main=0.9, 
        col=(SNPpk$frSNPp<0 )+1, # black dot if R>0, red if <0
        ylim=c(0,1))
    #mtext(paste(nrow(SNPpk), 'SNPps'), cex=0.8)
    mtext(paste(nrow(SNPpk), 'SNPps', '; (' , lowP, 'low P)'), cex=0.8)
    abline(v=0.05, lty=3, col='blue')
    
    SNPp_Test2 <- rbind(SNPp_Test2 , SNPpk) # JM
    }

    SNPp <-  SNPp_Test2 # JM

    #dev.copy2pdf(file='R2_x_p-val_14plots.pdf')


###### 4) Ajout diffSel_NS

    SNPp$diffSel_NS<-abs(SNPp$`Supporting_Reads_PCT_(REU-NS)`-SNPp$`Supporting_Reads_PCT_(REU-DeltaR)`)
    hist(SNPp$diffSel_NS)


    colnames(SNPp)



    # Filtrage
    # myTableFil<-rbind(myTableFil,pval_diffSel_NS_01_15=table(SNPp$dowstream_gene[SNPp$adj_p_val<0.01 & SNPp$diffSel_NS>15]))
    # myTableFil<-rbind(myTableFil,pval_diffSel_NS_01_25=table(SNPp$dowstream_gene[SNPp$adj_p_val<0.01 & SNPp$diffSel_NS>25]))
    # 
    # Rappel ls polymorph & cov sont filtrés 
    myTableFil<-rbind(myTableFil,pval_diffSel_NS_05_15=table(SNPp$dowstream_gene[SNPp$adj_p_val<0.05 & SNPp$diffSel_NS>15]))
    myTableFil<-rbind(myTableFil,pval_diffSel_NS_05_25=table(SNPp$dowstream_gene[SNPp$adj_p_val<0.05 & SNPp$diffSel_NS>25]))
    myTableFil<-rbind(myTableFil,pval_diffSel_NS_05=table(SNPp$dowstream_gene[SNPp$adj_p_val<0.05]))
    myTableFil<-rbind(myTableFil,pval_diffSel_NS_25=table(SNPp$dowstream_gene[SNPp$diffSel_NS>25]))

    myTableFil<-rbind(myTableFil,pval_diffSel_NS_05_25_negatif=table(SNPp$dowstream_gene[SNPp$adj_p_val<0.05 
                                                                                        & SNPp$diffSel_NS>25
                                                                                        & SNPp$frSNPp<0]))
    myTableFil<-rbind(myTableFil,pval_diffSel_NS_05_25_positif=table(SNPp$dowstream_gene[SNPp$adj_p_val<0.05 & SNPp$diffSel_NS>25 & 
                                                                                        SNPp$frSNPp>0]))



    # install.packages("ggplot2")
    # library(ggplot2)
    # 
    # ggplot(data.frame(SNPp), aes(x = SNPp$adj_p_val)) +
    #   geom_histogram( fill =  'blue') +
    #   labs(title = "Histogramme SNPp$adj_p_val",
    #        x = "Valeurs",
    #        y = "SNPp$adj_p_val")
    # 
    # 
    # ggplot(data.frame(SNPp), aes(x = SNPp$diffSel_NS)) +
    #   geom_histogram( fill =  'darkgreen',color='black') +
    #   labs(title = "Histogramme SNPp$diffSel_NS",
    #        x = "Valeurs",
    #        y = "SNPp$diffSel_NS")+
    #   theme_classic()




# II - plot à partir du filtrage des données

###### 1) biplots adj-pval ~ diffSel_NS avec les coloriages necessaires, x14 genes. 

Tous les SNPps polymorphes doivent apparaitre, seul la couleur va signaler les SNPps a retenir.

    SNPp_test3 <- rbind()
    g_labels <-  rbind()

    k=8
    for (k in 1:length(unique(SNPp$dowstream_gene))) {
    
    SNPpk <- SNPp[SNPp$dowstream_gene==unique(SNPp$dowstream_gene)[k],]
    
    head(SNPpk)
    dim(SNPpk)
    
    colors <- ifelse(
        SNPpk$adj_p_val <0.05,
        # 'red', 'black')
        (SNPpk$frSNPp<0) + 1, 'grey') # color will be black (=1) i slope positive, red (=2) if negative
    
    table(colors)
    
    SNPpk$couleur <-  colors
    
    txt <- paste(unique(SNPpk$gene_description), '\n','SNPps: polym=', nrow(SNPpk), '\n',
                'kept_p=', nrow(SNPpk[SNPpk$adj_p_val<0.01 & SNPpk$diffSel_NS>25 & SNPpk$frSNPp>0,]), '\n',
                'kept-m=', nrow(SNPpk[SNPpk$adj_p_val<0.01 & SNPpk$diffSel_NS>25 & SNPpk$frSNPp<0,])
    )
    
    ligne=c(downstream_gene =as.character(unique(SNPpk$dowstream_gene)), txt=txt)
    
    g_labels <-  rbind(g_labels , ligne)
    
    SNPp_test3 <- rbind(SNPp_test3, SNPpk)
    
    }

    SNPp <-  SNPp_test3
    rownames(g_labels) <-  1:14


    ##
    graphics.off()
    quartz(h=9,w=12)
    par(mfrow=c(3,5))
    k=7
    for (k in 1:length(unique(SNPp$dowstream_gene))) {
    
    SNPpk <- SNPp[SNPp$dowstream_gene==unique(SNPp$dowstream_gene)[k],]
    
    plot(x=SNPpk$diffSel_NS, 
        #y=SNPp_k$adj_p_val, 
        y=log10(SNPpk$adj_p_val), 
        #col=colors, 
        col=SNPpk$couleur, 
        main =paste(unique(SNPp$dowstream_gene)[k], "pval_0.05"),
        ylab = "log10(adj_pval)",xlab = "diffSel_NS (%)",
        #xlim = c(0,100), 
        xlim = c(0,60), 
        ylim = c(-5, 0)
    )
    # 
    abline(v=25, col='orange', lty=1) # p<0.01
    #abline(v=15, col='pink', lty=1) # p<0.01
    #abline(h=-2, col='skyblue', lty=2) # p<0.01
    abline(h=log10(0.05), col='darkgreen', lty=2) # p<0.01
    
    txt <- paste(unique(SNPpk$gene_description), 'SNPps: polym=', nrow(SNPpk), 'kept=', 
                nrow(SNPpk[SNPpk$adj_p_val<0.05 & SNPpk$diffSel_NS>25 & SNPpk$min_reads>20 & SNPpk$polymorphic,]))
    
    mtext(txt, cex=0.7)
    }

    #dev.copy2pdf(file='biplots_adj-pval~diffSel_NS.pdf')

    View(SNPpk[,c('Start', 'diffSel_NS','adj_p_val', 'p_val','frSNPp', 'couleur')])

Le code suivant permet de réaliser graphique du script 1 en ayant les coordonnées précises des posisitons d'un gèn présenté dans le papiere:

    #a)300356280#ok
    # Start           diffSel_NS adj_p_val     p_val   frSNPp couleur
    # 3339 300356280   1.090064 0.6441357 0.5642303 7.978482    grey#lui
    # 3340 300356280  10.578265 0.5681380 0.4875918 9.377783    grey

    #b)300355854#ok (2 coordonnées)
    # Start diffSel_NS adj_p_val      p_val   frSNPp couleur
    # 3241 300355854   27.73723 0.1382688 0.05080956 15.17851    grey

    #b.2)300356182
    # SNPpk[SNPpk$Start=='300356182',c('Start', 'diffSel_NS','adj_p_val', 'p_val','frSNPp', 'couleur')]
    # Start diffSel_NS  adj_p_val      p_val     frSNPp couleur
    # 3295 300356182   44.56621 0.32969994 0.21952173  -4.824144    grey
    # 3297 300356182    0.00000 0.07603197 0.01765718 -20.907493    grey

    #c)300357590#ok
    #          Start diffSel_NS  adj_p_val       p_val  frSNPp couleur
    # 3571 300357590   4.506937 0.01814986 0.002307104 63.7001       1
    # #c.2
    # Start diffSel_NS  adj_p_val       p_val   frSNPp couleur
    # 3174 300355627   14.02963 0.02570447 0.003579104 -30.5166       2

    #d)300356909#ok
    # SNPpk[SNPpk$Start=='300356909',c('Start', 'diffSel_NS','adj_p_val', 'p_val','frSNPp', 'couleur')]
    #       Start        diffSel_NS   adj_p_val        p_val   frSNPp couleur
    # 3499 300356909   39.51613 0.001868016 5.223672e-05 20.64789       1

    # # #verification rapide
    # tempo<-SNPpk[SNPpk$Start=='300356280'|SNPpk$Start=='300355854'|SNPpk$Start=='300357590'|SNPpk$Start=='300356909'
    #               ,c('Start', 'diffSel_NS','adj_p_val', 'couleur')]
    # tempo<-SNPpk[SNPpk$Start=='300356909',c('Start', 'diffSel_NS','adj_p_val', 'couleur')]
    # # tempo <-SNPpk[SNPpk$Start=='300355854',c('Start', 'diffSel_NS','adj_p_val', 'p_val','frSNPp', 'couleur')]
    # 
    # plot(x=tempo$diffSel_NS,
    #      #y=SNPp_k$adj_p_val,
    #      y=log10(tempo$adj_p_val),
    #      #col=colors,
    #      col=tempo$couleur,
    #      main =paste(unique(SNPp$dowstream_gene)[k], "pval_0.05"),
    #      ylab = "log10(adj_pval)",xlab = "diffSel_NS (%)",
    #      #xlim = c(0,100),
    #      xlim = c(0,60),
    #      ylim = c(-5, 0)
    # )
    # #
    # abline(v=25, col='orange', lty=1) # p<0.01
    # #abline(v=15, col='pink', lty=1) # p<0.01
    # #abline(h=-2, col='skyblue', lty=2) # p<0.01
    # abline(h=log10(0.05), col='darkgreen', lty=2) # p<0.01
    # 
    # txt <- paste(unique(tempo$gene_description), 'SNPps: polym=', nrow(tempo), 'kept=',
    #              nrow(tempo[tempo$adj_p_val<0.01 & tempo$diffSel_NS>25,]))
    # # #ok

###### 2) plots  genomiques frises 1/2

Calcul des couleurs et formes des boules

    g_pch_col <- function(i, SNPpk){
    
    pch_g <-  NA; col_g <-  NA
    ## pas de boule (pch=NA) si adj_pval>0.05 ou si diffSel_NS≤15%
    ## boule pleine (pch=19) si  0.01 > adj_pval  &  diffSel_NS > 15%, rouge ou noire selon slope
    ## boule grise vide (pch=1) si  0.01 <adj_pval<0.05 & diffSel_NS > 15%
    ## triangle vide (pch=2) si  0.01 <adj_pval<0.05   &   25% <diffSel_NS < 15%, gris
    ## triangle plein (pch=17) si  adj_pval<0.01   &   25% <diffSel_NS < 15%, gris
    
    if(!is.na(SNPpk$R2[i])){
        
        ## pas de boule (pch=NA) si adj_pval>0.05 ou si diffSel_NS≤15%
        if(SNPpk$adj_p_val[i] > 0.05 | SNPpk$diffSel_NS[i] < 15){
        pch_g <-  NA; col_g <-  NA
        }
        
        ## boule pleine (pch=19) si  0.01 > adj_pval  &  diffSel_NS > 15%, rouge ou noire selon slope
        #if(SNPpk$adj_p_val[i] < 0.01 | SNPpk$diffSel_NS[i] > 15){
        if(SNPpk$adj_p_val[i] < 0.01 & SNPpk$diffSel_NS[i] > 15){
        pch_g <-  19;  if(SNPpk$frSNPp[i]>0) col_g <- "black" else col_g='red' #col_g <-  SNPpk$couleur[i]
        }
        
        ## boule grise vide (pch=1) si  0.01 <adj_pval<0.05 & diffSel_NS > 15%
        if(SNPpk$adj_p_val[i] < 0.05 & SNPpk$adj_p_val[i] > 0.01 &
        SNPpk$diffSel_NS[i] > 15){
        pch_g <-  1; col_g <-  'grey'
        }
        
        ## triangle vide (pch=2) si  0.01 <adj_pval<0.05   &   25% <diffSel_NS < 15%, gris
        if( SNPpk$adj_p_val[i] < 0.05 & SNPpk$adj_p_val[i] > 0.01 &
            SNPpk$diffSel_NS[i] > 15 & SNPpk$diffSel_NS[i] < 25){
        
        pch_g <-  2; col_g <- 'grey'
        }
        
        
        ## triangle plein (pch=17) si  adj_pval<0.01   &   25% <diffSel_NS < 15%, gris
        if(SNPpk$adj_p_val[i] < 0.01 &
        SNPpk$diffSel_NS[i] > 15 & SNPpk$diffSel_NS[i] < 25){
        pch_g <-  17; col_g <- 'grey'
        }
        
    }
    
    ans <-  c(pch_g=pch_g, col_g=col_g)
    
    return(ans)
    
    }

    SNPp_Test4 <-  rbind()
    for (k in 1:length(unique(SNPp$dowstream_gene))) {
    SNPpk <- SNPp[SNPp$dowstream_gene==unique(SNPp$dowstream_gene)[k],]
    
    hans <- t(mapply(g_pch_col, i=1:nrow(SNPpk), MoreArgs = list (SNPpk=SNPpk)))

    SNPpk <- cbind(SNPpk, hans)
    SNPpk$pch_g <-  as.numeric(SNPpk$pch_g)
    
    SNPp_Test4 <-  rbind(SNPp_Test4, SNPpk)
    }
    SNPp <- SNPp_Test4


Calcul des frises Minuscules/majuscules: atgc=0; ATGC=1 :

    #setwd('./Sequences_TIGERISK/promoteursSeq_fasta/')#_JM
    setwd("~/Desktop/Documents-stage/Data_analyses_TIGERISK2/promoteursSeq_fasta")#_Sarah
    list.files()

    seqDF <- rbind()
    k=2
    for (k in 1:length(unique(SNPp$dowstream_gene))) {
    
    genek <-  unique(SNPp$dowstream_gene)[k]
    
    SNPpk <- SNPp[SNPp$dowstream_gene==genek,]
    
    uu <- read.table(paste0('~/Desktop/Documents-stage/Data_analyses_TIGERISK2/promoteursSeq_fasta/', genek, '.fasta'), skip = 1)
    uu
    dim(uu)
    nchar(uu[1,1])
    # convert table to string
    myseq=paste()
    for(i in 1:nrow(uu)){
        
        myseq=paste0(myseq, uu[i,1])
        nchar(myseq) # the seq length
        
        head(myseq)
        head(strsplit(myseq, "")) # WRONG
        head(strsplit(myseq, "")[[1]])
        # sequence numbering
        #theSeq <- unlist(strsplit(myseq, "")) # WRONG
        theSeq <-strsplit(myseq, "")[[1]]
        
        length(theSeq)
    head(theSeq, n=10)
        
        # lower vs uppercase ?
        masked <-  function(i, seq=theSeq)
        theSeq[i] %in% letters
        
        # masked(100)
        # masked(10)
        # mapply(masked, i=1:100)
        
        mask <- mapply(masked, i=1:length(theSeq))
        
        seqDFk <- cbind(position=1:length(theSeq), base=theSeq ,masked=mask)
        seqDFk <-  as.data.frame(seqDFk)
        
        head(seqDFk, n=5)
        seqDFk[45:65,] # case identified
        
        table(seqDFk$masked)
        
        }
    
    head(seqDFk, n=5)
    tail(seqDFk)
    
    #write.table(x = seqDFk,file = paste0('seqDF_', SNPpk$dowstream_gene[1]) )
    
    seqDFk$gene <- genek
    seqDFk <-  seqDFk[,c("gene", "position", "base", "masked" )]
    
    seqDF <- rbind(seqDF, seqDFk)
    }
    #write.table(x = seqDF,file = 'seqDF.csv')

    list.files()
    tail(SNPp_Test4)


###### 3) plots genomiques avec frise  2/2
    graphics.off()
    quartz(h=10, w=18)
    par(mfrow=c(3,5))

    for (k in 1:length(unique(SNPp$dowstream_gene))) {
    
    genek <-  unique(SNPp$dowstream_gene)[k]
    
    #SNPpk <- SNPp[SNPp$dowstream_gene==unique(SNPp$dowstream_gene)[k],]
    SNPpk <- SNPp[SNPp$dowstream_gene==genek,]
    
    #seqDFk <- read.table(file = paste('seqDFk', unique(SNPp$dowstream_gene)[k], sep='_'))
    seqDFk <- seqDF[seqDF$gene==genek,]
    
    
    plot(x=SNPpk$Distance_to_prom, 
        y=SNPpk$diffSel_NS, 
        xlim=c(-2500,0),
        ylim=c(0,100),
        type = "h",
        col=SNPpk$couleur,
        main =unique(SNPp$dowstream_gene)[k],
        ylab = "diffSel_NS", xlab = "Distance to promoter",
        # yaxs = "i",xaxs = "i",
    )
    
    # # Ajout des boules/triangles
    points(x=SNPpk$Distance_to_prom,
            y=SNPpk$diffSel_NS, xlim=c(-2500,0),
            col=SNPpk$col_g, pch=SNPpk$pch_g)
    
    # frise des minuscules/majuscules
    if (unique(SNPpk$gene_strand)=="-") {
        points(x=-as.numeric(seqDFk$position), y=80-as.numeric(as.logical(seqDFk$masked))*5, type='l', col='grey')
    }else{
        seqDFk$position<-rev(as.numeric(seqDFk$position))
        
        points(x=-as.numeric(seqDFk$position), y=80-as.numeric(as.logical(seqDFk$masked))*5, type='l', col='grey')
    }
    
    }

    # dev.copy2pdf(file="genomicPlots_frises.pdf")

    # # Pour les gènes chevauchants : LOC...259 et LOC...284 ####
    # 
    # SNPp_259 <- SNPp[SNPp$dowstream_gene==unique(SNPp$dowstream_gene)[7],]
    # 
    # SNPp_284 <- SNPp[SNPp$dowstream_gene==unique(SNPp$dowstream_gene)[8],]
    # SNPp_284$diffSel_NS <- -SNPp_284$diffSel_NS
    # 
    # SNPp_loc <- rbind(SNPp_259, SNPp_284)
    # 
    # ## Plot
    # graphics.off()
    # # quartz()
    # 
    # plot(x=SNPp_loc$Start, 
    #      y=SNPp_loc$diffSel_NS,
    #      #xlim=c(-2500,0),
    #      ylim=c(-70,70),
    #      type = "h",
    #      col=SNPpk$couleur,
    #      ylab = "diffSel_NS", xlab = "SNPp position",
    #      main =paste(unique(SNPp$dowstream_gene)[7],unique(SNPp$dowstream_gene)[8])
    #      )
    # 
    # ## Ajout des boules/triangles
    # points(x=SNPp_loc$Start,
    #        y=SNPp_loc$diffSel_NS, 
    #        #xlim=c(-2500,0),
    #        col=SNPp_loc$col_g, pch=SNPp_loc$pch_g)
    # 
    # ## frise des minuscules/majuscules
    # 
    # ### unique(SNPp_loc$gene_strand)
    # ### unique(SNPp_284$gene_strand)# -
    # 
    # #### le brin +
    # #seqDF_259 <- read.table(file = paste('seqDFk', unique(SNPp$dowstream_gene)[7], sep='_'))#+
    # seqDF_259 <- seqDF[seqDF$gene=='LOC109426259',]
    # 
    # 
    # # head(seqDF_259)
    # seqDF_259$Distance_to_prom <- rev(as.numeric(seqDF_259$position))
    # seqDF_259$Distance_to_prom <- -seqDF_259$Distance_to_prom
    # # head(seqDF_259)
    # # match(SNPp_259$Distance_to_prom,seqDF_259$Distance_to_prom,nomatch = NA_integer_, incomparables = NULL)
    # 
    # seqDF_259_test<-merge(x = SNPp_259, y = seqDF_259, 
    #                       by = 'Distance_to_prom', 
    #                       all.x = T, all.y = T 
    #                       )
    # # tail(seqDF_259_test)
    # 
    # # la frise minuscule/majuscule en bleu
    # points(x=as.numeric(seqDF_259_test$Start), 
    #        y=60-as.numeric(as.logical(seqDF_259_test$masked))*5, 
    #        type='b', 
    #        #type='c', 
    #        col='skyblue',
    #        cex= 0.23
    #        )
    # 
    # 
    # #dev.copy2pdf(file='./Graphics/UDPGT_cluster_genomic.pdf')


###### 4) Filtrage polyallelique :Reduction a une seule ligne des SNPs polyalleliques 

si une meme position retient plusieurs variants,
on veut pouvoir ne garder que le SNP ayant le diffSel_NS le plus elevé

On veut un tableau avec 2500 lignes, pas plus , apres merge !


    # head(SNP)
    # head(SNPp)
    dim(SNP)
    dim(SNPp)

    colnames(SNP)
    colnames(SNPp)

    uu <-  merge(SNP, SNPp[,setdiff(colnames(SNPp), colnames(SNP))], by=0, all.x = T, sort = F)

    dim(uu)
    colnames(uu)

    head(uu[,1:10])
    rownames(uu) <-  uu[,1]
    uu <-  uu[,-1]

    SNP <-  uu
    unique(SNP$gene_strand)
    length(unique(SNP$dowstream_gene))

    colnames(SNP)
    colnames(SNPp)

    # ne garder qu'une seule ligne quand une position est polyallelique


    # for (k in 1:length(unique(SNP$dowstream_gene))) {
    #   SNPk <- SNP[SNP$dowstream_gene==unique(SNP$dowstream_gene)[k],]
    #   
    #   outk <- rbind()
    #   
    #   for(d in 1:length(unique(SNPk$Distance_to_prom))){
    #     
    #     d <- unique(SNPk$Distance_to_prom)[d]
    #     
    #     vv <-  SNPk[SNPk$Distance_to_prom==d,] 
    #     
    #     if(nrow(vv)==1){
    #       polyallelic=FALSE
    #       
    #       vv <- cbind(vv,polyallelic=polyallelic )
    #       
    #       outk <- rbind(outk,vv)
    #     }
    #     
    #     if(nrow(vv)>1){
    #       polyallelic=TRUE
    #       
    #       vv <- cbind(vv,polyallelic=polyallelic )
    #       
    #       vv <- vv[vv$polymorphic==TRUE,]
    #       
    #       if(nrow(vv)>1){
    #         vv <- vv[order(vv$diffSel_NS, decreasing = T),]
    #         vv <-  vv[1,] 
    #       }
    #       outk <- rbind(outk,vv)
    #     }
    #     
    #     # outk <- rbind(outk,vv)
    #     
    #   }
    #   
    #   outk_all <- rbind(outk_all,outk)
    # }
    outk_all <- rbind()
    for (k in 1:length(unique(SNP$dowstream_gene))) {
    SNPk <- SNP[SNP$dowstream_gene==unique(SNP$dowstream_gene)[k],]
    
    outk <- rbind()
    
    for(d in 1:length(unique(SNPk$Distance_to_prom))){
        
        d <- unique(SNPk$Distance_to_prom)[d]
        
        vv <-  SNPk[SNPk$Distance_to_prom==d,] 
        
        if(nrow(vv)==1){
        polyallelic=FALSE
        
        vv <- cbind(vv,polyallelic=polyallelic )
        
        outk <- rbind(outk,vv)
        }
        
        if(nrow(vv)>1){
        polyallelic=TRUE
        #stop()
        
        vv <- cbind(vv,polyallelic=polyallelic )
        
        # vv <- vv[vv$polymorphic==TRUE,]
        
        # if(nrow(vv)>1){
        vv <- vv[order(vv$diffSel_NS, decreasing = T),]
        vv <-  vv[1,] 
        #}
        outk <- rbind(outk,vv)
        }
        
        # outk <- rbind(outk,vv)
        
    }
    
    outk_all <- rbind(outk_all,outk)
    }

    head(outk)
    head(SNPk)

    unique(outk_all$gene_strand)
    dim(outk_all)
    dim(SNP[SNP$polymorphic==T,])
    table(SNP$polymorphic)
    table(outk_all$polymorphic)
    table(outk_all$polyallelic)
    table(SNP$polyallelic)



    SNP <- outk_all

    # unique(SNP$dowstream_gene)


    # mettre pour polyallelique
    myTableFil<-rbind(myTableFil,Total_nbre_polyAllelique=table(SNP$dowstream_gene[SNP$polyallelic]))#==T
    myTableFil<-rbind(myTableFil,Total_without_polyAllelique=table(SNP$dowstream_gene[SNP$polyallelic==F]))
    ## avec les seuils ajouté (parmi ceux qu'on retient combien sont polyall)
    myTableFil<-rbind(myTableFil,Total_nbre_polyAllelique_pval_05_25=table(SNP$dowstream_gene[SNP$adj_p_val<0.05 & SNP$diffSel_NS>25 & SNP$polyallelic]))
    myTableFil<-rbind(myTableFil,Total_without_polyAllelique_pval_05_25=table(SNP$dowstream_gene[SNP$adj_p_val<0.05 & SNP$diffSel_NS>25 & SNP$polyallelic==F]))

    myTableFil<-rbind(myTableFil,Total_nbre_polyAllelique_pval_05_15=table(SNP$dowstream_gene[SNP$adj_p_val<0.05 & SNP$diffSel_NS>15 & SNP$polyallelic]))
    myTableFil<-rbind(myTableFil,Total_without_polyAllelique_pval_05_15=table(SNP$dowstream_gene[SNP$adj_p_val<0.05 & SNP$diffSel_NS>15 & SNP$polyallelic==F]))


    #View(myTableFil)
    myTableFil<-as.data.frame(myTableFil)

    # myTableFil$sum <- apply(myTableFil,1,sum)
    myTableFil<-cbind(myTableFil,sum=apply(myTableFil, 1, sum))


    # sum(SNPp$adj_p_val<0.01 & SNPp$diffSel_NS>25) # 63
    # sum(SNPp$adj_p_val<0.01 & SNPp$diffSel_NS>15) # 171
    # sum(SNPp$adj_p_val<0.05 & SNPp$diffSel_NS>25) # 96
    # sum(SNPp$adj_p_val<0.05 & SNPp$diffSel_NS>15) # 302

    write.csv(myTableFil, file = "~/Desktop/Documents-stage/Data_analyses_TIGERISK2/myTableFil_data2.csv")#Sarah
    # write.csv(myTableFil, file = "~/Desktop/Documents-stage/Data_analyses_TIGERISK2/myTableFil_data.csv")#Sarah
    #write.table(myTableFil, file = "./OUTPUT_data/resume_filtrage_SNPs.csv", sep='\t')#JeanM




###### 5) Tableau Génomique 

    unique(SNP$dowstream_gene)
    unique(SNP$dowstream_gene)[k]

    SNPp$dowstream_gene <- factor(as.character(SNPp$dowstream_gene),
                                #levels = levels(SNP$dowstream_gene)[c(1:7,9:14,8)]
                                levels = c( "LOC115259486", "LOC109412993", "LOC115257289" ,"LOC109413318",
                                            "LOC109408672",
                                            "LOC109410283" , "LOC115259235" ,"LOC109409164", "LOC109407008" ,
                                            "LOC115262073",
                                        "LOC109429691" ,  "LOC109426259" ,"LOC109421559",
                                        "LOC109410284" # en dernier car cause pb
                                        )
                                )
    
    levels(SNPp$dowstream_gene)
    SNPp <-  SNPp[order(SNPp$dowstream_gene),]

    SNP_mergin_all <- rbind() 
    # k=2
    # k=5
    k=8
    for (k in 1:length(unique(SNPp$dowstream_gene))) {
    
    genek <-  unique(SNPp$dowstream_gene)[k]
    
    # le tableau SNPpk à récupérer
    SNPk <- SNP[SNP$dowstream_gene==genek,]
    unique(SNPk$dowstream_gene)
    
    # le tableau seqDFk 
    seqDFk <- seqDF[seqDF$gene==genek,]
    
    dim(seqDFk)
    dim(SNPk)
    head(seqDFk, n=10)
    tail(seqDFk, n=10)
    
    head(SNPk, n=3)
    head(SNPk[,c('Reference', 'Variant_Type','Variant_Allele', 'Distance_to_prom')])
    seqDFk[seqDFk$Distance_to_prom %in% c(-2472,-2471,-2469,-2430),]
    head(seqDFk)
    
    # conversion seqDFk en dist_to_prom selon si c'est '+' ou '-'
    if (unique(SNPk$gene_strand)=="-") {

        seqDFk$Distance_to_prom <- as.numeric(seqDFk$position)
        seqDFk$Distance_to_prom <- -seqDFk$Distance_to_prom
        # head(seqDFk)

    }else{ # si c'est "+"

        seqDFk$Distance_to_prom <- rev(as.numeric(seqDFk$position))
        seqDFk$Distance_to_prom <- -seqDFk$Distance_to_prom

    }
    
    # On merge 

    SNP_mergin <- merge(y = SNPk, x = seqDFk, 
                        by = 'Distance_to_prom', 
                        all.x = T
    )
    
    uu <- SNP_mergin[
        which(SNP_mergin$Variant_Type=="Substitution" &
        nchar(SNP_mergin$Reference)==1)  ,]
    
    if(! all(toupper(uu$base)==toupper(uu$Reference)))
        stop("SNP Reference base differs from expected AALBF3 base")
    
    head(uu[! toupper(uu$base)==toupper(uu$Reference),])
    dim(uu[! toupper(uu$base)==toupper(uu$Reference),])
    
    # give gene attributes to positions with no SNP in SNP_mergin
    {
        head(SNP_mergin[,1:12])
        tail(SNP_mergin[,1:12])
        head(SNP_mergin[! is.na(SNP_mergin$p_val),1:12])
    colnames(SNP_mergin)[c(6,7,13:14, 27)]
    # "gene_description" "Chromosome"       "dowstream_gene"  
    #  "gene_strand"      "gene_class"
    
    SNP_mergin$gene_description <- as.character(SNP_mergin$gene_description)
    ugd <- unique(SNP_mergin$gene_description)
    ugd <-  ugd[!is.na(ugd)]
    ugd
    SNP_mergin$gene_description <- ugd
    
    chr <- unique(SNP_mergin$Chromosome)
    chr <-  chr[!is.na(chr)]
    SNP_mergin$Chromosome <-  chr
    
    SNP_mergin$dowstream_gene <- SNP_mergin$gene
    SNP_mergin <-  SNP_mergin[,- which(colnames(SNP_mergin)=='gene')]
    
    stra <- unique(SNP_mergin$gene_strand)
    stra <-  stra[! is.na(stra)]
    SNP_mergin$gene_strand <-  stra
    
    
    cl <- unique(SNP_mergin$gene_class)
    cl <-  cl[! is.na(cl)]
    SNP_mergin$gene_class <-  cl
    }
    
        # Pour sauvegarder le tout dans une data.frame
    SNP_mergin_all <- rbind(SNP_mergin_all, SNP_mergin)
    
    }


    # uu <- SNP_mergin_all[
    #   SNP_mergin_all$Variant_Type=="Substitution" &
    #   nchar(SNP_mergin_all$Reference)==1  ,]
    #                        
    # if(! all(toupper(uu$base)==toupper(uu$Reference)))
    # #if(toupper(prok2$Reference[red]) != toupper(prok2$base[red]))
    #   stop("SNP Reference base differs from expected AALBF3 base")


    #write.table(SNP_mergin_all, './OUTPUT_data/promoters_base_by_base.tab', sep='\t')#JeanM
    write.table(SNP_mergin_all, '~/Desktop/Documents-stage/Data_analyses_TIGERISK2/myTableFil_data.csv', sep='\t')#Sarah

    ## follow up: FrisesG_v3_cov.R, frises genomiques avec couvertures

    # Vérification du tableau final : ####
    # ## plots genomiques avec frise
    # 
    # SNPp <- SNP_mergin_all
    # 
    # #setwd("~/Desktop/Documents-stage/Data_analyses_TIGERISK2/promoteursSeq_fasta")#_Sarah
    # list.files()
    # graphics.off()
    # quartz(h=10, w=15)
    # par(mfrow=c(3,5))
    # 
    # k=1
    # 
    # unique(SNPp$dowstream_gene)
    # 
    # for (k in 1:length(unique(SNPp$dowstream_gene))) {
    #   SNPpk <- SNPp[SNPp$dowstream_gene==unique(SNPp$dowstream_gene)[k],]
    # 
    #   seqDFk <- read.table(file = paste('seqDFk', unique(SNPp$dowstream_gene)[k], sep='_'))
    # 
    #   plot(x=SNPpk$Distance_to_prom,
    #        y=SNPpk$diffSel_NS,
    #        xlim=c(-2500,0),
    #        ylim=c(0,100),
    #        type = "h",
    #        col=SNPpk$couleur,
    #        main =unique(SNPp$dowstream_gene)[k],
    #        ylab = "diffSel_NS", xlab = "Distance to promoter",
    #        # yaxs = "i",xaxs = "i",
    #   )
    # 
    #   # # Ajout des boules/triangles
    #   points(x=SNPpk$Distance_to_prom,
    #          y=SNPpk$diffSel_NS, xlim=c(-2500,0),
    #          col=SNPpk$col_g, pch=SNPpk$pch_g)
    # 
    #   # frise des minuscules/majuscules
    #   if (unique(SNPpk$gene_strand)=="-") {
    #     points(x=-as.numeric(seqDFk$position), y=80-as.numeric(as.logical(seqDFk$masked))*5, type='l', col='grey')
    #   }else{
    #     seqDFk$position<-rev(as.numeric(seqDFk$position))
    # 
    #     points(x=-as.numeric(seqDFk$position), y=80-as.numeric(as.logical(seqDFk$masked))*5, type='l', col='grey')
    #   }
    # 
    # }
    # 
