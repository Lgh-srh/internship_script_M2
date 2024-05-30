# La fréquence d'un SNP permet-il de prédire le niveau d'expression d'un gène ? (modèle linéaire)
Contact: sarah.loughani@univ-tlse3.fr

## Table of Contents
- **[Introduction ](#)**
  - Ne garder que les colonnes d'interet
  - Dedoublement LOC109410284 / LOC109426259
  - Distance SNP to transcription start 
  - Donner une couleur pour la catégorie de couleur
- **[I - Le modèle linéaire ](#)**
  - Ne garder que les colonnes d'interet
  - Réajustement du tableau de départ
  - Fonction du mooddèle linéaire expr=f(frSNP)
  - plot all 14 R2X p-val plots
  - SNP with top R2 per gene
  - R2 Boxplot for the 10 top SNP with highest freq range: 14 genes
- **[II -  Modèle linéaire & recherche de filtres à selectionner ](#)**
  - Filtrage couvreture minimale/maximale
  - Histograms minMax coveerag with quantiles
  - Histograms slopes with quantiles
  - Scatter plots R2 ~ slope
  - plots version2.0 : R2~pval & R2~minMaxDiff, & slope~minMaxDiff
  - Gene expression plots

# Introduction
Paramètrer les tableaux:
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

    list.files()

    # Nombre de reads RNA normalisés par gene

    library(readxl)

    normReads <- read_excel("Clean_ALL_Data_Assembly_filtered.xlsx",sheet =3)#_Sarah

    normReads <-  as.data.frame(normReads)


    # Les SNPs des promoteurs
    SNPs0 <- read_excel("Clean_ALL_Data_Assembly_filtered.xlsx",sheet =2) #_Sarah
    SNPs0 <- as.data.frame(SNPs0)
# Introduction
Paramètrer les tableaux:
###### 1) Ne garder que les colonnes d'interet

    SNP <- SNPs0[,c(1:9,
                    grep("Supporting", colnames(SNPs0)),
                    grep("Total Reads", colnames(SNPs0))
    )]

    # adapter les noms des colonnes a R
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

Principe: si c'est - alors on fait GeneStart= colEnd sinon GeneStart= colStart et donc on aura GeneStart-position_SNP

    dprom <-  function(i){
      SNP_start<-SNP[i,"Start"]
      
      if (SNP[i,'gene_strand']=='-'){
        # selectionne la position Gen end
        Gen_start<-as.numeric(region_list[region_list$`gene ID`==SNP[i,'dowstream_gene'],"Gene end" ])
        ans <-  Gen_start-SNP_start
      }  
      else if(SNP[i,'gene_strand']=='+'){
        # selectionne la position Gen end
        Gen_start<-as.numeric(region_list[region_list$`gene ID`==SNP[i,'dowstream_gene'],"Gene start" ])
        ans <-  SNP_start-Gen_start
      }
      else
        warning("Le marqueur génétique n'est ni '-' ni '+' pour la ligne ", i)
      
      return(ans)
    }


    dprom(4)
    SNP$Distance_to_prom <- mapply(dprom, i=1:nrow(SNP))

    hist(SNP$Distance_to_prom)  



###### 4) Donner une couleur pour la catégorie de couleur

    unique(SNP$gene_description)

    uu <-  data.frame(gene_description=unique(SNP$gene_description), gene_class=NA)
    uu$gene_class[grep('P450', uu$gene_description)] <-  'P450'
    uu$gene_class[grep('GST', uu$gene_description)] <-  'GST'
    uu$gene_class[grep('UDPGT', uu$gene_description)] <-  'UDPGT'
    uu$gene_class[grep('ABC', uu$gene_description)] <-  'ABC transporter'


    uu
    vv <-  data.frame(gene_class= unique(uu$gene_class), colG = NA)

    vv$colG <- rainbow(4)

    uu
    vv

    uu <- merge(uu, vv, by='gene_class')

    head(SNP)
    head(uu)

    SNP <-  merge(SNP, uu, by='gene_description')
    SNP <-  SNP[order(SNP$Chromosome, SNP$Start),]

# I- Le modèle linéaire

Le but étant de créer un graphique Reads vs SNPs.

###### 1) Réajustement du tableau de départ
Principe: 
Mettre en regard les frequences d'un SNP dans les 6 lignees et les valeurs d'expression du gene -> 4 replicats de mesure d'expression > 1 tableau  6 lignes x(4+1) colonnes

    colnames(SNP)[10:15]

    colnames(normReads)[2:25]

    matrix(colnames(normReads)[2:25], nrow=4)
    t(matrix(colnames(normReads)[2:25], nrow=4))

    corresp <- cbind(colnames(SNP)[10:15], t(matrix(colnames(normReads)[2:25], nrow=4)))

    colnames(corresp) <- c('frSNP', 'norm_expr1','norm_expr2','norm_expr3','norm_expr4')

    corresp <-  as.data.frame(corresp)
    corresp$line <- c('SRun', 'NS', 'Sel', 'SPau', 'StLou', 'Ssuz')
    corresp # voila le tableau a remplir pour chaque SNP

###### 2) Fonction du mooddèle linéaire expr=f(frSNP)
    i=1
    R2__maker <-  function(i,SNP=SNP, PLOT=FALSE, LEGEND=FALSE){ 
      SNP[i,]
      gene=SNP[i,'dowstream_gene'] 
      # #gene
      # 
      # 
      expr <- t(matrix(normReads[normReads$GeneID==SNP[i,'dowstream_gene'],2:25], nrow=4))
      #expr
      
      rownames(expr) <- corresp$line
      colnames(expr) <- colnames(corresp)[2:5]
      
      #expr # ressemble a un tableau
      #class(expr)
      
      # exporter le tableau de R pour le ré-importer
      write.table(expr, './temp')
      expr <- read.table('./temp', h=T)
      
      datai <- cbind( expr, frSNP=unlist(as.vector(SNP[i, 10:15])) )
      
      title <-  paste(SNP[i, 'gene_description'], SNP$dowstream_gene[i], 
                      #SNP$Start[i],
                      SNP$Distance_to_prom[i],
                      SNP$Reference[i], SNP$Variant_Allele[i], sep ='_')
      
      # Un tableau 24 lignes x 2 colonnes, necessaire pour un modele lineaire
      expr2 <- data.frame(frSNP=rep(datai$frSNP, each=4), 
                          norm_expr=as.vector(t(datai[,1:4])))
      rownames(expr2) <- paste(rep(rownames(datai), each=4), rep(1:4, 6), sep='_')
      
      
      modele1 <- lm(norm_expr~frSNP, data = expr2)
      if(is.na(modele1$coefficients[2])){
        p_val = NA
        R2 <-  NA
      } else{
        
        p_val<-summary(modele1)$coefficients['frSNP', 'Pr(>|t|)'] # p value for H0
        R2<-summary(modele1)$r.squared # R2
        
      }
      
      if(! is.na(modele1$coefficients[2]))
        if(PLOT){
          plot(norm_expr1~frSNP, data = datai, xlab='variant SNP (%)', 
              ylim=c(0, max(expr)), 
              xlim=c(0,100), #
              ylab='norm. expression', main=title, col=1:6, cex.axis=0.8, 
              cex.main=0.8)
          
          points(norm_expr2~frSNP, data = datai, col=1:6)
          points(norm_expr3~frSNP, data = datai, col=1:6)
          points(norm_expr4~frSNP, data = datai, col=1:6)
          
          
          abline(a = modele1$coefficients[1], b = modele1$coefficients[2], lty=2)
          
          # limod <- paste("R2=", round(summary(modele1)$r.squared, 2), ", P-value=", 
          #                round(summary(modele1)$coefficients['frSNP', 'Pr(>|t|)'],2))
          limod <- paste("R2=", round(summary(modele1)$r.squared, 2), ", Slope=", 
                        round(modele1$coefficients[2], 3),"p-value",round(summary(modele1)$coefficients['frSNP', 'Pr(>|t|)'],2))
          
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
    
    hans2 <- mapply(R2__maker, i=1:nrow(SNP), MoreArgs = list(PLOT=F, SNP=SNP))
    hans2 <- t(hans2)
    tail(hans2)
    dim(SNP)

    SNP_test <- cbind(SNP, hans2[,-1])

###### 3) plot all 14 R2X p-val plots

    SNP_test$R2 <- as.numeric(SNP_test$R2)
    SNP_test$p_val <- as.numeric(SNP_test$p_val)

    graphics.off()
    quartz(h=9,w=12)
    par(mfrow=c(3,5))

    for(k in 1: length(unique(SNP_test$dowstream_gene))){
      SNPk <- SNP_test[SNP_test$dowstream_gene==unique(SNP_test$dowstream_gene)[k],]
      SNPk$adj_p_val <- p.adjust(SNPk$p_val, method = 'fdr')
      
      lowP <- nrow(SNPk[SNPk$adj_p_val<0.05,])
      
      titre <- paste(SNPk$dowstream_gene[1],":", SNPk$gene_description[1])
      
      plot(R2~adj_p_val, data=SNPk, main=titre, cex=0.3, cex.main=0.9, 
          col=(SNPk$frSNP<0 )+1, # black dot if R>0, red if <0
          ylim=c(0,1))
      #mtext(paste(nrow(SNPk), 'SNPs'), cex=0.8)
      mtext(paste(nrow(SNPk), 'SNPs', '; (' , lowP, 'low P)'), cex=0.8)
      abline(v=0.05, lty=3, col='blue')
      
    }

    #dev.copy2pdf(file='R2_x_p-val_14plots.pdf')


>Attention, le code ci-dessous a été réalisé a posteriori afin de récupérer des données de coordonnées précises:

    #R2__maker(1, PLOT=T)
    # R2__maker(1, PLOT=T, SNP=SNP)
    # #R2__maker(97, PLOT=T)
    # R2__maker(97, SNP=SNP,PLOT=T, LEGEND=T)
    # R2__maker(1, SNP=SNP,PLOT=T, LEGEND=T)
    # R2__maker(100, SNP=SNP,PLOT=T, LEGEND=T)
    # dim(SNP)

    #pour un gene k = 7=> choisir à partir duscript V6i: LOC109426259
    #a)300356280#ok
    # Start           diffSel_NS adj_p_val     p_val   frSNPp couleur
    # 3339 300356280   1.090064 0.6441357 0.5642303 7.978482    grey#lui
    # 3340 300356280  10.578265 0.5681380 0.4875918 9.377783    grey

    #b)300355854#ok (2 coordonnées)
    # Start diffSel_NS adj_p_val      p_val   frSNPp couleur
    # 3241 300355854   27.73723 0.1382688 0.05080956 15.17851    grey

    #b.2)300356182
    # Start            diffSel_NS  adj_p_val      p_val     frSNPp couleur
    # 3295 300356182   44.56621 0.32969994 0.21952173  -4.824144    grey#lui
    # 3297 300356182    0.00000 0.07603197 0.01765718 -20.907493    grey

    #c)300357590#ok
    #          Start diffSel_NS  adj_p_val       p_val  frSNPp couleur
    # 3571 300357590   4.506937 0.01814986 0.002307104 63.7001       1

    #c.2) 300355627
    #             Start diffSel_NS  adj_p_val       p_val   frSNPp couleur
    # 3174 300355627   14.02963 0.02570447 0.003579104 -30.5166       2#lui

    #d)300356909#ok
    # SNPpk[SNPpk$Start=='300356909',c('Start', 'diffSel_NS','adj_p_val', 'p_val','frSNPp', 'couleur')]
    #       Start      diffSel_NS   adj_p_val        p_val   frSNPp couleur
    # 3499 300356909   39.51613 0.001868016.  5.223672e-05 20.64789       1#lui
    # 39.51613
    # 0.001868016
    # 5.223672e-05
    # 20.647893
    # 1


      k=1
      # ags <- SNP[SNP$dowstream_gene=='LOC109426259'
      #            & SNP$Start=='300356471'|SNP$Start=='300356182'|SNP$Start=='300357590'|SNP$Start=='300356909'
      #            ,]

      300356909

      for (q in c(#'300356280','300356182','300355627',
                  '300356909')) {
        ags <- SNP[SNP$dowstream_gene=='LOC109426259',]
        ags_2 <- ags[ags$Start==q,]
        mapply(R2__maker, i=1:nrow(ags_2), MoreArgs = list(PLOT=T, SNP=ags_2))
        
      }


###### 4) SNP with top R2 per gene

    # a subfile with the SNP with the highest R2 for each gene
    topSNPs <- rbind()
    #k=1
    for(k in 1: length(unique(SNP_test$dowstream_gene))){
      SNPk <- SNP_test[SNP_test$dowstream_gene==unique(SNP_test$dowstream_gene)[k],]
      
      SNPk$adj_p_val <- p.adjust(SNPk$p_val, method = 'fdr')
      
      SNPk <- SNPk[order(SNPk$R2, decreasing = T),]
      
      topSNPs <- rbind(topSNPs, SNPk[1,])
    }

    head(topSNPs, n=3) # a subfile with the SNP with the highest R2 for each gene


    graphics.off()
    quartz(h=9,w=12)
    par(mfrow=c(3,5))

    # hans2 <- mapply(R2__maker, k=1:nrow(SNP), MoreArgs = list(PLOT=F, SNP=SNP))
    #mapply(R2__maker, k=1:nrow(topSNPs), MoreArgs = list(PLOT=TRUE, SNP=topSNPs))
    mapply(R2__maker, i=1:nrow(topSNPs), MoreArgs = list(PLOT=TRUE, SNP=topSNPs))
    #
    # legende bidouillee
    plot(x = rep(0.5, 6), y = 2:7, col=1:6, ylim=c(0,8), yaxt = "n", xaxt='n',  xlab='', ylab='')
    #text(x = rep(0.56, 6), y = 2:7, labels = rownames(datai))
    text(x = rep(0.56, 6), y = 2:7, labels =corresp$line)

    #dev.copy2pdf(file='topSNP_x_14.pdf')

    ## plot Slope vs R2  (gene 2, GST_Theta ) ######
    #k=1 # do not use: possibly duplicated gene, no need for promoter mutations
    k=2

    #for(k in 1: length(unique(SNP_test$dowstream_gene))){
    SNPk <- SNP_test[SNP_test$dowstream_gene==unique(SNP_test$dowstream_gene)[k],]

    SNPk$adj_p_val <- p.adjust(SNPk$p_val, method = 'fdr')

    SNPk <- SNPk[order(SNPk$R2, decreasing = T),]

    #topSNPs <- rbind(topSNPs, SNPi[1,])
    #}

    head(SNPk,n=4)
    SNPk

    cl <- function(i){
      col='grey'
      if(! is.na(SNPk$adj_p_val[i]))
        if(SNPk$adj_p_val[i]<0.05)
          col <-  "black"

        return(col)
    }

    SNPk$couleur <- mapply(cl, i=1:nrow(SNPk))

    graphics.off()
    quartz(h=9,w=12)
    par(mfrow=c(3,5))

    plot(frSNP~R2, data=SNPk, type='h', main=SNPk$dowstream_gene[1], 
        ylab='Slope', cex.axis=0.8, col=SNPk$couleur)

    legend('topright',title='adj_p_val',pch=1,
          legend = c('≥0.05', '<0.05'), col=c('grey', 'black'))

    text(as.numeric(SNPk$R2[13]), col='red',
        y = as.numeric(SNPk$frSNP[13]) *2.8, labels = '13º',)


    mapply(R2__maker, i=1:13, MoreArgs = list(PLOT=T, SNP=SNPk))

    # legende bidouillee
    plot(x = rep(0.5, 6), y = 2:7, col=1:6, ylim=c(0,8), yaxt = "n", xaxt='n',  xlab='', ylab='')
    text(x = rep(0.56, 6), y = 2:7, labels =corresp$line)

    #dev.copy2pdf(file='Slope_vs_R2.pdf')

    #' need to plot not Slope but Slope X max(delta(variant SNP %))

###### 5) R2 Boxplot for the 10 top SNP with highest freq range: 14 genes

    #a boxplot of R2 values with, for each gene, for the 10 SNPs with the highest min-max diff in variant SNP %


    MINMAX <-  function(i, SNPk){
      
      #SNPk[i,grep('Supporting_Reads_PCT', colnames(SNPk))]
      MIN <- min(SNPk[i,grep('Supporting_Reads_PCT', colnames(SNPk))])
      MAX <- max(SNPk[i,grep('Supporting_Reads_PCT', colnames(SNPk))])
      ans <- c(MIN=MIN, MAX=MAX)
      
      return(ans)
    }

    MINMAX(3, SNPk)
    mapply(MINMAX, i=1:12,MoreArgs=list(SNPk=SNPk))
    hans <-  t(mapply(MINMAX, i=1:nrow(SNPk), MoreArgs=list(SNPk=SNPk)))
    SNPk <- cbind(SNPk, hans)
    head(SNPk, n=4)

    SNPk$minmaxDiff <- SNPk$MAX -SNPk$MIN

    SNPk <- SNPk[order(SNPk$minmaxDiff, decreasing = T),]
    topminmax <- SNPk[1:10,]
    topminmax$R2

    max(SNPk$R2)

    ###
    nbestDelta <- function(k, n=10){
      SNPk <- SNP_test[SNP_test$dowstream_gene==unique(SNP_test$dowstream_gene)[k],]
      
      SNPk$adj_p_val <- p.adjust(SNPk$p_val, method = 'fdr')
      
      SNPk <- SNPk[order(SNPk$R2, decreasing = T),]
      
      cl <- function(i){
        col='grey'
        if(! is.na(SNPk$adj_p_val[i]))
          if(SNPk$adj_p_val[i]<0.05)
            col <-  "black"
        
        #print(i)
        return(col)
      }
      
      SNPk$couleur <- mapply(cl, i=1:nrow(SNPk))
      
      hans <-  t(mapply(MINMAX, i=1:nrow(SNPk), MoreArgs=list(SNPk=SNPk)))
      SNPk <- cbind(SNPk, hans)
      
      SNPk$minmaxDiff <- SNPk$MAX -SNPk$MIN
      
      SNPk <- SNPk[order(SNPk$minmaxDiff, decreasing = T),]
      topminmax <- SNPk[1:10,]
      
      return(topminmax)
      
    }


    nbestDelta(6)

    topMm14 <- rbind()

    for(k in 1:14)
      topMm14 <- rbind(topMm14, nbestDelta(k))

    colnames(topMm14)
    topMm14$dowstream_gene <-  factor(topMm14$dowstream_gene, levels=unique(topMm14$dowstream_gene))
    #head(topMm14[,c(7,19,25)])
    head(topMm14[,c(7,21,27)])
    head(topMm14[,])

    library(ggplot2)

    graphics.off()
    quartz()
    #Basic box plot
    #ggplot(df, aes(x = "", y = y)) +
    g <-  ggplot(topMm14, aes(x = "", y = R2)) +
      geom_boxplot() +
      geom_jitter(colour=topMm14$colG)+
      ggtitle('ten_SNPs_with_highest varFreq ranges')+
      theme(plot.title = element_text(hjust = 0.5))

    g

    g + facet_wrap(.~ dowstream_gene, ncol=5) #+

    #dev.copy2pdf(file='ten_SNPs_with_highest varFreq ranges.pdf')










# II- Modèle linéaire & recherche de filtres à selectionner
###### 1) Filtrage couvreture minimale/maximale

On crée un tablau mettant en evidence la couvertur minimale et maximale pour les différentes lignées à chacune des poistions.


    # Histograms minMax. ##########
    SNP_mM <-  rbind()

    SNPmMmaker <- function(k){
      SNPk <- SNP_test[SNP_test$dowstream_gene==unique(SNP_test$dowstream_gene)[k],]
      
      SNPk$adj_p_val <- p.adjust(SNPk$p_val, method = 'fdr')
      
      SNPk <- SNPk[order(SNPk$R2, decreasing = T),]
      
      cl <- function(i){
        col='grey'
        if(! is.na(SNPk$adj_p_val[i]))
          if(SNPk$adj_p_val[i]<0.05)
            col <-  "black"
        
        #print(i)
        return(col)
      }
      
      SNPk$couleur <- mapply(cl, i=1:nrow(SNPk))
      
      hans <-  t(mapply(MINMAX, i=1:nrow(SNPk), MoreArgs=list(SNPk=SNPk)))
      SNPk <- cbind(SNPk, hans)
      
      SNPk$minmaxDiff <- SNPk$MAX -SNPk$MIN
      
      SNPk <- SNPk[order(SNPk$minmaxDiff, decreasing = T),]
      
      return(SNPk)
      
    }

    for(k in 1:length(unique(SNP_test$dowstream_gene)))
      SNP_mM <-  rbind(SNP_mM, SNPmMmaker(k))

    dim(SNP_mM)
    dim(SNP)
    colnames(SNP_mM)

    # l <- ggplot(SNP_mM, aes(x=minmaxDiff)) +
    #   geom_histogram()
    # 
    # l
    # 
    # l+ facet_wrap(.~ dowstream_gene, ncol=5) +
    #   ggtitle('Distribution of  maximal differences in SNP frequence ')+
    #   ylab('Number of SNPs')+ 
    #   xlab('Percent_max - Percent_min') +
    #   theme(plot.title = element_text(hjust = 0.5))


    SNP_mM$dowstream_gene <- as.character(SNP_mM$ dowstream_gene)
    SNP_mM$dowstream_gene <- factor(SNP_mM$ dowstream_gene, levels=unique(SNP_mM$ dowstream_gene))

    SNP_mM$gene_class <-  as.character(SNP_mM$gene_class)
    SNP_mM$gene_class <-  factor(SNP_mM$gene_class, levels=unique(SNP_mM$gene_class))

    g <- ggplot(SNP_mM, aes(x=minmaxDiff, colour=gene_class, fill=gene_class)) +
      geom_histogram()

    g

    g+ facet_wrap(.~ dowstream_gene, ncol=5) +
      ggtitle('Distribution of  maximal differences in SNP frequence ')+
      ylab('Number of SNPs')+ 
      xlab('Percent_max - Percent_min') +
      theme(plot.title = element_text(hjust = 0.5))

    #dev.copy2pdf(file='ten_SNPs_with_highest varFreq ranges.pdf')


###### 2) Histograms minMax coveerag with quantiles

adding quantile lines
see https://stackoverflow.com/questions/30568873/plot-quantiles-of-distribution-in-ggplot2-with-facets

    library (dplyr)
    ##
    colnames(SNP_mM)

    d2b <- SNP_mM %>%
      #group_by(model, type) %>%
      group_by(dowstream_gene) %>%
      summarize(q90= quantile(minmaxDiff, probs = .90),
                q95 = quantile(minmaxDiff, probs = .95))

    d2b

    l2 <- ggplot(SNP_mM, aes(x=minmaxDiff, colour=gene_class, fill=gene_class)) +
      geom_histogram()

    l2

    l2+ facet_wrap(.~ dowstream_gene, ncol=5) +
      ggtitle('Distribution of  maximal differences in SNP frequence ')+
      ylab('Number of SNPs')+ 
      xlab('max_% - min_%') +
      theme(plot.title = element_text(hjust = 0.5)) +
      geom_vline(data = d2b, aes(xintercept = q90)) +
      geom_vline(data = d2b, aes(xintercept = q95))

    ## Scatter plots R2 ~ minMax #####


    d2c <- SNP_mM %>%
      #group_by(model, type) %>%
      group_by(dowstream_gene) %>%
      summarize(q90= quantile(R2, probs = .90, na.rm=T),
                q95 = quantile(R2, probs = .95, na.rm=T))

    d2c
    graphics.off()
    quartz()

    h <- ggplot(SNP_mM, aes(x=minmaxDiff, y=R2, colour=gene_class, fill=gene_class)) +
      geom_point(size=0.4,)

    h
    #quartz(h=8, w=10)
    h + facet_wrap(.~ dowstream_gene, ncol=5)+
      ggtitle('R2_minMaxDiff biplots ')+
      theme(plot.title = element_text(hjust = 0.5)) +
      geom_vline(data = d2b, aes(xintercept = q90)) +
      geom_vline(data = d2b, aes(xintercept = q95)) +
      
      geom_hline(data = d2c, aes(yintercept = q90)) +
      geom_hline(data = d2c, aes(yintercept = q95))

    #dev.copy2pdf(file='R2_minMaxDiff_biplots.pdf')

###### 3) Histograms slopes  with quantiles 

adding quantile lines
see https://stackoverflow.com/questions/30568873/plot-quantiles-of-distribution-in-ggplot2-with-facets

    library (dplyr)
    ##
    colnames(SNP_mM)

    head(SNP_mM$frSNP)
    SNP_mM$frSNP <-  as.numeric(SNP_mM$frSNP)
    hist(abs(SNP_mM$frSNP))
    hist(log10(abs(SNP_mM$frSNP)))


    # #l3 <- ggplot(SNP_mM, aes(y=abs(frSNP), colour=gene_class, fill=gene_class)) +
    # l3 <- ggplot(SNP_mM, aes(y=log10(abs(frSNP)), colour=gene_class, fill=gene_class)) +
    # #geom_boxplot()
    # geom_histogram()
    # 
    # l3
    # 
    # quartz()
    # l3+ facet_wrap(.~ dowstream_gene, ncol=5) #+


    d2d <- SNP_mM %>%
      group_by(dowstream_gene) %>%
      summarize(q10= quantile(log10(abs(SNP_mM$frSNP)), probs = .10, na.rm=T),
                q05 = quantile(log10(abs(SNP_mM$frSNP)), probs = .05, na.rm=T))

    d2d # Wrong output for some obscure reason... all values equal whatever gene !?

    # fix d2d ! ...
    q105 <-  function(k){
      SNPk <- SNP_test[SNP_test$dowstream_gene==unique(SNP_test$dowstream_gene)[k],]
      unique(SNPk$dowstream_gene)
      head(SNPk$frSNP)
      SNPk$frSNP <- as.numeric(SNPk$frSNP)
      q10 <- quantile(log10(abs(SNPk$frSNP)), probs = .10, na.rm=T)
      q05 <- quantile(log10(abs(SNPk$frSNP)), probs = .05, na.rm=T)
      
      ans <- c(q05, q10)
      
      return(ans)
    }

    q105(4)
    mapply(q105, k=1:5)
    quant_log_slope <- t(mapply(q105, k=1:14))
    rownames(quant_log_slope) <- unique(SNP_test$dowstream_gene)
    colnames(quant_log_slope) <- c('q5', 'q10')
    d2c
    d2d
    quant_log_slope <- as_tibble(quant_log_slope)

    d2d$q10 <- quant_log_slope$q10
    d2d$q05 <- quant_log_slope$q5

    ##
    #l3 <- ggplot(SNP_mM, aes(y=log10(abs(frSNP)), colour=gene_class, fill=gene_class)) +
    l3 <- ggplot(SNP_mM, aes(x=log10(abs(frSNP)), colour=gene_class, fill=gene_class)) +
      geom_histogram()
    l3

    l3 <-  l3+ facet_wrap(.~ dowstream_gene, ncol=5)
    l3


    l3 + 
      geom_vline(data = d2d, aes(xintercept = q10)) +
      geom_vline(data = d2d, aes(xintercept = q05))+
      # geom_hline(data = d2d, aes(yintercept = q10)) +
      # geom_hline(data = d2d, aes(yintercept = q05))+
      ggtitle('Distribution of slopes')+
      theme(plot.title = element_text(hjust = 0.5)) +
      ylab('number of SNPs')

    #dev.copy2pdf(file='histogram_slopes.pdf')

###### 4) Scatter plots R2 ~ slope 

    #graphics.off()
    quartz()

    #h <- ggplot(SNP_mM, aes(x=minmaxDiff, y=R2, colour=gene_class, fill=gene_class)) +
    m <- ggplot(SNP_mM, aes(x=log10(abs(frSNP)), y=R2, colour=gene_class, fill=gene_class)) +
      geom_point(size=0.4,)

    m
    #quartz(h=8, w=10)
    m + facet_wrap(.~ dowstream_gene, ncol=5)+
      ggtitle('R2~Slope biplots ')+
      theme(plot.title = element_text(hjust = 0.5)) +
      geom_vline(data = d2d, aes(xintercept = q10)) +
      geom_vline(data = d2d, aes(xintercept = q05)) +
      
      geom_hline(data = d2c, aes(yintercept = q90)) +
      geom_hline(data = d2c, aes(yintercept = q95))

    #dev.copy2pdf(file='R2~slope_biplots.pdf')




    # 
    #dev.copy2pdf(file='MIN_MAX_histograms.pdf')




###### 5) plots version2.0 : R2~pval & R2~minMaxDiff, & slope~minMaxDiff

    k=2
    SNPk <- SNP_mM[SNP_mM$dowstream_gene==unique(SNP_mM$dowstream_gene)[k],]

    PVAL=0.05

    clr <- function(i, PVAL=0.05){
      #' couleur rouge si 
      #' adj pval signif & 
      #' slope <0 & 
      #' Sel in minN &
      #' NS or SRun in maxN
      #' 
      #' couleur noire si 
      #' adj pval signif & 
      #' slope >0 & 
      #' Sel in minN &
      #' NS or SRun in maxN
      
      MIN <- min(SNPk[i,grep('Supporting_Reads_PCT', colnames(SNPk))])
      MAX <- max(SNPk[i,grep('Supporting_Reads_PCT', colnames(SNPk))])
      
      suppR <- SNPk[i,grep('Supporting_Reads_PCT', colnames(SNPk))]
      
      minN <- names(suppR)[which(suppR==MIN)]
      minN <- sub('Supporting_Reads_PCT_\\(', '', minN)
      minN <- sub('\\)', '', minN)
      minN
      
      maxN <- names(suppR)[which(suppR==MAX)]
      maxN <- sub('Supporting_Reads_PCT_\\(', '', maxN)
      maxN <- sub('\\)', '', maxN)
      maxN
      
      
      col='grey'
      
      if(! is.na(SNPk$adj_p_val[i])){
        if(SNPk$adj_p_val[i]<PVAL)
          if(SNPk[i,'frSNP']>0)
            col <-  "black" else col <- 'red'
            
      }
      
      if(col != 'grey'){
        if(SNPk$frSNP[i] >0)
          condi <- "REU-DeltaR" %in% maxN &
            ("S-Run" %in% minN | "REU-NS" %in% minN)
        if( SNPk$frSNP[i] <0)
          condi <- "REU-DeltaR" %in% minN &
            ("S-Run" %in% maxN | "REU-NS" %in% maxN)
        
        if(!condi)
          col <-  'grey'
      }
      #print(i)
      return(col)
    }

    clr(8)
    #col <- cl(i, PVAL=PVAL)
    mapply(clr, i=1:50)
    SNPk$couleur <- mapply(clr, i=1:nrow(SNPk))


    uu <-  rbind()

    for(k in 1: length(unique(SNP_mM$dowstream_gene))){
      SNPk <- SNP_mM[SNP_mM$dowstream_gene==unique(SNP_mM$dowstream_gene)[k],]
      
      SNPk$couleur <- mapply(clr, i=1:nrow(SNPk), MoreArgs = list(PVAL=PVAL))
      
      uu <- rbind(uu, SNPk)
    }

    SNP_mM <-  uu

    ##
    #graphics.off()
    quartz(h=9,w=12)
    par(mfrow=c(3,5))

    for(k in 1: length(unique(SNP_mM$dowstream_gene))){
      SNPk <- SNP_mM[SNP_mM$dowstream_gene==unique(SNP_mM$dowstream_gene)[k],]
      
      titre <- paste(SNPk$dowstream_gene[1],":", SNPk$gene_description[1])
      
      plot(R2~adj_p_val, data=SNPk, main=titre, cex=0.4, cex.main=0.9,
          #col=(SNPk$frSNP<0 )+1, # black dot if R>0, red if <0
          col=couleur, 
          ylim=c(0,1), xlim=c(0, 0.3))
      
      abline(v=0.05, lty=3, col='blue')
      txt <- paste(nrow(SNPk), 'SNPs \n',
                  'black :' , table(SNPk$couleur)["black"], "\n",
                  'red :' , table(SNPk$couleur)["red"]
      )
      text(x=0.25, y=0.9, txt, cex=0.8)
    }

    #dev.copy2pdf(file='R2_x_p-val_14plotsB.pdf')

    ##
    quartz(h=9,w=12)
    par(mfrow=c(3,5))

    for(k in 1: length(unique(SNP_test$dowstream_gene))){
      
      SNPk <- SNP_mM[SNP_mM$dowstream_gene==unique(SNP_mM$dowstream_gene)[k],]
      
      plot(R2~minmaxDiff, data=SNPk, col=couleur, xlim=c(0,100), cex=0.5)
      
      mtext(unique(SNPk$dowstream_gene), col=unique(SNPk$colG), line=1)
      
      q90 <- quantile(SNPk$minmaxDiff, probs = .90)
      
      abline(v=q90, col='blue')
    }

    #dev.copy2pdf(file='R2_minMaxDiff_biplotsB.pdf')

    ##
    #graphics.off()
    quartz(h=9,w=12)
    par(mfrow=c(3,5))

    for(k in 1: length(unique(SNP_test$dowstream_gene))){
      
      SNPk <- SNP_mM[SNP_mM$dowstream_gene==unique(SNP_mM$dowstream_gene)[k],]
      
      
      plot(y=log10(abs(SNPk$frSNP)), x=SNPk$minmaxDiff, col=SNPk$couleur, cex=0.5, 
          xlab="max(freqVar) - min(freqVar)", 
          ylab='log10(abs(slope))')
      
      mtext(unique(SNPk$dowstream_gene), col=unique(SNPk$colG), line=1)
      
    }

    #dev.copy2pdf(file='slope~minMaxDiff.pdf')

###### 6) Gene expression plots 


    library(reshape2)


    exprim <-  function(k){
      SNPk <- SNP_mM[SNP_mM$dowstream_gene==unique(SNP_mM$dowstream_gene)[k],]
      
      gene=SNPk[1,'dowstream_gene'] 
      
      #expr <- t(matrix(normReads[normReads$GeneID==SNPk[i,'dowstream_gene'],2:25], nrow=4))
      expr <- t(matrix(normReads[normReads$GeneID==SNPk[1,'dowstream_gene'],2:25], nrow=4))
      #expr
      
      rownames(expr) <- corresp$line
      colnames(expr) <- colnames(corresp)[2:5]
      
      # exporter le tableau de R pour le ré-importer
      write.table(expr, './temp')
      expr <- read.table('./temp', h=T)
      
      expr3 <- cbind(site=rownames(expr),expr, col=as.factor(1:6), gene=gene)
      
      expr.m <-    melt(expr3)
      
      colnames(expr.m)[5] <- "norm_expr"
      expr.m <-  expr.m[,-4]
      
      expr.m$site <- factor(expr.m$site, 
                            levels=c( "SRun",  "NS", "StLou",  "SPau", "Ssuz",  "Sel"))
      
      return(expr.m)
    }

    # expr14.m <- rbind()

    exprim(4)

    expr14.m <- rbind(exprim(1), 
                      exprim(2),
                      exprim(3),
                      exprim(4),
                      exprim(5),
                      exprim(6),
                      exprim(7),
                      exprim(8),
                      exprim(9),
                      exprim(10),
                      exprim(11),
                      exprim(12),
                      exprim(13),
                      exprim(14)
    )
    head(expr14.m)

    unique(expr14.m$gene)

    quartz(h=9, w=14)

    b14 <-  ggplot(expr14.m , aes(x=site, y=norm_expr))+
      geom_boxplot()

    b14

    anno2 <- data.frame(x=6, y=0, hjust=0.5, size=0.1,
                        #lab=paste0("FC_NS_Sel=",round(shift$FC, 1)),
                        lab=paste0("(",round(shift$FC, 1), ')'),
                        gene=unique(expr14.m$gene))


    b14+
      geom_jitter(colour=expr14.m$col, width = 0.1) +
      geom_text(data = anno2, aes(x = x,  y = y, #size=size, # hjust=1, 
                                  label = lab))+
      facet_wrap(.~ gene, ncol=5, scales='free_y')


    dev.copy2pdf(file='./Graphics/expression_boxplots.pdf')

    # visualize single criteria selection methods ####
    head(SNPk)
    #k=3
    k=2

    SNPk <- SNP_mM[SNP_mM$dowstream_gene==unique(SNP_mM$dowstream_gene)[k],]


    graphics.off()
    quartz(h=9,w=12)
    par(mfrow=c(3,5))


    # view the 5 top SNPs by lowest pval= highest R2
    SNPk <- SNPk[order(SNPk$R2, decreasing = T),]
    mapply(R2__maker, i=1:5, MoreArgs = list(PLOT=T, SNP=SNPk))

    # view the 5 top SNPs by highest pval= lowest R2
    SNPk <- SNPk[order(SNPk$R2, decreasing = FALSE),]
    mapply(R2__maker, i=1:5, MoreArgs = list(PLOT=T, SNP=SNPk))
    head(SNPk$adj_p_val)
    # [1] 0.9789345 0.9738891 0.9686142 0.9686142 0.9686142
    # [6] 0.9686142

    # # view the 5 top SNPs by  delta freq range
    SNPk <- SNPk[order(SNPk$minmaxDiff, decreasing = TRUE),]
    mapply(R2__maker, i=1:5, MoreArgs = list(PLOT=T, SNP=SNPk))

    #dev.copy2pdf(file='GST_Theta_3methods.pdf')

    ## maximize the deltaR_NS freq diff and view the plots
    SNPk <- SNPk[order(abs(SNPk$`Supporting_Reads_PCT_(REU-DeltaR)`- SNPk$`Supporting_Reads_PCT_(REU-NS)`), decreasing = TRUE),]
    mapply(R2__maker, i=1:15, MoreArgs = list(PLOT=T, SNP=SNPk))

    #dev.copy2pdf(file='GST_Theta_topSelvsNS_SNPs.pdf')

    SNPk[ 1:15,c('downstream_gene', 'adj_p_val')]
    SNPk[ 1:15,c('downstream_gene')]

    cbind(adj_p_val=round(SNPk$adj_p_val[1:15],3), 
          minmaxDiff=SNPk$minmaxDiff[1:15],
          slope=round(as.numeric(SNPk$frSNP[1:15]),1))

