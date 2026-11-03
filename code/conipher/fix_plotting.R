#Some of th eplotting stuff colouring and ordering went wrong in original code. Fixed that here
library(CONIPHER)
library(igraph)

treebuilding_plot <- function(sample_pyclone_tree) {
  cat('\n Plotting inferred phylogenetic tree \n')
  require(mapplots)
  sampleID  <-  sample_pyclone_tree$parameters$sampleID
  prefix <- sample_pyclone_tree$parameters$prefix
  generalSave <- sample_pyclone_tree$parameters$generalSave
  ccf_buffer <- sample_pyclone_tree$parameters$ccf_buffer
  pval_cutoff <- sample_pyclone_tree$parameters$pval_cutoff
  use_boot <- sample_pyclone_tree$parameters$use_boot
  merge_clusters <- sample_pyclone_tree$parameters$merge_clusters
  correct_cpn_clusters <- sample_pyclone_tree$parameters$correct_cpn_clusters
  adjust_noisy_clusters <- sample_pyclone_tree$parameters$adjust_noisy_clusters
  adjust_noisy_clusters_prop <- sample_pyclone_tree$parameters$adjust_noisy_clusters_prop
  min_ccf <- sample_pyclone_tree$parameters$min_ccf
  min_cluster_size           <- sample_pyclone_tree$parameters$min_cluster_size

  nested_pyclone <- sample_pyclone_tree$nested_pyclone
  pyclone_tree   <- sample_pyclone_tree$graph_pyclone
  clonality_table <- sample_pyclone_tree$clonality_out$clonality_table_corrected
  clonality_out <- sample_pyclone_tree$clonality_out
  test_pyclone  <- sample_pyclone_tree$ccf_table_pyclone_clean
  cpn_removed_clusters <- sample_pyclone_tree$cpn_removed_clusters
  if(length(cpn_removed_clusters)==0)
  {x
    cpn_removed_clusters <- NA
  }

  merged_clusters <- sample_pyclone_tree$merged_clusters

  ### Plot trees -- AUTOMATIC
  date <- gsub('-', '', substr(Sys.time(), 1, 10))

  pdfname <- file.path(generalSave, 'corrected_pytree_and_bar.pdf')

  height.mult.factor <- ceiling(nrow(nested_pyclone$ccf_cluster_table)/25)
  width.mult.factor  <- ceiling(nrow(nested_pyclone$ccf_cluster_table)/25)

  #Hello world
  pdf(pdfname, width=22*width.mult.factor, height=12*height.mult.factor)
  {
    par(mar=c(0,0,0,0))
    layout(cbind(1:(nrow(nested_pyclone$ccf_cluster_table)+2),rep(nrow(nested_pyclone$ccf_cluster_table)+3,nrow(nested_pyclone$ccf_cluster_table)+2),rep(nrow(nested_pyclone$ccf_cluster_table)+3,nrow(nested_pyclone$ccf_cluster_table)+2)))
    require(beeswarm)

    tmp <- nested_pyclone$ccf_cluster_table
    main <- "" #paste(substr(colnames(tmp)[1], 1, 8), '\ Phylo CCF values', sep = '')
    
   
    colnames(tmp) <- c("Lymph node", "Primary")
    suppressPackageStartupMessages(require(gplots))
    plot.new()
    par(mar=c(2,2,2,2))
    title(main, cex = 2)

    colours.to.use <- color.tree(1:nrow(nested_pyclone$ccf_cluster_table))
    par(mar=c(0.1,5,0.1,2),lend=1)

    for (j in 1:nrow(nested_pyclone$ccf_cluster_table))
    {

      if(j==1)
      {
        border.col <- ifelse(clonality_table[j,]=='clonal','black','grey')
        bp <- barplot(nested_pyclone$ccf_cluster_table[j,],las=1,col=colours.to.use[j],border=border.col,names="",ylab=paste("Cl",rownames(nested_pyclone$ccf_cluster_table)[j],sep=" "),ylim=c(0,115),yaxt='n',cex.axis=1.5)

      }
      if(j!=1)
      {
        border.col <- ifelse(clonality_table[j,]=='clonal','black','grey')
        bp <- barplot(nested_pyclone$ccf_cluster_table[j,],las=1,col=colours.to.use[j],border=border.col,names="",ylab=paste("Cl",rownames(nested_pyclone$ccf_cluster_table)[j],sep=" "),ylim=c(0,115),yaxt='n',cex.axis=1.5)

      }
      axis(side = 2,at = c(0,50,100),labels=c(c(0,50,100)),las=1)
      if(j ==nrow(nested_pyclone$ccf_cluster_table))
      {
        axis(side=1,at=bp,labels=c("Primary", "Lymph node")
             ,tick=FALSE
             ,cex.axis=2)

      }
      abline(h=0)
      abline(h=100,lty='dashed')
      abline(h=50,lty='dashed')
      for (bar in 1:length(bp))
      {

        beeswarm(test_pyclone[test_pyclone[,'PycloneCluster']%in%rownames(nested_pyclone$ccf_cluster_table)[j],grep('PhyloCCF',colnames(test_pyclone))[bar]]*100
                 ,at=bp[bar]
                 ,add=TRUE
                 ,corralWidth = 0.5
                 ,method='swarm'
                 ,corral='wrap'
                 ,pch=21
                 ,col=colours.to.use[j]
                 ,bg='grey')
        segments(x0 = bp[bar],x1 = bp[bar],y0 = nested_pyclone$ccf_ci_lower[j,bar],y1 = nested_pyclone$ccf_ci_upper[j,bar],lwd=5)
        text(x=bp[bar],y=25,labels=nested_pyclone$ccf_cluster_table[j,bar],cex =1.5)
      }
    }

    plot.new()
    par(mar=c(2.1, 2.1, 4.1, 38), xpd=TRUE)

    g <- graph.data.frame(pyclone_tree$default_tree,directed = FALSE)
    indx <- V(g)$name
    vcol <- setNames(color.tree(pyclone_tree$edgelength), names(pyclone_tree$edgelength))[indx]
    l <- layout_as_tree(g, root = pyclone_tree$trunk)

    pie.size <- ncol(sample_pyclone_tree$nested_pyclone$ccf_cluster_table)
    node.shape <- setNames(rep('pie', length(vcol)), names(vcol))
    pie.slices <- lapply(1:length(vcol), function(x) rep(1, pie.size))
    empty.col = '#bdbdbd'#'white'

    node_size_factor <- log2(max(pyclone_tree$edgelength)) / 30
    node.size <- log2(pyclone_tree$edgelength) / node_size_factor
    node.size <- node.size[names(node.shape)]
    mat <- sample_pyclone_tree$nested_pyclone$ccf_cluster_table
    matchesss <- match(names(vcol), rownames(sample_pyclone_tree$nested_pyclone$ccf_cluster_table))

    base.pie.colours <- colours.to.use[matchesss]

    pie.colors <- sample_pyclone_tree$nested_pyclone$ccf_cluster_table[
    match(names(vcol), rownames(sample_pyclone_tree$nested_pyclone$ccf_cluster_table)),
    ,
    drop = FALSE
    ]

    pie.colors <- ifelse(pie.colors>=90,99,pie.colors)
    pie.colors <- ifelse(pie.colors<10&pie.colors>=1,10,pie.colors)

    pie.colors <- lapply(1:nrow(pie.colors), function(x) {
    if(!all(is.na(pie.colors[x,]))){
        tmp <- pie.colors[x,]
        tmp2 <- tmp
        colfunc <- colorRampPalette(c("white", base.pie.colours[x]))
        speccolours <- colfunc(100)
        tmp[tmp>0] <- speccolours[tmp]
        tmp[tmp2 == 0] <- empty.col
        tmp
    }
})


    g_dir <- graph.data.frame(pyclone_tree$default_tree,directed = TRUE)
    edges <- get.edgelist(g_dir)
    ecol <- setNames(rep('#bdbdbd', nrow(edges)),edges[,2])# baseline, set edge color to black
    ewidth <- rep(1,length(ecol))

    #label consensus edges in other colour
    ecol[paste(edges[,1],edges[,2],sep=":")%in%pyclone_tree$consensus_branches] <- '#000000'
    ewidth[paste(edges[,1],edges[,2],sep=":")%in%pyclone_tree$consensus_branches] <-150

    plot(g
         , layout=l
         , vertex.color = vcol[indx]
         , vertex.frame.color=vcol[indx]
         , vertex.shape = node.shape
         , vertex.lwd=5
         , vertex.pie.lwd=3
         , vertex.pie = pie.slices
         , vertex.pie.color = pie.colors
         , vertex.size = node.size
         , edge.color=ecol
         , edge.size=ewidth
         , vertex.label.cex=2
         , vertex.label.pos=2
         , vertex.label.dist=0
         , vertex.label.family='Helvetica'
         , vertex.label.font=2
         , vertex.label.color = 'black')

    title(main = sampleID, cex.main = 2)
    legend.pie(1,1,labels=gsub(paste0(substr(colnames(tmp)[1], 1, 8), "_"), "", colnames(tmp)), radius=0.2, bty="n", col='#bdbdbd',
               cex=1.5, label.dist=0.8
               ,border='white')

    snv_clusters <- sort(pyclone_tree$edgelength[indx], decreasing = T)

    snv_clusters_removed <- pyclone_tree$edgelength
    snv_clusters_removed <- sort(snv_clusters_removed[!names(snv_clusters_removed) %in% indx], decreasing = T)
    if(!is.na(cpn_removed_clusters[1]))
    {
      snv_clusters_removed <- c(snv_clusters_removed,table(sample_pyclone_tree$ccf_table_pyclone[,'PycloneCluster'])[cpn_removed_clusters])
    }

    tmp <- legend('topright', inset = c(-0.3, 0), legend = paste(names(snv_clusters), ' (', snv_clusters,' SNVs)', sep = ''), col = vcol[names(snv_clusters)], pch = 19, title = 'Clusters included:', bty = 'n', cex=2)  ## inset option controls how far from x and y margins
    if (length(snv_clusters_removed) > 0) {
      if(!is.na(cpn_removed_clusters[1]))
      {
        to_plot <- table(sample_pyclone_tree$ccf_table_pyclone[,'PycloneCluster'])[cpn_removed_clusters]
        legend(x=tmp$rect$left,y = 0,inset = c(-0.3, 0),legend = paste(names(to_plot), ' (', to_plot, ' SNVs)', sep = ''), col = vcol[names(to_plot)], pch = 19, title = 'Copy# clusters removed:', bty = 'n', cex=1.5)

      }
      legend('bottomright', inset = c(-0.3, 0), legend = paste(names(snv_clusters_removed), ' (', snv_clusters_removed, ' SNVs)', sep = ''), col = vcol[names(snv_clusters_removed)], pch = 19, title = 'Clusters removed:', bty = 'n', cex=1.5)
    }




  }
  dev.off()

  #next, plot all the possible trees
  trees.to.plot <- pyclone_tree$alt_trees
  if(length(trees.to.plot)==0)
  {
    #nothing to plot here.
  }
  if(length(trees.to.plot)!=0)
  {

    date <- gsub('-', '', substr(Sys.time(), 1, 10))

    pdfname <- file.path(generalSave, 'corrected_pytree_multipletrees.pdf')

    mult.factor <- ceiling(length(trees.to.plot)/50)

    pdf(pdfname, width=12*mult.factor, height=12*mult.factor)
    {
      nr.trees              <- length(trees.to.plot)
      columnnum             <- 1
      rownum                <- nr.trees/columnnum

      if(nr.trees<=50)
      {
        nr.to.use <- nr.trees
      }

      if(nr.trees>50)
      {
        nr.to.use <- signif(nr.trees+5,2)
      }

      for(i in 1: nr.to.use) {
        if((nr.to.use %% i) == 0) {
          if((i+(nr.to.use/i))<(columnnum+rownum))
          {
            columnnum <- i
            rownum    <- nr.to.use/columnnum
          }
        }
      }

      if(columnnum==1)
      {
        columnnum <- ceiling(columnnum*2)
        rownum    <- ceiling(rownum/2)
      }

      par(mfrow=c(rownum,columnnum),xpd=TRUE,mar=c(1, 1,1, 1))
      for (i in 1:nr.trees)
      {
        auto_tree    <- trees.to.plot[[i]]
        g <- graph.data.frame(auto_tree,directed = FALSE)
        indx <- V(g)$name
        vcol <- setNames(color.tree(pyclone_tree$edgelength), names(pyclone_tree$edgelength))[indx]

        l <- layout_as_tree(g, root = pyclone_tree$trunk)

        pie.size <- ncol(nested_pyclone$ccf_cluster_table)
        node.shape <- setNames(rep('pie', length(vcol)), names(vcol))
        pie.slices <- lapply(1:length(vcol), function(x) rep(1, pie.size))
        empty.col = 'gray85'

        node_size_factor <- log2(max(pyclone_tree$edgelength)) / 30
        node.size <- log2(pyclone_tree$edgelength) / node_size_factor
        node.size <- node.size[names(node.shape)]

        pie.colors <- nested_pyclone$ccf_cluster_table[match(names(vcol), rownames(nested_pyclone$ccf_cluster_table)),, drop = F]
        pie.colors <- lapply(1:nrow(pie.colors), function(x) {
          if(!all(is.na(pie.colors[x,]))){
            tmp <- pie.colors[x,]
            tmp[tmp > 0] <- vcol[rownames(pie.colors)[x]]
            tmp[tmp == 0] <- empty.col
            tmp
          }
        })

        g_dir <- graph.data.frame(auto_tree,directed = TRUE)
        edges <- get.edgelist(g_dir)
        ecol <- setNames(rep('#bdbdbd', nrow(edges)),edges[,2])# baseline, set edge color to black
        ewidth <- rep(1,length(ecol))

        #label consensus edges in other colour
        ecol[paste(edges[,1],edges[,2],sep=":")%in%pyclone_tree$consensus_branches] <- '#000000'
        ewidth[paste(edges[,1],edges[,2],sep=":")%in%pyclone_tree$consensus_branches] <- 2

        plot(g, main = sampleID, cex.main=2,
             , layout = l
            , vertex.color = vcol[indx]
            , vertex.shape = node.shape
             , vertex.pie = pie.slices
             , vertex.pie.color = pie.colors
             , vertex.pie.lty = 0
             , vertex.size = node.size
             , edge.width =ewidth
             , edge.color=ecol
             , arrow.size =0
             ,arrow.width=0
             ,arrow.mode=0
        )

      }

    }
    dev.off()
  }
}

color.tree <- function(edgelength, opacity = 255) {
  suppressPackageStartupMessages(require(RColorBrewer))
  ncols <- length(edgelength)
  max.cols <- 12
  cols     <- paste(brewer.pal(min(max.cols,ncols), name = "Paired"), sep = "")
  cols     <- rep(cols, ceiling(ncols / max.cols))[1:ncols]
  #Make colors with more than 12 clusters more transparent
  nopaq <- round(opacity / ceiling(ncols / 12))
  opacity_val <- c(sapply(1:ceiling(ncols / 12), function(x) rep(opacity - nopaq * (x - 1), 12)))[1:ncols]
  cols.opac <- paste(cols, as.hexmode(opacity_val), sep = "")
  return(cols.opac)
}



test <- readRDS("/DATA/s.wijnen/Projects/ClonalEvolution/clean_code/ClonalEvolution/output/CONIPHER_oncoanlyser/DONIMI_D14/DONIMI_D14.tree.RDS")

test$parameters$sampleID

treebuilding_plot(test)


annonymous_ids <- read.csv("/DATA/s.wijnen/Projects/ClonalEvolution/output/Anonymous_ptid_mapper.csv")

ff <- list.files("/DATA/s.wijnen/Projects/ClonalEvolution/clean_code/ClonalEvolution/output/CONIPHER_oncoanlyser/", pattern = ".tree.RDS", recursive=TRUE, full.names=TRUE)
for(f in ff) {
    tree <- readRDS(f)
    anonymous_id <- annonymous_ids %>% dplyr::filter(sampleID == tree$parameters$sampleID) %>% dplyr::pull(annonID)
    tree$parameters$sampleID <- anonymous_id
    treebuilding_plot(tree)
}

