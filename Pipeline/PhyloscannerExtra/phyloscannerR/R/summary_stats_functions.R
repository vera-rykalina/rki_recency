# Collects a variety of statistics about a single patient in a single tree

#' @keywords internal
#' @export calc.subtree.stats
#' @importFrom ape unroot

calc.subtree.stats <- function(host.id, tree.id, tree, tips.for.hosts, splits.table, tip.regex, no.read.counts, verbose = F){
  
  if(verbose) cat("Calculating statistics for host ",host.id,"\n", sep="")
  
  relevant.reads <- splits.table %>% filter(host == host.id)
  
  subgraphs <- length(unique(relevant.reads$subgraph))
  
  all.tips <- tips.for.hosts[[host.id]]
  
  subtree.all <- NULL

  if(length(all.tips)>1){
  
    subtree.all <- drop.tip(tree, tip=tree$tip.label[!(tree$tip.label %in% all.tips)])
    
    # just in case pruning changes the tip order
    
    if(!no.read.counts){
      reads.per.tip <- sapply(subtree.all$tip.label, function(x) as.numeric(read.count.from.label(x, tip.regex)))
    } else {
      reads.per.tip <- rep(1, length(subtree.all$tip.label))
    }
    
    names(reads.per.tip) <- subtree.all$tip.label
    
    overall.rtt <- calc.mean.root.to.tip(subtree.all, reads.per.tip)
    
    
    if(length(subtree.all$tip.label)>2){
      unrooted.subtree <- unroot(subtree.all)
      max.branch.length <- max(unrooted.subtree$edge.length)
    } else {
      # ape won't unroot two-tip trees
      max.branch.length <- sum(subtree.all$edge.length)
    }
    
    pat.distances <- cophenetic(subtree.all)
    max.pat.distance <- max(pat.distances)
    
    global.mean.pat.distance <- mean(pat.distances[upper.tri(pat.distances)])
    
    if(subgraphs==1){
      subgraph.mean.pat.distance <- global.mean.pat.distance
      largest.rtt <- overall.rtt
      
    } else {
      
      splits <- unique(relevant.reads$subgraph)
      
      reads.per.split <- sapply(splits, function(x) sum(splits.table$reads[which(splits.table$subgraph==x)] ) )
      
      winner <- splits[which(reads.per.split==max(reads.per.split))]
      
      if(length(winner)>1){
        #Rare but possible. Pick an arbitrary winner
        winner <- winner[1]
      }
      
      winner.tips <- relevant.reads$tip[which(relevant.reads$subgraph==winner)]
      if(length(winner.tips)==1){
        largest.rtt <- 0
        subgraph.mean.pat.distance <- NA
      } else {
        subtree <- drop.tip(tree, tip=tree$tip.label[!(tree$tip.label %in% winner.tips)])
        
        # just in case pruning changes the tip order
        
        if(!no.read.counts){
          reads.per.tip <- sapply(subtree$tip.label, function(x) as.numeric(read.count.from.label(x, tip.regex)))
        } else {
          reads.per.tip <- rep(1, length(subtree$tip.label))
        }
        
        names(reads.per.tip) <- subtree$tip.label
        
        largest.rtt <- calc.mean.root.to.tip(subtree, reads.per.tip)
        pat.distances <- cophenetic(subtree)
        subgraph.mean.pat.distance <- mean(pat.distances[upper.tri(pat.distances)])
      }
    }
    
  } else {
    
    # A patient with zero or one tips has no branches and patristic distances are meaningless.
    
    max.branch.length <- NA
    max.pat.distance <- NA
    global.mean.pat.distance <- NA
    subgraph.mean.pat.distance <- NA
    
    if(length(all.tips)==1){
      # A patient with one tip does have a root
      overall.rtt <- 0
      largest.rtt <- 0
    } else {
      # An absent patient has no root or tips
      overall.rtt <- NA
      largest.rtt <- NA
    }
  }
  
  return(list(overall.rtt = overall.rtt, largest.rtt = largest.rtt, max.branch.length = max.branch.length, max.pat.distance = max.pat.distance, 
              global.mean.pat.distance=global.mean.pat.distance, subgraph.mean.pat.distance = subgraph.mean.pat.distance))
}

# Calculates tip and read counts for all patients in a given window

#' @keywords internal
#' @export get.tip.and.read.counts

get.tip.and.read.counts <- function(ptree, hosts, tip.regex, has.read.counts, verbose = F){
  
  tree.id <- ptree$id
  
  if(verbose) cat("Calculating tip and read counts for tree ID ",ptree$id,"\n",sep="")
  
  tree <- ptree$tree
  blacklist <- ptree$blacklist
  
  # A vector of patients for each tip
  
  hosts.for.tips <- sapply(tree$tip.label, function(x) host.from.label(x, tip.regex))
  hosts.for.tips[blacklist] <- NA
  
  # If no patients from the input file are actually here, give a warning
  
  hosts.present <- intersect(hosts, unique(hosts.for.tips))
  if(length(hosts.present)==0){
    warning(paste("No listed hosts appear in tree ",tree.id,"\n",sep=""))
  }
  
  # A list of tips for each patient 
  
  tips.for.hosts <- map(setNames(hosts, hosts), function(x) tree$tip.label[which(hosts.for.tips==x)])
  
  # Make the tibble
  
  window.table <- tibble(host.id=hosts)
  window.table <- window.table %>% mutate(tree.id = tree.id) 
  window.table <- window.table %>% mutate(tips = map_int(host.id, function(x)  (length(tips.for.hosts[[x]]))))
  window.table <- window.table %>% mutate(reads = map_int(host.id, function(x){
    if(length(tips.for.hosts[[x]])==0){
      return(0L)
    } else {
      if(has.read.counts){
        return(sum(map_int(tips.for.hosts[[x]], function(y) read.count.from.label(y, tip.regex))))
      } else {
        return(as.integer(length(tips.for.hosts[[x]])))
      }
    }
  })) 
  
  window.table
}

#' @keywords internal
#' @export calc.all.stats.in.window

calc.all.stats.in.window <- function(ptree, hosts, tip.regex, has.read.counts, verbose = F){

  if(verbose) cat("Calculating host statistics for tree ID ",ptree$id,"\n",sep="")
  
  id <- ptree$id
  
  # todo this needs to go 
  
  if(is.null(id)){
    id <- ptree$suffix
  }
  
  tree <- ptree$tree
  blacklist <- ptree$blacklist
  
  # Get the splits
  
  splits.table <- ptree$splits.table
  
  # Find the clades
  
  clade.mrcas.by.host <- ptree$clade.mrcas.by.host
  clades.by.host <- ptree$clades.by.host
  
  # A vector of patients for each tip
  
  hosts.for.tips <- sapply(tree$tip.label, function(x) host.from.label(x, tip.regex))
  hosts.for.tips[blacklist] <- NA
  
  # If no patients from the input file are actually here, give a warning
  
  hosts.present <- intersect(hosts, unique(hosts.for.tips))
  
  if(length(hosts.present)==0){
    warning(paste("No listed hosts appear in tree ",ptree$id,"\n",sep=""))
  }
  
  # A list of tips for each patient 
  
  tips.for.hosts <- lapply(setNames(hosts, hosts), function(x) tree$tip.label[which(hosts.for.tips==x)])
  
  # Make the tibble
  
  window.table <- tibble(host.id = hosts, 
                         tree.id = id, 
                         xcoord = ptree$xcoord, 
                         tips = sapply(hosts, function(x) as.numeric(length(tips.for.hosts[[x]]))),
                         reads = sapply(hosts, function(x){
                           if(length(tips.for.hosts[[x]])==0){
                             return(0)
                           } else {
                             if(has.read.counts){
                               return(sum(sapply(tips.for.hosts[[x]], function(y) as.numeric(read.count.from.label(y, tip.regex)))))
                             } else {
                               return(as.numeric(length(tips.for.hosts[[x]])))
                             }
                           }
                         } 
                         ),
                         subgraphs =  sapply(hosts, function(x) length(unique(splits.table[which(splits.table$host==x),]$subgraph))),
                         clades = sapply(hosts, function(x) length(clades.by.host[[x]]))
  )
  
  
  
  new.cols <- map_dfr(hosts, function(x) calc.subtree.stats(x, id, tree, tips.for.hosts, splits.table, tip.regex, !has.read.counts, verbose))

  normalised.new.cols <- new.cols %>% transmute(across(.cols = 1:6, .fns = function(x) x/ptree$normalisation.constant))
  
  colnames(normalised.new.cols) <- paste0("normalised.", colnames(new.cols))
  
  window.table <- bind_cols(window.table, new.cols) 
  window.table <- bind_cols(window.table, normalised.new.cols) 
  
  recomb.file.name <- ptree$recombination.file.name
  
  
  
  if (!is.null(recomb.file.name)) {
    
    recomb.df <- read_csv(recomb.file.name, show_col_types = FALSE)
    col.names <- colnames(recomb.df)
    for (expected.col.name in c("Bam file","Recombination metric")) {
      if (! expected.col.name %in% col.names) {
        stop("Expected column name ", expected.col.name, " missing from ",
             recomb.file.name, ". Quitting.")
      }
    }
    missing.IDs <- setdiff(hosts, recomb.df$`Bam file`)
    extra.IDs <- setdiff(recomb.df$`Bam file`, hosts)
    if (length(missing.IDs) > 0) {
      stop(paste0("Error: the following bam file IDs are missing from ",
                  recomb.file.name, ": ", paste(missing.IDs, collapse = " "), "\nQuitting."))
    }
    if (length(extra.IDs) > 0) {
      stop(paste0("Error: the following bam file IDs were unexpected in ",
                  recomb.file.name, ": ", paste(extra.IDs, collapse = " "), "\nQuitting."))
    }
    recomb.df <- recomb.df %>% 
      select(`Bam file`, `Recombination metric`) %>% 
      rename(host.id = `Bam file`, recombination.metric = `Recombination metric`)

    window.table <- window.table %>% inner_join(recomb.df, by = "host.id")
  }
  
  # If you did dual detection, add that in
  
  if(!is.null(ptree$dual.detection.splits)){
    # Q: is there a tidyverse match?
    
    window.table$solo.dual.count <- ptree$dual.detection.splits$count[match(window.table$host.id, ptree$dual.detection.splits$host)]
    
  }
  

  window.table
}

# Get the proportion of reads from a given patient in each of its subgraphs in a single tree

#' @keywords internal
#' @export get.read.proportions

get.read.proportions <- function(host.id, tree.id, splits.table){
  
  this.pat.splits <- splits.table[which(splits.table$host==host.id),]
  
  if(nrow(this.pat.splits)>0){
    this.pat.reads.by.split <- aggregate(this.pat.splits$reads, by=list(Category=this.pat.splits$subgraph), sum)
    
    this.pat.reads.by.split <- this.pat.reads.by.split[order(this.pat.reads.by.split$x, decreasing = T),]
    return(this.pat.reads.by.split$x/sum(this.pat.reads.by.split$x))
    
  } else {
    return(NA)
  }
}

