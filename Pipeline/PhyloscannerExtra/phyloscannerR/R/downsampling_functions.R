# This returns the tips to get rid of (i.e. blacklist) per host for a given tree

#' @keywords internal
#' @importFrom extraDistr rmvhyper
#' @export downsample.host

downsample.host <- function(host, tree, number, tip.regex, host.ids, rename = F , exclude.underrepresented = F, no.read.counts = T, name = NA, verbose=F){
  if(verbose) cat("Downsampling tips for host ", host, "...\n", sep="")
  
  tips.from.host <- which(host.ids==host)
  labels.from.host <- tree$tip.label[tips.from.host]
  
  if(!no.read.counts){
    read.counts <- as.numeric(sapply(labels.from.host, function(tip) read.count.from.label(tip, tip.regex))) 
    
    if(any(is.na(read.counts))){
      stop("Some tips have no read counts despite matching a tip regex which requires them. This sometimes happens when the tip regex does not end with a $")
    }
    
  } else {
    read.counts <- rep(1, length(labels.from.host))
  }
  
  total.reads <- sum(read.counts)
  
  # if you have less reads than you are trying to sample, return them all
  
  if(total.reads < number){
    if(!exclude.underrepresented){
      # warning("Insufficient reads for downsampling host ",host," at a count of ",number," in tree ",name,", returning all.\n", sep="")
      return(list(blacklist=vector(), map=NULL))
    } else {
      # warning("Insufficient reads for downsampling host ",host," at a count of ",number," in tree ",name,", blacklisting this host.\n", sep="")
      
      if(rename){
        label.map <- lapply(1:length(labels.from.host), function(x){
          paste0(labels.from.host[x],"_X_UNDERREPRESENTED")
        })
        
        names(label.map) <- labels.from.host
        
        labels.from.host <- unlist(label.map[labels.from.host])
        
        return(list(blacklist=labels.from.host, map=label.map))
      } else {
        return(list(blacklist=labels.from.host, map=NULL))
      }
    }
  }
  
  # whether intentionally or not, rmvhyper does not work as hoped if n has length 1
  
  if(length(read.counts)==1){
    sample <- matrix(number, nrow=1, ncol=1)
  } else {
    sample <- rmvhyper(nn=1, n = read.counts, k = number)
  }
  
  if(rename){
    if(verbose) cat("Renaming tree tips with new read counts...\n", sep="")
    
    # need a new regex to do the find and replace. This is probably preferable to asking the user for one.
    
    # remove the first opening bracket
    new.tip.regex <- sub("\\(", "", tip.regex)
    # remove the second opening bracket
    new.tip.regex <- sub("\\(", "", new.tip.regex)
    # remove the first closing bracket
    new.tip.regex <- sub("\\)", "", new.tip.regex)
    # remove the second closing bracket
    new.tip.regex <- sub("\\)", "", new.tip.regex)
    # Start a new capture group start after the third one
    new.tip.regex <- sub("\\)", "\\)\\(", new.tip.regex)
    # End a new capture group start before the third one
    new.tip.regex <- sub("\\(", "\\)\\(", new.tip.regex)
    # Deal with the beginning
    if(substr(new.tip.regex,1,1)=="^"){
      new.tip.regex <- gsub("\\^", "\\^\\(", new.tip.regex)
    } else {
      new.tip.regex <- paste0("^(", new.tip.regex)
    }
    # Deal with the end
    if(substr(new.tip.regex,nchar(new.tip.regex),nchar(new.tip.regex))=="$"){
      new.tip.regex <- gsub("\\$", "\\)\\$", new.tip.regex)
    } else {
      new.tip.regex <- paste0(new.tip.regex,')$')
    }
    
    label.map <- lapply(1:length(labels.from.host), function(x){
      gsub(new.tip.regex, paste0("\\1", sample[x], "\\3"), labels.from.host[x])
    })
    
    never.selected.tips <- which(sample[1,]==0)
    
    label.map[never.selected.tips] <- paste0(label.map[never.selected.tips],"_X_DOWNSAMPLED")

    names(label.map) <- labels.from.host
    
    labels.from.host <- unlist(label.map[labels.from.host])
  } else {
    label.map <- NULL
  }
  
  sampled.names <- labels.from.host[which(sample > 0)]
  
  if(verbose){
    cat("Blacklisting ",length(setdiff(labels.from.host, sampled.names))," tips from host ",host,".\n",sep="")
  }

  return(list(blacklist = setdiff(labels.from.host, sampled.names), map=label.map))
}

# This returns the tips to get rid of (i.e. the blacklist) for a given tree

#' @keywords internal
#' @export downsample.tree

downsample.tree <- function(ptree, hosts.to.include, max.reads, rename = F, exclude.underrepresented = F, no.read.counts = T, tip.regex, seed=NA, verbose=F) {

  if(verbose){
    cat("Downsampling reads on tree ",ptree$id," to ",max.reads," per host...\n",sep = "")
  }


  if(!is.na(seed)){
    set.seed(seed)
  }

  if(is.null(ptree$tree)){
    
    if(verbose){
      cat("Loading tree...\n",sep = "")
    }
  
    tree <- read.tree(ptree$tree.file.name)

  } else {
    tree <- ptree$tree
  }
  
  
  if(!is.null(ptree$blacklist)){
    blacklist <- ptree$blacklist
  } else if(!is.null(ptree$blacklist.file.name)){
    blacklist <- vector()
  
    if(file.exists(ptree$blacklist.file.name)){
      if(verbose) cat("Reading existing blacklist file",ptree$blacklist.file.name,'\n')
      blacklisted.tips <- read.table(ptree$blacklist.file.name, sep=",", header=F, stringsAsFactors = F, col.names="read")
       
      if(nrow(blacklisted.tips)>0){
        blacklist <- c(blacklist, sapply(blacklisted.tips, get.tip.no, tree=tree))
      }
    } else {
      warning("File ",ptree$blacklist.file.name," does not exist; skipping.", sep="")
    }	
  } else {
    blacklist <- vector()
  }

  
  
  tree.1 <- drop.tip(tree, blacklist)

  host.ids      <- sapply(tree.1$tip.label, function(tip) host.from.label(tip, tip.regex))
  
  hosts.present <- host.ids[!is.na(host.ids)]
  hosts.present <- unique(hosts.present)
  
  if(!is.null(hosts.to.include)){
    if(length(intersect(hosts.to.include, hosts.present)) == 0){
      stop(paste("No entries in give vector of hosts are actually present in the tree", sep=""))
    }
    
    if(length(setdiff(hosts.to.include, hosts.present)) > 0){
      warning(paste("Not all hosts in file requested are present in the tree", sep=""))
    }
    hosts.to.include <- intersect(hosts.to.include, hosts.present)
    
  } else {
    hosts.to.include <- hosts.present
  }
  
  new.tip.labels <- tree$tip.label
  excluded <- unlist(lapply(hosts.to.include, function(x){
    result <- downsample.host(x, 
                              tree=tree.1, 
                              number = max.reads, 
                              tip.regex=tip.regex, 
                              host.ids=host.ids, 
                              rename, 
                              exclude.underrepresented, 
                              no.read.counts, 
                              name = ptree$id, 
                              verbose)
    if(rename & !is.null(result$map)){
      new.tip.labels <<- sapply(new.tip.labels, function(y){
        if(y %in% names(result$map)){
          result$map[[y]]
        } else {
          y
        }
      })
    }
    return(result$blacklist)
  }))
  

  tree$tip.label <- new.tip.labels
  ptree$tree <- tree
  
  excluded.nos  <- which(tree$tip.label %in% excluded) 
  
  if(length(excluded.nos)>0){
    ptree$bl.report$status[excluded.nos] <- "downsampled"
    ptree$bl.report$kept[excluded.nos]   <- F
  }
  
  new.blacklist <- c(blacklist, excluded.nos)
  ptree$blacklist <- new.blacklist
  
  ptree
}

