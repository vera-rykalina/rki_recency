\dontrun{
require(data.table)
require(tidyverse)
require(phyloscannerR)

#	specify path to phyloscanner_analyse_trees
prog.phyloscanner_analyse_trees <- '/Users/Oliver/git/phyloscanner/phyloscanner_analyse_trees.R'
#	specify out directory
outdir <- '/Users/Oliver/sandbox/DeepSeqProjects/RakaiPopSample_phsc_out190512'
#	specify valid input arguments to phyloscanner_analyse_trees
valid.input.args <- cmd.phyloscanner.analyse.trees.valid.args(prog.phyloscanner_analyse_trees)

#	set phyloscanner variables
#	arguments as used for the Rakai population-based analysis
control	<- list()
control$allow.mt <- TRUE				
control$alignment.file.directory = NULL 
control$alignment.file.regex = NULL
control$blacklist.underrepresented = FALSE	
control$count.reads.in.parsimony = TRUE
control$distance.threshold <- '0.025 0.05'
control$do.dual.blacklisting = FALSE					
control$duplicate.file.directory = NULL
control$duplicate.file.regex = NULL
control$file.name.regex = "^\\D*([0-9]+)_to_([0-9]+)\\D*$"
control$guess.multifurcation.threshold = FALSE
control$max.reads.per.host = 50
control$min.reads.per.host <- 30
control$min.tips.per.host <- 1	
control$multifurcation.threshold = 1e-5
control$multinomial= TRUE
control$norm.constants = NULL
control$norm.ref.file.name = system.file('HIV_DistanceNormalisationOverGenome.csv',package='phyloscannerR')
control$norm.standardise.gag.pol = TRUE
control$no.progress.bars = TRUE
control$outgroup.name = "REF_CPX_AF460972"
control$output.dir = outdir
control$parsimony.blacklist.k = 20
control$prune.blacklist = FALSE
control$post.hoc.count.blacklisting= TRUE
control$ratio.blacklist.threshold = 0 
control$raw.blacklist.threshold = 20					
control$recombination.file.directory = NULL
control$recombination.file.regex = NULL
control$relaxed.ancestry = TRUE
control$sankoff.k = 20
control$sankoff.unassigned.switch.threshold = 0
control$seed = 42
control$splits.rule = 's'
control$tip.regex = "^(.*)_fq[0-9]+_read_([0-9]+)_count_([0-9]+)$"
control$tree.file.regex = "^ptyr[0-9]+_InWindow_([0-9]+_to_[0-9]+)\\.tree$"
control$use.ff = FALSE
control$user.blacklist.directory = NULL 
control$user.blacklist.file.regex = NULL
control$verbosity = 1

#
#	Example 1: make bash for one file
#
tree.input <- system.file(file.path('extdata','Rakai_run192_trees.zip'),package='phyloscannerR')
control$output.string <- 'Rakai_run192'
cmd <- cmd.phyloscanner.analyse.trees(prog.phyloscanner_analyse_trees, 
		tree.input, 
		control,
		valid.input.args=valid.input.args)
cat(cmd)

#
#	Example 2: make bash for many files
#
#	download the phyloscanner tree of the Rakai population-based analysis  
tmp <- "https://datadryad.org/bitstream/handle/10255/dryad.208473/Dataset_S1.tar?sequence=1"
#	specify directory to untar public data
tree.dir <- "RakaiPopSample_deepseqtrees"
#	download and untar
download.file(tmp, destfile="Dataset_S1.tar", method="curl")
untar("Dataset_S1.tar", exdir=tree.dir, extras='-xvf')	
#	list zipped tree files. One zip file contains the viral trees of individuals in one putative transmission network. 
df <- tibble(F=list.files(tree.dir))
df <- df %>% 
		mutate(TYPE:= gsub('ptyr([0-9]+)_(.*)','\\2', F),
			RUN:= as.integer(gsub('ptyr([0-9]+)_(.*)','\\1', F))) %>%
		mutate(TYPE:= gsub('^([^\\.]+)\\.[a-z]+$','\\1',TYPE)) %>%
		spread(TYPE, F) %>%
		set_names(~ str_to_upper(.))
#	make one bash script for processing the viral trees of individuals in one putative transmission network.
cmds <- vector('list',nrow(df))
for(i in seq_len(nrow(df)))
{
	control$output.string <- paste0('ptyr',df$RUN[i])	
	tree.input <- file.path(indir, df$TREES_NEWICK[i])
	cmd <- cmd.phyloscanner.analyse.trees(prog.phyloscanner_analyse_trees, 
			tree.input, 
			control,
			valid.input.args=valid.input.args)
	cmds[[i]] <- cmd		
}	
#	output
cat(cmds[[100]])
}
