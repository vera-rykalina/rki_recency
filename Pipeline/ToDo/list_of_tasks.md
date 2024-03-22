# List of tasks and ideas to be tested

- Try a fresh installation of phyloscannerR to see if it works:
devtools::install on the phyloscannerR directory 

- Try to install R, phyloscannerR and all its dependencies, using conda -> checked: no conda installation for phyloscannerR

- Distribute the scripts and tool's scripts as intended by shiver & phyloscanner & HIV-PhyloTSI (now all scripts are merged within bin) -> checked: shouldn't be done so!

- Clarify if --x-iqtree is an option now 7-> checked: not yet

- Create custom HIV_DistanceNormalisationOverGenome.csv and compare it with an available file

- Compare IVA (reproducipility fails) with SPADES

- Double check MAF-associated sctipts -> done: Do not replace NaN with 0's.

- Optimise HPC resources for all process if needed -> mostly done

- Clarify excluded windows: 9470,9720,9480,9730,9490,9740 from windows250_VR_norms_oneline.txt (Tanya)

- Clarify the option --read-names-2. It does not exists anymore! -> done: the default setting in the current version

- Clafiry the differences in testpatstats.csv: Prop.gp.1 - Prop.gp.215 -> done: sample-specific differences

- Check on 4refs_HXB2_ABCD.fasta vs 2refs_HXB2_C.BW.fasta for --alignment-of-other-refs. Should I expect any dramatic difference? -> done: it is fine to use 2refs fasta

- Work on IVA's unexpected error (errorStrategy might help here): 
The following command failed with exit code 1
  smalt map -r 1  -n 16 -O -i 800 -y 0.5 iteration.4.10.1.map.map_index iteration.3.filtered.subiter.10.reads_1.fa iteration.3.filtered.subiter.10.reads_2.fa | samtools view -bS -T iteration.4.10.1.map.ref.fa  - > iteration.4.10.1.map.bam

- Check on --amplicons flag (compare: --amplicons True and --amplicons False)