% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/general_functions.R
\name{find.bam.and.references}
\alias{find.bam.and.references}
\title{Find bam and corresponding reference files}
\usage{
find.bam.and.references(
  data.dir,
  regex.person = "^([A-Z0-9]+-[A-Z0-9]+)-.*$",
  regex.bam = "^(.*)\\\\.bam$",
  regex.ref = "^(.*)_ref\\\\.fasta$",
  verbose = 1
)
}
\arguments{
\item{data.dir}{Full path of data directory}

\item{regex.person}{Regular expression with one set of round brackets, which identifies the person ID in the file name of bams and references}

\item{regex.bam}{Regular expression that identifies bam files, with one set of round brackets that identifies the sample ID.}

\item{regex.ref}{Regular expression that identifies ref files, with one set of round brackets that identifies the sample ID.}
}
\value{
tibble with rows 'IND' (individual identifier), 'SAMPLE' (sample identifier), 'BAM' (bam file), and 'REF' (reference file).
}
\description{
This function finds bam and corresponding reference files in a given directory,
and groups them by a common sample ID as well as by an individual ID.
}
