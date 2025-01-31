% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/PCIP_getBreakPoints.R
\name{PCIP_getBreakPoints}
\alias{PCIP_getBreakPoints}
\title{Get Read Breakpoints and ShearSites}
\usage{
PCIP_getBreakPoints(
  PAF = NULL,
  lengthTarget = NULL,
  targetName = NULL,
  gapAlignment = NA,
  distanceToLTR = NA,
  returnFILTEREDout = FALSE
)
}
\arguments{
\item{PAF}{tibble. Pairwise Alignment file (PAF) generated with minimap2, red with \code{readPairwiseAlignmentFile} and formated with \code{PCIP_filter}}

\item{lengthTarget}{numeric. Length of the TARGET chromosome in base pairs}

\item{targetName}{character. Name of the TARGET chromosome}

\item{gapAlignment}{numeric or NA. Maximum mean gap tolerated between substrings. For debugging purpose. NA deactivate this option.}

\item{distanceToLTR}{numeric or NA. Exclude reads with alignment gaps near the LTR edge > distanceToLTR. For debugging purpose.}

\item{returnFILTEREDout}{logical. Directly returns reads with large distanceToLTR. For debugging purpose.}
}
\value{
Returns a tibble. Field description:
\itemize{
\item readID (character) ID of the read.
\item seqnames.genome (character) HOST Chromosome name.
\item shearSite.genome (numeric) ShearSite position in the HOST genome.
\item strand.target (character) TARGET substring orientation.
\item strand.genome (character) GENOME substring orientation.
\item context (character) TARGET - HOST substrings orientations related to each other.
\item ligation (character) TARGET - HOST structure of the read as binary code (0 = HOST and 1 = TARGET).
\item minDistanceEdge (double) Minimal distance to one TARGET edge.
\item edge5_breakPoint (double) Position of the TARGET 5' edge in the HOST genome.
\item edge3_breakPoint (double) Position of the TARGET 3' edge in the HOST genome.
}
}
\description{
Extract TARGET-HOST junctions (= breakpoints) and DNA fragmentation site (= shearsite) created by sonication.
}
\author{
Vincent Hahaut
}
\keyword{PCIP}
