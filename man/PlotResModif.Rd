\name{PlotResModification}
\alias{PlotResModification}
\title{Plot the Residue Modification Position in Selected Ubiquitin Sequence}

\describe{
\item{Protein1}{Selet interested protein uniprot id}
\item{Protein2}{Selet interested protein uniprot id}
}

\usage{
PlotResModification(protein1, protein2)
}

\value{
Return a plot to show the residue modification positions in selected protein sequences.
}

\description{
The output is a plot which shows the amino acid modification positions in the selected proteins' sequences.
}

\examples{
PlotResModification("Q9UHB7", "Q9UKV5");
}
