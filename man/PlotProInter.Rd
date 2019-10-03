\name{PlotProteinInteractions}
\alias{PlotProteinInteractions}
\title{
Plot Protein Interactors for Selected Proteins.
}
\description{
Plot protein interactions and comparing two proteins whether have the same interactors.
}
\usage{
PlotProteinInteractions(protein1, protein2)
}
\value{
Return a plot to show the interactions between selected proteins and interactors.
}

\arguments{
  \item{protein1}{Selet interested protein uniprot id}
  \item{protein2}{Selet interested protein uniprot id}
}

\examples{
PlotProteinInteractions("Q9UKV5", "Q9UHB7")
}

