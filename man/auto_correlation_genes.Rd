% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Estimators.R
\name{auto_correlation_genes}
\alias{auto_correlation_genes}
\title{Get gene auto correlation}
\usage{
auto_correlation_genes(
  dt,
  dist = 6,
  by.codon = TRUE,
  codon.vs.nt = FALSE,
  genes,
  fun = acf,
  mean = FALSE
)
}
\arguments{
\item{dt}{a data.table of counts}

\item{dist}{numeric, default 6. Distance in nt or codons to check}

\item{by.codon}{logical, default TRUE. Else by nt}

\item{codon.vs.nt}{logical, default FALSE Else convert codon to nt space
(1,0,0)}

\item{genes}{numeric vector, genes, index of which gene this is from}

\item{fun}{which auto correlation function, default stats::acf}

\item{mean}{logical, default FALSE Else get only mean acf per position.}
}
\value{
a data.table of counts
}
\description{
Get gene auto correlation
}
