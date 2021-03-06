% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/src.r
\docType{methods}
\name{multirec}
\alias{multirec}
\title{Multi components reconstruction}
\usage{
multirec(
  object,
  imagedim,
  B = NULL,
  mask = NULL,
  midx = 1,
  comps = NULL,
  XY = c("X", "Y", "XY")[1],
  signflip = FALSE
)
}
\arguments{
\item{object}{an object of class \code{msma}. Usually, a result of a call to \code{msma}}

\item{imagedim}{a vector for original dimension.}

\item{B}{a list or a matrix.}

\item{mask}{a list or a matrix.}

\item{midx}{a block number.}

\item{comps}{a component sequence to be plotted.}

\item{XY}{a character, indicating "X" or "Y". "XY" for the scatter plots using X and Y scores from \code{msma}.}

\item{signflip}{a logical if the sign in the block is flipped to pose the super as possitive.}
}
\description{
This is a function that returns the weight vector of multiple components obtained by the \code{msma} function applied after dimension reduction by the radial basis function to the same dimension as the original image.
}
\details{
\code{multirec} requires the output result of \code{msma} function.
}
\examples{
\donttest{
data(baseimg)
data(diffimg)
data(mask)
img1 = simbrain(baseimg = baseimg, diffimg = diffimg, mask=mask)
B1 = rbfunc(imagedim=img1$imagedim, seppix=2, hispec=FALSE, mask=img1$brainpos)
SB1 = basisprod(img1$S, B1)
fit111 = msma(SB1, comp=2)
ws = multirec(fit111, imagedim=img1$imagedim, B=B1, mask=img1$brainpos)
}
}
