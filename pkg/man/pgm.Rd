\name{pgm}
\alias{pgm}
\title{Parse PGM Images}

\description{
  This function parses a PGM image file.
}

\usage{
  pgm(filename)
}

\arguments{
  \item{filename}{path to the pgm file}
}

\value{
  An instance of the class \sQuote{pgm} with the following items
  \item{width}{the image width in pixels}
  \item{height}{the image height in pixels}
  \item{maxval}{the maximum pixel value}
  \item{data}{a matrix of pixel values}
}

\details{
  This function parses a PGM image file. The parser only checks for the correct \sQuote{magic number}. No other error checking is performed.
}

\author{Matt Shotwell}
