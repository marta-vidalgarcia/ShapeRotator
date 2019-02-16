\name{match.datasets}
\title{match.datasets}
\alias{match.datasets}

\description{
Matches two arrays by their dimnames
}
\usage{
match.datasets(data.1, data.2)

}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data.1}{
an object of the class "array" that contains three-dimensional landmarks for structure #1 (e.g. skull)
}
  \item{data.2}{
an object of the class "array" that contains three-dimensional landmarks for structure #2 (e.g. mandible)
}
}

\value{
This function returns a list with the two translated arrays of dimensions (p, 3, n), in which dimnames associated with each specimen (n) are in the same order.
}

\Details{
Takes the two translated arrays and joined them by dimnames. Unmatched datasets (specimens missing in data.1 or data.2) are deleted.
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~


\examples{
# We need to match the dimnames of the translated datasets, but we have two extra specimens in data.1_t
length(dimnames(data.1_t)[[3]]) # n=12
length(dimnames(data.2_t)[[3]])  # n=10

matched <- match.datasets(data.1_t, data.2_t)

# When we matched the arrays, only matrices from specimens present in both datasets will be kept
length(dimnames(matched$matched1)[[3]]) # n=10
length(dimnames(matched$matched2)[[3]])  # n=10

}

\references{
M. Vidal-García, L. Bandara and J.S. Keogh. 2018. ShapeRotator: An R tool for standardized rigid rotations of articulated three-dimensional structures with application for geometric morphometrics. Ecology and Evolution. DOI: 10.1002/ece3.4018

}
\author{
Marta Vidal-Garcia (marta.vidalga@gmail.com)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
