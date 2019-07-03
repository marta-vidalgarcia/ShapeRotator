\name{plot.rotation.3D}
\title{plot.rotation.3D}
\alias{plot.rotation.3D}

\description{
Plots the two rotated datasets in 3D
}
\usage{
plot.rotation.3D(joined.data, data.1, data.2, specimen.num)

}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{joined.data}{
an object of the class "array" that contains the joined rotated arrays (data.1three-dimensional landmarks for structure #1 (e.g. skull). It can either be the original dataset, the translated array, or the rotated array
}
  \item{data.1}{
an object of the class "array" that contains three-dimensional landmarks for structure #1 (e.g. skull). It can either be the original dataset, the translated array, or the rotated array
}
  \item{data.2}{
an object of the class "array" that contains three-dimensional landmarks for structure #2 (e.g. mandible). It can either be the original dataset, the translated array, or the rotated array
}
  \item{specimen.num}{
an integer object that indicates which specimen will be plotted. The number of the specimen to be plotted corresponds to the position in the "joined.data" array
}
}

\value{
This function returns a 3D plot of the rotated structures (e.g. mandible and skull).
}

\details{
Plots the two rotated datasets in 3D for a specimen (specimen.num), and each structure is plotted with a different colour. Default is black and red. This uses the function plot3D from the R package plot3D (https://cran.r-project.org/web/packages/plot3D/) internally.
}

\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~


\examples{
plot.rotation.3D(rotation_45, data.1, data.2, specimen.num = 1)
}

\references{
http://www.rforscience.com/rpackages/visualisation/oceanview/

http://www.rforscience.com/rpackages/visualisation/plot3d/

M. Vidal-García, L. Bandara and J.S. Keogh. 2018. ShapeRotator: An R tool for standardized rigid rotations of articulated three-dimensional structures with application for geometric morphometrics. Ecology and Evolution. DOI: 10.1002/ece3.4018

}
\author{
Marta Vidal-Garcia (marta.vidalga@gmail.com)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
