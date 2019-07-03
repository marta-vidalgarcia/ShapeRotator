\name{join.arrays}
\title{join.arrays}
\alias{join.arrays}

\description{
Joins the two rotated arrays and returns a single array
}
\usage{
join.arrays(data.1, data.2)

}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data.1}{
an object of the class "array" that contains the rotated three-dimensional landmarks for structure #1 (e.g. skull)
}
  \item{data.2}{
an object of the class "array" that contains the rotated three-dimensional landmarks for structure #2 (e.g. mandible)
}

}

\value{
This function returns an array of the joined two arrays comprising the two rotated structures (e.g. mandible and skull).
}


\details{
This function join the two arrays that were rotated using double.rotation() or simple.rotation() and return a big array from joining the two smaller arrays. It is possible to translate one of the objects after the rotation, in the case of not wanting them in contact to one another, or in a slightly different location. The vector used in the translation might have to change depending on the angle used in the rotation function, in order to represent biologically-relevant positions of the two structures relative to one another.
}

\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~


\examples{
rotation_45_j <- join.arrays(rotation_45$rotated1, rotation_45$rotated2)

# The position of the skull relative to the mandible can be changed after the rotation process with the vector 'skull_translate' by including a translation within the function
skull_translate <- c(2.2,-0.1, 0)
rotation_45_t <- join.arrays(rotation_45$rotated1, translate(rotation_45$rotated2, land.e , skull_translate))

plot.rotation.3D(rotation_45_t, data.1, data.2, specimen.num = 1)
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
