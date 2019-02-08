\name{norm_3D}
\title{norm_3D}
\alias{norm_3D}

\description{
This is an internal function for the ShapeRotator tool. This function allows you to obtain the length of a vector u.
}
\usage{
norm_3D(u)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{u}{
 an object of the class "vector" in three-dimensional Cartesian coordinate system (x, y, z)
}
}

\value{
an object of the class "numeric" as calculated by √(x^2, y^2, z^2)

}

\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~


\examples{
vector_1 <- c(1, 2, 3)
vector_1_length <- norm_3D(vector_1)
>vector_1_length
[1] 3.741657

vector_2 <- c(1, 1, 1)
vector_2_length <- norm_3D(vector_2)
vector_2_length
[1] 1.732051
}

\references{
M. Vidal-García, L. Bandara and J.S. Keogh. ShapeRotator: An R tool for standardized rigid rotations of articulated three-dimensional structures with application for geometric morphometrics. Ecology and Evolution. In press. DOI: 10.1002/ece3.4018
}
\author{
Marta Vidal-Garcia (marta.vidalga@gmail.com)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
