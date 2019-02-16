\name{double.rotation}
\title{double.rotation}
\alias{double.rotation}

\description{
Rotates two arrays in a two-points articulation (e.g. skull and mandible)
}
\usage{
double.rotation(data.1, data.2, land.a, land.b, land.c, land.d, land.e, land.f, land.g, land.h, angle)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data.1}{
an object of the class "array" that contains translated three-dimensional landmarks for structure #1 (e.g. skull)
}
  \item{data.2}{
an object of the class "array" that contains translated three-dimensional landmarks for structure #2 (e.g. mandible)
}
}

\value{
Returns a list of two object where $rotated1 gives the first rotated array and $rotated2 gives the second
}

\Details{
This function allows the rigid rotation of two structures from a two-points articulation (e.g. skull and mandible). It allows to select the desired angle between the two structures so that we can perform the rigid rotation of each structure positioning them in the selected angle in relation to each other, by providing the desired angle in degrees (from 0° to 360°). Even though each of these rotations are calculated internally (only the four landmarks per structure and the desired angle between them need to be provided), it will be beneficial to choose landmarks that are spatially arranged in a way that facilitates the rotation process, and results in the rotating multi-structure being placed in a biologically-relevant angle between each sub-structure. ShapeRotator will give a warning message if the landmarks chosen are not optimal (i.e. the position and orientation of the two structures vary across the specimens due to morphological variability and not because of the rigid rotation).
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~


\examples{
# We can rigidly rotate each structure relative to each other in multiple angles (e.g. 0˚, 15˚, 45˚, 90˚, or 120˚)
rotation_0 = double.rotation(data.1, data.2, land.a, land.b, land.c, land.d, land.e, land.f, land.g, land.h, 0)

rotation_15 = double.rotation(data.1, data.2, land.a, land.b, land.c, land.d, land.e, land.f, land.g, land.h, 15)

rotation_45 = double.rotation(data.1, data.2, land.a, land.b, land.c, land.d, land.e, land.f, land.g, land.h, 45)

rotation_90 = double.rotation(data.1, data.2, land.a, land.b, land.c, land.d, land.e, land.f, land.g, land.h, 90)

rotation_120 = double.rotation(data.1, data.2, land.a, land.b, land.c, land.d, land.e, land.f, land.g, land.h, 120)

# The function returns a list of two objects that correspond to each rotated structure
rotation_0$rotated1 #rotated dataset 1 (e.g. skull)
rotation_0$rotated2 #rotated dataset 2 (e.g. mandible)
}

\references{
M. Vidal-García, L. Bandara and J.S. Keogh. 2018. ShapeRotator: An R tool for standardized rigid rotations of articulated three-dimensional structures with application for geometric morphometrics. Ecology and Evolution. DOI: 10.1002/ece3.4018

}
\author{
Marta Vidal-Garcia (marta.vidalga@gmail.com)
}

