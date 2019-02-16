### DOUBLE.ROTATION
# function rotation for a two-points articulation (2 different rotations for each translated dataset)
# Returns a list where $rotated1 gives the first object and $rotated2 gives the second
double.rotation <- function(data.1, data.2, land.a, land.b, land.c, land.d, land.e, land.f, land.g, land.h, angle)
{
  # First rotation.
  rot.array.1 <- array(data = NA, dim = c(dim(data.1)[1], dim(data.1)[2], dim(data.1)[3]),dimnames = list(NULL, NULL, c(dimnames(data.1)[[3]])));
  angle.v <-vector.angle(angle)


  for( i in 1:(dim(data.1)[3]) ){
    # Landmark A should be sitting on (0,0,0)

    # Untranslated specimen?
    if ( ! isTRUE(  all.equal((data.1[,,i])[land.a,] , c(0,0,0)) ) ) {
      warning (sprintf("Landmark A is not at the origin (data.1[,,%i])[%i,] = (%f,%f,%f)?", i, land.a, (data.1[,,i])[land.a,1], (data.1[,,i])[land.a,2], (data.1[,,i])[land.a,3]), immediate. = immediate_on)
    }

    # Rotate so that landmark C is sitting on the Z axis
    rot.array.1[,,i] <- rotate.orientation(data.1[,,i], land.c, c(0,0,1))


    #Sanity check: we need to have rot.array.1[,,i][land.c,] = (0, 0, y)
    if ( ! isTRUE(  all.equal( (rot.array.1[,,i])[land.c,1] , 0)) || ! isTRUE(  all.equal( (rot.array.1[,,i])[land.c,2] , 0)) ) {
      warning (sprintf("First rotation gone wrong in (rot.array.1[,,%i])[%i,] = (%f,%f,%f)?", i, land.c, (rot.array.1[,,i])[land.c,1], (rot.array.1[,,i])[land.c,2], (rot.array.1[,,i])[land.c,3]), immediate. = immediate_on)
    }

    # Now, we want to move landmark B to the X-Z axis.
    unit.landmark.b.pos.proj = unit_3D( c( (rot.array.1[,,i])[land.b,1], (rot.array.1[,,i])[land.b,2],0) )	# Unit vector with landmark.b.pos projected to the X-Y plane.
    rotmat <- rotmat_3D( cross_3D(unit.landmark.b.pos.proj, c(1,0,0)), angle_3D(unit.landmark.b.pos.proj, c(1,0,0)) )
    rot.array.1[,,i] <- rotveclist_3D ( rotmat, rot.array.1[,,i])

    # We want B_x > 0.
    if ( isTRUE((rot.array.1[,,i])[land.b,1] < 0) ) {
      rotmat <- rotmat_3D( c(0,0,1), pi)
      rot.array.1[,,1] <- rotveclist_3D ( rotmat, rot.array.1[,,i])
    }


    # Let us now check that the landmark C is in the right spot. That is, we want B_y < D_y, if not, we rotate by angle pi in the X axis.
    if ( isTRUE((rot.array.1[,,i])[land.d,2] < (rot.array.1[,,i])[land.b,2]) ) {
      rotmat <- rotmat_3D( c(1,0,0), pi)
      rot.array.1[,,i] <- rotveclist_3D ( rotmat, rot.array.1[,,i])
    }

    # Landmark B has to be 0 on the Y axis, and we want to make sure that the orientation is right
    if ( ! isTRUE(  all.equal( (rot.array.1[,,i])[land.b,2] , 0)) || isTRUE((rot.array.1[,,i])[land.d,2] < (rot.array.1[,,i])[land.b,2]) || isTRUE(rot.array.1[,,i][land.b,1] < 0) ) {
      warning (sprintf("Second rotation gone wrong in (rot.array.1[,,%i])[%i,] = (%f,%f,%f)?",
                       i, land.b, (rot.array.1[,,i])[land.b,1], (rot.array.1[,,i])[land.b,2], (rot.array.1[,,i])[land.b,3]), immediate. = immediate_on)
    }

    # Rotate body to be the correct angle in the X-Y plane (so around the Z axis)
    rotmat <- rotmat_3D( c(0,0,1), angle * pi/180)
    rot.array.1[,,i] <- rotveclist_3D ( rotmat, rot.array.1[,,i])
  }


  # Establishing the rotation axes for data.2
  rot.array.2 <- array(data = NA, dim = c(dim(data.2)[1], dim(data.2)[2], dim(data.2)[3]),dimnames = list(NULL, NULL, c(dimnames(data.2)[[3]])));


  #Rotating for data.2 (ROTATION number 1)
  for( i in 1:(dim(data.2)[3]) ){
    # Untranslated specimen?
    if ( ! isTRUE(  all.equal((data.2[,,i])[land.e,] , c(0,0,0)) ) ) {
      warning (sprintf("Landmark E is not at the origin (data.2[,,%i])[%i,] = (%f,%f,%f)?", i, land.e, (data.2[,,i])[land.e,1], (data.2[,,i])[land.e,2], (data.2[,,i])[land.e,3]), immediate. = immediate_on)
    }

    # Rotate so that landmark F is sitting on the Z axis in the same direction as landmark C
    rot.array.2[,,i] <- rotate.orientation(data.2[,,i], land.f, unit_3D(rot.array.1[,,i][land.c,]) )

    #Sanity check: we need to have rot.array.1[,,i][land.c,] in the same direction as rot.array.2[,,i][land.f,]
    if ( ! isTRUE(  all.equal(unit_3D(rot.array.1[,,i][land.c,]), unit_3D(rot.array.2[,,i][land.f,]))) ) {
      warning (sprintf("Rotation of landmark F gone wrong?  (rot.array.1[,,%i])[%i,] = (%f,%f,%f), (rot.array.2[,,%i])[%i,] = (%f,%f,%f)",
                       i, land.c, (rot.array.1[,,i])[land.c,1], (rot.array.1[,,i])[land.c,2], (rot.array.1[,,i])[land.c,3],
                       i, land.f, (rot.array.2[,,i])[land.f,1], (rot.array.2[,,i])[land.f,2], (rot.array.2[,,i])[land.f,3]),
               immediate. = immediate_on)
    }


    # Now, we want to move landmark G to the X-Z axis.
    unit.landmark.g.pos.proj = unit_3D( c( (rot.array.2[,,i])[land.g,1], (rot.array.2[,,i])[land.g,2],0) )	# Unit vector with landmark.g.pos projected to the X-Y plane.
    rotmat <- rotmat_3D( cross_3D(unit.landmark.g.pos.proj, c(1,0,0)), angle_3D(unit.landmark.g.pos.proj, c(1,0,0)) )
    rot.array.2[,,i] <- rotveclist_3D ( rotmat, rot.array.2[,,i])

    # We want G_x > 0, because it should be in the same direction as B_x.
    if ( isTRUE((rot.array.2[,,i])[land.g,1] < 0) ) {
      rotmat <- rotmat_3D( c(0,0,1), pi)
      rot.array.2[,,1] <- rotveclist_3D ( rotmat, rot.array.2[,,i])
    }

    # This completely constrains the problem, we should now have G_y > H_y (because correcting for this would mean that landmark C ~ landmark F).
    # Nevertheless, we should check and issue a warning.

    if (  isTRUE(  rot.array.2[,,i][land.g,2] < rot.array.2[,,i][land.h,2]) ){
      warning (sprintf("Unsatisfied constraint in second rotation? (rot.array.2[,,%i])[%i,] = (%f,%f,%f), (rot.array.2[,,%i])[%i,] = (%f,%f,%f)",
                       i, land.g, (rot.array.2[,,i])[land.g,1], (rot.array.2[,,i])[land.g,2], (rot.array.2[,,i])[land.g,3],
                       i, land.h, (rot.array.2[,,i])[land.h,1], (rot.array.2[,,i])[land.h,2], (rot.array.2[,,i])[land.h,3]),
               immediate. = immediate_on)
    }
  }

  return ( list("rotated1" = rot.array.1, "rotated2" = rot.array.2))

}

