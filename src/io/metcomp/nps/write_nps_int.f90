! This module contains a subroutine to write unformatted Fortran output compatible
! with the intermediate file format used by WPS and NPS.
! It was taken from Ghassan Alaka's "unpack_netcdf" program.

subroutine write_2d_slab(ofile, version, hdate, xfcst, &
                map_source, field, units, desc, xlvl, &
                iproj, startloc, startlat, startlon, dx, dy, &
                deltalat, deltalon, earth_radius, truelat1, &
                truelat2, xlonc, nlats, is_wind_grid_rel, &
                slab, nx, ny)
    implicit none
    ! args
    character(len=300), intent(in)  :: ofile
    integer, intent(in)             :: version
    character(len=24), intent(in)   :: hdate
    real*4, intent(in)              :: xfcst    
    character(len=32), intent(in)   :: map_source
    character(len=9), intent(in)    :: field
    character(len=25), intent(in)   :: units
    character(len=46), intent(in)   :: desc
    real*4, intent(in)              :: xlvl
    integer, intent(in)             :: iproj    
    character(len=8), intent(in)    :: startloc
    real*4, intent(in)              :: startlat, startlon, dx, dy, deltalat, deltalon
    real*4, intent(in)              :: earth_radius, truelat1, truelat2
    real*4, intent(in)              :: xlonc, nlats
    logical, intent(in)     :: is_wind_grid_rel
    real*4, dimension (1:nx,1:ny), intent(in) :: slab
    integer, intent(in)             :: nx, ny    
    
    ! local variables
    integer, parameter  :: ounit = 10           ! Outfile unit

    open (ounit, FILE=ofile, form='unformatted',access='append')

    !  1) WRITE FORMAT VERSION
    write(unit=ounit) version
    
    !print *, 'version, hdate, xfcst, map_source, field  = ', version,  hdate, xfcst, map_source, field
    !print *, 'units, desc, xlvl, nx, ny, iproj = ', units, desc, xlvl, nx, ny,  iproj
    !print *, 'startloc, startlat, startlon = ', startloc, startlat, startlon
    !print *, 'deltalat, deltalon, earth_radius = ', deltalat, deltalon, earth_radius
    
    !  2) WRITE METADATA
    ! Cylindrical equidistant
    if (iproj == 0) then
          write(unit=ounit) hdate, xfcst, map_source, field, &
                units, desc, xlvl, nx, ny, iproj
          write(unit=ounit) startloc, startlat, startlon, &
                deltalat, deltalon, earth_radius
     
    ! Mercator
    else if (iproj == 1) then
          write(unit=ounit) hdate, xfcst, map_source, field, &
                units, desc, xlvl, nx, ny, iproj
          write(unit=ounit) startloc, startlat, startlon, dx, dy, &
                truelat1, earth_radius
     
    ! Lambert conformal
    else if (iproj == 3) then
          write(unit=ounit) hdate, xfcst, map_source, field, &
                units, desc, xlvl, nx, ny, iproj
          write(unit=ounit) startloc, startlat, startlon, dx, dy, &
                xlonc, truelat1, truelat2, earth_radius

    ! Gaussian
    else if (iproj == 4) then
          write(unit=ounit) hdate, xfcst, map_source, field, &
                units, desc, xlvl, nx, ny, iproj
          write(unit=ounit) startloc, startlat, startlon, &
                nlats, deltalon, earth_radius
     
    ! Polar stereographic
    else if (iproj == 5) then
          write(unit=ounit) hdate, xfcst, map_source, field, &
                units, desc, xlvl, nx, ny, iproj
          write(unit=ounit) startloc, startlat, startlon, dx, dy, &
                xlonc, truelat1, earth_radius
        
    end if
     
    !  3) WRITE WIND ROTATION FLAG
    write(unit=ounit) is_wind_grid_rel
     
    !  4) WRITE 2-D ARRAY OF DATA
    write(unit=ounit) slab

    ! 5) Close the file
    close (ounit)

end subroutine write_2d_slab
