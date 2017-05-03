! FILE: hpost_interp.F90



SUBROUTINE linear_interp_to_z(data_in,km,data_out,im,jm,kmz,z_in,z_out,undef,increases_up)
    !!
    ! Linearly interpolate a 3-d field from its input Z-levels to a given set of
    ! output Z-levels. 
    ! This was taken from HPost and modified so that (a) the input and output
    ! data are arranged such that the Z dimension is the outtermost index and
    ! (b) to account for round off errors in the given levels by having a 
    ! threshold instead of searching for equality.
    !
    ! This code was taken from HPost
    !
    implicit none
    Real(4), intent(in) :: data_in(km,im,jm) ! Input array to interpolate
    Real(4), intent(out) :: data_out(kmz,im,jm) ! Output array
    Integer, intent(in) :: km !number of input levels  :: im,jm,km,kmz,i,j,k,k0,m
    integer, intent(in) :: im ! grid points in i direction
    integer, intent(in) :: jm ! grid points in j direction
    integer, intent(in) :: kmz ! number of output levels
    real(4), intent(in) :: z_in(km,im,jm) ! pressure values corresponding to data_in
    real(4), intent(in) :: z_out(kmz) ! list of pressure levels to interpolate to
    real(4), intent(in) :: undef ! value used for missing values
    logical, intent(in) :: increases_up ! True if value increases upward (e.g. temperature), False otherwise (e.g. pressure)
    ! Local variables
    Real     :: aaa,bbb,xxx
    Integer :: i,j,k,k0,m
    Integer :: THRESHOLD

    THRESHOLD = 0.0001
    data_out = undef
    if( increases_up ) then
        DO j = 1, jm
            DO i = 1, im ; k0 = 1
                DO m = 1, kmz 
                    !xxx = z_out(m) - THRESHOLD
                    xxx = z_out(m)
                    if ( z_in(1,i,j) <= xxx .AND. xxx < z_in(km,i,j) ) then
                        DO k = k0, km-1 
                            IF ( xxx < z_in(k+1,i,j) ) EXIT 
                        ENDDO
                        aaa = (z_in(k+1,i,j)-z_out(m)) / (z_in(k+1,i,j)-z_in(k,i,j)) 
                        bbb = 1 - aaa
                        if(data_in(k,i,j) .ne. undef .and. &
                           data_in(k+1,i,j) .ne. undef) then
                            data_out(m,i,j) = aaa*data_in(k,i,j) +     &
                                            bbb*data_in(k+1,i,j)
                            k0 = k
                        endif
                    endif
                ENDDO
            ENDDO
        ENDDO
    else ! increases down (e.g. pressure)
        DO j = 1, jm
            DO i = 1, im ; k0 = 1
                DO m = 1, kmz 
                    !xxx = z_out(m) - THRESHOLD
                    xxx = z_out(m)
                    if ( z_in(1,i,j) >= xxx .AND. xxx > z_in(km,i,j) ) then
                        DO k = k0, km-1 
                            IF ( xxx > z_in(k+1,i,j) ) EXIT 
                        ENDDO
                        aaa = (z_out(m)-z_in(k+1,i,j))/(z_in(i,j,k)-z_in(k+1,i,j)) 
                        bbb = 1 - aaa
                        if(data_in(k,i,j) .ne. undef .and. &
                           data_in(k+1,i,j) .ne. undef) then
                            data_out(m,i,j) = aaa*data_in(k,i,j) + bbb*data_in(k+1,i,j)
                            k0 = k
                        endif
                    endif
                ENDDO
            ENDDO
        ENDDO
    
    endif

END SUBROUTINE


