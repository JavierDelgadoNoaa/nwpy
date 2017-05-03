!jdelgado mods:
!   - Removed the internal_time/internal_time_loop stuff. I don't know what
!     the purpose of that is, but it skips non-boundary points which I don't
!     care to do, so I don't need it.

Subroutine interp_press2press_lin (press_in, press_out, data_in, &
    data_out, generic, extrapolate, ignore_lowest, tfield, IDS, IDE, JDS, &
    JDE, KDS, KDE, IMS, IME, JMS, JME, KMS, KME, ITS, ITE, JTS, JTE, KTS, &
    KTE )

    ! Interpolates data from one set of pressure surfaces to
    ! another set of pressures
!
      Implicit None
      Integer, Intent(in) :: generic
      Integer, intent(in) :: IDS, IDE, JDS, JDE, KDS, KDE
      Integer, intent(in) :: IMS, IME, JMS, JME, KMS, KME
      Integer, intent(in) :: ITS, ITE, JTS, JTE, KTS, KTE 
      Real, Intent (In) :: press_in (IMS:IME, JMS:JME, generic)
      Real, Intent (In) :: press_out (IMS:IME, JMS:JME, KDS:KDE-1)
      Real, Intent (In) :: data_in (IMS:IME, JMS:JME, generic)
      Real, Intent (Out) :: data_out (IMS:IME, JMS:JME, KMS:KME)
      Logical, Intent (In) :: extrapolate, ignore_lowest
      Logical, Intent (In) :: tfield
      !Integer :: internal_time  !jdelgado
!
!    REAL, INTENT(IN)                   :: press_in(IMS:IME,generic,JMS:JME)
!    REAL, INTENT(IN)                   :: data_in(IMS:IME,generic,JMS:JME)
      Logical :: col_smooth
!
      Integer :: i, j
      Integer :: k, kk
      Real :: desired_press
      Real :: dvaldlnp, dlnp, tadiabat, tiso
!
      Real, Parameter :: ADIAFAC = 9.81 / 1004.
      Real, Parameter :: TSTEXTRAPFAC = .0065
!
    !jdelgado declare implicit variables
      Integer :: LMIN, max_smooth

      !jdelgado dbg
      integer :: ilook, jlook, klook
      ilook = 100
      jlook = 100
      klook = 1

      print *, 'jza33 ims=',ims, 'ime=', ime
      print *, 'jza33 jms=',jms, 'jme=', jme
      print *, 'jza33 kms=',kms, 'kme=', kme
      !jdelgado end

      Do k = KMS, KME
         Do j = JMS, JME
            Do i = IMS, IME
               data_out (i, j, k) = - 99999.9
            End Do
         End Do
      End Do
!
      If (ignore_lowest) Then
         LMIN = 2
      Else
         LMIN = 1
      End If
!
      Do j = JTS, Min (JTE, JDE-1)
         test_i: Do i = ITS, Min (ITE, IDE-1)
        !
             !jdelgado removed all this
             !IF (internal_time_loop .gt. 1) THEN
             !   IF (J .ne. JDS .and. J .ne. JDE-1 .and. &
             !     I .ne. IDS .and. I .ne. IDE-1 ) THEN
                    !! not on boundary
             !     CYCLE test_i
             !   ENDIF
             !ENDIF
        !
        !
            col_smooth = .False.
!
            output_loop: Do k = KDS, KDE - 1

               !print *, 'jza11 k = ', k
               desired_press = press_out (i, j, k)
!
               If (k .Gt. KDS) Then
                  If (TFIELD .And. col_smooth .And. desired_press .Le. &
                 & press_in(i, j, LMIN) .And. press_out(i, j, k-1) .Ge. &
                 & press_in(i, j, LMIN)) Then
                     max_smooth = k
!	  write(message,*) 'I,J, MAX_SMOOTH: ', I,J, MAX_SMOOTH
!         CALL wrf_debug(100,message)
                  End If
               End If
!
! keep track of where the extrapolation begins
!
               If (desired_press .Gt. press_in(i, j, LMIN)) Then
                  If (TFIELD .And. k .Eq. 1 .And. &
                 & (desired_press-press_in(i, j, LMIN)) .Gt. 3000.) &
                 & Then
                     col_smooth = .True. ! due to large extrapolation distance
                  End If
	
!
                  If ((desired_press-press_in(i, j, LMIN)) .Lt. 50.) &
                 & Then ! 0.5 mb
                     data_out (i, j, k) = data_in (i, j, LMIN)
                  Else
                     If (extrapolate) Then
                ! Extrapolate downward because desired P level is below
                ! the lowest level in our input data.  Extrapolate using simple
                ! 1st derivative of value with respect to ln P for the bottom 2
                ! input layers.
!
                ! Add a check to make sure we are not using the gradient of
                ! a very thin layer
!
                        print *, 'jza00 Extrapolating for desired p=',desired_press
                        If (TFIELD) Then
                           tiso = 0.5 * (data_in(i, j, 1)+data_in(i, j, &
                          & 2))
                        End If
!
                        print *, 'jza01 Extrapolating for desired p=',desired_press
!
                        If ((press_in(i, j, LMIN)-press_in(i, j, &
                       & LMIN+1)) .Gt. 500.) Then ! likely isobaric data
                           dlnp = Log (press_in(i, j, LMIN)) - Log &
                          & (press_in(i, j, LMIN+1))
                           dvaldlnp = (data_in(i, j, LMIN)-data_in(i, &
                          & j, LMIN+1)) / dlnp
                        Else ! assume terrain following
                           dlnp = Log (press_in(i, j, LMIN)) - Log &
                          & (press_in(i, j, LMIN+5))
                           dvaldlnp = (data_in(i, j, LMIN)-data_in(i, &
                          & j, LMIN+5)) / dlnp
                        End If

                        print *, 'jza02 Extrapolating for desired p=',desired_press

                        data_out (i, j, k) = data_in (i, j, LMIN) + &
                       & dvaldlnp * (Log(desired_press)-Log(press_in(i, &
                       & j, LMIN)))
!
                        print *, 'jza02 Extrapolating for desired p=',desired_press
                        
                        If (TFIELD .And. data_out(i, j, k) .Lt. &
                       & tiso-0.2) Then
!
! restrict slope to -1K/10 hPa
                           dvaldlnp = Max &
                          & (dvaldlnp,-1.0/Log(press_in(i, j, &
                          & LMIN)/(press_in(i, j, LMIN)-1000.)))
!
                           data_out (i, j, k) = data_in (i, j, LMIN) + &
                          & dvaldlnp * &
                          & (Log(desired_press)-Log(press_in(i, j, &
                          & LMIN)))
!
                        Else If (TFIELD .And. data_out(i, j, k) .Gt. &
                       & tiso+0.2) Then
!
! restrict slope to +0.8K/10 hPa
                           dvaldlnp = Min (dvaldlnp, &
                          & 0.8/Log(press_in(i, j, LMIN)/(press_in(i, &
                          & j, LMIN)-1000.)))
!
                           data_out (i, j, k) = data_in (i, j, LMIN) + &
                          & dvaldlnp * &
                          & (Log(desired_press)-Log(press_in(i, j, &
                          & LMIN)))
!
                        End If
                        
                        print *, 'jza02 Extrapolating for desired p=',desired_press
!
                     Else
                        data_out (i, j, k) = data_in (i, j, LMIN)
                     End If
                  End If

                   print *, 'jza04 Extrapolating for desired p=',desired_press
                   print *, 'jza041 pres_in(i,j,generic)=', press_in(i, j, generic)
               Else If (desired_press .Lt. press_in(i, j, generic)) then
                   print *, 'jza043 Extrapolating for desired p=',desired_press
                   If ((press_in(i, j, generic)-desired_press) .Lt. 10.) then
                     data_out (i, j, k) = data_in (i, j, generic)
                  Else
                     If (extrapolate) Then
                        print *, 'jza05 Extrapolating for desired p=',desired_press

                ! Extrapolate upward
                        If ((press_in(i, j, generic-1)-press_in(i, j, &
                                       generic)) .Gt. 50.) Then
                           dlnp = Log (press_in(i, j, generic)) - Log &
                                       (press_in(i, j, generic-1))
                           dvaldlnp = (data_in(i, j, &
                                       generic)-data_in(i, j, generic-1)) / dlnp
                        Else
                           dlnp = Log (press_in(i, j, generic)) - Log &
                          & (press_in(i, j, generic-2))
                           dvaldlnp = (data_in(i, j, &
                          & generic)-data_in(i, j, generic-2)) / dlnp
                        End If
                        data_out (i, j, k) = data_in (i, j, generic) + &
                       & dvaldlnp * (Log(desired_press)-Log(press_in(i, &
                       & j, generic)))
                     Else
                        data_out (i, j, k) = data_in (i, j, generic)
                     End If
                  End If
               Else
                  !We can trap between two levels and linearly interpolate
!
                  !print *, 'jza323 Interpolating for desired p=',desired_press

                  input_loop: Do kk = LMIN, generic - 1
                     If (desired_press .Eq. press_in(i, j, kk)) Then
                        data_out (i, j, k) = data_in (i, j, kk)
                        Exit input_loop
                     Else If ((desired_press .Lt. press_in(i, j, kk)) &
                    & .And. (desired_press .Gt. press_in(i, j, kk+1))) &
                    & Then
!
!       do trapped in lnp
!
                        dlnp = Log (press_in(i, j, kk)) - Log &
                               (press_in(i, j, kk+1))
                        dvaldlnp = (data_in(i, j, kk)-data_in(i, j, &
                           kk+1)) / dlnp
                        data_out (i, j, k) = data_in (i, j, kk+1) + &
                             dvaldlnp * (Log(desired_press)-Log(press_in(i, &
                       & j, kk+1)))
!
                        Exit input_loop
                     End If
                
                !jdelgado start dbg
                if( i .eq. ilook) then
                    if( j .eq. jlook) then
                        if( k .eq. klook ) then
                            print *, 'jza22 desired_pres = ', desired_press
                            print *, 'jza22 pres_in = ', press_in(i,j,k)
                            print *, 'jza22 data_out = ', data_out(i,j,k)
                        endif
                    endif
                endif
!               !jdelgado end

                  End Do input_loop ! Iteration thru KK (for interpolation)
               End If
            End Do output_loop ! iteration thru k (for extrap/interp)
!
            !print *, 'jza99 Doing column smooth - looks for data_out(i,j,k+1)!!'
            If (col_smooth) Then
               Do k = Max (KDS, max_smooth-4), max_smooth + 4
                  data_out (i, j, k) = 0.5 * (data_out(i, j, &
                 & k)+data_out(i, j, k+1))
               End Do
            End If
!
         End Do test_i
      End Do
End Subroutine interp_press2press_lin
!
Subroutine wind_adjust (press_in, press_out, U_in, V_in, U_out, V_out, &
& generic, depth_replace, IDS, IDE, JDS, JDE, KDS, KDE, IMS, IME, JMS, &
& JME, KMS, KME, ITS, ITE, JTS, JTE, KTS, KTE)
!
    !implicit none
      Integer :: IDS, IDE, JDS, JDE, KDS, KDE
      Integer :: IMS, IME, JMS, JME, KMS, KME
      Integer :: ITS, ITE, JTS, JTE, KTS, KTE, generic
      Integer :: MAXLIN, MAXLOUT
!
      Real, Intent (In) :: press_in (IMS:IME, JMS:JME, generic)
      Real, Intent (In) :: press_out (IMS:IME, JMS:JME, KDS:KDE-1)
      Real, Intent (In) :: U_in (IMS:IME, JMS:JME, generic)
      Real, Intent (In) :: V_in (IMS:IME, JMS:JME, generic)
      Real, Intent (Inout) :: U_out (IMS:IME, KMS:KME, JMS:JME)
      Real, Intent (Inout) :: V_out (IMS:IME, KMS:KME, JMS:JME)
      Real :: p1d_in (generic)
      Real :: p1d_out (KDS:KDE-1)
!
    !jdelgado declare variables that were implicit
      Integer :: LMIN
      Integer :: zmax
!
      Do j = JTS, Min (JTE, JDE-1)
         Do i = ITS, Min (ITE, IDE-1)
!
!        IF (press_out(I,J,1) .lt. press_in(I,J,2)) then
            If ((press_in(i, j, 2)-press_out(i, j, 1)) .Gt. 200.) Then
!
               U_out (i, 1, j) = U_in (i, j, 2)
               V_out (i, 1, j) = V_in (i, j, 2)
!
               INLOOP: Do L = 2, generic
                  p1d_in (L) = - 9999.
                  If ((press_in(i, j, 2)-press_in(i, j, L)) .Lt. &
                 & depth_replace) Then
                     p1d_in (L) = (press_in(i, j, 2)-press_in(i, j, L))
                     MAXLIN = L
                  Else
                     p1d_in (L) = (press_in(i, j, 2)-press_in(i, j, L))
                     Exit INLOOP
                  End If
               End Do INLOOP
!
               OUTLOOP: Do L = KDS, KDE - 1
                  p1d_out (L) = - 9999.
                  If ((press_out(i, j, 1)-press_out(i, j, L)) .Lt. &
                 & depth_replace) Then
                     p1d_out (L) = (press_out(i, j, 1)-press_out(i, j, &
                    & L))
                     MAXLOUT = L
                  Else
                     Exit OUTLOOP
                  End If
               End Do OUTLOOP
!
               Do L = 1, MAXLOUT
                  ptarg = p1d_out (L)
!
                  FINDLOOP: Do LL = 2, MAXLIN
!
                     If (p1d_in(LL) .Lt. ptarg .And. p1d_in(LL+1) .Gt. &
                    & ptarg) Then
!
                        dlnp = Log (p1d_in(LL)) - Log (p1d_in(LL+1))
                        dudlnp = (U_in(i, j, LL)-U_in(i, j, LL+1)) / &
                       & dlnp
                        dvdlnp = (V_in(i, j, LL)-V_in(i, j, LL+1)) / &
                       & dlnp
                        U_out (i, L, j) = U_in (i, j, LL) + dudlnp * &
                       & (Log(ptarg)-Log(p1d_in(LL)))
                        V_out (i, L, j) = V_in (i, j, LL) + dvdlnp * &
                       & (Log(ptarg)-Log(p1d_in(LL)))
!
                        Exit FINDLOOP
                     End If
!
                  End Do FINDLOOP
               End Do ! MAXLOUT loop
!
!
            End If
!
         End Do
      End Do
!
!
End Subroutine wind_adjust
!--------------------------------------------------------------------
!
Subroutine interp_press2press_log (press_in, press_out, data_in, &
& data_out, generic, extrapolate, ignore_lowest, IDS, IDE, JDS, JDE, &
& KDS, KDE, IMS, IME, JMS, JME, KMS, KME, ITS, ITE, JTS, JTE, KTS, KTE)
!& internal_time)  !jdelgado
!
    ! Interpolates ln(data) from one set of pressure surfaces to
    ! another set of pressures
!
      Implicit None
      Integer :: IDS, IDE, JDS, JDE, KDS, KDE
      Integer :: IMS, IME, JMS, JME, KMS, KME
      Integer :: ITS, ITE, JTS, JTE, KTS, KTE, generic
      !Integer :: internal_time !jdelgado
!
!    REAL, INTENT(IN)                   :: press_in(IMS:IME,generic,JMS:JME)
      Real, Intent (In) :: press_in (IMS:IME, JMS:JME, generic)
      Real, Intent (In) :: press_out (IMS:IME, JMS:JME, KDS:KDE-1)
!    REAL, INTENT(IN)                   :: data_in(IMS:IME,generic,JMS:JME)
!    REAL, INTENT(IN)                   :: data_in(IMS:IME,JMS:JME,generic)
      Real :: data_in (IMS:IME, JMS:JME, generic)
      Real, Intent (Out) :: data_out (IMS:IME, JMS:JME, KMS:KME)
      Logical, Intent (In) :: extrapolate, ignore_lowest
!
      Integer :: i, j
      Integer :: k, kk
      Real :: desired_press
      Real :: dlnvaldlnp, dlnp
!
    !jdelgado declare implicit variables
      Integer :: LMIN, max_smooth
!
      Do k = 1, generic
         Do j = JTS, Min (JTE, JDE-1)
            Do i = ITS, Min (ITE, IDE-1)
               data_in (i, j, k) = Max (data_in(i, j, k), 1.e-13)
            End Do
         End Do
      End Do
!
      Do k = KMS, KME
         Do j = JMS, JME
            Do i = IMS, IME
               data_out (i, j, k) = - 99999.9
            End Do
         End Do
      End Do
!
      If (ignore_lowest) Then
         LMIN = 2
      Else
         LMIN = 1
      End If
!
      Do j = JTS, Min (JTE, JDE-1)
         test_i: Do i = ITS, Min (ITE, IDE-1)
            
            !jdelgado start- commented this out
            !If (internal_time .Gt. 1) Then
            !   If (j .Ne. JDS .And. j .Ne. JDE-1 .And. i .Ne. IDS .And. &
            !  & i .Ne. IDE-1) Then
            !! not on boundary
            !!      Cycle test_i
            !!   End If
            !!End If
            !jdelgado end
!
!
            output_loop: Do k = KDS, KDE - 1
!
               desired_press = press_out (i, j, k)
!
               If (desired_press .Gt. press_in(i, j, LMIN)) Then
!
                  If ((desired_press-press_in(i, j, LMIN)) .Lt. 10.) &
                 & Then ! 0.1 mb
                     data_out (i, j, k) = data_in (i, j, LMIN)
                  Else
                     If (extrapolate) Then
                ! Extrapolate downward because desired P level is below
                ! the lowest level in our input data.  Extrapolate using simple
                ! 1st derivative of value with respect to ln P for the bottom 2
                ! input layers.
!
                ! Add a check to make sure we are not using the gradient of
                ! a very thin layer
!
                        If ((press_in(i, j, LMIN)-press_in(i, j, &
                       & LMIN+1)) .Gt. 100.) Then
                           dlnp = Log (press_in(i, j, LMIN)) - Log &
                          & (press_in(i, j, LMIN+1))
                           dlnvaldlnp = (Log(data_in(i, j, &
                          & LMIN))-Log(data_in(i, j, LMIN+1))) / dlnp
!
                        Else
!
                           dlnp = Log (press_in(i, j, LMIN)) - Log &
                          & (press_in(i, j, LMIN+2))
                           dlnvaldlnp = (Log(data_in(i, j, &
                          & LMIN))-Log(data_in(i, j, LMIN+2))) / dlnp
!
                        End If
!
                        data_out (i, j, k) = Exp (Log(data_in(i, j, &
                       & LMIN))+dlnvaldlnp*(Log(desired_press)-&
                       & Log(press_in(i, j, LMIN))))
                     Else
                        data_out (i, j, k) = data_in (i, j, LMIN)
                     End If
                  End If
               Else If (desired_press .Lt. press_in(i, j, generic)) &
              & Then
                  If ((press_in(i, j, generic)-desired_press) .Lt. 10.) &
                 & Then
                     data_out (i, j, k) = data_in (i, j, generic)
                  Else
                     If (extrapolate) Then
                ! Extrapolate upward
                        If ((press_in(i, j, generic-1)-press_in(i, j, &
                       & generic)) .Gt. 50.) Then
                           dlnp = Log (press_in(i, j, generic)) - Log &
                          & (press_in(i, j, generic-1))
                           dlnvaldlnp = (Log(data_in(i, j, &
                          & generic))-Log(data_in(i, j, generic-1))) / &
                          & dlnp
                        Else
                           dlnp = Log (press_in(i, j, generic)) - Log &
                          & (press_in(i, j, generic-2))
                           dlnvaldlnp = (Log(data_in(i, j, &
                          & generic))-Log(data_in(i, j, generic-2))) / &
                          & dlnp
                        End If
                        data_out (i, j, k) = Exp (Log(data_in(i, j, &
                       & generic))+dlnvaldlnp*(Log(desired_press)-&
                       & Log(press_in(i, j, generic))))
                     Else
                        data_out (i, j, k) = data_in (i, j, generic)
                     End If
                  End If
               Else
            ! We can trap between two levels and linearly interpolate
!
                  input_loop: Do kk = LMIN, generic - 1
                     If (desired_press .Eq. press_in(i, j, kk)) Then
                        data_out (i, j, k) = data_in (i, j, kk)
                        Exit input_loop
                     Else If ((desired_press .Lt. press_in(i, j, kk)) &
                    & .And. (desired_press .Gt. press_in(i, j, kk+1))) &
                    & Then
!
!       do trapped in lnp
!
                        dlnp = Log (press_in(i, j, kk)) - Log &
                       & (press_in(i, j, kk+1))
                        dlnvaldlnp = (Log(data_in(i, j, &
                       & kk))-Log(data_in(i, j, kk+1))) / dlnp
                        data_out (i, j, k) = Exp (Log(data_in(i, j, kk+&
                       & 1))+dlnvaldlnp*(Log(desired_press)-&
                       & Log(press_in(i, j, kk+1))))
!
                        Exit input_loop
!
                     End If
!
                  End Do input_loop
               End If
            End Do output_loop
         End Do test_i
      End Do
End Subroutine interp_press2press_log
