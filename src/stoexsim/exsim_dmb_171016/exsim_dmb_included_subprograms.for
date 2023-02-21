! --------------------------- BEGIN TRIM_C -----------------------
      subroutine trim_c(cstr, nchar)

! strips leading and trailing blanks from cstr, returning the
! result in cstr, which is now nchar characters in length

! Strip off tabs also.

! Here is a sample use in constructing a column header, filled out with 
! periods:

!* Read idtag:
!        idtag = ' '
!        read(nu_in, '(1x,a)') idtag
!        call trim_c(idtag, nc_id)
!* Set up the column headings:
!        colhead = ' '
!        colhead = idtag(1:nc_id)//'......' ! nc_id + 6 > length of colhead

! Dates: 12/23/97 - written by D. Boore
!        12/08/00 - pad with trailing blanks.  Otherwise some comparisons
!                   of the trimmed character string with other strings
!                   can be in error because the original characters are left
!                   behind after shifting.  For example, here is a string
!                   before and after shifting, using the old version:
!                      col:12345
!                           MTWH  before
!                          MTWHH  after (but nc = 4).
!        03/21/01 - Check for a zero length input string
!        11/09/01 - Change check for zero length string to check for all blanks
!        10/19/09 - Strip off tabs
!        04/28/15 - Replaced comment characters * or C with ! (The Fortran 95 standard)

      character cstr*(*)

      if(cstr .eq. ' ') then
        nchar = 0
        return
      end if

      nend = len(cstr)

! Replace tabs with blanks:

      do i = 1, nend
        if(ichar(cstr(i:i)) .eq. 9) then
           cstr(i:i) = ' '
        end if
      end do



!      if(nend .eq. 0) then
!        nchar = 0
!        return
!      end if

      do i = nend, 1, -1
        if (cstr(i:i) .ne. ' ') then
           nchar2 = i
           goto 10
        end if
      end do

10    continue

      do j = 1, nchar2
        if (cstr(j:j) .ne. ' ') then
          nchar1 = j
          goto 20
        end if
      end do

20    continue
   
      nchar = nchar2 - nchar1 + 1
      cstr(1:nchar) = cstr(nchar1: nchar2)
      if (nchar .lt. nend) then
        do i = nchar+1, nend
          cstr(i:i) = ' '
        end do
      end if

      return
      end
! --------------------------- END TRIM_C -----------------------

! --------------------- BEGIN UPSTR ----------------------------------
      Subroutine UPSTR ( text )
! Converts character string in TEXT to uppercase
! Dates: 03/12/96 - Written by Larry Baker
!        04/28/15 - Replaced comment characters * or C with ! (The Fortran 95 standard)

!
      Implicit   None
!
      Character  text*(*)
!
      Integer    j
      Character  ch
!
      Do 1000 j = 1,LEN(text)
         ch = text(j:j)
         If ( LGE(ch,'a') .and. LLE(ch,'z') ) Then
            text(j:j) = CHAR ( ICHAR(ch) - ICHAR('a') + ICHAR('A') )
         End If
 1000    Continue
!
      Return
      End
! --------------------- END UPSTR ----------------------------------
! ---------------------- BEGIN SKIP -------------------
      subroutine SKIP(lunit, nlines)
        if (nlines .lt. 1) then
          return
        else
          do i = 1, nlines
             read(lunit, *)
          end do
          return
        end if
      end
! ---------------------- END SKIP -------------------
! ------------------------------------------------------------------ skipcmnt
      subroutine skipcmnt(nu, comment, ncomments)

! Skip text comments in the file attached to unit nu, but save skipped 
! comments in character array comment.  Skip at least one line, and more as 
! long as the lines are preceded by "|" or "!".

! Dates: 04/16/01 - Written by D. Boore
!        12/07/01 - Added check for eof
!        11/04/03 - Use trim_c to trim possible leading blank
!        02/03/07 - Initialize comments to blank
!        04/28/15 - Replaced comment characters * or C with ! (The Fortran 95 standard)

      character comment(*)*(*), buf*80

      ncomments = 0
100   buf = ' '
      read (nu,'(a)',end=999) buf
      call trim_c(buf,nc_buf)
      if (buf(1:1) .eq.'!' .or. buf(1:1) .eq.'|' .or. 
     :                     ncomments + 1 .eq. 1) then
        ncomments = ncomments + 1
        comment(ncomments) = ' '
        comment(ncomments) = buf(1:nc_buf)
        goto 100
      else 
        backspace nu
      end if

999   continue
 
      return
      end
! ------------------------------------------------------------------ skipcmnt

! --------------------------- BEGIN GET_LUN ----------------
      subroutine get_lun(lun)

! Finds a logical unit number not in use; returns
! -1 if it cannot find one.

! Dates -- 05/19/98 - Written by D. Boore, following
!                     Larry Baker's suggestion
!        04/28/15 - Replaced comment characters * or C with ! (The Fortran 95 standard)

      logical isopen
      do i = 99,10,-1
        inquire (unit=i, opened=isopen)
        if(.not.isopen) then
          lun = i
          return
        end if
      end do
      lun = -1

      return
      end
! --------------------------- END GET_LUN ----------------
     

!----------------- BEGIN AccSqInt -----------------------------
      subroutine accsqint(acc, npts, dt, rmv_trnd, a_sq_int)


! Form integral of acceleration squared, assuming that the acceleration
! is represented by straight lines connecting the digitized values.  This
! routine can be used to compute Arias intensity, defined as

!            Ixx = (pi/2g)*int(acc^2*dt), integrating from 0.0 to the total
!  duration of the record.  The units of Ixx are 
!  velocity [ l^(-1)t^2*(l/t^2)^2*t ] =  l^(-1+2)*t^(2-4+1) = l*t^(-1) = l/t

! Be sure to use consistent units for the acceleration of gravity (g) and acc.
! I am not sure what is conventionally used, but Wilson (USGS OFR 93-556) 
! uses m/sec.

! Dates: 01/13/99 - Written by D.M. Boore
!        04/28/15 - Replaced comment characters * or C with ! (The Fortran 95 standard)

      real a_sq_int(*), acc(*)
      logical rmv_trnd
      double precision cum, a1, a2, ddt_3

      if (rmv_trnd) then      
! remove trend first
        call dcdt(acc, dt, npts, 1, npts, .false., .true.)
      end if

! compute integral of squared acceleration (assume a_sq_int = 0 for first point)

      ddt_3 = dble(dt/3)

      cum = 0.0

      a_sq_int(1) = sngl(cum)
      do j=2,npts
        a1 = acc(j-1)
        a2 = acc(j)
        cum = cum + (a1**2+a1*a2+a2**2)*ddt_3
        a_sq_int(j) = sngl(cum)
      end do

! high pass filter the velocity (may want to allow this in a future version;
! as it is, the acceleration time series can be filtered, so there is no need
! to do it again).

      return
      end
!----------------- END AccSqInt -----------------------------


! ------------------------ begin dcdt -------------------
      subroutine dcdt (y,dt,npts,indx1,indx2,ldc,ldt)
!+
!  dcdt - fits dc or trend between indices indx1 and indx2.
!         then removes dc or detrends whole trace.
!         y is real, dt = delta t.
!         if remove dc, ldc = .true.
!         if detrend, ldt = .true.
!-

! Dates: 12/14/00 - Cleaned up formatting of original program
!        04/28/15 - Replaced comment characters * or C with ! (The Fortran 95 standard)

      real y(*)
      logical ldc,ldt

      if (.not. ldc .and. .not. ldt) then
        return
      end if

!
!...fit dc and trend between indices indx1 and indx2.
      nsum = indx2-indx1+1
      sumx = 0.0
      sumx2 = 0.0
      sumy = 0.0
      sumxy = 0.0
      do i=indx1,indx2
         xsubi = (i-1)*dt
         sumxy = sumxy+xsubi*y(i)
         sumx = sumx+xsubi
         sumx2 = sumx2+xsubi*xsubi
         sumy = sumy+y(i)
      end do
!
!... remove dc.
      if (ldc) then
        avy = sumy/nsum
        do i=1,npts
          y(i) = y(i)-avy
        end do
! Debug
        write(*,'(a)') ' indx1, indx2, avy'
        write(*, *)      indx1, indx2, avy
! Debug



        return
      endif
!
!... detrend. see draper and smith, p. 10.
      if (ldt) then
        bxy = (sumxy-sumx*sumy/nsum)/(sumx2-sumx*sumx/nsum)
        axy = (sumy-bxy*sumx)/nsum
        qxy = dt*bxy
        do i=1,npts
          y(i) = y(i)-(axy+(i-1)*qxy)
        end do
        return
      endif
!
      return
      end
! ------------------------ end dcdt -------------------


!* --------------------- BEGIN LOCATE -----------------
      SUBROUTINE locate(xx,n,x,j)
      
! Comments added by D. Boore on 26feb2010:
!  finds j such that xx(j) < x <= xx(j+1)
!  EXCEPT if x = xx(1), then j = 1 (logically it would be 0 from
!  the above relation, but this is the same returned value of j
!  for a point out of range).
!  Also, if x = xx(n), j = n-1, which is OK
!  Note that j = 0 or j = n indicates that x is out of range.
!
! See the program test_locate.for to test this routine.

! Dates: 04/28/15 - Replaced comment characters * or C with ! (The Fortran 95 standard)

      INTEGER j,n
      REAL x,xx(n)
      INTEGER jl,jm,ju
      jl=0
      ju=n+1
10    if(ju-jl.gt.1)then
        jm=(ju+jl)/2
        if((xx(n).ge.xx(1)).eqv.(x.ge.xx(jm)))then
          jl=jm
        else
          ju=jm
        endif
      goto 10
      endif
      if(x.eq.xx(1))then
        j=1
      else if(x.eq.xx(n))then
        j=n-1
      else
        j=jl
      endif
      return
      END
!* --------------------- END LOCATE -----------------

! --------------------- BEGIN LOCATE_D -----------------
      SUBROUTINE locate_d(xx,n,x,j)
! Dates: 04/28/15 - Replaced comment characters * or C with ! (The Fortran 95 standard)
      INTEGER j,n
      double precision x,xx(n)
      INTEGER jl,jm,ju
      jl=0
      ju=n+1
10    if(ju-jl.gt.1)then
        jm=(ju+jl)/2
        if((xx(n).ge.xx(1)).eqv.(x.ge.xx(jm)))then
          jl=jm
        else
          ju=jm
        endif
      goto 10
      endif
      if(x.eq.xx(1))then
        j=1
      else if(x.eq.xx(n))then
        j=n-1
      else
        j=jl
      endif
      return
      END
! --------------------- END LOCATE_D -----------------
! --------------------- BEGIN ZBRENT -----------------
      FUNCTION zbrent(func,x1,x2,tol)
      
! Dates: 04/28/15 - Replaced "pause" statements with write statements

      INTEGER ITMAX
      REAL zbrent,tol,x1,x2,func,EPS
      EXTERNAL func
      PARAMETER (ITMAX=100,EPS=3.e-8)
      INTEGER iter
      REAL a,b,c,d,e,fa,fb,fc,p,q,r,s,tol1,xm
      a=x1
      b=x2
      fa=func(a)
      fb=func(b)
!      if((fa.gt.0..and.fb.gt.0.).or.(fa.lt.0..and.fb.lt.0.))pause
!     *'root must be bracketed for zbrent'
      if((fa.gt.0..and.fb.gt.0.).or.(fa.lt.0..and.fb.lt.0.)) then
        write(*,*) 
     :    ' root must be bracketed for zbrent; STOPPING EXECUTION'
        STOP
      end if
      c=b
      fc=fb
      do 11 iter=1,ITMAX
        if((fb.gt.0..and.fc.gt.0.).or.(fb.lt.0..and.fc.lt.0.))then
          c=a
          fc=fa
          d=b-a
          e=d
        endif
        if(abs(fc).lt.abs(fb)) then
          a=b
          b=c
          c=a
          fa=fb
          fb=fc
          fc=fa
        endif
        tol1=2.*EPS*abs(b)+0.5*tol
        xm=.5*(c-b)
        if(abs(xm).le.tol1 .or. fb.eq.0.)then
          zbrent=b
          return
        endif
        if(abs(e).ge.tol1 .and. abs(fa).gt.abs(fb)) then
          s=fb/fa
          if(a.eq.c) then
            p=2.*xm*s
            q=1.-s
          else
            q=fa/fc
            r=fb/fc
            p=s*(2.*xm*q*(q-r)-(b-a)*(r-1.))
            q=(q-1.)*(r-1.)*(s-1.)
          endif
          if(p.gt.0.) q=-q
          p=abs(p)
          if(2.*p .lt. min(3.*xm*q-abs(tol1*q),abs(e*q))) then
            e=d
            d=p/q
          else
            d=xm
            e=d
          endif
        else
          d=xm
          e=d
        endif
        a=b
        fa=fb
        if(abs(d) .gt. tol1) then
          b=b+d
        else
          b=b+sign(tol1,xm)
        endif
        fb=func(b)
11    continue
      write(*,*) ' zbrent exceeding maximum iterations'
!      pause 'zbrent exceeding maximum iterations'
      zbrent=b
      return
      END
! --------------------- END ZBRENT -----------------
! --------------------------- BEGIN DIST_3DF ---------------------
      subroutine dist_3df(alat_sta, along_sta, 
     :   alat_ref, along_ref, h_ref, strike_f, dip_f, w1, w2, s1, s2, 
     :   h_min_c, 
     :   d_jb, az_jb, 
     :   d_cd2f, az_cd2f, d_c, az_c,
     :   d_sta_n, d_sta_e, icase_cd2f, icase_c, icase_jb)

! Computes various distance measures from a station to a fault.  The 
! orientation of the fault is that used by Spudich et al., Yucca Mt. project.
! The fault is assumed to be a rectangle whose upper and lower edges are
! horizontal.

! Input:
!      alat_sta, along_sta:  latitude and longitude of station, in degrees,
!                            west longitude is negative
!      alat_ref, along_ref:  as above, for reference point used in defining
!                            fault orientation
!      h_ref:                depth to reference point
!      strike_f, dip_f:      azimuth and dip of fault, in degrees.  strike_f is
!                            measured positive clockwise from north; dip_f is
!                            measured from the horizontal.  When looking in 
!                            the direction strike_f, a positive dip is down to
!                            the right.
!      w1, w2, s1, s2:       distances from the reference point to the edges of
!                            the fault.  s1 and s2 are the distances along
!                            strike to the near and far edges; w1 and w2 are
!                            the distances along the dip direction to the upper
!                            and lower edges of the fault.
!      h_min_c:              minimum depth for computing Campbell's distance
!                            (usually 3.0 km)

! Output:
!      d_jb, d_cd2f, d_c:    Joyner & Boore, closest distance to fault surface,
!                            and Campbell distance, respectively.
!      az_jb, az_cd2f, az_c: as above, for azimuths (NOT YET IMPLEMENTED IN THIS
!                            SUBROUTINE)
!      d_sta_n, d_sta_e:     north and east components of station location
!                            relative to the reference point
!      irgn_cd2f, etc:       region in fault-plane coordinates used to 
!                            compute distances.  I could include a sketch here,
!                            but I will not take the time.  These output 
!                            variables were included mainly to help me check
!                            the subroutine.




! Dates: 12/06/98 - Written by D. Boore
!        09/17/00 - Bring in subroutines via include statement;
!                   renamed faz and az_f to fstrike, strike_f
!        11/12/01 - Changed specification in headers above to indicate that
!                   west longitude is negative (consistent with revision to
!                   subroutine deg2km_f)
!        09/16/10 - Compute az_jb (following Kaklamanos et al definitions).
!                 - Changed "rgn" to "case" (Kaklamanos et al., terminology and definition).
!        04/28/15 - Replaced comment characters * or C with ! (The Fortran 95 standard)

      real dist_sta(3)
      real ix(3), iy(3), iz(3)

      pi = 4.0*atan(1.0)
      dtor = pi/ 180.0

! set up unit vectors in fault coordinates in terms of north, east, down 
! unit vectors.

! Convert angles to radians:
      fstrike =  dtor * strike_f
      fdip = dtor * dip_f

! Initialize arrays:

      do i = 1, 3
         ix(i) = 0.0
         iy(i) = 0.0
         iz(i) = 0.0
      end do

! Compute unit vectors:        ! 1, 2, 3 correspond to n, e, d

      ix(1) = cos(fstrike)
      ix(2) = sin(fstrike)
      ix(3) = 0.0

      iy(1) = -sin(fstrike)*sin(fdip)
      iy(2) =  cos(fstrike)*sin(fdip)
      iy(3) =          -cos(fdip)

      iz(1) = -sin(fstrike)*cos(fdip)
      iz(2) =  cos(fstrike)*cos(fdip)
      iz(3) =           sin(fdip)

! Convert station lat, long into distance north and east:
      call deg2km_f(alat_sta, along_sta, alat_ref, along_ref,
     :              dist_sta(1), dist_sta(2))
      dist_sta(3) = -h_ref    ! note minus sign
      
      d_sta_n = dist_sta(1)
      d_sta_e = dist_sta(2)

! Convert coordinates of reference-to-station vector from n,e,d coordinates
! into fault coordinates:
 
      rx = 0.0
      ry = 0.0
      rz = 0.0
  
      do i = 1, 3
        rx = rx + dist_sta(i) * ix(i)
        ry = ry + dist_sta(i) * iy(i)
        rz = rz + dist_sta(i) * iz(i)
      end do

! Find region and closest distance to fault in the fault plane coordinates:

      call find_h(rx, rz, w1, w2, s1, s2, h_cd2f,
     :            icase_cd2f)                         ! cd2f = Closest Distance
                                                      !        to Fault

! Now do it for Campbell:

! Define w1 for Campbell (I assume that w2 does not need defining; in other
! words, not all of the fault plane is above the Campbell depth)

      d2top_c = h_min_c
      d2top = h_ref + w1 * iz(3)        ! iz(3) = sin(fdip)
      if ( d2top .lt. d2top_c .and. iz(3) .ne. 0.0) then
        w1_c = (d2top_c - h_ref)/ iz(3)
      else
        w1_c = w1
      end if

      call find_h(rx, rz, w1_c, w2, s1, s2, h_c,      ! c = Campbell
     :            icase_c)

! (Work on azimuths later)
      az_cd2f = 999.9
      az_c    = 999.9

! Now do it for Joyner-Boore:

! Need to find rx, ry, rz, w1, w2, s1, s2 in terms of coordinates
! of the fault plane projected onto the surface:

      s1_jb = s1
      s2_jb = s2
      w1_jb = w1 * cos(fdip)
      w2_jb = w2 * cos(fdip)

      rx_jb = rx
      rz_jb = -sin(fstrike) * dist_sta(1) + cos(fstrike) * dist_sta(2)

! Then find the region and distance in the plane to the fault surface
      call find_h(rx_jb, rz_jb, 
     :            w1_jb, w2_jb, s1_jb, s2_jb, h_jb,
     :            icase_jb)

! Now compute the distances:

      d_cd2f = sqrt(h_cd2f**2 + ry   **2)
      d_c    = sqrt(h_c   **2 + ry   **2)
      d_jb   = h_jb

      if (icase_jb == 1 .or. icase_jb == 2 .or. icase_jb == 3 ) then
        az_jb = atan2(rz_jb - w1_jb, rx_jb - s2_jb)/dtor
      else if (icase_jb == 4) then
        az_jb = -90.0
      else if (icase_jb == 5 .or. icase_jb == 6) then
        az_jb = +90.0
      else if (icase_jb == 7 .or. icase_jb == 8 .or. icase_jb == 9 )
     :                                                             then
        az_jb = atan2(rz_jb - w1_jb, rx_jb - s1_jb)/dtor
      else
        write(*,*) 
     :   ' Computing az_jb in dist_3df, icase_jb /= 1-9; QUIT!'
        stop        
      end if
      
      return

      end
! --------------------------- END DIST_3DF ---------------------

!-------------------- BEGIN FIND_H ----------------------------
      subroutine find_h(rx, rz, w1, w2, s1, s2, h, icase)

! 09/16/10 - Redefine regions to be those of Kaklamanos et al
!        04/28/15 - Replaced comment characters * or C with ! (The Fortran 95 standard)

! Here is a conversion table:

!    My notes (old definition of region)   Kaklamanos et al case
!           1                                     7
!           2                                     4
!           3                                     1
!           4                                     2
!           5                                     3
!           6                                     6
!           7                                     9
!           8                                     8
!           9                                     5

! Now it is easy to see where the station lies with respect to the fault;
! there are 9 possibilities:  

      if (   rx .le. s1 .and. rz .le. w1 ) then
! old region 1 (see notes)
!        iregion = 1
        icase = 7
        h = sqrt( (s1-rx)**2 + (w1-rz)**2 )
       
      else if (   rz .le. w1 .and. rx .ge. s1 .and. rx .le. s2 ) then
! old region 2 (see notes)
!        iregion = 2
        icase = 4
        h = w1 - rz

      else if (   rx .ge. s2 .and. rz .le. w1 ) then
! old region 3 (see notes)
!        iregion = 3
        icase = 1
        h = sqrt( (rx-s2)**2 + (w1-rz)**2 )
        
      else if (   rx .ge. s2 .and. rz .ge. w1 .and. rz .le. w2 ) then
! old region 4 (see notes)
!        iregion = 4
        icase = 2
        h = rx - s2

      else if (   rx .ge. s2 .and. rz .ge. w2 ) then
! old region 5 (see notes)
!        iregion = 5
        icase = 3
        h = sqrt( (rx-s2)**2 + (rz-w2)**2 )
        
      else if (   rz .ge. w2 .and. rx .ge. s1 .and. rx .le. s2 ) then
! old region 6 (see notes)
!        iregion = 6
        icase = 6
        h = rz - w2

      else if (   rz .ge. w2 .and. rx .le. s1 ) then
! old region 7 (see notes)
!        iregion = 7
        icase = 9
        h = sqrt( (s1-rx)**2 + (rz-w2)**2 )
        
      else if (   rx .le. s1 .and. rz .ge. w1 .and. rz .le. w2 ) then
! old region 8 (see notes)
!        iregion = 8
        icase = 8
        h = s1 - rx

      else if (      rx .ge. s1 .and. rx .le. s2 
     :         .and. rz .ge. w1 .and. rz .le. w2 ) then
! old region 9 (see notes)
!        iregion = 9
        icase = 5
        h = 0.0

      else
! reaching this is an error
        write(*,'(a)') ' ERROR: Region not found in find_h'

      end if

      return
      end
!-------------------- END FIND_H ----------------------------

!      include '\forprogs\deg2km_f.for'
!      include '\forprogs\locate.for'

!-------------------- BEGIN KM2DEG_F ----------------------------
      subroutine km2deg_f( vn_in, ve_in, alat_ref_in, along_ref_in, 
     :               vlat_out, vlong_out )
        
! convert km north and east from a reference point into lat, long

! assumes positive latitude between 0 and 70 degrees
! assumes east longitude is positive
! assumes angles in degrees

! WARNING: NEEDS DOUBLE PRECISION VERSION OF LOCATE (ATTACHED HERE)

! Dates:  10/01/95 - written by D. Boore
!         05/27/98 - Name changed to km2deg_f
!         06/01/98 - Changed to double precision
!         02/14/09 - Changed input to single precision
!        04/28/15 - Replaced comment characters * or C with ! (The Fortran 95 standard)
               
      double precision alat_tbl(71), b_tbl(71), adcoslat_tbl(71)
      double precision vn, ve, alat_ref, along_ref, vlat, vlong
      Data alat_tbl /
     : 0.000000, 1.000000, 2.000000, 3.000000, 4.000000, 5.000000,
     : 6.000000, 7.000000, 8.000000, 9.000000,10.000000,11.000000,
     :12.000000,13.000000,14.000000,15.000000,16.000000,17.000000,
     :18.000000,19.000000,20.000000,21.000000,22.000000,23.000000,
     :24.000000,25.000000,26.000000,27.000000,28.000000,29.000000,
     :30.000000,31.000000,32.000000,33.000000,34.000000,35.000000,
     :36.000000,37.000000,38.000000,39.000000,40.000000,41.000000,
     :42.000000,43.000000,44.000000,45.000000,46.000000,47.000000,
     :48.000000,49.000000,50.000000,51.000000,52.000000,53.000000,
     :54.000000,55.000000,56.000000,57.000000,58.000000,59.000000,
     :60.000000,61.000000,62.000000,63.000000,64.000000,65.000000,
     :66.000000,67.000000,68.000000,69.000000,70.000000
     :/
      Data b_tbl /
     : 1.842808, 1.842813, 1.842830, 1.842858, 1.842898, 1.842950,
     : 1.843011, 1.843085, 1.843170, 1.843265, 1.843372, 1.843488,
     : 1.843617, 1.843755, 1.843903, 1.844062, 1.844230, 1.844408,
     : 1.844595, 1.844792, 1.844998, 1.845213, 1.845437, 1.845668,
     : 1.845907, 1.846153, 1.846408, 1.846670, 1.846938, 1.847213,
     : 1.847495, 1.847781, 1.848073, 1.848372, 1.848673, 1.848980,
     : 1.849290, 1.849605, 1.849992, 1.850242, 1.850565, 1.850890,
     : 1.851217, 1.851543, 1.851873, 1.852202, 1.852531, 1.852860,
     : 1.853188, 1.853515, 1.853842, 1.854165, 1.854487, 1.854805,
     : 1.855122, 1.855433, 1.855742, 1.856045, 1.856345, 1.856640,
     : 1.856928, 1.857212, 1.857490, 1.857762, 1.858025, 1.858283,
     : 1.858533, 1.858775, 1.859008, 1.859235, 1.859452
     :/
      Data adcoslat_tbl /
     : 1.855365, 1.855369, 1.855374, 1.855383, 1.855396, 1.855414,
     : 1.855434, 1.855458, 1.855487, 1.855520, 1.855555, 1.855595,
     : 1.855638, 1.855683, 1.855733, 1.855786, 1.855842, 1.855902,
     : 1.855966, 1.856031, 1.856100, 1.856173, 1.856248, 1.856325,
     : 1.856404, 1.856488, 1.856573, 1.856661, 1.856750, 1.856843,
     : 1.856937, 1.857033, 1.857132, 1.857231, 1.857331, 1.857435,
     : 1.857538, 1.857643, 1.857750, 1.857858, 1.857964, 1.858074,
     : 1.858184, 1.858294, 1.858403, 1.858512, 1.858623, 1.858734,
     : 1.858842, 1.858951, 1.859061, 1.859170, 1.859276, 1.859384,
     : 1.859488, 1.859592, 1.859695, 1.859798, 1.859896, 1.859995,
     : 1.860094, 1.860187, 1.860279, 1.860369, 1.860459, 1.860544,
     : 1.860627, 1.860709, 1.860787, 1.860861, 1.860934
     :/

      pi = 4.0*atan(1.0)
      d2r = pi/ 180.0

      vn = dble(vn_in)
      ve = dble(ve_in)
      alat_ref =  dble(alat_ref_in)
      along_ref =  dble(along_ref_in) 

! interpolate to find proper arc distance:

      call locate_d( alat_tbl, 71, alat_ref, j)
      b = b_tbl(j) + (alat_ref-alat_tbl(j))*
     :  (b_tbl(j+1)-b_tbl(j))/
     :  (alat_tbl(j+1)-alat_tbl(j))

      adcoslat = adcoslat_tbl(j) + (alat_ref-alat_tbl(j))*
     :  (adcoslat_tbl(j+1)-adcoslat_tbl(j))/
     :  (alat_tbl(j+1)-alat_tbl(j))

      a = adcoslat * cos(d2r*alat_ref)

      dlambda = +ve/a ! version with minus used if assume west long is +
!      dlambda = -ve/a ! minus; positve ve corresponds to decrease in long
      dphi    =  vn/b

! convert from minutes of arc to degrees:
      dlambda = dlambda / 60.0
      dphi    = dphi    / 60.0

      vlat  = alat_ref  + dphi
      vlong = along_ref + dlambda

! Consider using the simpler sphere approximation:
!      vlat = alat_ref + vn/(6371.0 * d2r)
!      vlong = along_ref + ve/(6371.0 * d2r * 
!     :        cos(0.5 * (alat_ref + vlat) * d2r))

      vlat_out = sngl(vlat)
      vlong_out = sngl(vlong)
      
      return
      end
!-------------------- END KM2DEG_F ----------------------------

!-------------------- BEGIN DEG2KM_F ----------------------------
      subroutine deg2km_f( alat_sta, along_sta, alat_ref, along_ref, 
     :                       d_sta_n, d_sta_e   )
        
! convert lat, long into km north and east from a reference point

! assumes latitude between 0 and 70 degrees
! assumes west longitude is negative
! assumes angles in degrees

! Dates:  12/06/98 - written by D. Boore, based on km2deg_f
!         12/18/98 - modified to allow for negative latitudes
!         09/16/00 - Removed double precision
!        04/28/15 - Replaced comment characters * or C with ! (The Fortran 95 standard)
      
      real alat_tbl(71), b_tbl(71), adcoslat_tbl(71)
      real a, b, dphi, dlambda
      Data alat_tbl /
     : 0.000000, 1.000000, 2.000000, 3.000000, 4.000000, 5.000000,
     : 6.000000, 7.000000, 8.000000, 9.000000,10.000000,11.000000,
     :12.000000,13.000000,14.000000,15.000000,16.000000,17.000000,
     :18.000000,19.000000,20.000000,21.000000,22.000000,23.000000,
     :24.000000,25.000000,26.000000,27.000000,28.000000,29.000000,
     :30.000000,31.000000,32.000000,33.000000,34.000000,35.000000,
     :36.000000,37.000000,38.000000,39.000000,40.000000,41.000000,
     :42.000000,43.000000,44.000000,45.000000,46.000000,47.000000,
     :48.000000,49.000000,50.000000,51.000000,52.000000,53.000000,
     :54.000000,55.000000,56.000000,57.000000,58.000000,59.000000,
     :60.000000,61.000000,62.000000,63.000000,64.000000,65.000000,
     :66.000000,67.000000,68.000000,69.000000,70.000000
     :/
      Data b_tbl /
     : 1.842808, 1.842813, 1.842830, 1.842858, 1.842898, 1.842950,
     : 1.843011, 1.843085, 1.843170, 1.843265, 1.843372, 1.843488,
     : 1.843617, 1.843755, 1.843903, 1.844062, 1.844230, 1.844408,
     : 1.844595, 1.844792, 1.844998, 1.845213, 1.845437, 1.845668,
     : 1.845907, 1.846153, 1.846408, 1.846670, 1.846938, 1.847213,
     : 1.847495, 1.847781, 1.848073, 1.848372, 1.848673, 1.848980,
     : 1.849290, 1.849605, 1.849992, 1.850242, 1.850565, 1.850890,
     : 1.851217, 1.851543, 1.851873, 1.852202, 1.852531, 1.852860,
     : 1.853188, 1.853515, 1.853842, 1.854165, 1.854487, 1.854805,
     : 1.855122, 1.855433, 1.855742, 1.856045, 1.856345, 1.856640,
     : 1.856928, 1.857212, 1.857490, 1.857762, 1.858025, 1.858283,
     : 1.858533, 1.858775, 1.859008, 1.859235, 1.859452
     :/
      Data adcoslat_tbl /
     : 1.855365, 1.855369, 1.855374, 1.855383, 1.855396, 1.855414,
     : 1.855434, 1.855458, 1.855487, 1.855520, 1.855555, 1.855595,
     : 1.855638, 1.855683, 1.855733, 1.855786, 1.855842, 1.855902,
     : 1.855966, 1.856031, 1.856100, 1.856173, 1.856248, 1.856325,
     : 1.856404, 1.856488, 1.856573, 1.856661, 1.856750, 1.856843,
     : 1.856937, 1.857033, 1.857132, 1.857231, 1.857331, 1.857435,
     : 1.857538, 1.857643, 1.857750, 1.857858, 1.857964, 1.858074,
     : 1.858184, 1.858294, 1.858403, 1.858512, 1.858623, 1.858734,
     : 1.858842, 1.858951, 1.859061, 1.859170, 1.859276, 1.859384,
     : 1.859488, 1.859592, 1.859695, 1.859798, 1.859896, 1.859995,
     : 1.860094, 1.860187, 1.860279, 1.860369, 1.860459, 1.860544,
     : 1.860627, 1.860709, 1.860787, 1.860861, 1.860934
     :/

      pi = 4.0*atan(1.0)
      d2r = pi/ 180.0

! interpolate to find proper arc distance:

      call locate( alat_tbl, 71, abs(alat_ref), j)

      b = b_tbl(j) + (abs(alat_ref)-alat_tbl(j))*
     :  (b_tbl(j+1)-b_tbl(j))/
     :  (alat_tbl(j+1)-alat_tbl(j))

      adcoslat = adcoslat_tbl(j) + (abs(alat_ref)-alat_tbl(j))*
     :  (adcoslat_tbl(j+1)-adcoslat_tbl(j))/
     :  (alat_tbl(j+1)-alat_tbl(j))

      a = adcoslat * cos(d2r*abs(alat_ref))

! compute lat,long relative to reference:
      dphi = alat_sta - alat_ref
      dlambda = along_sta - along_ref

! convert from degrees to minutes of arc:
      dlambda = dlambda * 60.0
      dphi    = dphi    * 60.0

! compute distances (positive ve corresponds to increase in longitude:
!                    vn positive to the north, ve positive to the east)
      d_sta_e =  a * dlambda 
      d_sta_n =  b * dphi

! Consider replacing the above computation with the following simple
! computation based on assuming that the Earth is a perfect sphere:
!      vn = (alat_sta - alat_ref)*d2r*6371.0
!      ve = (along_sta - along_ref)*d2r*
!     :      cos(0.5*(alat_sta+alat_ref)*d2r)*6371.0

      return
      end
!-------------------- END DEG2KM_F ----------------------------

! --------------------- BEGIN YINTRF ------------------------------------
      function yintrf( x, xin, yin, n)
!
! returns an interpolated value (yintrf) based on straight line
! interpolation of the data in xin and yin.

! Needs Numerical recipe routine locate

!
! dates:  3/14/85 - written
!        11/30/95 - substituted LOCATE instead of starting from beginning
!                   each time
!        03/13/96 - added code to deal with xin increasing or decreasing
!        12/12/00 - Stripped off "locate.for"
!        04/28/15 - Replaced comment characters * or C with ! (The Fortran 95 standard)

      dimension xin(1), yin(1)
      logical incrs

! Is xin increasing or decreasing?
      incrs = .true.
      if (xin(n) .lt. xin(1)) incrs = .false.

! Set value if x is outside the range of xin:
      if (incrs) then
        if ( x .le. xin(1) ) then
            yintrf = yin(1)
            return
        end if
        if ( x .ge. xin(n) ) then
            yintrf = yin(n)
            return
        end if  
      else
        if ( x .ge. xin(1) ) then
            yintrf = yin(1)
            return
        end if
        if ( x .le. xin(n) ) then
            yintrf = yin(n)
            return
        end if  
      end if

! Locate the proper cell and interpolate:
      call locate(xin, n, x, j)
      yintrf = yin(j) + (x-xin(j))*(yin(j+1) - yin(j))/
     * (xin(j+1)-xin(j))

      return
      end
! --------------------- END YINTRF ------------------------------------


! ----------------------------- BEGIN FORK --------------------------
      SUBROUTINE FORK(LX,CX,SIGNI)
! FAST FOURIER                                  2/15/69
!                          LX
!    CX(K) = SQRT(1.0/LX)* SUM (CX(J)*EXP(2*PI*SIGNI*I*(J-1)*(K-1)/LX))
!                          J=1                        FOR K=1,2,...,LX
!
!  THE SCALING BETWEEN FFT AND EQUIVALENT CONTINUUM OUTPUTS
!  IS AS FOLLOWS.
!
!
!     GOING FROM TIME TO FREQUENCY:
!             F(W)=DT*SQRT(LX)*CX(K)
!
!                  WHERE W(K)=2.0*PI*(K-1)*DF

!                  and    DF = 1/(LX*DT)
!
!
!     GOING FROM FREQUENCY TO TIME, WHERE THE FREQUENCY
!     SPECTRUM IS GIVEN BY THE DIGITIZED CONTINUUM SPECTRUM:
!
!             F(T)=DF*SQRT(LX)*CX(K)
!
!                  WHERE T(K)=(K-1)*DT
!
!
!  THE RESULT OF THE SEQUENCE...TIME TO FREQUENCY,POSSIBLE MODIFICATIONS
!  OF THE SPECTRUM (FOR FILTERING,ETC.), BACK TO TIME...
!  REQUIRES NO SCALING.
!
!
!  THIS VERSION HAS A SLIGHT MODIFICATION TO SAVE SOME TIME...
!  IT TAKES THE FACTOR 3.1415926*SIGNI/L OUTSIDE A DO LOOP (D.BOORE 12/8
!  FOLLOWING A SUGGESTION BY HENRY SWANGER).
!

! Some brief notes on usage:

! "signi" is a real variable and should be called either with the value "+1.0"
! of "-1.0".  The particular value used depends on the conventions being used
! in the application (e.g., see Aki and Richards, 1980, Box 5.2, pp. 129--130).

! Time to frequency:
! In calling routine,
! 
!       do i = 1, lx
!         cx(i) = CMPLX(y(i), 0.0)
!       end do
!  where y(i) is the time series and lx is a power of 2
! 
!  After calling Fork with the complex array specified above, the following 
! symmetries exist:
! 
!        cx(1)        = dc value (f = 0 * df, where df = 1.0/(lx*dt))
!        cx(lx/2 + 1) = value at Nyquist (f = (lx/2+1-1)*df = 1.0/(2*dt))
!        cx(lx)       = CONJG(cx(2))
!        cx(lx-1)     = CONJG(cx(3))
!         |           =      |
!        cx(lx-i+2)   = CONJG(cx(i))
!         |           =      |
!        cx(lx/2+2)   = CONJG(cx(lx/2))
! 
! where "CONJG" is the Fortran complex conjugate intrinsic function
! 
! This symmetry MUST be preserved if modifications are made in the frequency 
! domain and another call to Fork (with a different sign for signi) is used
! to go back to the time domain.  If the symmetry is not preserved, then the
! time domain array will have nonzero imaginary components.  There is one case
! where advantage can be taken of this, and that is to find the Hilbert 
! transform and the window of a time series with only two calls to Fork (there 
! is a short note in BSSA {GET REFERENCE} discussing this trick, which amounts 
! to zeroing out the last half of the array and multiplying all but the dc and 
! Nyquist values by 2.0; in the time domain, REAL(cx(i)) and AIMAG(cx(i)) 
! contain the filtered (if a filter was applied) and Hilbert transform of the 
! filtered time series, respectively, while CABS(cx(i)) and ATAN2(AIMAG(cx(i)), 
! REAL(cx(i))) are the window and instantaneous phase of the filtered time 
! series, respectively.

! Some references:

! Farnbach, J.S. (1975). The complex envelope in seismic signal analysis, 
! BSSA 65, 951--962. 
! He states that the factor of 2 is applied for i = 2...npw2/2 (his indices 
! start at 0, I've added 1), which is different than the next reference:

! Mitra, S.K. (2001). Digital Signal Processing, McGraw-Hill, New York.
! He gives an algorithm on p. 794 (eq. 11.81), in which the factor of 2 is 
! applied from 0 frequency to just less than Nyquist.

! 
! The easiest way to ensure the proper symmetry is to zero out the
! last half of the array (as discussed above), but the following is what
! I usually use:  
! modify (filter) only half
! of the cx array:
! 
!       do i = 1, lx/2
!         cx(i) = filter(i)*cx(i)
!       end do
! 
! where "filter(i)" is a possibly complex filter function (and recall that 
! the frequency corresponding to i is f = float(i-1)*df).  After this, fill out
! the last half of the array using
!       
!       do i = lx/2+2, lx
!         cx(i) = CONJG(cx(lx+2-j))
!       end do
! 
! Note that nothing is done with the Nyquist value.  I assume (but am not sure!)
! that this value should be 0.0
! 
! Dates: xx/xx/xx - Written by Norm Brenner(?), Jon Claerbout(?)
!        12/21/00 - Replaced hardwired value of pi with pi evaluated here,
!                     and added comments regarding usage.  Also deleted
!                     dimension specification of cx(lx) and replace it with
!                     cx(*) in the type specification statement.  I also
!                     cleaned up the formatting of the subroutine.
!        08/28/01 - Added comment about variable "signi" being real, and 
!                   added "float" in equations for "sc" and "temp", although 
!                   not strictly required.
!        06/19/02 - Added some comments about computing envelopes and
!                   instantaneous frequencies
!        01/19/15 - Modernize code (get rid of go to statements)
             
      complex cx(*),carg,cexp,cw,ctemp

      pi = 4.0*atan(1.0)

      j=1
      sc=sqrt(1./real(lx))

      do i=1,lx
      
        if (i <= j) then
          ctemp=cx(j)*sc
          cx(j)=cx(i)*sc
          cx(i)=ctemp
        end if
        
        m=lx/2        
        
        DO
          if (j <= m) EXIT
          j=j-m
          m=m/2
          if (m < 1) EXIT
        END DO
        
        j = j + m
        
      end do

      l=1
      DO WHILE (l < lx)
        istep=2*l
        temp= pi * signi/real(l)

        do m=1,l
          carg=(0.,1.)*temp*(m-1)
          cw=cexp(carg)
          do i=m,lx,istep
            ctemp=cw*cx(i+l)
            cx(i+l)=cx(i)-ctemp
            cx(i)=cx(i)+ctemp
          end do
        end do

        l=istep
      END DO
 
! Previous code (as of 18 July 2010):
!      pi = 4.0*atan(1.0)
!
!      j=1
!      sc=sqrt(1./float(lx))
!
!      do i=1,lx
!        if(i.gt.j) go to 2
!        ctemp=cx(j)*sc
!        cx(j)=cx(i)*sc
!        cx(i)=ctemp
!2       m=lx/2
!3       if(j.le.m) go to 5
!        j=j-m
!        m=m/2
!        if(m.ge.1) go to 3
!5       j=j+m
!      end do
!
!      l=1
!6     istep=2*l
!      temp= pi * signi/float(l)
!
!      do m=1,l
!        carg=(0.,1.)*temp*(m-1)
!        cw=cexp(carg)
!        do i=m,lx,istep
!          ctemp=cw*cx(i+l)
!          cx(i+l)=cx(i)-ctemp
!          cx(i)=cx(i)+ctemp
!        end do
!      end do
!
!      l=istep
!      if(l.lt.lx) go to 6

      return
      end
! ----------------------------- END FORK --------------------------


! ------------------------------------------------------------- Get_NPW2
      subroutine get_npw2(npts,signnpw2,npw2)

! Find npw2 (less than npts if signnpw2 < 0)

! Dates: 12/12/00 - Written by D. Boore
!        04/04/10 - Correct error that it does not return
!                   npw2 = npts, if npts is a power of 2.
!        04/28/15 - Replaced comment characters * or C with ! (The Fortran 95 standard)

      npw2_exp = int( alog(float(npts))/alog(2.0) )
      if (signnpw2 < 0.0) then
        npw2 = 2.0**npw2_exp
      else 
        npw2_temp = 2.0**npw2_exp
        if (npw2_temp == npts) then 
          npw2 = npts
        else
          npw2 = 2.0**(npw2_exp+1)
        end if
      end if

      return
      end
! ------------------------------------------------------------- Get_NPW2


! ----------------------------------------------------------------- begin rscalc_interp_acc
      subroutine rscalc_interp_acc(acc, na, omega, damp_in, dt_in,
     :                             rd, rv, aa)
        
!-----------------------------------------------------------------------
! This version does not return response time series.

! Dates: 03/04/10 - Program rsp obtained from B. Chiou.  cmpmax written by
!                   I. Idriss; ucmpmx by R. Youngs.  D. Boore changed the input 
!                   parameters to be equivalent to rdrvaa.for
!        03/05/10 - Substituted rdrvaa for cmpmax
!        03/06/10 - Renamed from rsp_rdrvaa_ucmpmx to rscalc_interp_acc.
!                 - Renamed subroutine ucmpmx to icmpmx ("i" for "interpolate
!                   acceleration) and modified icmpmx.
!        08/13/12 - Norm Abrahamson suggests doing a more exact interpolation 
!                   when period < n*dt (n=10 here).   He suggests interpolating
!                   by adding zeros in the frequency domain, figuring
!                   out at the beginning what will be the shortest period desired,
!                   interpolating the input time series accordingly, and then feeding 
!                   this into rdrvaa.    This requires a major restructuring of this 
!                   subroutine.  I will not do this yet; I just wanted to write
!                   down the suggested revision someplace.
!        10/09/12 - Following up on the last comment, the interpolation is in the driver
!                   smc_interpolate_time_series_using_fork.for, and it uses 
!                   interpolate_time_series_using_fork.for.  I have not incorporated
!                   it yet into programs such as smc2rs, smc2rs2, or blpadflt.
!        04/28/15 - Replaced comment characters * or C with ! (The Fortran 95 standard)
!        05/16/15 - Use Implicit None and replace real*4 and double precision with modern equivalents
!                 - Get rvrdaa.for via an include statement rather than including it in this file
!                   in order to get the most recent version (I had done this before, but for
!                   some reason had copied rvrdaa.for to the bottom of this file, probably
!                   to simplify the distribution of the program).

      implicit none

      real(4) :: acc(*), omega, damp_in, dt_in, 
     :       rd, rv, aa, d0, v0 
      
      integer :: na, kg, kug, npr


      real(8) :: ug(:), pr, damp, dt, z(3)

      real(8) :: w, twopi
      
      allocatable :: ug
      
      integer :: i, nn

      dt = dble(dt_in)
      damp = dble(damp_in)
      
      kg = na
      
      allocate(ug(na))

          do i=1,kg
            ug(i) = dble(acc(i))
          enddo
!...
!... Compute response spectra
!...
           kug=kg-1
           
      w = dble(omega)
      
      twopi = 4.0d0*dasin(1.0d0)
      
      pr = twopi/w

      if(dt == 0.0d0 .or. pr < 10.0d0*dt) then
        call icmpmx(kug, ug, dt, pr, w, damp, z)
        rd = sngl(z(1))
        rv = sngl(z(2))
        aa = sngl(z(3))
      else
        d0 = 0.0
        v0 = 0.0
        call rdrvaa(acc,kg,omega,damp_in,dt_in,rd,rv,aa, d0, v0)
      endif
            
      deallocate(ug)

      return
      end

!-----------------------------------------------------------------------
      subroutine icmpmx(kug, ug, dt_in, pr, w, d, z)
      
! z(1) = SD
! z(2) = RV
! z(3) = AA

! Dates: 03/06/10 - Original program ucmpmx (written by Bob Youngs/I. Idriss)
!                   renamed icmpmx and modified to assume equal time spacing
!                   if the original and interpolated acceleration time
!                   series (the "u" in the original name referred to 
!                   "u"nequal spacing).
!        05/16/15 - Use Implicit None and replace real*4 and double precision with modern equivalents

      implicit none

! Input
      integer :: kug
      real(8) :: ug(*), pr, w, d, dt_in

! Output
      real(8) :: z(*)

! Working variables
      real(8) :: t(3), c(3), x(2,3)
      real(8) :: f1, f2, f3, f4, f5, f6, wd, w2, w3
      real(8) :: dt, e, g1, g2, h1, h2, dug, g, z1, z2, z3, z4
      real(8) :: a, b
      integer :: nn, i, k, ns, is, j
!
      nn=1
      wd=sqrt(1.-d*d)*w
      w2=w*w
      w3=w2*w
      DO i=1,3
        x(1,i)=0.
        z(i)=0.
      END DO
      
      f2=1./w2
      f3=d*w
      f4=1./wd
      f5=f3*f4
      f6=2.*f3
      
      ns= int(10.*dt_in/pr-0.01)+1   !! 05/05/2008
      dt=dt_in/real(ns)
      
      DO k=1,kug
      
        f1=2.*d/w3/dt
        e=dexp(-f3*dt)
        g1=e*dsin(wd*dt)
        g2=e*dcos(wd*dt)
        h1=wd*g2-f3*g1
        h2=wd*g1+f3*g2
        dug=(ug(k+1)-ug(k))/real(ns)
        g=ug(k)
        z1=f2*dug
        z3=f1*dug
        z4=z1/dt
        
        DO is=1,ns
          z2=f2*g
          b=x(1,1)+z2-z3
          a=f4*x(1,2)+f5*b+f4*z4
          x(2,1)=a*g1+b*g2+z3-z2-z1
          x(2,2)=a*h1-b*h2-z4
          x(2,3)=-f6*x(2,2)-w2*x(2,1)
          nn = nn + 1
          DO j=1,3
            c(j)=abs(x(2,j))
            IF (c(j) > z(j)) THEN
              z(j)=c(j)
            END IF
            x(1,j)=x(2,j)
          END DO
          
          g=g+dug
          
        END DO
      END DO
  
      RETURN
      END

!------------------------------
 
!      include '\forprogs\rdrvaa.for'

! ----------------------------------------------------------------- end rscalc_interp_acc

! NOTE: This program was incorporated into rscalc_interp_acc.for.
! Note that rdrvaa includes initial conditions in the argument
! list, while rscalc_interp_acc does not.

!----------------- BEGIN RDRVAA -----------------------------
      subroutine rdrvaa(acc,na,omega,damp,dt,rd,rv,aa, d0, v0)
! This is a modified version of "Quake.For", originally
! written by J.M. Roesset in 1971 and modified by
! Stavros A. Anagnostopoulos, Oct. 1986.  The formulation is that of
! Nigam and Jennings (BSSA, v. 59, 909-922, 1969).  

!   acc = acceleration time series
!    na = length of time series
! omega = 2*pi/per
!  damp = fractional damping (e.g., 0.05)
!    dt = time spacing of input
!    rd = relative displacement of oscillator
!    rv = relative velocity of oscillator
!    aa = absolute acceleration of oscillator
! d0,v0 = initial displacement and velocity (usually set to 0.0)

! Dates: 02/11/00 - Modified by David M. Boore, based on RD_CALC
!        03/11/01 - Double precision version
!        03/14/01 - Added d0, v0 (note on 05 March 2010: I recommend 
!                   that they not be used, by setting them to 0.0.  
!                   I've kept them as arguments in the subroutine call
!                   so that I do not have to modify programs that use
!                   this subroutine).                   
!        03/14/01 - Changed name back to rdrvaa
!        01/31/03 - Moved implicit statement before the type declarations
!        10/10/07 - Initial variable assignments and iteration loop modified 
!                   to double-precision (Chris Stephens)
!        03/05/10 - Delete old (single precision) lines of code
!        12/22/10 - Remove minus sign in front of the initialization of y, ydot. 
!                   The minus sign was a remnant of an earlier version where I did not
!                   understand the meaning of y and ydot.
!        04/28/15 - Replaced comment characters * or C with ! (The Fortran 95 standard)
!        05/16/15 - Use Implicit None and replace real*4 and double precision with modern equivalents

      real(4) :: acc(*), omega, damp, dt, rd, rv, aa, d0, v0
      
      real(8) :: d2, bom, d3, omd, om2, c1, omt, omdt, c2, c3, c4, ss, 
     :           cc, bomt, ee, s1, s2, s3, a11, a12, a21, a22,
     :           s4, s5, b11, b12, b21, b22, y, ydot, y1, z, z1, z2, ra
     
      integer :: i, na
      
      
      d2=1.d0-dble(damp)*dble(damp)
      d2=dsqrt(d2)
      bom=dble(damp)*dble(omega)
      d3 = 2.d0*bom                 ! for aa
      omd=dble(omega)*d2
      om2=dble(omega)*dble(omega)
      c1=1.d0/om2
      
      omt=dble(omega)*dble(dt)
      omdt=omd*dble(dt)
      c2=2.d0*dble(damp)/(om2*omt)
      c3=c1+c2
      c4=1.d0/(dble(omega)*omt)
      ss=dsin(omdt)
      cc=dcos(omdt)
      bomt=dble(damp)*omt
      ee=dexp(-bomt)
      ss=ss*ee
      cc=cc*ee
      s1=ss/omd
      s2=s1*bom
      s3=s2+cc
      a11=s3
      a12=s1
      a21=-om2*s1
      a22=cc-s2
      s4=c4*(1.d0-s3)
      s5=s1*c4+c2
      b11=s3*c3-s5
      b12=-c2*s3+s5-c1
      b21=-s1+s4
      b22=-s4
      
      rd=0.
      rv = 0.                           ! for rv
      aa = 0.                           ! for aa
      
      y=    dble(d0)
      ydot= dble(v0)    
!      y=0.
!      ydot=0.

      do i=1, na-1
      
        y1=a11*y+a12*ydot+b11*dble(acc(i))+b12*dble(acc(i+1))
        ydot=a21*y+a22*ydot+b21*dble(acc(i))+b22*dble(acc(i+1))
        y=y1    ! y is the oscillator output at time corresponding to index i
        z=dabs(y)
        if (z > rd) rd=z
        z1 = dabs(ydot)                   ! for rv
        if (z1 > rv) rv = z1            ! for rv
        ra = -d3*ydot -om2*y1            ! for aa
        z2 = dabs(ra)                     ! for aa
        if (z2 > aa) aa = z2            ! for aa
        
      end do
      
      return
      end
!----------------- END RDRVAA -----------------------------



! --------------- BEGIN TIME_DIFF ---------------------------------
      subroutine time_diff(time_start, time_stop, time_elapsed)

! Dates: 02/18/09 - Written by D.M. Boore
!        04/28/15 - Replaced comment characters * or C with ! (The Fortran 95 standard)
! To be used with
! Standard Fortran 90 intrinsic Subroutine DATE_AND_TIME
!      character datx*8, timx*10
!      call DATE_AND_TIME( datx, timx )
! Date is returned as 'CCYYMMDD'
! Time is returned as 'hhmmss.sss'

      implicit none
      character, intent(in) :: time_start*(*), time_stop*(*)
      real, intent(out) :: time_elapsed
      real ::  secb, sece 
      integer :: ihb, imb, ihe, ime

      read(time_start(1:10),'(i2,i2,f6.3)') 
     :                       ihb, imb, secb
      read(time_stop(1:10),'(i2,i2,f6.3)') 
     :                       ihe, ime, sece
      time_elapsed = 
     :  3600.0*float(ihe-ihb) + 60.0*float(ime-imb) + sece-secb 

      end subroutine time_diff
! --------------- END TIME_DIFF ---------------------------------


!  ------------------- BEGIN BUTTRLCF -------------------------------
      function buttrlcf( f, fcut, norder)
!
! Computes the response of an norder, bidirectional
! high-pass Butterworth filter.  This is the filter
! used by the AGRAM processing (the equation was
! provided by April Converse).

! Modification: 3/27/89 - created by modifying HiPassF
!        04/28/15 - Replaced comment characters * or C with ! (The Fortran 95 standard)

      buttrlcf = 1.0
      if ( fcut.eq.0.0 ) return

      buttrlcf = 0.0

      if ( f .eq. 0.0) return

      buttrlcf = 1.0/ (1.0+(fcut/f)**(2.0*norder))

      return
      end
!  ------------------- END BUTTRLCF -------------------------------
!----------------- BEGIN Acc2V -----------------------------
      subroutine acc2v(acc, npts, dt, rmv_trnd, vel)


! Compute velocity time series from acceleration,
! assuming that the acceleration
! is represented by straight lines connecting the digitized values.

! Dates: 01/06/09 - Written by D.M. Boore, patterned after smc2vd
!        04/28/15 - Replaced comment characters * or C with ! (The Fortran 95 standard)
!        05/15/15 - Remove ddt, which was set but not used

      real acc(*), vel(*) 
      logical rmv_trnd
      double precision cumv, a1, a2,
     : ddt_2

      if (rmv_trnd) then      
! remove trend first (straight line between first and last points)
! Note: acc is replaced with detrended time series
!        call dcdt(acc, dt, npts, 1, npts, .false., .true.)  ! old routine,
!                                                         ! gives steps at ends

         call rmvtrend(acc, npts)
      end if

! compute velocity and displacement, using analytical formulas based
! on representing the acceleration as a series of straightline segments.

      ddt_2   = 0.5d0 * dble(dt)

      cumv = 0.0d0

      vel(1) = sngl(cumv)
      do j=2,npts
        a1 = dble(acc(j-1))
        a2 = dble(acc(j))
        cumv = cumv + (a1 + a2)*ddt_2
        vel(j) = sngl(cumv)
      end do

      return
      end
!----------------- END Acc2V -----------------------------

! ----------------------------- BEGIN RMVTREND ----------------
      subroutine rmvtrend(y, n)

! Removes a straightline fit to first and last points, replacing
! the input array with the detrended array

! Dates: 02/09/99 - written by D. Boore
!        04/28/15 - Replaced comment characters * or C with ! (The Fortran 95 standard)


      real y(*)

      y1 = y(1)
      y2 = y(n)
      slope = (y2 - y1)/float(n-1)

      do i = 1, n
        y(i) = y(i) - (y1 + slope*float(i-1))
      end do

      return
      end
! ----------------------------- END RMVTREND ----------------

! *************** Begin ran1
      FUNCTION ran1(idum)
! Dates: 04/28/15 - Replaced comment characters * or C with ! (The Fortran 95 standard)
      INTEGER idum,IA,IM,IQ,IR,NTAB,NDIV
      REAL ran1,AM,EPS,RNMX
      PARAMETER (IA=16807,IM=2147483647,AM=1./IM,IQ=127773,IR=2836,
     *NTAB=32,NDIV=1+(IM-1)/NTAB,EPS=1.2e-7,RNMX=1.-EPS)
      INTEGER j,k,iv(NTAB),iy
      SAVE iv,iy
      DATA iv /NTAB*0/, iy /0/
      if (idum.le.0.or.iy.eq.0) then
        idum=max(-idum,1)
        do 11 j=NTAB+8,1,-1
          k=idum/IQ
          idum=IA*(idum-k*IQ)-IR*k
          if (idum.lt.0) idum=idum+IM
          if (j.le.NTAB) iv(j)=idum
11      continue
        iy=iv(1)
      endif
      k=idum/IQ
      idum=IA*(idum-k*IQ)-IR*k
      if (idum.lt.0) idum=idum+IM
      j=1+iy/NDIV
      iy=iv(j)
      iv(j)=idum
      ran1=min(AM*iy,RNMX)
      return
      END
! *************** End ran1
      FUNCTION gasdev(idum)
      INTEGER idum
      REAL gasdev
CU    USES ran1
      INTEGER iset
      REAL fac,gset,rsq,v1,v2,ran1
      SAVE iset,gset
      DATA iset/0/
      if (idum.lt.0) iset=0
      if (iset.eq.0) then
1       v1=2.*ran1(idum)-1.
        v2=2.*ran1(idum)-1.
        rsq=v1**2+v2**2
        if(rsq.ge.1..or.rsq.eq.0.)goto 1
        fac=sqrt(-2.*log(rsq)/rsq)
        gset=v1*fac
        gasdev=v2*fac
        iset=1
      else
        gasdev=gset
        iset=0
      endif
      return
      END
!----------------- BEGIN GSPRD_F -----------------------------

! I added entry points so that the program could be called using a single 
! argument, without passing the other arguments through common.  Using
! a single argument is necessary when called by sme Numerical Recipes programs.

! Use:

! Call the setup entry point:
!      dummy =  gsprd_f_setup(r_ref,nsprd_segs,rlow,
!     :                     a_s,b_s,m_s,
!     :                     amag)

! Call as a function:
!            gsprd_n = gsprd_f(rn)

! Deallocate arrays:
!      dummy =  gsprd_f_deallocate()





! Dates: 06/07/95 - Written by D.M. Boore
!        07/02/99 - Added magnitude-dependent "depth" from Atkinson
!                   and Silva, which required adding some parameters to
!                   the passed arguments
!        06/05/00 - Added some explanation of r
!        06/08/00 - Make gsprd nondimensional through the use of r_ref, which 
!                   now appears in the definition of variable const
!                   in const_am0_gsprd
!        01/27/02 - Following Silva, parameters added to allow magnitude
!                   dependent slope (to capture finite fault effects)
!        11/27/05 - Remove deff for Atkinson (2005) source
!        04/24/07 - Put "rmod = r" in the if statement
!        11/13/08 - Redo the function so that it can be used in Numerical
!                   Recipes routines, which assume calls to function(x).
!                   Added entry points rather than using common blocks
!                   to do this (using Larry Baker's help).
!        04/08/11 - Removed using AS00 deff to compute rmod, because the 
!                   application is for small size faults for which the
!                   finite-fault effect approximated by deff is not relevant
!                   (of course, amag will be small, so this would reduce 
!                   the impact of using deff).
!                 - Removed numsource from argument list.
!        04/28/15 - Replaced comment characters * or C with ! (The Fortran 95 standard)

      function gsprd_f(r)
     
      save
      
      real rlow_init(*), a_s_init(*), 
     :                            b_s_init(*), 
     :                            m_s_init(*)
      real, allocatable :: rlow(:), a_s(:), b_s(:), m_s(:)
      real geff(10)
      
! Note that generally r = hypocentral distance.  For Atkinson and Silva 
! (BSSA 90, 255--274) r is the closest distance to the fault plane ("d" in 
! their paper; below their eq. 4), so that rmod is, in effect, accounting
! source depth twice.  See comments in AS00 section of subroutine
! spect_scale

      
!      if (numsource .eq. 9 ) then ! Atkinson and Silva (2000)                                                         
!        deff = 10.0**(-0.05 + 0.15 * amag)
!        rmod = sqrt(r**2 + deff**2)        
!      else      
!        rmod = r      
!      end if

      rmod = r
      
      geff(1) = r_ref/rlow(1)  ! usually set r_ref = 1.0 km.  Be careful
                               ! if a different value or different units are
                               ! used.  In particular, using different units
                               ! will require a change in the scaling factor
                               ! of 1.e-20 used in the definition of const in
                               ! const_am0_gsprd

      do i = 2, nsprd_segs
        slope = a_s(i-1) + b_s(i-1)*(amag - m_s(i-1))
        geff(i) = geff(i-1)*(rlow(i)/rlow(i-1))**slope
      end do
      if (rmod .le. rlow(1)) then
        j = 1
      else if (rmod .ge. rlow(nsprd_segs)) then
        j = nsprd_segs
      else
        call locate(rlow, nsprd_segs, rmod, j)
      end if
      slope = a_s(j) + b_s(j)*(amag - m_s(j))

      gsprd_f = (geff(j)) * (rmod/rlow(j))**slope
      
      return

      entry gsprd_f_setup(r_ref_init,nsprd_segs_init,rlow_init,
     :                  a_s_init,b_s_init,m_s_init,
     :                  amag_init)
     
      allocate(rlow(nsprd_segs_init), 
     :                                a_s(nsprd_segs_init), 
     :                                b_s(nsprd_segs_init),  
     :                                m_s(nsprd_segs_init))
      r_ref                    = r_ref_init
      nsprd_segs               = nsprd_segs_init
      rlow       = rlow_init(1:nsprd_segs_init) 
      a_s        = a_s_init(1:nsprd_segs_init) 
      b_s        = b_s_init(1:nsprd_segs_init) 
      m_s        = m_s_init(1:nsprd_segs_init)
      amag       = amag_init
      
      return
      
      entry gsprd_f_deallocate

      deallocate(rlow, a_s, b_s, m_s)
      
      gsprd_f_deallocate = 1.0
      
      return
     
      
      end
!----------------- END GSPRD_F -----------------------------
 
