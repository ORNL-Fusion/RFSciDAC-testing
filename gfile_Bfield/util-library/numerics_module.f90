!-----------------------------------------------------------------------------
!+ Math/geometric routines
! 
!  Subroutine quicksort
!  Subroutine int_curve_curve
!  Subroutine int_line_curve
!  Subroutine int_two_lines
!  Subroutine move_L_on_C
!  Subroutine linear_interp
!  Subroutine locate_bisect
!  Function rlinspace
!  Function inside_poly
!
!-----------------------------------------------------------------------------
Module math_geo_module
Implicit None
Private
Public :: linear_interp
Public :: rlinspace
Contains

!------------------------------------------------------------------------------
!+ Linear interpolation
!------------------------------------------------------------------------------
Subroutine linear_interp(xarr,yarr,narr,xin,yout,ierr)
Use kind_mod, Only: real64, int32
Implicit None
Integer(int32), Intent(In) :: narr
Real(real64), Intent(In) :: xarr(narr), yarr(narr), xin
Integer(int32), Intent(Out) :: ierr
Real(real64), Intent(Out) :: yout

Integer(int32) :: j

If (xin .lt. Minval(xarr)) Then
  Write(*,*) 'Error from linear_interp, x out of bounds lower (x,Minval(xarr))=',xin,Minval(xarr)
  ierr = 1
  Return
!  Stop "quitting"
Endif
If (xin .gt. Maxval(xarr)) Then
  Write(*,*) 'Error from linear_interp, x out of bounds upper (x,Maxval(xarr))=',xin,Maxval(xarr)
  ierr = 1
  Return
!  Stop "quitting"
Endif
Call locate_bisect(xarr,narr,xin,j,ierr)

If (j .eq. 0) Then
  If (abs(xin - xarr(1)) < 1.e-10) Then
    yout = yarr(1)
    ierr = 0
    Return
  Endif
Endif
If (j .eq. narr) Then
  If (abs(xin - xarr(narr)) < 1.e-10) Then
    yout = yarr(narr)
    ierr = 0
    Return
  Endif
Endif

If (ierr .eq. 1) Then
  Write(*,*) "xarr([1,end])",xarr(1),xarr(narr)
  Write(*,*) "xin",xin
  Write(*,*) 'xin - xarr(1)',xin-xarr(1)
  Stop "error linear interp!, ierr set from locate_bisect"  
  yout = 0.d0
  Return
Endif

yout = yarr(j) + ((xin-xarr(j))/(xarr(j+1)-xarr(j)))*(yarr(j+1)-yarr(j))

End Subroutine linear_interp



!------------------------------------------------------------------------------
!+ Search an ordered list of xarr for x, returns j where xarr(j) <= x <= xarr(j+1)
!------------------------------------------------------------------------------
Subroutine locate_bisect(xarr,narr,x,j,ierr)
!
! If x is off the table (including exactly equal to end values!) then
!  j is set to 0 (off left end) or narr (off right end)
!
Use kind_mod, Only: real64, int32
Implicit None
Integer(int32), Intent(In) :: narr
Real(real64), Intent(In) :: xarr(narr), x
Integer(int32), Intent(Out) :: j, ierr
Logical :: test1
Integer(int32) :: jlo, jup, jmid
ierr = 0
jlo = 0
jup = narr + 1
test1 = ( xarr(narr) .gt. xarr(1) )

!If ( ( x .gt. Maxval(xarr) ) .OR. ( x .lt. Minval(xarr) ) ) Then
!  ierr = 1
!  Stop "x out of range in locate_bisect"
!Endif

Do While (jup - jlo .gt. 1) 
  jmid = (jup + jlo)/2  ! compute midpoint
  If ( test1 .eqv. (x .gt. xarr(jmid)) ) Then
    jlo = jmid
  Else
    jup = jmid
  Endif
Enddo
j = jlo

If (( j .eq. 0 ) .OR. (j .eq. narr)) ierr = 1

End Subroutine locate_bisect

!------------------------------------------------------------------------------
!+ Returns a linearly spaced real vector given endpoints and number of elements
!------------------------------------------------------------------------------
Function rlinspace(xstart,xend,numel)  & 
Result(rlinvec)
!
! Description: 
!   This function returns a real vector of length(numel) with linearly spaced
!   values from xstart to xend.  Similar to the Matlab function.
!
! Inputs:
!  xstart,xend: Values of the first and last points of the array [real]
!  numel: Number of elements in the array [integer]
! Outputs:
!  rlinvec: The linearly spaced array
!         
! History:
!  Version   Date      Comment
!  -------   ----      -------
!  1.0     07/22/2009  Original Code.  JL
!  1.1     09/01/2010  Updated for PENTA3. JL
!  1.2     07/18/2011  Ported to w7 routines. JL
!
! Author(s): J. Lore 07/2009 - 7/18/2011
!
! Modules used:
Use kind_mod, Only: real64, int32

Implicit None

! Input/output                      !See above for descriptions
Real(real64),    Intent(in) :: xstart  
Real(real64),    Intent(in) :: xend
Integer(int32), Intent(in) :: numel
Real(real64)                :: rlinvec(numel)

! Local scalars
Integer(int32)   ::  ii
!- End of header -------------------------------------------------------------

If (numel .eq. 1) Then
  rlinvec(1) = xstart
  Return
Endif

Do ii = 1,numel
  rlinvec(ii) = ( Real(ii,real64) - 1._real64 ) * ( xend - xstart ) &
       / ( numel - 1._real64 ) + xstart
Enddo

EndFunction rlinspace


End Module math_geo_module
