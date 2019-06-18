Program test_ezfit
  Use Iso_fortran_env , Only: real64
  Implicit None

  ! gfile name and string length (trimmed in routine)
  Integer, Parameter :: fil_len = 128
  Character(len=fil_len) :: gfilename = '../gfiles/g148712.04101'

  ! Psi read from gfile is scaled by this
  ! No default is used, so set_eqd_psi_factor must be called!
  Real(real64), Parameter :: eqd_psi_factor = 1._real64 

  ! Parameters used when initializing splines
  ! eqd_g_tag = 0 means gfile format
  ! reverse_psi = 1 required if psi has a maximum on axis,
  !               which is not the case for this example
  Integer, Parameter :: eqd_g_tag = 0
  Integer, Parameter :: reverse_psi = 0

  ! Evaluation parameters
  Integer, Parameter :: Ntest = 20
  Real(real64) :: R(Ntest), Z(Ntest) ,Bvec_tmp(3), Bvec(Ntest,3)
  Real(real64), Parameter :: Rstart = 1.2_real64, Rend = 2.2_real64
  Real(real64), Parameter :: Ztest = 0.1_real64
  Integer :: ierr, ii

    
  Write(*,'(/a)') '-----------------------------------'
  Write(*,*) 'Reading: ',Trim(gfilename)

  ! Must set this before read!
  Call set_eqd_psi_factor(1._real64)
  
  Call readgfile(gfilename,fil_len)

  Call init_ez_spline(reverse_psi,eqd_g_tag)

  Write(*,'(a/)') '-----------------------------------'

  
  Do ii = 1, Ntest
     R(ii) = (Real(ii,real64)-1._real64)*(Rend-Rstart)  &
          / (Real(Ntest,real64)-1._real64) + Rstart
  End Do
  Z(:) = Ztest
  Bvec = 0._real64
  
  ! Returns (Br, Bz, Bphi)
  Do ii = 1, Ntest
     Call eval_B_vec(R(ii),Z(ii),Bvec_tmp,ierr)
     Bvec(ii,:) = Bvec_tmp
  End Do

  Write(*,*) '   R (m)             Z (m)             Br (T)' &
       // '            Bz (T)            Bphi (T)'
  Do ii = 1, Ntest
     Write(*,'(5f18.12)') R(ii),Z(ii),Bvec(ii,:)
  End Do
  
End Program test_ezfit
