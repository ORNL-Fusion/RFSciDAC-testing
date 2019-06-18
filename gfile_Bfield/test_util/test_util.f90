Program test_util
  Use Iso_fortran_env , Only: real64  
  Use g3d_module, Only : readg_g3d, bfield_geq_bicub, g_type
  Implicit None

  ! gfile name and type
  Character(len=128) :: gfilename = '../gfiles/g148712.04101'
  Type(g_type) :: g
  
  ! Evaluation parameters
  Integer, Parameter :: Ntest = 20
  Real(real64) :: R(Ntest), Z(Ntest),  Bvec(Ntest,3)
  Real(real64), Parameter :: Rstart = 1.2_real64, Rend = 2.2_real64
  Real(real64), Parameter :: Ztest = 0.1_real64
  Integer :: ierr, ii

  Write(*,'(/a)') '-----------------------------------'

  Call readg_g3d(gfilename,g)

  Write(*,'(a/)') '-----------------------------------'

  Do ii = 1, Ntest
     R(ii) = (Real(ii,real64)-1._real64)*(Rend-Rstart)  &
          / (Real(Ntest,real64)-1._real64) + Rstart
  End Do
  Z(:) = Ztest
  Bvec = 0._real64

  ! Returns Bvec(:,[Br, Bz, Bphi])  
  Call bfield_geq_bicub(g,R,Z,Ntest,Bvec,ierr)


  Write(*,*) '   R (m)             Z (m)             Br (T)' &
       // '            Bz (T)            Bphi (T)'
  Do ii = 1, Ntest
     Write(*,'(5f18.12)') R(ii),Z(ii),Bvec(ii,:)
  End Do

  
End program test_util

