cmake_minimum_required(VERSION 2.8)
project ( pspline_lib )
enable_language (Fortran)

include_directories(../ezcdf)

# Compiler flags
if (CMAKE_COMPILER_IS_GNUCC)
  # gfortran
  Message ("Using gfortran flags")
  set (CMAKE_Fortran_FLAGS_RELEASE "-O3")
else ()
  Message ("Did not recognize compiler -- no flags set!")
endif ()

# Set up directories and options
message ("CMAKE_CURRENT_SOURCE_DIR:" ${CMAKE_CURRENT_SOURCE_DIR})

add_library( pspline STATIC
	${CMAKE_CURRENT_SOURCE_DIR}/akherm1.f
	${CMAKE_CURRENT_SOURCE_DIR}/akherm2.f
	${CMAKE_CURRENT_SOURCE_DIR}/akherm3.f
	${CMAKE_CURRENT_SOURCE_DIR}/bcspeval.f
	${CMAKE_CURRENT_SOURCE_DIR}/bcspgrid.f
	${CMAKE_CURRENT_SOURCE_DIR}/bcspline.f
	${CMAKE_CURRENT_SOURCE_DIR}/bcspvec.f
	${CMAKE_CURRENT_SOURCE_DIR}/bpsplinb.f
	${CMAKE_CURRENT_SOURCE_DIR}/bpspline.f
	${CMAKE_CURRENT_SOURCE_DIR}/cspeval.f
	${CMAKE_CURRENT_SOURCE_DIR}/cspline.f
	${CMAKE_CURRENT_SOURCE_DIR}/cubsplb.f
	${CMAKE_CURRENT_SOURCE_DIR}/dnherm1.f
	${CMAKE_CURRENT_SOURCE_DIR}/dnherm2.f
	${CMAKE_CURRENT_SOURCE_DIR}/dnherm3.f
	${CMAKE_CURRENT_SOURCE_DIR}/evbicub.f
	${CMAKE_CURRENT_SOURCE_DIR}/evintrp2d.f
	${CMAKE_CURRENT_SOURCE_DIR}/evintrp3d.f
	${CMAKE_CURRENT_SOURCE_DIR}/evspline.f
	${CMAKE_CURRENT_SOURCE_DIR}/evtricub.f
	${CMAKE_CURRENT_SOURCE_DIR}/ezhybrid_init.f90
	${CMAKE_CURRENT_SOURCE_DIR}/ezhybrid_initx.f90
	${CMAKE_CURRENT_SOURCE_DIR}/ezlinear_init.f90
	${CMAKE_CURRENT_SOURCE_DIR}/ezspline_2netcdf.f90
	${CMAKE_CURRENT_SOURCE_DIR}/ezspline_cdfget3.f90
	${CMAKE_CURRENT_SOURCE_DIR}/ezspline_cdfput3.f90
	${CMAKE_CURRENT_SOURCE_DIR}/ezspline_cinterp.f90
	${CMAKE_CURRENT_SOURCE_DIR}/ezspline_derivative.f90
	${CMAKE_CURRENT_SOURCE_DIR}/ezspline_free.f90
	${CMAKE_CURRENT_SOURCE_DIR}/ezspline_gradient.f90
	${CMAKE_CURRENT_SOURCE_DIR}/ezspline_init.f90
	${CMAKE_CURRENT_SOURCE_DIR}/ezspline_interp.f90
	${CMAKE_CURRENT_SOURCE_DIR}/ezspline_isgridregular.f90
	${CMAKE_CURRENT_SOURCE_DIR}/ezspline_isindomain.f90
	${CMAKE_CURRENT_SOURCE_DIR}/ezspline_load.f90
	${CMAKE_CURRENT_SOURCE_DIR}/ezspline_modulo.f90
	${CMAKE_CURRENT_SOURCE_DIR}/ezspline_obj.f90
	${CMAKE_CURRENT_SOURCE_DIR}/ezspline_save.f90
	${CMAKE_CURRENT_SOURCE_DIR}/ezspline_setup.f90
	${CMAKE_CURRENT_SOURCE_DIR}/ezspline_setupx.f90
	${CMAKE_CURRENT_SOURCE_DIR}/ezspline_type.f90
	${CMAKE_CURRENT_SOURCE_DIR}/ezspline.f90
	${CMAKE_CURRENT_SOURCE_DIR}/f2test.f
	${CMAKE_CURRENT_SOURCE_DIR}/f3test.f
	${CMAKE_CURRENT_SOURCE_DIR}/genxpkg.f
	${CMAKE_CURRENT_SOURCE_DIR}/gridbicub.f
	${CMAKE_CURRENT_SOURCE_DIR}/gridherm1.f
	${CMAKE_CURRENT_SOURCE_DIR}/gridherm2.f
	${CMAKE_CURRENT_SOURCE_DIR}/gridherm3.f
	${CMAKE_CURRENT_SOURCE_DIR}/gridintrp2d.f
	${CMAKE_CURRENT_SOURCE_DIR}/gridintrp3d.f
	${CMAKE_CURRENT_SOURCE_DIR}/gridpc1.f
	${CMAKE_CURRENT_SOURCE_DIR}/gridpc2.f
	${CMAKE_CURRENT_SOURCE_DIR}/gridpc3.f
	${CMAKE_CURRENT_SOURCE_DIR}/gridspline.f
	${CMAKE_CURRENT_SOURCE_DIR}/gridtricub.f
	${CMAKE_CURRENT_SOURCE_DIR}/herm1ev.f
	${CMAKE_CURRENT_SOURCE_DIR}/herm2ev.f
	${CMAKE_CURRENT_SOURCE_DIR}/herm3ev.f
	${CMAKE_CURRENT_SOURCE_DIR}/ibc_ck.f
	${CMAKE_CURRENT_SOURCE_DIR}/mkbicub.f
	${CMAKE_CURRENT_SOURCE_DIR}/mkbicubw.f
	${CMAKE_CURRENT_SOURCE_DIR}/mkherm1.f
	${CMAKE_CURRENT_SOURCE_DIR}/mkherm2.f
	${CMAKE_CURRENT_SOURCE_DIR}/mkherm3.f
	${CMAKE_CURRENT_SOURCE_DIR}/mkintrp2d.f
	${CMAKE_CURRENT_SOURCE_DIR}/mkintrp3d.f
	${CMAKE_CURRENT_SOURCE_DIR}/mkspl2p.f
	${CMAKE_CURRENT_SOURCE_DIR}/mkspl2pb.f
	${CMAKE_CURRENT_SOURCE_DIR}/mkspl2z.f
	${CMAKE_CURRENT_SOURCE_DIR}/mkspl2zb.f
	${CMAKE_CURRENT_SOURCE_DIR}/mkspl3pb.f
	${CMAKE_CURRENT_SOURCE_DIR}/mkspl3zb.f
	${CMAKE_CURRENT_SOURCE_DIR}/mkspline.f
	${CMAKE_CURRENT_SOURCE_DIR}/mktricub.f
	${CMAKE_CURRENT_SOURCE_DIR}/mktricubw.f
	${CMAKE_CURRENT_SOURCE_DIR}/nspline.f
	${CMAKE_CURRENT_SOURCE_DIR}/pc1ev.f
	${CMAKE_CURRENT_SOURCE_DIR}/pc2ev.f
	${CMAKE_CURRENT_SOURCE_DIR}/pc3ev.f
	${CMAKE_CURRENT_SOURCE_DIR}/psp_tolsum.f
	${CMAKE_CURRENT_SOURCE_DIR}/pspline_calls.f
	${CMAKE_CURRENT_SOURCE_DIR}/pspline.f
	${CMAKE_CURRENT_SOURCE_DIR}/r8akherm1.f
	${CMAKE_CURRENT_SOURCE_DIR}/r8akherm2.f
	${CMAKE_CURRENT_SOURCE_DIR}/r8akherm3.f
	${CMAKE_CURRENT_SOURCE_DIR}/r8bcspeval.f
	${CMAKE_CURRENT_SOURCE_DIR}/r8bcspgrid.f
	${CMAKE_CURRENT_SOURCE_DIR}/r8bcspline.f
	${CMAKE_CURRENT_SOURCE_DIR}/r8bcspvec.f
	${CMAKE_CURRENT_SOURCE_DIR}/r8bpsplinb.f
	${CMAKE_CURRENT_SOURCE_DIR}/r8bpspline.f
	${CMAKE_CURRENT_SOURCE_DIR}/r8cspeval.f
	${CMAKE_CURRENT_SOURCE_DIR}/r8cspline.f
	${CMAKE_CURRENT_SOURCE_DIR}/r8cubsplb.f
	${CMAKE_CURRENT_SOURCE_DIR}/r8dnherm1.f
	${CMAKE_CURRENT_SOURCE_DIR}/r8dnherm2.f
	${CMAKE_CURRENT_SOURCE_DIR}/r8dnherm3.f
	${CMAKE_CURRENT_SOURCE_DIR}/r8evbicub.f
	${CMAKE_CURRENT_SOURCE_DIR}/r8evintrp2d.f
	${CMAKE_CURRENT_SOURCE_DIR}/r8evintrp3d.f
	${CMAKE_CURRENT_SOURCE_DIR}/r8evspline.f
	${CMAKE_CURRENT_SOURCE_DIR}/r8evtricub.f
	${CMAKE_CURRENT_SOURCE_DIR}/r8genxpkg.f
	${CMAKE_CURRENT_SOURCE_DIR}/r8gridbicub.f
	${CMAKE_CURRENT_SOURCE_DIR}/r8gridherm1.f
	${CMAKE_CURRENT_SOURCE_DIR}/r8gridherm2.f
	${CMAKE_CURRENT_SOURCE_DIR}/r8gridherm3.f
	${CMAKE_CURRENT_SOURCE_DIR}/r8gridintrp2d.f
	${CMAKE_CURRENT_SOURCE_DIR}/r8gridintrp3d.f
	${CMAKE_CURRENT_SOURCE_DIR}/r8gridpc1.f
	${CMAKE_CURRENT_SOURCE_DIR}/r8gridpc2.f
	${CMAKE_CURRENT_SOURCE_DIR}/r8gridpc3.f
	${CMAKE_CURRENT_SOURCE_DIR}/r8gridspline.f
	${CMAKE_CURRENT_SOURCE_DIR}/r8gridtricub.f
	${CMAKE_CURRENT_SOURCE_DIR}/r8herm1ev.f
	${CMAKE_CURRENT_SOURCE_DIR}/r8herm2ev.f
	${CMAKE_CURRENT_SOURCE_DIR}/r8herm3ev.f
	${CMAKE_CURRENT_SOURCE_DIR}/r8mkbicub.f
	${CMAKE_CURRENT_SOURCE_DIR}/r8mkbicubw.f
	${CMAKE_CURRENT_SOURCE_DIR}/r8mkherm1.f
	${CMAKE_CURRENT_SOURCE_DIR}/r8mkherm2.f
	${CMAKE_CURRENT_SOURCE_DIR}/r8mkherm3.f
	${CMAKE_CURRENT_SOURCE_DIR}/r8mkintrp2d.f
	${CMAKE_CURRENT_SOURCE_DIR}/r8mkintrp3d.f
	${CMAKE_CURRENT_SOURCE_DIR}/r8mkspl2p.f
	${CMAKE_CURRENT_SOURCE_DIR}/r8mkspl2pb.f
	${CMAKE_CURRENT_SOURCE_DIR}/r8mkspl2z.f
	${CMAKE_CURRENT_SOURCE_DIR}/r8mkspl2zb.f
	${CMAKE_CURRENT_SOURCE_DIR}/r8mkspl3pb.f
	${CMAKE_CURRENT_SOURCE_DIR}/r8mkspl3zb.f
	${CMAKE_CURRENT_SOURCE_DIR}/r8mkspline.f
	${CMAKE_CURRENT_SOURCE_DIR}/r8mktricub.f
	${CMAKE_CURRENT_SOURCE_DIR}/r8mktricubw.f
	${CMAKE_CURRENT_SOURCE_DIR}/r8nspline.f
	${CMAKE_CURRENT_SOURCE_DIR}/r8pc1ev.f
	${CMAKE_CURRENT_SOURCE_DIR}/r8pc2ev.f
	${CMAKE_CURRENT_SOURCE_DIR}/r8pc3ev.f
	${CMAKE_CURRENT_SOURCE_DIR}/r8psp_tolsum.f
	${CMAKE_CURRENT_SOURCE_DIR}/r8pspline.f
	${CMAKE_CURRENT_SOURCE_DIR}/r8seval.f
	${CMAKE_CURRENT_SOURCE_DIR}/r8seval2.f
	${CMAKE_CURRENT_SOURCE_DIR}/r8seval3.f
	${CMAKE_CURRENT_SOURCE_DIR}/r8speval.f
	${CMAKE_CURRENT_SOURCE_DIR}/r8spgrid.f
	${CMAKE_CURRENT_SOURCE_DIR}/r8splaan.f
	${CMAKE_CURRENT_SOURCE_DIR}/r8splbrk.f
	${CMAKE_CURRENT_SOURCE_DIR}/r8spleen.f
	${CMAKE_CURRENT_SOURCE_DIR}/r8splinck.f
	${CMAKE_CURRENT_SOURCE_DIR}/r8spline.f
	${CMAKE_CURRENT_SOURCE_DIR}/r8spvec.f
	${CMAKE_CURRENT_SOURCE_DIR}/r8tcspeval.f
	${CMAKE_CURRENT_SOURCE_DIR}/r8tcspgrid.f
	${CMAKE_CURRENT_SOURCE_DIR}/r8tcspline.f
	${CMAKE_CURRENT_SOURCE_DIR}/r8tcspvec.f
	${CMAKE_CURRENT_SOURCE_DIR}/r8tpsplinb.f
	${CMAKE_CURRENT_SOURCE_DIR}/r8tpspline.f
	${CMAKE_CURRENT_SOURCE_DIR}/r8util_bcherm1.f
	${CMAKE_CURRENT_SOURCE_DIR}/r8util_bcherm2.f
	${CMAKE_CURRENT_SOURCE_DIR}/r8util_bcherm3.f
	${CMAKE_CURRENT_SOURCE_DIR}/r8v_spline.f
	${CMAKE_CURRENT_SOURCE_DIR}/r8vecbicub.f
	${CMAKE_CURRENT_SOURCE_DIR}/r8vecherm1.f
	${CMAKE_CURRENT_SOURCE_DIR}/r8vecherm2.f
	${CMAKE_CURRENT_SOURCE_DIR}/r8vecherm3.f
	${CMAKE_CURRENT_SOURCE_DIR}/r8vecintrp2d.f
	${CMAKE_CURRENT_SOURCE_DIR}/r8vecintrp3d.f
	${CMAKE_CURRENT_SOURCE_DIR}/r8vecpc1.f
	${CMAKE_CURRENT_SOURCE_DIR}/r8vecpc2.f
	${CMAKE_CURRENT_SOURCE_DIR}/r8vecpc3.f
	${CMAKE_CURRENT_SOURCE_DIR}/r8vecspline.f
	${CMAKE_CURRENT_SOURCE_DIR}/r8vectricub.f
	${CMAKE_CURRENT_SOURCE_DIR}/r8xlookup.f
	${CMAKE_CURRENT_SOURCE_DIR}/r8zonfind.f
	${CMAKE_CURRENT_SOURCE_DIR}/seval.f
	${CMAKE_CURRENT_SOURCE_DIR}/seval2.f
	${CMAKE_CURRENT_SOURCE_DIR}/seval3.f
	${CMAKE_CURRENT_SOURCE_DIR}/speval.f
	${CMAKE_CURRENT_SOURCE_DIR}/spgrid.f
	${CMAKE_CURRENT_SOURCE_DIR}/splaan.f
	${CMAKE_CURRENT_SOURCE_DIR}/splbrk.f
	${CMAKE_CURRENT_SOURCE_DIR}/spleen.f
	${CMAKE_CURRENT_SOURCE_DIR}/splinck.f
	${CMAKE_CURRENT_SOURCE_DIR}/spline_ez1.f
	${CMAKE_CURRENT_SOURCE_DIR}/spvec.f
	${CMAKE_CURRENT_SOURCE_DIR}/tcspeval.f
	${CMAKE_CURRENT_SOURCE_DIR}/tcspgrid.f
	${CMAKE_CURRENT_SOURCE_DIR}/tcspline.f
	${CMAKE_CURRENT_SOURCE_DIR}/tcspvec.f
	${CMAKE_CURRENT_SOURCE_DIR}/tpsplinb.f
	${CMAKE_CURRENT_SOURCE_DIR}/tpspline.f
	${CMAKE_CURRENT_SOURCE_DIR}/util_bcherm1.f
	${CMAKE_CURRENT_SOURCE_DIR}/util_bcherm2.f
	${CMAKE_CURRENT_SOURCE_DIR}/util_bcherm3.f
	${CMAKE_CURRENT_SOURCE_DIR}/v_spline.f
	${CMAKE_CURRENT_SOURCE_DIR}/vecbicub.f
	${CMAKE_CURRENT_SOURCE_DIR}/vecherm1.f
	${CMAKE_CURRENT_SOURCE_DIR}/vecherm2.f
	${CMAKE_CURRENT_SOURCE_DIR}/vecherm3.f
	${CMAKE_CURRENT_SOURCE_DIR}/vecin2d_argchk.f
	${CMAKE_CURRENT_SOURCE_DIR}/vecin3d_argchk.f
	${CMAKE_CURRENT_SOURCE_DIR}/vecintrp2d.f
	${CMAKE_CURRENT_SOURCE_DIR}/vecintrp3d.f
	${CMAKE_CURRENT_SOURCE_DIR}/vecpc1.f
	${CMAKE_CURRENT_SOURCE_DIR}/vecpc2.f
	${CMAKE_CURRENT_SOURCE_DIR}/vecpc3.f
	${CMAKE_CURRENT_SOURCE_DIR}/vecspline.f
	${CMAKE_CURRENT_SOURCE_DIR}/vectricub.f
	${CMAKE_CURRENT_SOURCE_DIR}/xlookup.f
	${CMAKE_CURRENT_SOURCE_DIR}/zonfind.f
  )
