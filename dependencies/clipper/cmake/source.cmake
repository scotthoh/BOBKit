set( clipper-core_sources 
${WRK_DIR}/checkouts/clipper/clipper/core/atomsf.cpp
${WRK_DIR}/checkouts/clipper/clipper/core/coords.cpp
${WRK_DIR}/checkouts/clipper/clipper/core/nxmap_operator.cpp
${WRK_DIR}/checkouts/clipper/clipper/core/cell.cpp
${WRK_DIR}/checkouts/clipper/clipper/core/derivs.cpp
${WRK_DIR}/checkouts/clipper/clipper/core/ramachandran.cpp
${WRK_DIR}/checkouts/clipper/clipper/core/clipper_instance.cpp
${WRK_DIR}/checkouts/clipper/clipper/core/fftmap.cpp
${WRK_DIR}/checkouts/clipper/clipper/core/resol_basisfn.cpp
${WRK_DIR}/checkouts/clipper/clipper/core/clipper_memory.cpp
${WRK_DIR}/checkouts/clipper/clipper/core/fftmap_sparse.cpp
${WRK_DIR}/checkouts/clipper/clipper/core/resol_fn.cpp
${WRK_DIR}/checkouts/clipper/clipper/core/clipper_message.cpp
${WRK_DIR}/checkouts/clipper/clipper/core/hkl_compute.cpp
${WRK_DIR}/checkouts/clipper/clipper/core/resol_targetfn.cpp
${WRK_DIR}/checkouts/clipper/clipper/core/clipper_stats.cpp
${WRK_DIR}/checkouts/clipper/clipper/core/hkl_data.cpp
${WRK_DIR}/checkouts/clipper/clipper/core/rotation.cpp
${WRK_DIR}/checkouts/clipper/clipper/core/clipper_test.cpp
${WRK_DIR}/checkouts/clipper/clipper/core/hkl_datatypes.cpp
${WRK_DIR}/checkouts/clipper/clipper/core/spacegroup.cpp
${WRK_DIR}/checkouts/clipper/clipper/core/clipper_types.cpp
${WRK_DIR}/checkouts/clipper/clipper/core/hkl_info.cpp
${WRK_DIR}/checkouts/clipper/clipper/core/spacegroup_data.cpp
${WRK_DIR}/checkouts/clipper/clipper/core/clipper_util.cpp
${WRK_DIR}/checkouts/clipper/clipper/core/hkl_lookup.cpp
${WRK_DIR}/checkouts/clipper/clipper/core/symop.cpp
${WRK_DIR}/checkouts/clipper/clipper/core/container.cpp
${WRK_DIR}/checkouts/clipper/clipper/core/hkl_operators.cpp
${WRK_DIR}/checkouts/clipper/clipper/core/container_hkl.cpp
${WRK_DIR}/checkouts/clipper/clipper/core/map_interp.cpp
${WRK_DIR}/checkouts/clipper/clipper/core/container_map.cpp
${WRK_DIR}/checkouts/clipper/clipper/core/map_utils.cpp
${WRK_DIR}/checkouts/clipper/clipper/core/xmap.cpp
${WRK_DIR}/checkouts/clipper/clipper/core/container_types.cpp
${WRK_DIR}/checkouts/clipper/clipper/core/nxmap.cpp
${WRK_DIR}/checkouts/clipper/clipper/core/clipper_thread.cpp
${WRK_DIR}/checkouts/clipper/clipper/core/test_core.cpp
${WRK_DIR}/checkouts/clipper/clipper/core/test_data.cpp
)
set(clipper-core_headers 
${WRK_DIR}/checkouts/clipper/clipper/core/atomsf.h
${WRK_DIR}/checkouts/clipper/clipper/core/coords.h
${WRK_DIR}/checkouts/clipper/clipper/core/nxmap_operator.h
${WRK_DIR}/checkouts/clipper/clipper/core/cell.h
${WRK_DIR}/checkouts/clipper/clipper/core/derivs.h
${WRK_DIR}/checkouts/clipper/clipper/core/ramachandran.h
${WRK_DIR}/checkouts/clipper/clipper/core/clipper_instance.h
${WRK_DIR}/checkouts/clipper/clipper/core/fftmap.h
${WRK_DIR}/checkouts/clipper/clipper/core/resol_basisfn.h
${WRK_DIR}/checkouts/clipper/clipper/core/clipper_memory.h
${WRK_DIR}/checkouts/clipper/clipper/core/fftmap_sparse.h
${WRK_DIR}/checkouts/clipper/clipper/core/clipper_precision.h
${WRK_DIR}/checkouts/clipper/clipper/core/resol_fn.h
${WRK_DIR}/checkouts/clipper/clipper/core/clipper_message.h
${WRK_DIR}/checkouts/clipper/clipper/core/hkl_compute.h
${WRK_DIR}/checkouts/clipper/clipper/core/resol_targetfn.h
${WRK_DIR}/checkouts/clipper/clipper/core/clipper_stats.h
${WRK_DIR}/checkouts/clipper/clipper/core/hkl_data.h
${WRK_DIR}/checkouts/clipper/clipper/core/rotation.h
${WRK_DIR}/checkouts/clipper/clipper/core/clipper_test.h
${WRK_DIR}/checkouts/clipper/clipper/core/hkl_datatypes.h
${WRK_DIR}/checkouts/clipper/clipper/core/spacegroup.h
${WRK_DIR}/checkouts/clipper/clipper/core/clipper_sysdep.h
${WRK_DIR}/checkouts/clipper/clipper/core/clipper_types.h
${WRK_DIR}/checkouts/clipper/clipper/core/hkl_info.h
${WRK_DIR}/checkouts/clipper/clipper/core/spacegroup_data.h
${WRK_DIR}/checkouts/clipper/clipper/core/clipper_util.h
${WRK_DIR}/checkouts/clipper/clipper/core/hkl_lookup.h
${WRK_DIR}/checkouts/clipper/clipper/core/symop.h
${WRK_DIR}/checkouts/clipper/clipper/core/container.h
${WRK_DIR}/checkouts/clipper/clipper/core/hkl_operators.h
${WRK_DIR}/checkouts/clipper/clipper/core/container_hkl.h
${WRK_DIR}/checkouts/clipper/clipper/core/map_interp.h
${WRK_DIR}/checkouts/clipper/clipper/core/container_map.h
${WRK_DIR}/checkouts/clipper/clipper/core/map_utils.h
${WRK_DIR}/checkouts/clipper/clipper/core/xmap.h
${WRK_DIR}/checkouts/clipper/clipper/core/container_types.h
${WRK_DIR}/checkouts/clipper/clipper/core/nxmap.h
${WRK_DIR}/checkouts/clipper/clipper/core/clipper_thread.h
${WRK_DIR}/checkouts/clipper/clipper/core/test_core.h
${WRK_DIR}/checkouts/clipper/clipper/core/test_data.h
)

set(clipper-contrib_sources
${WRK_DIR}/checkouts/clipper/clipper/contrib/convolution_search.cpp
${WRK_DIR}/checkouts/clipper/clipper/contrib/sfcalc.cpp
${WRK_DIR}/checkouts/clipper/clipper/contrib/edcalc.cpp
${WRK_DIR}/checkouts/clipper/clipper/contrib/sfcalc_obs.cpp
${WRK_DIR}/checkouts/clipper/clipper/contrib/fffear.cpp
${WRK_DIR}/checkouts/clipper/clipper/contrib/sfscale.cpp
${WRK_DIR}/checkouts/clipper/clipper/contrib/function_object_bases.cpp
${WRK_DIR}/checkouts/clipper/clipper/contrib/sfweight.cpp
${WRK_DIR}/checkouts/clipper/clipper/contrib/mapfilter.cpp
${WRK_DIR}/checkouts/clipper/clipper/contrib/skeleton.cpp
${WRK_DIR}/checkouts/clipper/clipper/contrib/originmatch.cpp
${WRK_DIR}/checkouts/clipper/clipper/contrib/test_contrib.cpp
)

set(clipper-contrib_headers
${WRK_DIR}/checkouts/clipper/clipper/contrib/convolution_search.h
${WRK_DIR}/checkouts/clipper/clipper/contrib/sfcalc.h
${WRK_DIR}/checkouts/clipper/clipper/contrib/edcalc.h
${WRK_DIR}/checkouts/clipper/clipper/contrib/sfcalc_obs.h
${WRK_DIR}/checkouts/clipper/clipper/contrib/fffear.h
${WRK_DIR}/checkouts/clipper/clipper/contrib/sfscale.h
${WRK_DIR}/checkouts/clipper/clipper/contrib/function_object_bases.h
${WRK_DIR}/checkouts/clipper/clipper/contrib/sfweight.h
${WRK_DIR}/checkouts/clipper/clipper/contrib/mapfilter.h
${WRK_DIR}/checkouts/clipper/clipper/contrib/skeleton.h
${WRK_DIR}/checkouts/clipper/clipper/contrib/originmatch.h
${WRK_DIR}/checkouts/clipper/clipper/contrib/test_contrib.h
)

set(clipper-phs_sources
${WRK_DIR}/checkouts/clipper/clipper/phs/phs_io.cpp
)

set(clipper-phs_headers
${WRK_DIR}/checkouts/clipper/clipper/phs/phs_io.h
)

set(clipper-cns_sources
${WRK_DIR}/checkouts/clipper/clipper/cns/cns_hkl_io.cpp
${WRK_DIR}/checkouts/clipper/clipper/cns/cns_map_io.cpp
)

set(clipper-cns_headers
${WRK_DIR}/checkouts/clipper/clipper/cns/cns_hkl_io.h
${WRK_DIR}/checkouts/clipper/clipper/cns/cns_map_io.h
)

#set(clipper-mmdb_sources
#${WRK_DIR}/checkouts/clipper/clipper/mmdb/clipper_mmdb.cpp
#)

#set(clipper-mmdb_headers
#${WRK_DIR}/checkouts/clipper/clipper/mmdb/clipper_mmdb.h
#)

set(clipper-cif_sources
${WRK_DIR}/checkouts/clipper/clipper/cif/cif_data_io.cpp
)

set(clipper-cif_headers
${WRK_DIR}/checkouts/clipper/clipper/cif/cif_data_io.h
)

set(clipper-ccp4_sources
${WRK_DIR}/checkouts/clipper/clipper/ccp4/ccp4_map_io.cpp
${WRK_DIR}/checkouts/clipper/clipper/ccp4/ccp4_mtz_io.cpp
${WRK_DIR}/checkouts/clipper/clipper/ccp4/ccp4_mtz_types.cpp
${WRK_DIR}/checkouts/clipper/clipper/ccp4/ccp4_utils.cpp
)

set(clipper-ccp4_headers
${WRK_DIR}/checkouts/clipper/clipper/ccp4/ccp4_map_io.h
${WRK_DIR}/checkouts/clipper/clipper/ccp4/ccp4_mtz_io.h
${WRK_DIR}/checkouts/clipper/clipper/ccp4/ccp4_mtz_types.h
${WRK_DIR}/checkouts/clipper/clipper/ccp4/ccp4_utils.h
)

# added clipper-gemmi files
set(clipper-gemmi_sources
${WRK_DIR}/checkouts/clipper/clipper/gemmi/clipper_gemmi.cpp
${WRK_DIR}/checkouts/clipper/clipper/gemmi/clipper_gemmi_model.cpp
)

set(clipper-gemmi_headers
${WRK_DIR}/checkouts/clipper/clipper/gemmi/clipper_gemmi.h
${WRK_DIR}/checkouts/clipper/clipper/gemmi/clipper_gemmi_model.h
)

# changed minimol_io to minimol_io_seq, minimol_io_mmdb, minimol_io_gemmi
set(clipper-minimol_sources
${WRK_DIR}/checkouts/clipper/clipper/minimol/container_minimol.cpp
${WRK_DIR}/checkouts/clipper/clipper/minimol/minimol.cpp	
${WRK_DIR}/checkouts/clipper/clipper/minimol/minimol_data.cpp
${WRK_DIR}/checkouts/clipper/clipper/minimol/minimol_io_seq.cpp
${WRK_DIR}/checkouts/clipper/clipper/minimol/minimol_io_gemmi.cpp
#${WRK_DIR}/checkouts/clipper/clipper/minimol/minimol_io_mmdb.cpp
${WRK_DIR}/checkouts/clipper/clipper/minimol/minimol_seq.cpp
${WRK_DIR}/checkouts/clipper/clipper/minimol/minimol_utils.cpp
${WRK_DIR}/checkouts/clipper/clipper/minimol/test_minimol_gemmi.cpp
)

set(clipper-minimol_headers
${WRK_DIR}/checkouts/clipper/clipper/minimol/container_minimol.h
${WRK_DIR}/checkouts/clipper/clipper/minimol/minimol.h
${WRK_DIR}/checkouts/clipper/clipper/minimol/minimol_data.h
${WRK_DIR}/checkouts/clipper/clipper/minimol/minimol_io_gemmi.h
#${WRK_DIR}/checkouts/clipper/clipper/minimol/minimol_io_mmdb.h
${WRK_DIR}/checkouts/clipper/clipper/minimol/minimol_io_seq.h
${WRK_DIR}/checkouts/clipper/clipper/minimol/minimol_seq.h
${WRK_DIR}/checkouts/clipper/clipper/minimol/minimol_utils.h
${WRK_DIR}/checkouts/clipper/clipper/minimol/test_minimol_gemmi.h
)

