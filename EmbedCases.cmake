# -*- mode: cmake; tab-width: 2; indent-tabs-mode: t; truncate-lines: t; compile-command: "cmake -Wdev" -*-
# vim: set filetype=cmake autoindent tabstop=2 shiftwidth=2 noexpandtab softtabstop=2 nowrap:

macro (opm_pack_file case_dir case_name suffix in_ext out_ext)
  # put the output in the same relative path in the output; we then
  # get no change if we build in-source
  set (rel_file "${case_dir}/${case_name}${suffix}")
  set (input_file "${PROJECT_SOURCE_DIR}/${rel_file}${in_ext}")
  set (output_file "${PROJECT_BINARY_DIR}/${rel_file}${out_ext}")

  # make sure that the output directory exists
  get_filename_component (output_dir "${output_file}" PATH)
  file (MAKE_DIRECTORY "${output_dir}")

  # run the shell script to encode the file
  set (pack_script "${PROJECT_SOURCE_DIR}/benchmarks/input/create_hex_data_file.sh")
  add_custom_command (
	OUTPUT "${output_file}"
	COMMAND "${pack_script}"
	ARGS "${input_file}" "${output_file}"
	DEPENDS "${input_file}" "${pack_script}"
	COMMENT "Creating packed binary of ${rel_file}"
	)

  # cannot add files to targets other than in add_custom_target,
  # and that command can only run once, so we must return a list
  # of dependencies that is added
  list (APPEND ${case_name}_DEPENDS "${output_file}")
endmacro (opm_pack_file)

# each case consists of a .grdecl file and a .data file
macro (opm_pack_case test_exe case_name)
  opm_pack_file ("benchmarks/input" "${case_name}" "_grid" ".grdecl" ".grdecl.gz.hex")
  opm_pack_file ("benchmarks/input" "${case_name}" "_upscaled_relperm" ".out" ".out.gz.hex")

  # we cannot add files directly (sic) but must wrap in a target
  add_custom_target (${case_name} DEPENDS ${${case_name}_DEPENDS})
  add_dependencies ("${test_exe}" "${case_name}")
	list(APPEND OPM_BENCHMARKS ${test_exe})
endmacro (opm_pack_case)

# rel.perm curve is packed separately because it is common for all cases
macro (opm_pack_stone test_exe)
	opm_pack_file ("benchmarks/input" "stonefile" "_benchmark" ".txt" ".txt.gz.hex")
	add_custom_target (stonefile ALL DEPENDS ${stonefile_DEPENDS})
	add_dependencies ("${test_exe}" "stonefile")
endmacro (opm_pack_stone)

# pack these cases which are alternatives in the code
if(Boost_IOSTREAMS_FOUND)
  opm_pack_stone (upscale_relperm_benchmark)
  opm_pack_case (upscale_relperm_benchmark benchmark20)

  if(INSTALL_BENCHMARKS)
    add_custom_target(benchmarks ALL DEPENDS upscale_relperm_benchmark)
    set_target_properties(upscale_relperm_benchmark PROPERTIES EXCLUDE_FROM_ALL 0)
  else()
    add_custom_target(benchmarks DEPENDS upscale_relperm_benchmark)
  endif()
endif()
