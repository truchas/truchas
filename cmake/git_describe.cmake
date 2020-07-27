set(TRUCHAS_VER 3.1.0-alpha)

message(DEBUG "GIT_FOUND: ${GIT_FOUND}")
message(DEBUG "GIT: ${GIT}")
message(DEBUG "INFILE: ${INFILE}")
message(DEBUG "OUTFILE: ${OUTFILE}")

if(GIT_FOUND)
  execute_process(
    COMMAND ${GIT} describe --tags --dirty
    RESULT_VARIABLE result
    OUTPUT_VARIABLE TRUCHAS_VER
    OUTPUT_STRIP_TRAILING_WHITESPACE
    )
endif()
configure_file(${INFILE} ${OUTFILE} @ONLY)
