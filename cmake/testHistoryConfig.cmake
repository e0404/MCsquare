if(NOT DEFINED src OR NOT DEFINED dst)
  message(FATAL_ERROR "Need -Dsrc and -Ddst")
endif()

# Read, replace 'Num_Primaries: 10000000' -> 'Num_Primaries: 1000', write
file(READ "${src}" _c)
string(REGEX REPLACE "Num_Primaries[ \t]*10000000"
                     "Num_Primaries \t1000"
                     _c "${_c}")
file(WRITE "${dst}" "${_c}")