cmake_minimum_required (VERSION 3.1)

project(BaseVarC)

set(CMAKE_CXX_FLAGS "-g -Wall -O3 -pthread")
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
# Don't use e.g. GNU extension (like -std=gnu++11) for portability
set(CMAKE_CXX_EXTENSIONS OFF)
set(EXECUTABLE_OUTPUT_PATH ${CMAKE_SOURCE_DIR}/bin)


add_custom_target(
  seqlib_build ALL
  COMMAND ./configure && make
  WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}/SeqLib
  COMMENT "build seqlib"
  )

# add_custom_command(
#   OUTPUT ${CMAKE_SOURCE_DIR}/SeqLib
#   COMMAND ./configure && make
#   WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}/SeqLib
#   )

# add_custom_target(
#   seqlib_build ALL
#   DEPENDS ${CMAKE_SOURCE_DIR}/SeqLib
#   )

add_subdirectory(fmt EXCLUDE_FROM_ALL)
add_subdirectory(src)
