set(SRC_LIST BaseVarC.cpp BamProcess.cpp BaseType.cpp algorithm.c)

include_directories("${CMAKE_SOURCE_DIR}/include" "${CMAKE_SOURCE_DIR}/SeqLib" "${CMAKE_SOURCE_DIR}/SeqLib/htslib")
# link_directories won't search the directory recursively
link_directories("${CMAKE_SOURCE_DIR}/SeqLib/src" "${CMAKE_SOURCE_DIR}/SeqLib/htslib")

add_executable(BaseVarC ${SRC_LIST})
target_link_libraries(BaseVarC PRIVATE fmt-header-only)
target_link_libraries(BaseVarC PRIVATE seqlib)
target_link_libraries(BaseVarC PRIVATE hts)

add_dependencies(BaseVarC seqlib_build)
