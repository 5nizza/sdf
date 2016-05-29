cmake_minimum_required(VERSION 2.6)
project(sdf)
set(PROJECT sdf)
set(VERSION 0.0.0)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -std=c++11 -pthread")


set(SPDLOG_INCLUDE_DIR    "$ENV{SDF_TP}/spdlog/")
set(AIGER_INCLUDE_DIR    "$ENV{SDF_TP}/aiger-1.9.4/")
set(AIGER_LIBRARIES      "$ENV{SDF_TP}/aiger-1.9.4/aiger.o")
#set(CUDD_INCLUDE_DIR   "$ENV{SDF_TP}/cudd/cudd/" "$ENV{SDF_TP}/cudd/obj/" "$ENV{SDF_TP}/cudd/mtr/" "$ENV{SDF_TP}/cudd/epd/" "$ENV{SDF_TP}/cudd/util/")
#set(CUDD_LIBRARIES    "$ENV{SDF_TP}/cudd/obj/libobj.a" "$ENV{SDF_TP}/cudd/cudd/libcudd.a" "$ENV{SDF_TP}/cudd/mtr/libmtr.a" "$ENV{SDF_TP}/cudd/epd/libepd.a" "$ENV{SDF_TP}/cudd/util/libutil.a" "$ENV{SDF_TP}/cudd/st/libst.a")
set(CUDD_INCLUDE_DIR   $ENV{SDF_TP}/cudd-3.0.0/cudd/ $ENV{SDF_TP}/cudd-3.0.0/cplusplus/ $ENV{SDF_TP}/cudd-3.0.0/mtr/)
set(CUDD_LIBRARIES   $ENV{SDF_TP}/cudd-3.0.0/cudd/.libs/libcudd.a)

add_subdirectory(src)

enable_testing()
add_subdirectory(func_tests)

## how are shared libraries found ?    ## AK: should it be in the beginning?
set(LD_LIB_VAR "LD_LIBRARY_PATH")
