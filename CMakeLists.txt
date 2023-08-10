cmake_minimum_required (VERSION 3.12)
project(psdetsim LANGUAGES CXX)

cmake_policy(SET CMP0086 NEW)
cmake_policy(SET CMP0078 NEW)
set (CMAKE_CXX_STANDARD 11)


### import modules
find_package(SWIG)
if (SWIG_FOUND)
  include(${SWIG_USE_FILE})
  include_directories(${CMAKE_SOURCE_DIR} /usr/include/python3.6m)
endif()


list(APPEND CMAKE_PREFIX_PATH $ENV{ROOTSYS})


#---Locate the ROOT package and defines a number of variables (e.g. ROOT_INCLUDE_DIRS)
find_package(ROOT REQUIRED COMPONENTS Core RIO RooFit RooFitCore RooStats Hist Tree Matrix Physics MathCore Net)

#---Define useful ROOT functions and macros (e.g. ROOT_GENERATE_DICTIONARY)
include(${ROOT_USE_FILE})

include_directories(${CMAKE_SOURCE_DIR} ${ROOT_INCLUDE_DIRS})
add_definitions(${ROOT_CXX_FLAGS})

include_directories(${PROJECT_SOURCE_DIR}/include)


ROOT_GENERATE_DICTIONARY(G__decaypath  ${PROJECT_SOURCE_DIR}/include/decaypath.hh  ${PROJECT_SOURCE_DIR}/include/common.hh LINKDEF  ${PROJECT_SOURCE_DIR}/include/decaypathLinkDef.hh)
add_library(decaypath SHARED  ${PROJECT_SOURCE_DIR}/include/decaypath.hh  ${PROJECT_SOURCE_DIR}/src/decaypath.cc G__decaypath.cxx)
target_link_libraries(decaypath ${ROOT_LIBRARIES})

add_executable(simulationC
    ${PROJECT_SOURCE_DIR}/src/simulation.cc ${PROJECT_SOURCE_DIR}/include/simulation.hh
    mainsimulation.cc)
target_link_libraries(simulationC ${ROOT_LIBRARIES})
target_link_libraries(simulationC decaypath)

if (SWIG_FOUND)
  set_property(SOURCE ${PROJECT_SOURCE_DIR}/include/simulation.i PROPERTY CPLUSPLUS ON)
  swig_add_library(simulation
    TYPE SHARED
    LANGUAGE python
    SOURCES ${PROJECT_SOURCE_DIR}/include/simulation.i  ${PROJECT_SOURCE_DIR}/src/simulation.cc
    )
  swig_link_libraries(simulation decaypath)
  target_include_directories(simulation PUBLIC "")
endif()

#----------------------------------------------------------------------------
# Copy all resources including scripts, examples, text files, etc... to the build directory, i.e. the directory in which we
# build project. This is so that we can us the resources because it
# relies on these scripts being in the current working directory.
#
set(PROJ_DATA
  simparmsex.txt
  parmsex.txt
  FRDM-QRPA12-halflife.txt
  stable.csv
  hx.root
  hy.root
  himpx.root
  himpy.root
  hdx.root
  hdy.root
  tdiffdataex.root
  )


foreach(_script ${PROJ_DATA})
  configure_file(
    ${PROJECT_SOURCE_DIR}/${_script}
    ${PROJECT_BINARY_DIR}/${_script}
    COPYONLY
    )
endforeach()

