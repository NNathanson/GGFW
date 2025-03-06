cmake_minimum_required(VERSION 3.0 FATAL_ERROR)
project(GGFW LANGUAGES CXX)
list (APPEND CMAKE_PREFIX_PATH $ENV{ROOTSYS})

find_package(ROOT 6.16 CONFIG REQUIRED)


include ("${ROOT_USE_FILE}")
set(CMAKE_CXX_FLAGS "${ROOT_CXX_FLAGS}")

include_directories(${ROOT_INCLUDE_DIRS} ${CMAKE_CURRENT_SOURCE_DIR})
message("${ROOT_INCLUDE_DIRS}")
message("${CMAKE_CURRENT_SOURCE_DIR}")

ROOT_GENERATE_DICTIONARY(G__GGFW
                BootstrapProfile.h
                Calculator.h
                FlowContainer.h
                FlowPtContainer.h
                GFW.h
                GFWBinningCuts.h
                GFWCorrConfigs.h
                GFWCumulant.h
                GFWPowerArray.h
                GFWRegions.h
                OutputHandler.h
                ProfileSubset.h
                StatHandler.h
                LINKDEF GGFWLinkDef.h) 

add_library(GGFW SHARED
        G__GGFW.cxx
        BootstrapProfile.cxx
        Calculator.cxx
        FlowContainer.cxx
        FlowPtContainer.cxx
        GFW.cxx
        GFWCumulant.cxx
        GFWPowerArray.cxx
        OutputHandler.cxx
        ProfileSubset.cxx
        StatHandler.cxx
        )

set_target_properties(GGFW PROPERTIES ENABLE_EXPORTS 1)
target_include_directories(GGFW PUBLIC ${CMAKE_SOURCE_DIR})
target_link_libraries(GGFW PUBLIC ${ROOT_LIBRARIES})

