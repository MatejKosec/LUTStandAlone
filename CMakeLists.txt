##########################################
#CMake file for LuT standalone
#########################################
cmake_minimum_required(VERSION 2.8)
project (LuT)
file(GLOB_RECURSE LuT_SOURCES "src/*.cpp")
file(GLOB_RECURSE LuT_HEADERS "src/*.h")
set (LuT_INCLUDE_DIRS "")
foreach (_headerFile ${LuT_HEADERS})
    get_filename_component(_dir ${_headerFile} PATH)
    list (APPEND LuT_INCLUDE_DIRS ${_dir})
endforeach()
list(REMOVE_DUPLICATES LuT_INCLUDE_DIRS)
include_directories(${LuT_INCLUDE_DIRS})
add_executable (LuTExe ${LuT_SOURCES})
