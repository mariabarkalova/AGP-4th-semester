cmake_minimum_required(VERSION 3.26)

set(PROJECT_NAME_STR AGP-4th-semester)
set(CMAKE_CXX_STANDARD 17)

project (${PROJECT_NAME_STR})

find_package(OpenMP REQUIRED)

include_directories(${CMAKE_CURRENT_SOURCE_DIR}/src/ ${CMAKE_CURRENT_SOURCE_DIR}/include/ ${CMAKE_CURRENT_SOURCE_DIR}/libs/sample_src ${CMAKE_CURRENT_SOURCE_DIR}/libs/src)

file(GLOB SRC_FILES 
                    ${CMAKE_CURRENT_SOURCE_DIR}/libs/sample_src/*.h
                    ${CMAKE_CURRENT_SOURCE_DIR}/libs/sample_src/*.hpp
                    ${CMAKE_CURRENT_SOURCE_DIR}/libs/sample_src/GKLS/*.hpp
                    ${CMAKE_CURRENT_SOURCE_DIR}/libs/sample_src/GKLS/*.cpp
                    ${CMAKE_CURRENT_SOURCE_DIR}/libs/sample_src/Grishagin/*.hpp
                    ${CMAKE_CURRENT_SOURCE_DIR}/libs/sample_src/Grishagin/*.cpp
                    ${CMAKE_CURRENT_SOURCE_DIR}/libs/sample_src/Hill/*.cpp
                    ${CMAKE_CURRENT_SOURCE_DIR}/libs/sample_src/Shekel/*.cpp
                    ${CMAKE_CURRENT_SOURCE_DIR}/libs/src/*.cpp
                    ${CMAKE_CURRENT_SOURCE_DIR}/libs/src/*.h
                    ${CMAKE_CURRENT_SOURCE_DIR}/libs/src/*.hpp
		    ${CMAKE_CURRENT_SOURCE_DIR}/evolvent/*.cpp
   	            ${CMAKE_CURRENT_SOURCE_DIR}/evolvent/*.hpp
                    ${CMAKE_CURRENT_SOURCE_DIR}/src/*.cpp
		    ${CMAKE_CURRENT_SOURCE_DIR}/src/*.h
    	            ${CMAKE_CURRENT_SOURCE_DIR}/src/*.hpp)

add_executable(${PROJECT_NAME_STR} ${SRC_FILES})

target_link_libraries(${PROJECT_NAME_STR} PRIVATE OpenMP::OpenMP_CXX)