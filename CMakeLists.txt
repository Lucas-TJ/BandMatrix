# CMakeLists.txt for BandMatrix plugin
cmake_minimum_required(VERSION 3.12)
project(BandMatrix VERSION 1.0 LANGUAGES CXX)

find_package(Sofa.Core REQUIRED)

set(HEADER_FILES
    src/BandMatrix/config.h.in
    src/BandMatrix/init.h
    src/BandMatrix/BandMatrix.h
    src/BandMatrix/BandMatrix.inl
    

)

set(SOURCE_FILES
    src/BandMatrix/init.cpp
    src/BandMatrix/BandMatrix.cpp
    
)

add_library(${PROJECT_NAME} SHARED ${HEADER_FILES} ${SOURCE_FILES})
target_link_libraries(${PROJECT_NAME} PUBLIC Sofa.Core)

sofa_create_package_with_targets(
    PACKAGE_NAME ${PROJECT_NAME}
    PACKAGE_VERSION ${PROJECT_VERSION}
    TARGETS ${PROJECT_NAME} AUTO_SET_TARGET_PROPERTIES
    INCLUDE_SOURCE_DIR "src"
    INCLUDE_INSTALL_DIR ${PROJECT_NAME}
    RELOCATABLE "plugins"
    )
    
cmake_dependent_option(BANDMATRIX_BUILD_TESTS "Compile the automatic tests" ON "SOFA_BUILD_TESTS OR NOT DEFINED SOFA_BUILD_TESTS" OFF)
if(BANDMATRIX_BUILD_TESTS)
    enable_testing()
    add_subdirectory(tests)
endif()
