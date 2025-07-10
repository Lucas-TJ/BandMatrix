cmake_minimum_required(VERSION 3.12)
project(BandMatrix_test VERSION 1.0)
find_package(Sofa.Testing REQUIRED)

set(HEADER_FILES
    BandMatrix_test.h
)

set(SOURCE_FILES
    test.cpp
    BandMatrix_test.cpp
)

add_executable(${PROJECT_NAME} ${SOURCE_FILES})

target_link_libraries(${PROJECT_NAME} Sofa.Testing BandMatrix Sofa.LinearAlgebra.Testing )

add_test(NAME ${PROJECT_NAME} COMMAND ${PROJECT_NAME})
