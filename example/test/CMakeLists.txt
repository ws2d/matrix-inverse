cmake_minimum_required(VERSION 3.13)

## We must specify the C and C++ compilers before the project command
# set(CMAKE_C_COMPILER "/usr/bin/clang")
# set(CMAKE_CXX_COMPILER "/usr/bin/clang++")
set(CMAKE_CXX_FLAGS "-Wall")
set(CMAKE_CXX_STANDARD 17)

project(test_example)

# Import Microsoft SEAL
find_package(SEAL 4.1.2 EXACT REQUIRED
              PATHS /usr/local/include/) 

# 指定 NTL 的头文件路径
include_directories(${CMAKE_SOURCE_DIR}/../build/include)

# 指定 NTL 的库路径
link_directories(${CMAKE_SOURCE_DIR}/../build/lib)
# # Import Microsoft SEAL
# find_package(SEAL 4.1.1 EXACT REQUIRED
#               PATHS /home/hs/seal/SEAL-main) 



#   找到seal版本
# find_package(SEAL 4.1 REQUIRED)
# set(SEAL_LIBRARY "SEAL::seal")

# Import OpenMP
# FIND_PACKAGE( OpenMP REQUIRED)
# if(OPENMP_FOUND)
#     message("OPENMP FOUND")
#     set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
#     set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
#     set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
# endif()


set(REAL_SOURCE_DIR ${PROJECT_SOURCE_DIR}/../src)
set(TEST_SOURCE_DIR ${PROJECT_SOURCE_DIR})



add_executable(test_example
               ${TEST_SOURCE_DIR}/test_example.cpp
            #    ${TEST_SOURCE_DIR}/../src/utils.cpp
               ${TEST_SOURCE_DIR}/../src/matrix.cpp
               ${TEST_SOURCE_DIR}/../src/Bootstrapper.cpp
            # #    ${TEST_SOURCE_DIR}/../src/enc_mat_mul.cpp
            #    ${TEST_SOURCE_DIR}/../src/long_enc_cipher.h
            #    ${TEST_SOURCE_DIR}/../src/jiang_enc_cipher.h
            #    ${TEST_SOURCE_DIR}/../src/strassen_block.h
            #    ${TEST_SOURCE_DIR}/../src/R_T_enc_cipher.h
            #    ${TEST_SOURCE_DIR}/../src/coprime_enc_cipher.h
            #    ${TEST_SOURCE_DIR}/../src/block_long_cipher.h
               )
# 链接 NTL 库
target_link_libraries(test_example ntl gmp m)        


set_target_properties(test_example PROPERTIES LINK_FLAGS "-static")
                    
if(TARGET SEAL::seal)
    target_link_libraries(test_example SEAL::seal)
elseif(TARGET SEAL::seal_shared)
    target_link_libraries(test_example SEAL::seal_shared)
else()
    message(FATAL_ERROR "Cannot find target SEAL::seal or SEAL::seal_shared")
endif()
