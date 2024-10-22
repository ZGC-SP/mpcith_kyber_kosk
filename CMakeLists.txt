cmake_minimum_required(VERSION 3.16)
project(PSS)

set(CMAKE_CXX_STANDARD 14)

# link the GMP and NTL libraries for algebraic operations
INCLUDE_DIRECTORIES(${PROJECT_SOURCE_DIR}/lib/gmp/include)
LINK_DIRECTORIES(${PROJECT_SOURCE_DIR}/lib/gmp/lib)

INCLUDE_DIRECTORIES(${PROJECT_SOURCE_DIR}/lib/ntl/include)
LINK_DIRECTORIES(${PROJECT_SOURCE_DIR}/lib/ntl/lib)

INCLUDE_DIRECTORIES(${PROJECT_SOURCE_DIR}/lib/openssl/include)
LINK_DIRECTORIES(${PROJECT_SOURCE_DIR}/lib/openssl/lib)

add_executable(PSS main.cpp src/packed_secret_sharing.cpp src/packed_secret_sharing.h src/utils.cpp src/utils.h src/plain_LWE_prover.cpp src/plain_LWE_prover.h benchmark.h src/frodo_prover.cpp src/frodo_prover.h test.h src/kyber_prover.cpp src/kyber_prover.h)

add_compile_options(-O3)

target_link_libraries(PSS gmp ntl m pthread ssl crypto)