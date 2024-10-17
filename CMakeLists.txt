cmake_minimum_required(VERSION 3.29)
project(mpcith_kyber_kosk)

set(CMAKE_CXX_STANDARD 11)

set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -O3 -g")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -g")

INCLUDE_DIRECTORIES(${PROJECT_SOURCE_DIR}/ntl/include)
LINK_DIRECTORIES(${PROJECT_SOURCE_DIR}/ntl/lib)

add_executable(mpcith_kyber_kosk main.cpp
        kyber/api.h
        kyber/cbd.c
        kyber/cbd.h
        kyber/fips202.c
        kyber/fips202.h
        kyber/indcpa.c
        kyber/indcpa.h
        kyber/kem.c
        kyber/kem.h
        kyber/ntt.c
        kyber/ntt.h
        kyber/params.h
        kyber/poly.c
        kyber/poly.h
        kyber/polyvec.c
        kyber/polyvec.h
        kyber/randombytes.c
        kyber/randombytes.h
        kyber/reduce.c
        kyber/reduce.h
        kyber/symmetric.h
        kyber/symmetric-shake.c
        kyber/verify.c
        kyber/verify.h
        utils/gf3329.c
        utils/gf3329.h
        utils/precomputed_kyber.c
        utils/precomputed_kyber.h
        ss.cpp
        ss.hpp
        params.hpp
        kosk.cpp
        kosk.hpp
        mlwe_prover.cpp
        mlwe_prover.hpp
        mlwe_verifier.cpp
        mlwe_verifier.hpp
)

target_link_libraries(mpcith_kyber_kosk ntl m pthread)