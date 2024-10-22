#include <iostream>

#include "src/utils.h"
#include "benchmark.h"
#include "test.h"
#include <thread>

int main() {
//    PlainTest(32, true);
//    PlainTest(64, false);
//    TestAddEvaluation(128, true);
//    TestAddEvaluation(128, false);
//    TestMulEvaluation(32);

    long prime = GenPrime_long(61);
    ZZ_p::init(ZZ(prime));
//    ZZ_p::init(ZZ(3329));

//    std::cout << "Multiplication Gate: " << std::endl;
//    PackedSecretSharingMultBenchmark();

    std::cout << "Linear Transformation Gate: " << std::endl;
    LinearTransformBenchmark();

    std::cout << "Pre-processing Gate: " << std::endl;
    PreprocessingBenchmark();

//    std::cout << "Plain LWE Batch Prove: " << std::endl;
//    PlainLweBatchProveBenchmark();

//    std::cout << "Frodo Prover: " << std::endl;
//    FrodoProverBenchmark();

//    std::cout << "Kyber Prover: " << std::endl;
//    KyberProverBenchmark();

    return 0;
}
