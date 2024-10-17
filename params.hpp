
#ifndef MPCITH_KYBER_KOSK_PARAMS_HPP
#define MPCITH_KYBER_KOSK_PARAMS_HPP
extern "C" {
#include "kyber/params.h"
}

#ifndef KYBER_K
#define KYBER_K 2	/* Change this for different security strengths */
#endif

#if KYBER_K == 2
#define MPCITH_N 1454
#define MPCITH_T 150
#define MPCITH_L KYBER_N
#define MPCITH_K 70
#define MPCITH_V (KYBER_K * 2)

// TODO: confirm the parameters
#elif KYBER_K == 3
#define MPCITH_N 1454
#define MPCITH_T 150
#define MPCITH_L KYBER_N
#define MPCITH_K 70
#define MPCITH_V (KYBER_K * 2)

#elif KYBER_K == 4
#define MPCITH_N 1454
#define MPCITH_T 150
#define MPCITH_L KYBER_N
#define MPCITH_K 70
#define MPCITH_V (KYBER_K * 2)
#else
#error "just support KYBER_K = 2/3/4 level currently"
// more parameter settings for high security levels are needed to be computed
#endif

#endif //MPCITH_KYBER_KOSK_PARAMS_HPP
