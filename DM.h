#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>
#include <gmp.h>

void square_multiply(mpz_t r, mpz_t a, mpz_t n, mpz_t h);

int testFermat(mpz_t n, mpz_t k);

void my_pow_2(mpz_t dest, mpz_t exp);

void decompo(mpz_t t, mpz_t s, mpz_t n);

int testMillerRabin(mpz_t n, mpz_t k);


