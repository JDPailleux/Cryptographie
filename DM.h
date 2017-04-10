#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <gmp.h>

void testFermat(mpz_t n, int k);

void testMillerRabin(mpz_t n, int k);

void square_multiply(mpz_t r, mpz_t a, mpz_t n, mpz_t h);
