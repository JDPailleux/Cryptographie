#include "DM.h"


void square_multiply(mpz_t r, mpz_t a, mpz_t n, mpz_t h){ // r : resultat de la fonction, n : modulo, h : exposant 
	mpz_t zero, un, deux, modu;	// "modu" va contenir le reste de la division euclidienne de h par 2 à chaque itération.	
	mpz_init_set_str(zero,"0",10); // zero = 0 
	mpz_init_set_str(un,"1",10); // un = 1
	mpz_init_set_str(deux,"2",10); // deux = 2 
	mpz_add(r, zero ,a);	    // r = 0 + a => r = a

	while (mpz_cmp(h, zero) > 0){
		
		mpz_mul(r, r, r); // r = r*r
		mpz_mod(r, r, n); // r = r*r modulo n
		mpz_fdiv_qr(h, modu, h, deux); // modu = reste de la div et h = h/2 pour la prochaine itération
		if ( mpz_cmp(modu,un) == 0){ // On test h%2 == 1
			mpz_mul(r, r, a); // r = r*a
			mpz_mod(r, r, n); // r = r*a modulo n
		}
	}
}

void testFermat(mpz_t n, int k){	// n: entier et k nombre de répétition
	int i, res =1;
	mpz_t val, a, un, n_1;
	mpz_init_set_str(un,"1",10);   // un = 1 
	mpz_init(a); 
	mpz_init(val); 
	mpz_init(n_1);
	mpz_sub(n_1, n, un);	//n_1 = n-1	

	// Valeur aléatoire 
	gmp_randstate_t state;
	gmp_randinit_mt(state);
	gmp_randseed_ui(state, time(NULL));

	for (i= 1; i <= k; i++){
		mpz_urandomm(a, state, n); 
		square_multiply(val, a, n, n_1);
		if ( mpz_cmp(val, un) == 0) res = 0; 	// n est composé
	}
	
	if ( res == 1) gmp_printf("Le nombre %Zd est supposé premier \n",n);
	else  gmp_printf("Le nombre %Zd est composé \n \n",n);
}

/*void testMillerRabin(mpz_t n, int k){


}*/

int main(){ //int argc, char** argv
	
	mpz_t n;
	mpz_init_set_str(n,"11",10); 
	// Test de la fonction de Fermat:
	printf("#### TEST DE FERMAT #### \n n = 11, rep = 4\n");
	
	testFermat(n,4);

	// Test de la bibliothèque GMP:
	mpz_t rop, base, exp, mod;
	mpz_init(rop);
	mpz_init_set_str(base,"2",10);
	mpz_init_set_str(exp, "5", 10);
	mpz_init_set_str(mod, "15",10);
	square_multiply(rop, base, mod, exp);
	gmp_printf(" %Zd ^ %Zd mod %Zd = %Zd\n",base,exp, mod, rop);
	return 0;
}
