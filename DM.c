#include "DM.h"


void square_multiply(mpz_t r, mpz_t a, mpz_t n, mpz_t h){ // r : resultat de la fonction, n : modulo, h : exposant 
	mpz_t modu;	// "modu" va contenir le reste de la division euclidienne de h par 2 à chaque itération.	
	mpz_init(modu);
	mpz_set(r,a);
	while (mpz_cmp_ui(h, 0) > 0){
		
		mpz_mul(r, r, r); // r = r*r
		mpz_mod(r, r, n); // r = r*r modulo n
		mpz_fdiv_qr_ui(h, modu, h, 2); // modu = reste de la div et h = h/2 pour la prochaine itération
		if ( mpz_cmp_ui(modu,1) == 0){ // On test h%2 == 1
			mpz_mul(r, r, a); // r = r*a
			mpz_mod(r, r, n); // r = r*a modulo n
		}
	}
}

void testFermat(mpz_t n, int k){	// n: entier et k nombre de répétition
	int i, res =1;
	mpz_t val, a, n_1, n_3; 		 // n_1 : correspond à n-1, n_3 : correspond à n-3,
	mpz_inits(a, val, n_1, n_3, (void*)0);   // a: Correspond à la valeur aléatoire, val : stock la valeur de retour de suare_multiply
	
	mpz_sub_ui(n_1, n, 1);	//n_1 = n-1	
	mpz_sub_ui(n_3, n, 3);	//n_1 = n-3

	// Valeur aléatoire 
	gmp_randstate_t state;
	gmp_randinit_mt(state);
	gmp_randseed_ui(state, time(NULL));

	for (i= 1; i <= k; i++){ 
		mpz_urandomm(a, state, n_3); // On a : 0 <= a <= n-4
		mpz_add_ui(a,a,2); // a = a+2. On doit alors avoir 1 < a < n-1
		square_multiply(val, a, n, n_1); 
		
		if ( mpz_cmp_ui(val, 1) == 0) res = 0; 	// n est composé
	}
	
	if ( res == 1) gmp_printf("Le nombre %Zd est premier \n",n);
	else  gmp_printf("Le nombre %Zd est composé \n \n",n);
}
/*
void decompo(mpz_t t, mpz_t s, mpz_t n){ // Fonction qui décompose n-1. (n vaudra n-1)
	int boolp = 0; // Entier qui va nous permettre de sortir des deux boucles.	
	mpz_t tmp;
	mpz_init(tmp);

	// On initialise t et s.
	mpz_init(s);
	mpz_init(t);
	// Deux boucles imbriquées qui vont trouver t et s tel que 2^s * t = n-1
	for ( mpz_sqrt(s, n); mpz_cmp_ui(s,1) >= 0; mpz_sub_ui(s,s,1)){
		mpz_set_ui(t,1); // t = 1
		
		while (1){

		}
		mpz_add_ui(t,t,2); // t += 2 (car t doit être impair)
	}
}*/

void testMillerRabin(mpz_t n, mpz_t k){ // n : nombre à tester, k : nombre de répétition

	if ( mpz_cmp(n,2) == 0) gmp_printf("premier\n");
	else if ( mpz_mod(n,2) == 0 ) gmp_printf("composé\n"); // cas nombre paire != 2
	else {

	int booleen =1;
	mpz_t t, i, a, y, j, s, n1; // n1 => n+1
	mpz_inits(t, i, a, y, j, s, (void*) 0);

	for (mpz_set_ui(i,1); mpz_cmp_ui(i,k) <= 0; mpz_add_ui(i, 1)){ //for (i= 1; i <= k; i++)

		mpz_urandomm(a, state, n1);  // a appartient à [0,n]
		square_multiply(y, a, n, t);
		if ( (mpz_cmp_ui(y, 1) != 1) && (mpz_cmp_ui(y, 1) != -1){
			for (mpz_set_ui(i,1); mpz_cmp_ui(i,k) <= 0; mpz_add_ui(i, 1)){ //for (j= 1; j <= s-1; j++)
				
				mpz_mul(y, y, y); // y = y*y
				mpz_mod(y, y, n); // y = (y*y) % n
				if ( mpz_cmp_ui(y,1) == 0 ) return "composé";
				if ( mpz_cmp_ui(y,-1) == 0 ) break;
			}
		}
	}
	
	if ( booleen == 1) gmp_printf("Le nombre %Zd est premier \n",n);
	else  gmp_printf("Le nombre %Zd est composé \n \n",n);

}

int main(){ //int argc, char** argv
	
	mpz_t n;
	mpz_init_set_str(n,"11",10); 
	// Test de la fonction de Fermat:
	printf("#### TEST DE FERMAT #### \n n = 11, rep = 4\n");
	
	testFermat(n,15);

	/*// Test de la bibliothèque GMP:
	mpz_t rop, base, exp, mod;
	mpz_init(rop);
	mpz_init_set_str(base,"2",10);
	mpz_init_set_str(exp, "5", 10);
	mpz_init_set_str(mod, "15",10);
	square_multiply(rop, base, mod, exp); */
	return 0;
}
