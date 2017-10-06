#include "DM.h"

void square_multiply(mpz_t r, mpz_t a, mpz_t n, mpz_t h);

int testFermat(mpz_t n, mpz_t k){	// n: entier et k nombre de répétition
	if ( (mpz_cmp_ui(n,3) == 0) || (mpz_cmp_ui(n,2) == 0)) return 1;
	if ( mpz_cmp_ui(n,1) == 0 ) return 0;
	int res =1;
	mpz_t i, val, a, n_1, n_3; 		 // n_1 : correspond à n-1, n_3 : correspond à n-3,
	mpz_inits(i, a, val, n_1, n_3, (void*)0);// a: Correspond à la valeur aléatoire, val : stock la valeur de retour de suare_multiply
	
	mpz_sub_ui(n_1, n, 1);	//n_1 = n-1	
	mpz_sub_ui(n_3, n, 3);	//n_3 = n-3

	// Valeur aléatoire 
	gmp_randstate_t state;
	gmp_randinit_mt(state);
	gmp_randseed_ui(state, time(NULL));
	
	for (mpz_set_ui(i,1); mpz_cmp(i,k) <= 0; mpz_add_ui(i, i, 1)){ //for (i= 1; i <= k; i++)
		mpz_urandomm(a, state, n_3); // On a : 0 <= a <= n-4 
		mpz_add_ui(a,a,2); // a = a+2. On doit alors avoir 1 < a < n-1
		square_multiply(val, a, n, n_1); 
		mpz_sub_ui(n_1, n, 1);	//n_1 = n-1 car il est mis à zéro dans la fonction précédente

		if ( mpz_cmp_ui(val, 1) != 0) res = 0; 	// n est composé
	}
	
	if ( res == 1) return 1;
	else  return 0;	
	
	// Libération mémoire
	mpz_clear(i); mpz_clear(a); mpz_clear(val); mpz_clear(n_3); mpz_clear(n_1); gmp_randclear(state);
}

void my_pow_2(mpz_t dest, mpz_t exp){ // calcul 2^exp et stock le résultat dans dest.
	
	mpz_t i;
	mpz_init(i);
	mpz_set_ui(dest, 1); // Dest = 1 

	for (mpz_set_ui(i,1); mpz_cmp(i,exp) <= 0; mpz_add_ui(i, i, 1)){ //for (i= 1; i <= exp; i++)
		mpz_mul_ui(dest,dest,2); // dest = dest*2
	}
	mpz_clear(i);
}

void decompo(mpz_t t, mpz_t s, mpz_t n){ // Fonction qui décompose n-1. (n vaudra n-1)
	int boolp = 0; // Entier qui va nous permettre de sortir des deux boucles.	
	mpz_t tmp;

	// On initialise t et s.
	mpz_inits(s, t, tmp, (void*)0);

	// Deux boucles imbriquées qui vont trouver t et s tel que 2^s * t = n-1
	for ( mpz_sqrt(s, n); mpz_cmp_ui(s,1) >= 0; mpz_sub_ui(s,s,1)){
		mpz_set_ui(t,1); // t = 1
		do{
			my_pow_2(tmp,s);
			mpz_mul(tmp, tmp, t);
			if (mpz_cmp(tmp,n) == 0) { boolp = 1; break;} // On a trouvé t et s.
			mpz_add_ui(t,t,2); // t += 2 (car t doit être impair)
		}while (mpz_cmp(tmp,n) <= 0);

		if(boolp == 1) break; // On sort aussi de la deuxième boucle.
	}
}

int testMillerRabin(mpz_t n, mpz_t k){ // n : nombre à tester, k : nombre de répétition  (renvoie 1 : premier, renvoie 0 : composé)
	if ( mpz_cmp_ui(n,1) == 0) return 0;
	else if ( mpz_cmp_ui(n,2) == 0) return 1;
	else {
	int next;
	mpz_t tmp_t, t, i, a, y, j, s, n_1; // n_1 => n-1
	// tmp_t sera utilisé à la place de t dans square_multiply car la fonction va le mettre à zéro (donc perte de la valeur de t)
	mpz_inits(tmp_t, t, i, a, y, j, s, n_1,(void*) 0); 
	mpz_sub_ui(n_1,n,1);
	decompo(t,s,n_1); // On décompose n-1 + initilisation des valeurs de t et s
	mpz_set(tmp_t,t);

	// Valeur aléatoire 
	gmp_randstate_t state;
	gmp_randinit_mt(state);
	gmp_randseed_ui(state, time(NULL));

	for (mpz_set_ui(i,1); mpz_cmp(i,k) <= 0; mpz_add_ui(i, i, 1)){ //for (i= 1; i <= k; i++)
		next = 0;
		mpz_urandomm(a, state, n_1);  // on tire a entre 0 et n-2
		mpz_add_ui(a,a,1); // On a : 0 < a < n
		square_multiply(y, a, n, tmp_t);
		mpz_set(tmp_t,t); // tmp_t = t car il a été mis à zéro par square_multiply
		
		if ( (mpz_cmp_ui(y, 1) != 0) && (mpz_cmp(y, n_1) != 0)){

				for (mpz_set_ui(j,1); mpz_cmp(j,s) < 0; mpz_add_ui(j, j, 1)){ //for (j= 1; j <= s-1; j++)
					
				mpz_mul(y, y, y); // y = y*y
				mpz_mod(y, y, n); // y = (y*y) % n
				if ( mpz_cmp_ui(y,1) == 0 ){ 

					// Libération mémoire
					mpz_clear(t); mpz_clear(i); mpz_clear(a); mpz_clear(y); mpz_clear(j); mpz_clear(s);
					mpz_clear(n_1);gmp_randclear(state);
					return 0; // On arrête la fonction.
				}
				if ( mpz_cmp(y,n_1) == 0 ) {
					next = 1;
					break; // On sort de la boucle j et on continue i
				}
			}
		if ( next == 1) break; // On continue 
		else {	// Libération mémoire
			mpz_clear(t); mpz_clear(i); mpz_clear(a); mpz_clear(y); mpz_clear(j); mpz_clear(s); mpz_clear(n_1);
			gmp_randclear(state);
			return 0;
			}
		}
	}

	// Libération mémoire
	mpz_clear(t); mpz_clear(i); mpz_clear(a); mpz_clear(y); mpz_clear(j); mpz_clear(s);
	mpz_clear(n_1);gmp_randclear(state);
	return 1;
	}
}


void square_multiply(mpz_t r, mpz_t a, mpz_t n, mpz_t h){ // r : resultat de la fonction, n : modulo, h : exposant 
	mpz_set(r,a);	//r = a

	char* str = mpz_get_str(NULL, 2, h); // On converti en base 2 l'exposant h
	long long int i, len = (long long int)strlen(str); // On stocke dans len la taille de la chaine de caractères contenant l'exposant

	for (i = 1; i < len; i++){
		mpz_mul(r, r, r); // r = r*r
		mpz_mod(r, r, n); // r = r*r modulo n
		if ( str[i] == '1'){ // Bit à 1
			mpz_mul(r, r, a); // r = r*a
			mpz_mod(r, r, n); // r = r*a modulo n
		}
	}
}

int main(){ 
	
	int res = 0;
	mpz_t n, k;
	mpz_init(n);
	mpz_init(k);
	int x = 1;

	while(x){
		printf("\nMerci de choisir un nombre n à tester : ");
		gmp_scanf("%Zd", &n);
		gmp_printf("Vous avez choisi n : %Zd \n", &n);
		printf("\nMerci de choisir k le nombre de répétition souhaité : ");
		gmp_scanf("%Zd", &k);
		gmp_printf("Vous avez choisi k : %Zd \n", &k);
    
		// Test de la fonction de Fermat:
		printf("\n### TEST DE FERMAT ###\n");	
		res = testFermat(n,k);
		if ( res == 1) gmp_printf("Le nombre %Zd est premier \n \n",n);
		else gmp_printf("Le nombre %Zd est composé \n \n",n);

		// Test de la fonction de miller rabin:
		printf("\n### TEST DE MILLER RABIN ###\n");
		res = testMillerRabin(n,k);
		if ( res == 1) gmp_printf("Le nombre %Zd est premier \n \n",n);
		else gmp_printf("Le nombre %Zd est composé \n \n",n);
		
		printf("Voulez vous continuer (Oui): 1 (Non): 0 : ");
		scanf("%d", &x);
	}
	

	return 0;
}

