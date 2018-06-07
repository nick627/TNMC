#include "stdafx.h"

#define DEBUG
//#define RELEASE

void PrintTime(std::ofstream& fout, clock_t t)
{
	fout << "Time: " << t << " clicks (" << ((float)t) / CLOCKS_PER_SEC << " seconds)" << std::endl << std::endl;
	fout << std::endl << std::endl;
	fout.close();

	std::cout << "Time: " << t << " clicks (" << ((float)t) / CLOCKS_PER_SEC << " seconds)" << std::endl << std::endl;
	std::cout << std::endl << std::endl;
}

//-----------------matrix_t-------------------------------------------------
//------------------------------------------------------------------------------
typedef struct
{
	mpz_t *MATRIX;
	mpz_t *IDENTITY;
	uint64_t rows;
	uint64_t cols;
	uint64_t next_free_row;	// used to insert rows in the matrix, point to the next free row to be inserted position
} matrix_t;

// allocates space for m rows * n columns matrix for MATRIX and IDENTITY 
void init_matrix(matrix_t * matrix, uint64_t m, uint64_t n)
{
	matrix->MATRIX = (mpz_t *)calloc(m, sizeof(mpz_t));
	matrix->IDENTITY = (mpz_t *)calloc(m, sizeof(mpz_t));
	matrix->rows = m;
	matrix->cols = n;
	matrix->next_free_row = 0;
}

void push_row(matrix_t * matrix, mpz_t row)
{
	mpz_set(matrix->MATRIX[matrix->next_free_row], row);

	mpz_init2(matrix->IDENTITY[matrix->next_free_row], matrix->cols); // initializes a n bit vector all set to 0
	mpz_setbit(matrix->IDENTITY[matrix->next_free_row], matrix->next_free_row); // set the next_free_row bit to 1
	matrix->next_free_row++;
}

void print_matrix_matrix(std::ofstream& fout, matrix_t * matrix)
{
	uint64_t i, j;

	printf("\nMATRIX\n");
	for (i = 0; i < matrix->rows; i++)
	{
		printf("[");
		for (j = 0; j < matrix->cols; j++)
		{
			printf("%2d", mpz_tstbit(matrix->MATRIX[i], j));
		}
		printf(" ]\n");
	}
	printf("\n");
}
void print_matrix_identity(std::ofstream& fout, matrix_t * matrix)
{
	uint64_t i, j;

	fout << std::endl << "IDENTITY" << std::endl;
	// std::cout << std::endl << "IDENTITY" << std::endl;
	for (i = 0; i < matrix->rows; i++)
	{
		fout << "[ ";
		// std::cout << "[ ";
		for (j = 0; j < matrix->cols; j++)
		{
			fout << mpz_tstbit(matrix->IDENTITY[i], j);
			// std::cout << mpz_tstbit(matrix->IDENTITY[i], j);
		}
		fout << " ]" << std::endl;
		// std::cout << " ]" << std::endl;
	}
	
	fout << std::endl;
	// std::cout << std::endl;
}

void free_matrix(matrix_t *matrix)
{
	free(matrix->MATRIX);
	free(matrix->IDENTITY);
}

#define MIN(a,b) ((a) < (b) ? (a) : (b))

// performs a Gauss elimination on matrix->MATRIX, result (linear dependence) will be in the matrix->IDENTITY
void gauss_elimination(matrix_t * matrix)
{
	// printf("\nPerforming Gauss elimination..\n");

	mpz_t *m = matrix->MATRIX;
	mpz_t *I = matrix->IDENTITY;

	uint64_t col, row, next_row, next_pivot;
	for (next_row = 0, col = 0; col < MIN(matrix->cols, matrix->rows); col++) // for all rows
	{
		next_pivot = -1;
		for (row = next_row; row < matrix->rows; row++) // search for the next pivot
		{
			if (mpz_tstbit(m[row], col))
			{
				next_pivot = row; // row contains the next pivot
				next_row++;
				break;
			}
		}

		if (next_pivot == -1)
			continue;

		if (next_pivot != next_row - 1) // current row is not the pivot, switch rows
		{
			mpz_swap(m[next_pivot], m[next_row - 1]);
			mpz_swap(I[next_pivot], I[next_row - 1]);
		}

		for (row = next_row; row < matrix->rows; row++)
		{
			if (mpz_tstbit(m[row], col))
			{
				mpz_xor(m[row], m[row], m[next_row - 1]); // XOR the rows to eliminate the 1 in position (row, next_row-1)
				mpz_xor(I[row], I[row], I[next_row - 1]);
			}
		}
	}
}

// does not check for bounds, the caller must
void get_matrix_row(mpz_t rop, matrix_t * matrix, uint64_t row_index)
{
	mpz_set(rop, matrix->MATRIX[row_index]);
}

void get_identity_row(mpz_t rop, matrix_t * matrix, uint64_t row_index)
{
	mpz_set(rop, matrix->IDENTITY[row_index]);
}
//------------------------------------------------------------------------------
//-----------------End matrix_t-------------------------------------------------


#define GET_BIT_AT(index) ((numbers[index>>3] & (1<<(index&7))) >> (index&7))
#define SET_BIT_AT(index) (numbers[index>>3] |= (1 << (index&7)))
#define CLEAR_BIT_AT(index) (numbers[index>>3] &= ~(1 << (index&7)))

char * numbers;
int64_t base_ref;

// Can handle up to base 2^31 * 8 * 8
int64_t sieve_primes_up_to(int64_t base)
{
	int64_t num_primes = 0;
	int64_t i;

	// Find primes, base included
	base++;

	numbers = (char *)calloc(base / 64 + 1, sizeof(uint64_t));
	base_ref = base;

	for (i = 0; i < base; i++)
		SET_BIT_AT(i);

	int64_t p = 2;
	num_primes++;
	int64_t offset;

	while (1)
	{
		offset = 2 * p;

		while (offset < base)
		{
			CLEAR_BIT_AT(offset);
			offset += p;
		}

		offset = p + 1;
		while (offset < base && (GET_BIT_AT(offset) == 0))
		{
			offset++;
		}

		if (offset == base)
			break;

		p = offset;
		num_primes++;
	}

	return num_primes;
}

// Fill the array with only primes where n is a quadratic residue: x^2 = n (mod p)
int fill_primes_with_quadratic_residue(uint64_t * primes_array, mpz_t n)
 {
	int64_t j, i;
	mpz_t b;
	mpz_init(b);

	primes_array[0] = 2;
	for (j = 1, i = 3; i < base_ref; i++)
	{
		mpz_set_ui(b, (unsigned long)i);
		if ((GET_BIT_AT(i)) == 1 && mpz_jacobi(n, b) == 1)
		{
			primes_array[j] = i;
			j++;
		}
	}

	free(numbers);
	return j;
}

// Solve the modular equatioon x^2 = n (mod p) using the Shanks-Tonelli
// algorithm. x will be placed in q and 1 returned if the algorithm is
// successful. Otherwise 0 is returned (currently in case n is not a quadratic
// residue mod p). A check is done if p = 3 (mod 4), in which case the root is
// calculated as n ^ ((p+1) / 4) (mod p).
//
// Note that currently mpz_legendre is called to make sure that n really is a
// quadratic residue. The check can be skipped, at the price of going into an
// eternal loop if called with a non-residue.
int mpz_sqrtm(mpz_t q, const mpz_t n, const mpz_t p)
{
	mpz_t w, n_inv, y;
	unsigned int i, s;
	//TMP_DECL;
	//TMP_MARK;

	if (mpz_divisible_p(n, p)) // Is n a multiple of p?
	{       
		mpz_set_ui(q, 0); // Yes, then the square root is 0.
		return 1;         // (special case, not caught
	}                     // otherwise)
						  // if(mpz_legendre(n, p) != 1) // Not a quadratic residue? 
						  // return 0; // No, so return error       
	if (mpz_tstbit(p, 1) == 1) // p = 3 (mod 4) ?
	{
		mpz_set(q, p);
		mpz_add_ui(q, q, 1);
		mpz_fdiv_q_2exp(q, q, 2);
		mpz_powm(q, n, q, p); // q = n ^ ((p+1) / 4) (mod p)
		return 1;
	}

	//MPZ_TMP_INIT(y, 2*SIZ(p));
	//MPZ_TMP_INIT(w, 2*SIZ(p));
	//MPZ_TMP_INIT(n_inv, 2*SIZ(p));

	mpz_init(y);
	mpz_init(w);
	mpz_init(n_inv);

	mpz_set(q, p);
	mpz_sub_ui(q, q, 1);                // q = p-1                          
	s = 0;                              // Factor out 2^s from q            
	while (mpz_tstbit(q, s) == 0) s++;	//
	mpz_fdiv_q_2exp(q, q, s);           // q = q / 2^s                      
	mpz_set_ui(w, 2);                   // Search for a non-residue mod p   
	while (mpz_legendre(w, p) != -1)    // by picking the first w such that
		mpz_add_ui(w, w, 1);            // (w/p) is -1                      
	mpz_powm(w, w, q, p);               // w = w^q (mod p)                  
	mpz_add_ui(q, q, 1);				//
	mpz_fdiv_q_2exp(q, q, 1);           // q = (q+1) / 2                    
	mpz_powm(q, n, q, p);               // q = n^q (mod p)                  
	mpz_invert(n_inv, n, p);			//
	for (;;) {							//
		mpz_powm_ui(y, q, 2, p);        // y = q^2 (mod p)                  
		mpz_mul(y, y, n_inv);
		mpz_mod(y, y, p);               // y = y * n^-1 (mod p) 
		i = 0;
		while (mpz_cmp_ui(y, 1) != 0)
		{
			i++;
			mpz_powm_ui(y, y, 2, p);    //  y = y ^ 2 (mod p)
		}
		if (i == 0) 
		{                    /* q^2 * n^-1 = 1 (mod p), return   */
										 //TMP_FREE;
			return 1;
		}
		if (s - i == 1)			// In case the exponent to w is 1,
		{					    // Don't bother exponentiating   
			mpz_mul(q, q, w);           
		}
		else 
		{
			mpz_powm_ui(y, w, 1 << (s - i - 1), p);
			mpz_mul(q, q, y);
		}
		mpz_mod(q, q, p);               // r = r * w^(2^(s-i-1)) (mod p)
	}

	mpz_clear(w); mpz_clear(n_inv); mpz_clear(y);
	return 0;
}

typedef struct modular_root
{
	unsigned long root1;
	unsigned long root2;
} modular_root_t;

typedef struct
{
	mpz_t value_x;
	mpz_t value_x_squared;
	// uint64_t *factors_exp;      
	// this can be used to keep track for the full factorization of
	// an element x?-n. For the use of this version, refer to the
	// qs_exponent_vector.c version

	mpz_t factors_vect;
	// this nb_primes_in_base bit vector is a substitute of the
	// factors_exp array, all we actually need for the elimination
	// in the Linear algebra step is the exponents modulo 2.
	// This saves huge space within the sieving step
} smooth_number_t;

mpz_t N; // number to factorize
matrix_t matrix;

uint64_t nb_smooth_numbers_found = 0;
uint64_t nb_qr_primes; // number of primes p where N is a quadratic residue mod p
uint64_t *primes; // array holding the primes of the smoothness base

smooth_number_t *smooth_numbers;
int NB_VECTORS_OFFSET = 5; // number of additional rows in the matrix, to make sure that a linear relation exists 

//-----------------------------------------------------------
// base <- exp((1/2) sqrt(ln(n) ln(ln(n))))
//-----------------------------------------------------------
void get_smoothness_base(mpz_t base, mpz_t n)
{
	mpfr_t fN, lnN, lnlnN;
	mpfr_init(fN), mpfr_init(lnN), mpfr_init(lnlnN);

	mpfr_set_z(fN, n, MPFR_RNDU);
	mpfr_log(lnN, fN, MPFR_RNDU);
	mpfr_log(lnlnN, lnN, MPFR_RNDU);

	mpfr_mul(fN, lnN, lnlnN, MPFR_RNDU);
	mpfr_sqrt(fN, fN, MPFR_RNDU);
	mpfr_div_ui(fN, fN, 2, MPFR_RNDU);
	mpfr_exp(fN, fN, MPFR_RNDU);

	mpfr_get_z(base, fN, MPFR_RNDU);

	mpfr_clears(fN, lnN, lnlnN, NULL);
}

// returns the index of the first element to start sieving from
// (first multiple of root that is directly greater than start)
// res = p * t + root >= start
void get_sieving_start_index(mpz_t res, mpz_t start, mpz_t p, unsigned long root) 
{
	mpz_t q, r;
	mpz_init(q);
	mpz_init(r);

	mpz_sub_ui(start, start, root);
	mpz_fdiv_qr(q, r, start, p);

	if (mpz_cmp_ui(r, 0) != 0)
		mpz_add_ui(q, q, 1);

	mpz_mul(q, q, p); // next element p*q+root that is directly >= start
	mpz_add_ui(q, q, root);
	mpz_set(res, q);
	mpz_clear(q);
	mpz_clear(r);
}

// given an array of exponents of factors on the base [primes], reconstructs the number
void reconstruct_mpz(mpz_t rop, uint64_t *factors_exp)
{
	uint64_t i;
	mpz_t t;
	mpz_init(t);
	mpz_t p_exp;
	mpz_init(p_exp);

	mpz_set_ui(t, 1);
	for (i = 0; i < nb_qr_primes; i++)
	{
		mpz_set_ui(p_exp, primes[i]);
		mpz_pow_ui(p_exp, p_exp, factors_exp[i]);
		mpz_mul(t, t, p_exp);
	}

	mpz_set(rop, t);
}

// save the sooth number n to the smooth_numbers array, and at the same time its exponents vector to the matrix
mpz_t tmp_matrix_row;
void save_smooth_number(smooth_number_t n)
{
	mpz_clear(tmp_matrix_row); // tmp_matrix_row must be initialized already 
	mpz_init2(tmp_matrix_row, nb_qr_primes); // init a vector of *exactly* nb_qr_primes bits

	if (nb_smooth_numbers_found > nb_qr_primes + NB_VECTORS_OFFSET - 1) // if we have sufficient smooth numbers, skip saving
		return;

	smooth_number_t tmp;
	mpz_init(tmp.value_x);
	mpz_init(tmp.value_x_squared);

	mpz_set(tmp.value_x, n.value_x);
	mpz_pow_ui(tmp.value_x_squared, n.value_x, 2);
	mpz_sub(tmp.value_x_squared, tmp.value_x_squared, N); 
	// x^2 -N. Saving this will enable us to not go through exponents
	// and reconstruct the original number

	// otherwise we can reconstruct value_x_squared from the exponents vector, this is useful in the factoring step
	// to calculate the square modulo N from the factors directly.
	// It takes a lot of space, doesn't woth it anyways
	
	// tmp.factors_exp = calloc(nb_qr_primes, sizeof(uint64_t));
	// memcpy(tmp.factors_exp, n.factors_exp, nb_qr_primes * sizeof(uint64_t));
	// reconstruct_mpz(tmp.value_x_squared, tmp.factors_exp);
	smooth_numbers[nb_smooth_numbers_found++] = tmp;

	// reconstruct and saves the smooth number to the GF2 matrix ***//* already done in the sieving step
	 /* 
	 uint64_t i;
	 for(i=0; i<nb_qr_primes; i++)
	 {
	 if(n.factors_exp[i]&1)
	 mpz_setbit(tmp_matrix_row, i);
	 }
	 //*/
	
	// the coefficient vector in GF2 has already been constructed
	push_row(&matrix, n.factors_vect);
}


int main(int argc, char **argv)
{
	std::ofstream fout(argv[1]);

	clock_t my_time;

	mpz_init(N);
	mpz_t B;
	mpz_init(B);

	unsigned long int uBase;
	int64_t nb_primes;
	modular_root_t * modular_roots;

	uint64_t i, j;

	if (mpz_init_set_str(N, argv[2], 10) == -1)
	{
		fout << "Cannot load N "<< argv[2] << std::endl;
		std::cout << "Cannot load N " << argv[2] << std::endl;
		exit(2);
	}

	mpz_t sqrtN, rem;
	mpz_init(sqrtN);
	mpz_init(rem);
	mpz_sqrtrem(sqrtN, rem, N);
	
	if (mpz_cmp_ui(rem, 0) != 0) // if not perfect square, calculate the ceiling 
	{
		mpz_add_ui(sqrtN, sqrtN, 1);
	}
	else // N is a perfect square, factored!
	{
		fout << std::endl << "<<<[FACTOR]>>> " << mpz_get_str(NULL, 10, sqrtN) << std::endl;
		std::cout << std::endl << "<<<[FACTOR]>>> " << mpz_get_str(NULL, 10, sqrtN) << std::endl;
		return 0;
	}

	if (mpz_probab_prime_p(N, 10) > 0) // don't bother factoring
	{
		fout << "N: " << mpz_get_str(NULL, 10, N) << " is prime" << std::endl;
		std::cout << "N: " << mpz_get_str(NULL, 10, N) << " is prime" << std::endl;
		exit(0);
	}


	//--------------------------------------------------------
	//  calculate the smoothness base for the given N
	//--------------------------------------------------------

	get_smoothness_base(B, N); // if N is too small, the program will surely fail, please consider a pen and paper instead
	uBase = mpz_get_ui(B);

	fout << "n: " << mpz_get_str(NULL, 10, N) << " \tBase ( exp((1/2) sqrt(ln(n) ln(ln(n)))) ): " << mpz_get_str(NULL, 10, B) << std::endl;
	std::cout << "n: " << mpz_get_str(NULL, 10, N) << " \tBase ( exp((1/2) sqrt(ln(n) ln(ln(n)))) ): " << mpz_get_str(NULL, 10, B) << std::endl;

	//--------------------------------------------------------
	// sieve primes that are less than the smoothness base using Eratosthenes sieve
	//--------------------------------------------------------
	
	my_time = clock();

	nb_primes = sieve_primes_up_to((int64_t)(uBase));

#ifdef DEBUG
	fout << std::endl << "Primes found " << nb_primes << " [Smoothness Base " << uBase << "]" << std::endl;
	std::cout << std::endl << "Primes found " << nb_primes << " [Smoothness Base " << uBase << "]" << std::endl;
#endif

	//STOP_TIMER_PRINT_TIME("\tEratosthenes Sieving done");

	//--------------------------------------------------------
	// fill the primes array with primes to which n is a quadratic residue
	//--------------------------------------------------------

	//START_TIMER();

	primes = (uint64_t *)calloc(nb_primes, sizeof(int64_t));
	nb_qr_primes = fill_primes_with_quadratic_residue(primes, N);

#ifdef DEBUG
	for (i = 0; i < nb_qr_primes; i++)
	{
		fout << primes[i] << ", ";
		std::cout<< primes[i] << ", ";
	}
	fout << std::endl << std::endl << "N-Quadratic primes found " << i << std::endl;
	std::cout << std::endl << std::endl << "N-Quadratic primes found " << i << std::endl;
#endif

	//STOP_TIMER_PRINT_TIME("\tQuadratic prime filtering done");
	
	//--------------------------------------------------------
	// calculate modular roots
	//--------------------------------------------------------
	
	//START_TIMER();
	
	modular_roots = (modular_root_t *)calloc(nb_qr_primes, sizeof(modular_root_t));
	mpz_t tmp, r1, r2;
	mpz_init(tmp);
	mpz_init(r1);
	mpz_init(r2);

	for (i = 0; i < nb_qr_primes; i++)
	{
		mpz_set_ui(tmp, (unsigned long)primes[i]);
		mpz_sqrtm(r1, N, tmp); // calculate the modular root
		mpz_neg(r2, r1); // -q mod n
		mpz_mod(r2, r2, tmp);

		modular_roots[i].root1 = mpz_get_ui(r1);
		modular_roots[i].root2 = mpz_get_ui(r2);
	}
	mpz_clear(tmp);
	mpz_clear(r1);
	mpz_clear(r2);

	//STOP_TIMER_PRINT_TIME("\nModular roots calculation done");

#ifdef DEBUG
	for (i = 0; i < nb_qr_primes; i++)
	{
		fout << "p_" << i << " [" << primes[i] << " -> roots: " << modular_roots[i].root1 << " - " << modular_roots[i].root2 << "]" << std::endl;
		// std::cout << "p_" << i << " [" << primes[i] << " -> roots: " << modular_roots[i].root1 << " - " << modular_roots[i].root2 << "]" << std::endl;
	}
#endif


	//--------------------------------------------------------
	//         ***** initialize the matrix *****
	//--------------------------------------------------------

	//START_TIMER();
	init_matrix(&matrix, nb_qr_primes + NB_VECTORS_OFFSET, nb_qr_primes);
	mpz_init2(tmp_matrix_row, nb_qr_primes);

	//STOP_TIMER_PRINT_TIME("\nMatrix initialized");

	//--------------------------------------------------------
	// [Sieving]
	//--------------------------------------------------------

	//START_TIMER();

	mpz_t x, sieving_index, next_sieving_index;
	unsigned long ui_index, SIEVING_STEP = 50000; // we sieve for 50000 elements at each loop
	uint64_t p_pow;
	smooth_number_t *x_squared;

	x_squared = (smooth_number_t *)calloc(SIEVING_STEP, sizeof(smooth_number_t));
	smooth_numbers = (smooth_number_t *)calloc(nb_qr_primes + NB_VECTORS_OFFSET,
		sizeof(smooth_number_t));

	mpz_init_set(x, sqrtN);
	mpz_init_set(sieving_index, x);
	mpz_init_set(next_sieving_index, x);

	mpz_t p;
	mpz_init(p);
	mpz_t str;
	mpz_init_set(str, sieving_index);

#ifdef DEBUG
	fout << std::endl << "Sieving..." << std::endl << std::endl;
	std::cout << std::endl << "Sieving..." << std::endl << std::endl;
#endif
	
	//--------------------------------------------------------
	// Init before sieving
	//--------------------------------------------------------
	
	for (i = 0; i < SIEVING_STEP; i++)
	{
		mpz_init(x_squared[i].value_x);
		mpz_init(x_squared[i].value_x_squared);

		// the factors_exp array is used to keep track of exponents
		// x_squared[i].factors_exp = calloc(nb_qr_primes, sizeof(uint64_t));
		// we use directly the exponents vector modulo 2 to preserve space
		mpz_init2(x_squared[i].factors_vect, nb_qr_primes);
		mpz_add_ui(x, x, 1);
	}

	int nb_smooth_per_round = 0;
	char s[512];

	//--------------------------------------------------------
	// WHILE smooth numbers found less than the primes in the smooth base + NB_VECTORS_OFFSET
	//--------------------------------------------------------

	while (nb_smooth_numbers_found < nb_qr_primes + NB_VECTORS_OFFSET)
	{
		nb_smooth_per_round = 0;
		mpz_set(x, next_sieving_index); // sieve numbers from sieving_index to sieving_index + sieving_step
		mpz_set(sieving_index, next_sieving_index);

#ifdef DEBUG
		fout << "Sieving at: " << mpz_get_str(NULL, 10, sieving_index) << " <--> Smooth numbers found: " << nb_smooth_numbers_found << " " << nb_qr_primes << std::endl;
		// std::cout << "Sieving at: " << mpz_get_str(NULL, 10, sieving_index) << " <--> Smooth numbers found: " << nb_smooth_numbers_found << " " << nb_qr_primes << std::endl;
#endif 

		for (i = 0; i < SIEVING_STEP; i++)
		{
			mpz_set(x_squared[i].value_x, x);

			mpz_pow_ui(x_squared[i].value_x_squared, x, 2); // calculate value_x_squared <- x^2 -n 
			mpz_sub(x_squared[i].value_x_squared, x_squared[i].value_x_squared,
				N);

			mpz_clear(x_squared[i].factors_vect);
			mpz_init2(x_squared[i].factors_vect, nb_qr_primes); // reconstruct a new fresh 0ed vector of size nb_qr_primes bits

			mpz_add_ui(x, x, 1);
		}
		mpz_set(next_sieving_index, x);

		//--------------------------------------------------------
		// eliminate factors in the x_squared array, those who are 'destructed' to 1 are smooth
		//--------------------------------------------------------

		for (i = 0; i < nb_qr_primes; i++)
		{
			mpz_set_ui(p, (unsigned long)primes[i]);
			mpz_set(x, sieving_index);

			// get the first multiple of p that is directly larger that sieving_index
			// Quadratic SIEVING: all elements from this number and in positions multiples of root1 and root2
			// are also multiples of p
			get_sieving_start_index(x, x, p, modular_roots[i].root1);
			mpz_set(str, x);
			mpz_sub(x, x, sieving_index); // x contains index of first number that is divisible by p

			for (j = mpz_get_ui(x); j < SIEVING_STEP; j += primes[i])
			{
				p_pow = mpz_remove(x_squared[j].value_x_squared,
					x_squared[j].value_x_squared, p); // eliminate all factors of p

				if (p_pow & 1) // mark bit if odd power of p exists in this x_squared[j]
				{
					mpz_setbit(x_squared[j].factors_vect, i);
				}

				if (mpz_cmp_ui(x_squared[j].value_x_squared, 1) == 0)
				{
					save_smooth_number(x_squared[j]);
					nb_smooth_per_round++;
				}
				// sieve next element located p steps from here
			}

			// same goes for root2
			if (modular_roots[i].root2 == modular_roots[i].root1)
				continue;

			mpz_set(x, sieving_index);

			get_sieving_start_index(x, x, p, modular_roots[i].root2);
			mpz_set(str, x);
			mpz_sub(x, x, sieving_index);

			for (j = mpz_get_ui(x); j < SIEVING_STEP; j += primes[i])
			{
				p_pow = mpz_remove(x_squared[j].value_x_squared,
					x_squared[j].value_x_squared, p);

				if (p_pow & 1)
				{
					mpz_setbit(x_squared[j].factors_vect, i);
				}

				if (mpz_cmp_ui(x_squared[j].value_x_squared, 1) == 0)
				{
					save_smooth_number(x_squared[j]);
					nb_smooth_per_round++;
				}
			}
		}

#ifdef DEBUG
		fout << "Smooth numbers found " << nb_smooth_numbers_found << std::endl;
		// std::cout << "Smooth numbers found " << nb_smooth_numbers_found << std::endl;
#endif

		/*
		sprintf(s, "[start: %s - end: %s - step: %" PRId64 "] nb_smooth_per_round: %d",
		mpz_get_str(NULL, 10, sieving_index),
		mpz_get_str(NULL, 10, next_sieving_index),
		SIEVING_STEP,
		nb_smooth_per_round);
		APPEND_TO_LOG_FILE(s);
		//*/
	}

	//STOP_TIMER_PRINT_TIME("\nSieving DONE");

	uint64_t t = 0;

	//--------------------------------------------------------
	//the matrix ready, start Gauss elimination. The Matrix is filled on the call of save_smooth_number()
	//--------------------------------------------------------

	//START_TIMER();

	gauss_elimination(&matrix);

	//STOP_TIMER_PRINT_TIME("\nGauss elimination done");
	
#ifdef DEBUG
	//print_matrix_matrix(fout, &matrix);
	print_matrix_identity(fout, &matrix);
#endif

	uint64_t row_index = nb_qr_primes + NB_VECTORS_OFFSET - 1; // last row in the matrix
	int nb_linear_relations = 0;
	mpz_t linear_relation_z, solution_z;
	mpz_init(linear_relation_z);
	mpz_init(solution_z);

	get_matrix_row(linear_relation_z, &matrix, row_index--); // get the last few rows in the Gauss eliminated matrix
	while (mpz_cmp_ui(linear_relation_z, 0) == 0)
	{
		nb_linear_relations++;
		get_matrix_row(linear_relation_z, &matrix, row_index--);
	}

#ifdef DEBUG
	fout << "Linear dependent relations found : " << nb_linear_relations << std::endl;
	std::cout << "Linear dependent relations found : " << nb_linear_relations << std::endl;
#endif

	//--------------------------------------------------------
	// Factor
	//--------------------------------------------------------

	//We use the last linear relation to reconstruct our solution

	//START_TIMER();

#ifdef DEBUG
	fout << std::endl << "Factorizing..." << std::endl;
	std::cout << std::endl << "Factorizing..." << std::endl;
		
	fout << "s\tt" << std::endl;
	// std::cout << "s\tt" << std::endl;
#endif

	mpz_t solution_X, solution_Y;
	mpz_init(solution_X);
	mpz_init(solution_Y);

	// we start testing from the first linear relation encountered in the matrix
	for (j = nb_linear_relations; j > 0; j--)
	{
#ifdef DEBUG
		fout << nb_linear_relations - j + 1;
		// std::cout << "Trying " << nb_linear_relations - j + 1 << "..." << std::endl;
#endif
		mpz_set_ui(solution_X, 1);
		mpz_set_ui(solution_Y, 1);

		get_identity_row(solution_z, &matrix,
			nb_qr_primes + NB_VECTORS_OFFSET - j + 1);

		for (i = 0; i < nb_qr_primes; i++)
		{
			/*
#ifdef DEBUG
			fout << "vectx: " << solution_z << std::endl;
			// std::cout << "vectx: " << solution_z << std::endl;
#endif
			//*/

			if (mpz_tstbit(solution_z, i))
			{
				/*
#ifdef DEBUG
			fout << "vect: "<<smooth_numbers[i].value_x << std::endl;
			// std::cout << "vect: " << smooth_numbers[i].value_x << std::endl;
#endif
			//*/


				mpz_mul(solution_X, solution_X, smooth_numbers[i].value_x);
				mpz_mod(solution_X, solution_X, N); // reduce x to modulo N

				mpz_mul(solution_Y, solution_Y,
					smooth_numbers[i].value_x_squared);
				// TODO: handling huge stuff here, there is no modulo N like in the solution_X case!
				// eliminate squares as long as you go
			}
		}

		mpz_sqrt(solution_Y, solution_Y);
		mpz_mod(solution_Y, solution_Y, N); // y = sqrt(MUL(xi^2 - n)) mod N

#ifdef DEBUG
		fout << "\t" << solution_X << "\t" << solution_Y << std::endl;
		// std::cout << solution_X << "\t" << solution_Y << std::endl;
#endif

		mpz_sub(solution_X, solution_X, solution_Y);

		mpz_gcd(solution_X, solution_X, N);

		if (mpz_cmp(solution_X, N) != 0 && mpz_cmp_ui(solution_X, 1) != 0) // factor can be 1 or N, try another relation
			break;
	}
	mpz_cdiv_q(solution_Y, N, solution_X);


	my_time = clock() - my_time;

	fout << std::endl << ">>>>>> FACTORED " << mpz_get_str(NULL, 10, N) << " = "
		<< mpz_get_str(NULL, 10, solution_X) << " * " << mpz_get_str(NULL, 10, solution_Y) << std::endl;
	std::cout << std::endl << ">>>>>> FACTORED " << mpz_get_str(NULL, 10, N) << " = "
		<< mpz_get_str(NULL, 10, solution_X) << " * " << mpz_get_str(NULL, 10, solution_Y) << std::endl;
	
	PrintTime(fout, my_time);

	//STOP_TIMER_PRINT_TIME("\nFactorizing done");

	printf("\nCleaning memory..\n");

	//----------------------clear the x_squared array----------------------
	for (i = 0; i < SIEVING_STEP; i++)
	{
		mpz_clear(x_squared[i].value_x);
		mpz_clear(x_squared[i].value_x_squared);
		//free(x_squared[i].factors_exp);
		mpz_clear(x_squared[i].factors_vect);
	}
	free(x_squared);
	//---------------------clear the x_squared array ----------------------

	free(modular_roots);
	//***************** clear the smooth_numbers array -------------------
	for (i = 0; i < nb_qr_primes + NB_VECTORS_OFFSET; i++)
	{
		mpz_clear(smooth_numbers[i].value_x);
		mpz_clear(smooth_numbers[i].value_x_squared);
		//free(smooth_numbers[i].factors_exp);
	}
	free(smooth_numbers);
	//-------------------clear the smooth_numbers array--------------------

	free(primes);
	//------------------------clear mpz _t--------------------------
	mpz_clear(B);
	mpz_clear(N);
	sqrtN, rem;
	mpz_clear(x);
	mpz_clear(sieving_index);
	mpz_clear(next_sieving_index);
	mpz_clear(p);
	mpz_clear(str);
	//-------------------------- clear mpz _t ---------------------------

	free_matrix(&matrix);

	// getchar();

	return 0;
}