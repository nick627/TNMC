#include "stdafx.h"

#define DEBUG
//#define RELEASE

clock_t t_1 = 0, t_2 = 0, t_3 = 0;

gmp_randstate_t rstate;

void PrintTime(std::ofstream& fout, clock_t t)
{
#ifdef DEBUG
	fout << "Time: " << t << " clicks (" << ((float)t) / CLOCKS_PER_SEC << " seconds)" << std::endl << std::endl;
	fout << std::endl << std::endl;
	fout.close();

	std::cout << "Time: " << t << " clicks (" << ((float)t) / CLOCKS_PER_SEC << " seconds)" << std::endl << std::endl;
	std::cout << std::endl << std::endl;
#endif
}

void PrintProbablePrime(std::ofstream& fout, unsigned long long * i, mpz_t * a)
{
#ifdef DEBUG
	fout << *i << "\t" << *a << "\tProbable prime" << std::endl;
	std::cout << *i << "\t" << *a << "\tProbable prime" << std::endl;
#endif
}

void PrintComposite(std::ofstream& fout, unsigned long long * i, mpz_t * a)
{
#ifdef DEBUG
	fout << *i << "\t" << *a << "\tCompisite" << std::endl;
	std::cout << *i << "\t" << *a << "\tCompisite" << std::endl;
#endif
}

void Ferma(std::string file, const char * ch_n, const char * ch_k) // Расширенный алгоритм Евклида
{
	clock_t t;
	unsigned long long i = 1;

#ifdef DEBUG
	std::ofstream fout(file, std::ofstream::app);
	fout << "n = " << ch_n << std::endl;
	fout << "Ferma" << std::endl;
	fout << "i\tBase\tResult" << std::endl;

	std::cout << "n = " << ch_n << std::endl;
	std::cout << "Ferma" << std::endl;
#endif

	mpz_t unit;
	mpz_init_set_str(unit, "1", 10);

	mpz_t temp;
	mpz_init_set_str(temp, "2", 10);

	mpz_t n;
	mpz_init_set_str(n, ch_n, 10);

	mpz_t k;
	mpz_init_set_str(k, ch_k, 10);

	mpz_t a;
	mpz_init(a);

	mpz_t n_1;
	mpz_init(n_1);
	mpz_sub(n_1, n, unit); // n_1 = n - 1;

	if (mpz_cmp(n, temp) == 0) // if(n == 2)
	{
#ifdef DEBUG
		fout << "Prime" << std::endl;
		std::cout << "Prime" << std::endl;
#endif
		return; // return 1;
	}

	t = clock();

	while (mpz_sgn(k) == 1) // while(k > 0)
	{
		mpz_urandomm(a, rstate, n); // a = rand();

		mpz_powm(temp, a, n_1, n); // temp = a^(n - 1) (mod n)

		mpz_sub(temp, temp, unit); /// temp -= 1;
		if (mpz_sgn(temp) == 1) // if(temp > 1)
		{
#ifdef DEBUG
			PrintComposite(fout, &i, &a);
#endif

			mpz_sub(k, k, unit); // k -= 1
			i++;
			continue; // break; // return 0;
		}

		mpz_sub(k, k, unit); // k -= 1
		i++;
#ifdef DEBUG
		PrintProbablePrime(fout, &i, &a);
#endif

	}

	t = clock() - t;
	t_1 += t;
#ifdef DEBUG
	PrintTime(fout, t);
#endif
}

void Solovey_Shtrassen(std::string file, const char * ch_n, const char * ch_k)
{
	clock_t t;
	unsigned long long i = 1;

#ifdef DEBUG
	std::ofstream fout(file, std::ofstream::app);
	fout << "n = " << ch_n << std::endl;
	fout << "Solovey- Shtrassen" << std::endl;
	fout << "i\tBase\tResult" << std::endl;

	std::cout << "n = " << ch_n << std::endl;
	std::cout << "Solovey- Shtrassen" << std::endl;
#endif

	mpz_t unit;
	mpz_init_set_str(unit, "1", 10);

	mpz_t temp;
	mpz_init(temp);

	mpz_t n;
	mpz_init_set_str(n, ch_n, 10);

	mpz_t k;
	mpz_init_set_str(k, ch_k, 10);

	mpz_t a;
	mpz_init(a);

	mpz_t n_1_div_2;
	mpz_init(n_1_div_2);
	mpz_sub(n_1_div_2, n, unit);
	mpz_tdiv_q_2exp(n_1_div_2, n_1_div_2, 1); // n_1_div_2 = (n - 1) / 2;

	mpz_t r;
	mpz_init(r);

	signed long int s_signed_long_int;
	mpz_t s_mpz_t; // mpz_jacobi
	mpz_init(s_mpz_t);

	t = clock();

	while (mpz_sgn(k) == 1) // while(k > 0)
	{
		mpz_urandomm(a, rstate, n); // a = rand();

		mpz_gcd(temp, a, n);
		mpz_sub(temp, temp, unit); /// temp -= 1;
		if (mpz_sgn(temp) == 1) // if(gcd(a, n) > 1)
		{
#ifdef DEBUG
			PrintComposite(fout, &i, &a);
#endif

			mpz_sub(k, k, unit); // k -= 1
			i++;
			continue; // break; // return 0;
		}

		mpz_powm(r, a, n_1_div_2, n); // r = a^((n - 1) / 2) (mod n)

		s_signed_long_int = mpz_jacobi(a, n);
		mpz_set_si(s_mpz_t, s_signed_long_int); // s = jac(a, n)
		mpz_mod(s_mpz_t, s_mpz_t, n); // jac(a, n) (mod n)

		if (mpz_cmp(r, s_mpz_t) != 0) // if(r != (s % n))
		{
#ifdef DEBUG
			PrintComposite(fout, &i, &a);
#endif

			mpz_sub(k, k, unit); // k -= 1
			i++;
			continue; // break; // return 0;
		}

		mpz_sub(k, k, unit); // k -= 1
		i++;
#ifdef DEBUG
		PrintProbablePrime(fout, &i, &a);
#endif
	}

	t = clock() - t;
	t_2 += t;
#ifdef DEBUG
	PrintTime(fout, t);
#endif
}

int check_right_bit(mpz_t * res, mpz_t * bit, mpz_t * a)
{
	mpz_and(*res, *a, *bit); // res = a & bit

	return (mpz_sgn(*res) == 0) ? 1 : 0;
}

void Rabin_Miller(std::string file, const char * ch_n, const char * ch_k)
{
	clock_t time_clock_t;
	unsigned long long i = 1;

#ifdef DEBUG
	std::ofstream fout(file, std::ofstream::app);
	fout << "n = " << ch_n << std::endl;
	fout << "Rabin- Miller" << std::endl;
	fout << "i\tBase\tResult" << std::endl;

	std::cout << "n = " << ch_n << std::endl;
	std::cout << "Rabin- Miller" << std::endl;
#endif

	mpz_t unit;
	mpz_init_set_str(unit, "1", 10);

	mpz_t temp;
	mpz_init(temp);
	mpz_t x;
	mpz_init(x);

	mpz_t n;
	mpz_init_set_str(n, ch_n, 10);

	mpz_t k;
	mpz_init_set_str(k, ch_k, 10);

	mpz_t a;
	mpz_init(a);

	mpz_t t;
	mpz_init(t);
	mpz_sub(t, n, unit); // t = n - 1

	long int s = 0;

	mpz_t bit_1, for_check_res;
	mpz_init_set_str(bit_1, "1", 10); // for (n & 1)
	mpz_init(for_check_res);

	mpz_t two;
	mpz_init_set_str(two, "2", 10);

	time_clock_t = clock();

	// Представляем (n - 1) в виде (2^s * t), где t нечетно
	while (check_right_bit(&for_check_res, &bit_1, &t)) // while ( (t & 1) == 0)
	{
		mpz_tdiv_q_2exp(t, t, 1); // t /= 2
		s++;
	}

	while (mpz_sgn(k) == 1) // while(k > 0)
	{
		mpz_urandomm(a, rstate, n); // a = rand();

		mpz_powm(x, a, t, n); // x = a^t (mod n)

		mpz_set(temp, x); // temp = x
		mpz_sub(temp, temp, unit); // temp -= 1;
		if (mpz_sgn(temp) == 0) // if (x == 1)
		{
#ifdef DEBUG
			PrintProbablePrime(fout, &i, &a);
#endif

			mpz_sub(k, k, unit); // k -= 1
			i++;
			continue;
		}

		mpz_set(temp, x); // temp = x
		mpz_sub(temp, temp, n); mpz_add(temp, temp, unit); // temp -= (n - 1)
		if (mpz_sgn(temp) == 0) // if (x == (n - 1))
		{
#ifdef DEBUG
			PrintProbablePrime(fout, &i, &a);
#endif

			mpz_sub(k, k, unit); // k -= 1
			i++;
			continue;
		}

		long int j = 0;

		while (j < (s - 1))
		{
			mpz_powm(x, x, two, n); // x = x^2 (mod n)

			mpz_set(temp, x); // tenp = x
			mpz_sub(temp, temp, unit); // temp -= 1
			if (mpz_sgn(temp) == 0) // if(x == 1)
			{
				// PrintComposite(fout, &i, &a);
				break; // return; // return 0;
			}

			mpz_set(temp, x); // temp = x
			mpz_sub(temp, temp, n); mpz_add(temp, temp, unit); // temp -= (n - 1)
			if (mpz_sgn(temp) == 0) // if(x == (n - 1))
			{
				break; // return to while(k > 0)
			}

			j++;
		}

		mpz_set(temp, x); // temp = x
		mpz_sub(temp, temp, n); mpz_add(temp, temp, unit); // temp -= (n - 1)
		if (mpz_sgn(temp) != 0) // if(x != (n - 1))
		{
#ifdef DEBUG
			PrintComposite(fout, &i, &a);
#endif

			mpz_sub(k, k, unit); // k -= 1
			i++;
			continue; // break; // return 0;
		}

		mpz_sub(k, k, unit); // k -= 1
		i++;
#ifdef DEBUG
		PrintProbablePrime(fout, &i, &a);
#endif
	}

	time_clock_t = clock() - time_clock_t;
	t_3 += time_clock_t;
#ifdef DEBUG
	PrintTime(fout, time_clock_t);
#endif
}

void for_task(char *file, const char *a, const char *b)
{
#ifdef RELEASE
	for (int i = 0; i < 10000; i++)
#endif
	{
		Ferma(file, a, b);
		Solovey_Shtrassen(file, a, b);
		Rabin_Miller(file, a, b);
	}
	std::cout << a << std::endl;
	std::cout << "Time t_1: " << t_1 << " clicks (" << ((float)t_1) / CLOCKS_PER_SEC << " seconds)" << std::endl << std::endl;
	std::cout << "Time t_2: " << t_2 << " clicks (" << ((float)t_2) / CLOCKS_PER_SEC << " seconds)" << std::endl << std::endl;
	std::cout << "Time t_3: " << t_3 << " clicks (" << ((float)t_3) / CLOCKS_PER_SEC << " seconds)" << std::endl << std::endl;
}

int main(int argc, char *argv[])
{
	gmp_randinit_default(rstate);
	gmp_randseed_ui(rstate, time(0));


	for_task("5.txt", "561", "5");
	
	t_1 = t_2 = t_3 = 0;
	for_task("1.txt",
		"15006150061500615007",
		"5");
	
	t_1 = t_2 = t_3 = 0;
	for_task("2.txt",
		"2074863650801173054457570649619278292063",
		"5");

	t_1 = t_2 = t_3 = 0; 
	for_task("3.txt", 
		"264852565359818589896673892356334055903",
		"5");

	t_1 = t_2 = t_3 = 0;
	for_task("4.txt",
		"17819352767119642433001928618050874258632267419750360825444427478626936230345477",
		"5");

	if (argc >= 4)
	{
		for_task(argv[1], argv[2], argv[3]);
	}

	getchar();
	return 0;
}