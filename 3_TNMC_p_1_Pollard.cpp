#include "stdafx.h"

#define DEBUG
//#define RELEASE

#define _MPFR_RNDN MPFR_RNDU

clock_t t_1 = 0, t_2 = 0, t_3 = 0;

gmp_randstate_t g_rstate;

mpz_t g_one;
mpz_t g_two;

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

void end_func(std::ofstream& fout, clock_t t, mpz_t temp, mpz_t n, mpz_t d, mpz_t a)
{

	mpz_cdiv_q(temp, n, d);
	std::cout << std::endl << n << " = " << d << " * " << temp << std::endl;

#ifdef DEBUG
	fout << std::endl << "a = " << a << std::endl;

	fout << std::endl << n << " = " << d << " * " << temp << std::endl;
	PrintTime(fout, t);
#endif
}

void p_1_method_Pollard(std::string file, const char * ch_n)
{	// 215 Machovenko
	clock_t t;
	unsigned long long i = 0;

	//matrix_t matrix;
	mpfr_set_default_prec(100);


#ifdef DEBUG
	std::ofstream fout(file);
	fout << "n = " << ch_n << std::endl;
	fout << "p_1_method_Pollard" << std::endl;
	fout << "i\ts\ta\tl_i" << std::endl;

	std::cout << "n = " << ch_n << std::endl;
	std::cout << "p_1_method_Pollard" << std::endl;
#endif

	mpz_t n;
	mpz_init(n);
	mpz_set_str(n, ch_n, 10);

	std::cout << n << std::endl;

	mpfr_t ln_n;
	mpfr_init(ln_n);
	mpfr_set_z(ln_n, n, _MPFR_RNDN); // (mpfr_t)ln_n = (mpz_t)n
	mpfr_log(ln_n, ln_n, _MPFR_RNDN); // ln_n = ln(n)
	// mpfr_printf("ln(n) = %Rf\n", ln_n);

	mpz_t p_i;
	mpz_init(p_i);

	mpz_t l;
	mpz_init(l);

	mpz_t d; // d = 1;
	mpz_init_set_str(d, "1", 10);

	mpz_t a; // a = rand();
	mpz_init(a);

	mpz_t temp;
	mpz_init_set_str(temp, "2", 10);

	mpfr_t mpfrRes;
	mpfr_init(mpfrRes);

	long long s = 11;
	long long s_i = 0;

	mpz_t B[11];
	for (s_i = 0; s_i < sizeof(B) / sizeof(B[0]); s_i++)
	{
		mpz_init(B[i]);
	}
	mpz_init_set_str(B[0], "2", 10);
	mpz_init_set_str(B[1], "3", 10);
	mpz_init_set_str(B[2], "5", 10);
	mpz_init_set_str(B[3], "7", 10);
	mpz_init_set_str(B[4], "11", 10);
	mpz_init_set_str(B[5], "13", 10);
	mpz_init_set_str(B[6], "17", 10);
	mpz_init_set_str(B[7], "19", 10);
	mpz_init_set_str(B[8], "23", 10);
	mpz_init_set_str(B[9], "29", 10);
	mpz_init_set_str(B[10], "31", 10);

	/*
	mpz_set_str(B[0], "2", 10);
	mpz_set_str(B[1], "5", 10);
	mpz_set_str(B[2], "23", 10);
	mpz_set_str(B[3], "441751", 10);
	mpz_set_str(B[4], "68564555489", 10);
	//*/

	//mpz_urandomm(a, g_rstate, n); // a = rand()

	t = clock();

	while (1) // while(d == 1)
	{
		i++;

		mpz_init_set_str(p_i, "100000", 10); // 10^5 < p_i < 10^6

		mpz_urandomm(a, g_rstate, n); // a = rand()

		mpz_gcd(d, a, n); // d = gcd(a, n)

		if (mpz_cmp(d, g_two) >= 0)
		{
			t = clock() - t;
			t_1 += t;
			end_func(fout, t, temp, n, d, a);
			break;
		}

		// for (s_i = 0; s_i < 5; s_i++)
		// for (s_i = 0; s_i < sizeof(B)/ sizeof(B[0]); s_i++)
		for (s_i = 0; s_i < s; s_i++)
		{
			//mpz_set(p_i, B[s_i]); // p_i = B[s_i]
			mpz_next_prime_candidate(p_i, p_i, g_rstate); // p_i <--> prime

			mpfr_set_z(mpfrRes, p_i, _MPFR_RNDN); // (mpfr_t)mpfrRes = (mpz_t)p_i
			mpfr_log(mpfrRes, mpfrRes, _MPFR_RNDN); // mpfrRes = ln(mpfrRes) // ln(p_i)
			mpfr_div(mpfrRes, ln_n, mpfrRes, _MPFR_RNDN); // mpfrRes = (ln(n) / ln(p_i))
			// mpfr_printf("%Rf\n", mpfrRes);
			mpfr_get_z(l, mpfrRes, _MPFR_RNDN); // (mpz_t)l = (mpfr_t)[(ln(n) / ln(p_i))]

			mpir_ui mpir_l = mpz_get_ui(l);
			mpz_t mpz_t_l;
			mpz_init_set_ui(mpz_t_l, mpir_l);
#ifdef DEBUG
			// std::cout << i << "\t" << mpir_l << std::endl;
			fout << i << "\t" << s << "\t" << a << "\t" << mpir_l << std::endl;
#endif
			mpz_powm(temp, p_i, mpz_t_l, n); // p_i^l (mod n)
			mpz_powm(a, a, temp, n); // a = a^(p_i^l) (mod n)


			mpz_set(temp, a); /// temp = a
			mpz_sub(temp, temp, g_one); // temp = a - 1

			mpz_gcd(d, temp, n); // d = gcd((a - 1), n)

			if (mpz_cmp(d, g_one) != 0 && mpz_cmp(d, n) != 0)
			{
				// break;
				t = clock() - t;
				t_1 += t;

				mpfr_clear(mpfrRes); // free(mpfrRes)

				end_func(fout, t, temp, n, d, a);

				/*
				mpz_cdiv_q(temp, n, d);
				std::cout << std::endl << n << " = " << d << " * " << temp << std::endl;

#ifdef DEBUG
				fout << std::endl << "a = " << a << std::endl;

				fout << std::endl << n << " = " << d << " * " << temp << std::endl;
				PrintTime(fout, t);
#endif
				//*/

				return;
			}
			else
			{
				i++;
			}
			/*
			if (mpz_cmp(d, g_one) == 0)
			{
				// s++; // Base++
				continue;
			}
			else if (mpz_cmp(d, n) == 0)
			{
				// s--; // Base--
				continue;
			}
			else
			{
				break;
			}
			//*/
		}
		s++;
	}

}

void for_task(char * n_file, const char * a)
{
	p_1_method_Pollard(n_file, a);

	std::cout << a << std::endl;
	std::cout << "Time t_1: " << t_1 << " clicks (" << ((float)t_1) / CLOCKS_PER_SEC << " seconds)" << std::endl << std::endl;
}

int main(int argc, char * argv[])
{
	gmp_randinit_default(g_rstate);
	gmp_randseed_ui(g_rstate, time(0)); // for rand()

	mpz_init_set_str(g_one, "1", 10);
	mpz_init_set_str(g_two, "2", 10);

	std::cout << "\"File number\"" << std::endl;

	
	t_1 = t_2 = t_3 = 0;
	if (argc >= 3)
		for_task(argv[1], argv[2]);

//	for_task("test.txt", "1359331");
//	for_task("1.txt",
		//"6966346018918884971");
		//"17481286671975058943");

//	getchar();
	return 0;
}