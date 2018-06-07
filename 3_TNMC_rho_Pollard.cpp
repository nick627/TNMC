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

void Rho_method_Pollard(std::string file, const char * ch_n)
{
	clock_t t;
	unsigned long long i = 0;

	mpz_t n;
	mpz_init_set_str(n, ch_n, 10);

	mpz_t x, y; // x = y = 2
	mpz_init_set_str(x, "2", 10);
	mpz_init_set_str(y, "2", 10);

	mpz_t d; // d = 1;
	mpz_init_set_str(d, "2", 10);

	mpz_t c; // c = rand();
	mpz_init(c);
	mpz_urandomm(c, g_rstate, n);

#ifdef DEBUG
	std::ofstream fout(file);
	fout << "n = " << ch_n << std::endl;
	fout << "Rho_method_Pollard" << std::endl;
	fout << "f(x) = x^2 + " << c << std::endl;
	fout << "i\ta\tb\tgcd(a - b, n)" << std::endl;

	std::cout << "n = " << ch_n << std::endl;
	std::cout << "f(x) = x^2 + " << c << std::endl;
	std::cout << "Rho_method_Pollard" << std::endl;
#endif

	mpz_t temp;
	mpz_init_set_str(temp, "2", 10);

	t = clock();

	while (1) // while(d == 1)
	{
		i++;

		mpz_powm(x, x, g_two, n); /// x = x^2 (mod n)
		mpz_add(x, x, c); /// x += c
		mpz_add(x, x, n); /// x += n
		mpz_mod(x, x, n); // x(i + 1) = f(x)

		mpz_powm(y, y, g_two, n); /// y = x^2 (mod n)
		mpz_add(y, y, c); /// y += c
		mpz_add(y, y, n); /// y += n
		mpz_mod(y, y, n); /// y(i + 1) = f(y)

		mpz_powm(y, y, g_two, n); /// y = x^2 (mod n)
		mpz_add(y, y, c); /// y += c
		mpz_add(y, y, n); /// y += n
		mpz_mod(y, y, n); // y(i + 1) = f(f(y))

		mpz_sub(temp, x, y); /// temp = x - y
		mpz_gcd(d, temp, n); // d = gcd(( x - y), n)

#ifdef DEBUG
							 // std::cout << i << "\t" << x << "\t" << y << "\t" << d << std::endl;
		fout << i << "\t" << x << "\t" << y << "\t" << d << std::endl;
#endif

		if (mpz_cmp(d, g_one) != 0) // if (d != 1)
		{
			break;
		}
	}

	t = clock() - t;
	t_1 += t;


	mpz_cdiv_q(temp, n, d); // temp = n / d

	std::cout << std::endl << n << " = " << d << " * " << temp << std::endl;

#ifdef DEBUG
	fout << std::endl << n << " = " << d << " * " << temp << std::endl;
	PrintTime(fout, t);
#endif
}

void for_task(char * n_file, const char * a)
{
	Rho_method_Pollard(n_file, a);

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
		
	//for_task("1.txt",
		//"6966346018918884971");
	//	"17481286671975058943");


//	getchar();
	return 0;
}