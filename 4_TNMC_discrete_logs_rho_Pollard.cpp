#include "stdafx.h"

#define DEBUG
//#define RELEASE

clock_t t_1 = 0;

//gmp_randstate_t g_rstate;

mpz_t g_one, g_two;

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

mpz_t a, b, p, q;

mpz_t g_temp;

void my_func(mpz_t num)
{
	if (mpz_cmp(num, q) < 0)
	{
		mpz_mul(g_temp, a, num);
		mpz_mod(g_temp, g_temp, p);
		return; // return (a * num) % p
	}

	mpz_mul(g_temp, b, num);
	mpz_mod(g_temp, g_temp, p); // return (b * num) % p
}

// Элемент a имеет порядок q по модулю p. Найти дискретный логарифм x - такое целое число 1<x<q, что a^x = b(mod p)
void Rho_method_Pollard(std::string file, const char * ch_a, const char * ch_b, const char * ch_p, const char * ch_q)
{
	clock_t t;
	unsigned long long i = 0;

#ifdef DEBUG
	std::ofstream fout(file);
	fout << "a = " << ch_a << std::endl;
	fout << "b = " << ch_b << std::endl;
	fout << "p = " << ch_p << std::endl;
	fout << "q = " << ch_q << std::endl;
	fout << "Discrete logs rho method Pollard" << std::endl;
	fout << "i\tc\tlog_a_(c)\td\tlog_a_(d)" << std::endl;

	std::cout << "Discrete logs rho method Pollard" << std::endl;
#endif


	mpz_init_set_str(a, ch_a, 10);
	mpz_init_set_str(b, ch_b, 10);
	mpz_init_set_str(p, ch_p, 10);
	mpz_init_set_str(q, ch_q, 10);

	mpz_t u, v; // u = v = 2
	mpz_init_set_str(u, "2", 10);
	mpz_init_set_str(v, "2", 10);

	mpz_t log_a_c, log_a_d; // log_a_c = log_a_d = 2
	mpz_init_set_str(log_a_c, "2", 10);
	mpz_init_set_str(log_a_d, "2", 10);

	mpz_t count_log_a_c, cout_log_a_d; // count_log_a_c = cout_log_a_d = 2
	mpz_init_set_str(count_log_a_c, "2", 10);
	mpz_init_set_str(cout_log_a_d, "2", 10);

	mpz_t temp1, temp2;
	mpz_init(temp1);
	mpz_init(temp2);

	mpz_t c, d;
	mpz_init(c);
	mpz_init(d);

	mpz_t j;
	mpz_init_set_str(j, "1", 10);

	t = clock();

	mpz_powm(temp1, a, u, p); /// temp1 = a^u
	mpz_powm(temp2, b, v, p); /// temp2 = b^v
	mpz_mul(c, temp1, temp2); /// c = temp1 * temp2
	mpz_mod(c, c, p); // c = (a^u * b^v) (mod p)

	mpz_set(d, c); // d = c

	while (1)
	{
		i++;

		if (mpz_cmp(c, q) < 0)
		{
			mpz_add(log_a_c, log_a_c, g_one); // log_a_c += 1
		}
		else
		{
			mpz_add(count_log_a_c, count_log_a_c, g_one); // count_log_a_c += 1
		}

		my_func(c);
		mpz_set(c, g_temp); // c = my_func(c)

		if (mpz_cmp(d, q) < 0)
		{
			mpz_add(log_a_d, log_a_d, g_one); // log_a_d += 1
		}
		else
		{
			mpz_add(cout_log_a_d, cout_log_a_d, g_one); // cout_log_a_d += 1
		}

		my_func(d);
		mpz_set(d, g_temp); // d = my_func(d)

		if (mpz_cmp(d, q) < 0)
		{
			mpz_add(log_a_d, log_a_d, g_one); // log_a_d += 1
		}
		else
		{
			mpz_add(cout_log_a_d, cout_log_a_d, g_one); // cout_log_a_d += 1
		}
		my_func(d);
		mpz_set(d, g_temp); // d = my_func(d)


#ifdef DEBUG 
		fout << i << "\t" << c << "\t" << log_a_c << " + " << count_log_a_c << "*x\t" << d << "\t" << log_a_d << " + " << cout_log_a_d << "*x" << std::endl;
		// std::cout << i << "\t" << c << "\t" << log_a_c << " + " << count_log_a_c << "*x\t" << d << "\t" << log_a_d << " + " << cout_log_a_d << "*x" << std::endl;
#endif

		if (mpz_cmp(c, d) == 0) // c == d
		{
			break;
		}
	}


	while (mpz_cmp(j, q) <= 0) // while (j <= q)
	{
		mpz_mul(temp1, count_log_a_c, j); mpz_add(temp1, temp1, log_a_c); mpz_mod(temp1, temp1, q); // (log_a_c + count_log_a_c * j) % q

		mpz_mul(temp2, cout_log_a_d, j); mpz_add(temp2, temp2, log_a_d); mpz_mod(temp2, temp2, q); // (log_a_d + cout_log_a_d * j) % q

		if (mpz_cmp(temp1, temp2) == 0)
		{
			break;
		}

		mpz_add(j, j, g_one);
	}

	t = clock() - t;
	t_1 += t;

	mpz_mod(temp1, j, q); // temp1 = j % q

#ifdef DEBUG
	fout << "x = " << temp1 << " (mod " << q << ")" << std::endl;
	std::cout << "x = " << temp1 << " (mod " << q << ")" << std::endl;
	PrintTime(fout, t);
#endif
}

void for_task(char * n_file, const char * a, const char * b, const char * p, const char * q)
{
	Rho_method_Pollard(n_file, a, b, p, q);

	std::cout << a << std::endl;
	std::cout << "Time t_1: " << t_1 << " clicks (" << ((float)t_1) / CLOCKS_PER_SEC << " seconds)" << std::endl << std::endl;
}

int main(int argc, char * argv[])
{
//	gmp_randinit_default(g_rstate);
//	gmp_randseed_ui(g_rstate, time(0)); // for rand()

	mpz_init_set_str(g_one, "1", 10);
	mpz_init_set_str(g_two, "2", 10);
	mpz_init(g_temp);

	std::cout << "\"Format: name_file\ta\tb\tp\tq\"" << std::endl;

	t_1 = 0;
	if (argc >= 5)
		for_task(argv[1], argv[2], argv[3], argv[4], argv[5]);

	// for_task("test.txt", "10", "64","107","53");
	//	for_task("test.txt", "1359331");
		
	//for_task("1.txt",
		//"6966346018918884971");
	//	"17481286671975058943");


//	getchar();
	return 0;
}