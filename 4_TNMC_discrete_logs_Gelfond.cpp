#include "stdafx.h"

#define DEBUG
//#define RELEASE

clock_t t_1 = 0;

//gmp_randstate_t g_rstate;

mpz_t g_zero, g_one, g_two;

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

void func_with_egcd(mpz_t b, mpz_t n)
{
	mpz_t g, x, y;
	mpz_init(g);
	mpz_init(x);
	mpz_init(y);

	// egcd
	mpz_t a_egcd, b_egcd;
	mpz_init_set(a_egcd, b);
	mpz_init_set(b_egcd, n);

	mpz_t x2, x1, y2, y1;
	mpz_init_set_str(x1, "1", 10); // x2 = 1
	mpz_init_set_str(x2, "0", 10); // x1 = 0
	mpz_init_set_str(y1, "0", 10); // y2 = 0
	mpz_init_set_str(y2, "1", 10); // y1 = 1

	mpz_t q_egcd, r_egcd;
	mpz_init(q_egcd);
	mpz_init(r_egcd);

	while (mpz_cmp(b_egcd, g_zero) != 0)
	{
		mpz_cdiv_qr(q_egcd, r_egcd, a_egcd, b_egcd);

		//std::cout << "q = " << q_egcd << ", r = " << r_egcd << std::endl;
		mpz_set(a_egcd, b_egcd); mpz_set(b_egcd, r_egcd); // a = b; b = r

		mpz_mul(g_temp, q_egcd, x2); mpz_sub(x, x1, g_temp); // x = x2 - q * x1
		mpz_mul(g_temp, q_egcd, y2); mpz_sub(y, y1, g_temp); // y = y2 - q * y1

		mpz_set(x1, x2); // x2 = x1
		mpz_set(x2, x); // x1 = x
		mpz_set(y1, y2); // y2 = y1
		mpz_set(y2, y); // y1 = y

		//std::cout << "x0 = " << x1 << " y0 = " << y1 << std::endl;
		//std::cout << "x1 = " << x2 << " y1 = " << y2 << std::endl<<std::endl;
	}

	mpz_set(g, a_egcd);
	mpz_set(x, x1); mpz_mul_si(x, x, -1);
	mpz_set(y, y1); mpz_mul_si(y, y, -1);
	// end egcd

	mpz_set(g_temp, g);

	if (mpz_cmpabs(g, g_one) == 0)
	{
		mpz_mod(g_temp, x, n);

		// std::cout << g_temp;
	}
}

// Элемент a имеет порядок q по модулю p. Найти дискретный логарифм x - такое целое число 1<x<q, что a^x = b(mod p)
void p_1_discrete_logs_Gelfond(std::string file, const char * ch_a, const char * ch_b, const char * ch_p, const char * ch_q)
{
	clock_t t;
	unsigned long long i = 0;

#ifdef DEBUG
	std::ofstream fout(file);
	fout << "a = " << ch_a << std::endl;
	fout << "b = " << ch_b << std::endl;
	fout << "p = " << ch_p << std::endl;
	fout << "q = " << ch_q << std::endl;
	fout << "Discrete logs method Gelfond" << std::endl;
	fout << "k\t(b * a^(-k s) (mod p))" << std::endl;

	std::cout << "Discrete logs method Gelfond" << std::endl;
#endif
	
	mpz_init_set_str(a, ch_a, 10);
	mpz_init_set_str(b, ch_b, 10);
	mpz_init_set_str(p, ch_p, 10);
	mpz_init_set_str(q, ch_q, 10);

	mpz_t u, v;
	mpz_init(u); mpz_init(v);

	mpz_t temp1, temp2;
	mpz_init(temp1);
	mpz_init(temp2);

	mpz_t s;
	mpz_init(s); 
	
	mpz_sqrt(s, q); mpz_add(s, s, g_one); // s = [sqrt(q)] + 1

	/*
	mpz_t _s;
	mpz_init(_s);
	mpz_powm(_s, a, s, p);
	func_with_egcd(_s, p);
	mpz_set(g_temp, _s);
	std::cout << "_s = " << _s << std::endl; // _s = 1157715000490578439
	//*/

	mpir_ui size_ba_arr = mpz_get_ui(s);

	mpz_t * ba_arr = (mpz_t *)malloc(sizeof(mpz_t) * size_ba_arr); // b, ba^(-1*s), ba^(-2*s), ..., ba^-(s-1)s (mod p)
	for (mpir_ui i = 0; i < size_ba_arr; i++)	mpz_init(ba_arr[i]);
	
	mpir_si my_i = 0;


	t = clock();

	while (mpz_cmp_si(s, my_i) != 0) // while (i < s)
	{
		mpz_sub(temp1, q, s); mpz_mul_si(temp1, temp1, my_i); 
		mpz_powm(temp1, a, temp1, p); // temp = a ^ (i * (q - s)) (mod p)
		mpz_mul(temp1, b, temp1); mpz_mod(temp1, temp1, p); // temp1 = b * temp1 (mod p)

		mpz_set(ba_arr[my_i], temp1);

#ifdef DEBUG
		fout << my_i << "\t" << ba_arr[my_i] << std::endl;
		// std::cout << my_i << "\t" << ba_arr[my_i] << std::endl;
#endif

		my_i++;
	}

	fout << std::endl;

	my_i = 0;

	while (mpz_cmp_si(s, my_i) != 0) // while (i < s)
	{
		mpz_powm_ui(temp2, a, my_i, p);
		// mpz_set(a_arr[my_i], temp2);

#ifdef DEBUG
		fout << "a^" << my_i << " = " << temp2 << std::endl;
		// std::cout << "a^" << my_i << " = " << temp2 << std::endl;
#endif

		for (mpir_si k = 0; k < size_ba_arr; k++)
		{
			// if (mpz_cmp(a_arr[my_i], ba_arr[k]) == 0)
			if (mpz_cmp(temp2, ba_arr[k]) == 0)
			{
				t = clock() - t;
				
				mpz_set_si(u, k);
				mpz_set_si(v, my_i);

				mpz_mul(temp1, u, s); mpz_add(temp1, temp1, v); mpz_mod(temp1, temp1, q);
#ifdef DEBUG
				// fout << std::endl << "a^" << my_i << " = " << a_arr[my_i] << "; (" << ba_arr[k] << ", " << k << ")" << std::endl;
				fout << std::endl << "a^" << my_i << " = " << temp2 << "; (" << ba_arr[k] << ", " << k << ")" << std::endl;
				// std::cout << std::endl << "a^" << my_i << " = " << a_arr[my_i] << "; (" << ba_arr[k] << ", " << k << ")" << std::endl;
				std::cout << std::endl << "a^" << my_i << " = " << temp2 << "; (" << ba_arr[k] << ", " << k << ")" << std::endl;

				fout << "x = " << temp1 << std::endl;
				fout << "x = " << u << " * " << s << " + " << v << " (mod " << q << ")" << std::endl;
				std::cout << "x = " << temp1 << std::endl;
				std::cout << "x = " << u << " * " << s << " + " << v << " (mod " << q << ")" << std::endl;
				PrintTime(fout, t);
#endif

				free(ba_arr);
				return;
			}
		}

		my_i++;
	}
}

void for_task(char * n_file, const char * a, const char * b, const char * p, const char * q)
{
	p_1_discrete_logs_Gelfond(n_file, a, b, p, q);
}

int main(int argc, char * argv[])
{
//	gmp_randinit_default(g_rstate);
//	gmp_randseed_ui(g_rstate, time(0)); // for rand()

	mpz_init_set_str(g_zero, "0", 10);
	mpz_init_set_str(g_one, "1", 10);
	mpz_init_set_str(g_two, "2", 10);
	mpz_init(g_temp);

	std::cout << "\"Format: name_file\ta\tb\tp\tq\"" << std::endl;

	t_1 = 0;
	if (argc >= 5)
		;//	for_task(argv[1], argv[2], argv[3], argv[4], argv[5]);

	for_task("test.txt", "7", "167","587","293");
	//	for_task("test.txt", "1359331");
		
	//for_task("1.txt",
		//"6966346018918884971");
	//	"17481286671975058943");


//	getchar();
	return 0;
}
