#include "stdafx.h"

#define DEBUG
//#define RELEASE

clock_t t_1 = 0, t_2 = 0, t_3 = 0;

void first(std::string file, const char * ch_a, const char * ch_b) // Расширенный алгоритм Евклида
{
	clock_t t;

#ifdef DEBUG
	std::ofstream fout(file);
#endif

	mpz_t a, b, d, x, y;

	mpz_init_set_str(a, ch_a, 10);
	mpz_init_set_str(b, ch_b, 10);
	mpz_init(d);
	mpz_init(x);
	mpz_init(y);

#ifdef DEBUG
	std::cout << "Extended Evklid" << std::endl;
	std::cout << "a = " << a << std::endl;
	std::cout << "b = " << b << std::endl << std::endl;

	fout << "a = " << a << std::endl;
	fout << "b = " << b << std::endl;
#endif

	mpz_t q, r, x1, x2, y1, y2;

	mpz_init(q);
	mpz_init(r);
	mpz_init(x1);
	mpz_init(x2);
	mpz_init(y1);
	mpz_init(y2);

	unsigned long long i = 0;

	// add var
	mpz_t var;
	mpz_init(var);

	t = clock();

	if (mpz_sgn(b) == 0) // if (b == 0)
	{
		mpz_set(d, a); // d = a
		mpz_set_str(x, "1", 10); //	x = 1
		mpz_set_str(y, "0", 10); // y = 0
	}
	else
	{
		mpz_set_str(x2, "1", 10); // x2 = 1
		mpz_set_str(x1, "0", 10); // x1 = 0
		mpz_set_str(y2, "0", 10); // y2 = 0
		mpz_set_str(y1, "1", 10); // y1 = 1
	}

#ifdef DEBUG
	fout << "i\t q\t r\t x\t y\t a\t b\t x1\t x2\t y1\t y2" << std::endl;
#endif

	while (mpz_sgn(b) > 0) // while (b>0)
	{
		mpz_tdiv_q(q, a, b); // q = a / b
		mpz_mul(var, q, b);	mpz_sub(r, a, var); // r = a - q * b
		mpz_mul(var, q, x1); mpz_sub(x, x2, var); // x = x2 - q * x1
		mpz_mul(var, q, y1); mpz_sub(y, y2, var); // y = y2 - q * y1
		mpz_set(a, b); // a = b
		mpz_set(b, r); // b = r
		mpz_set(x2, x1); // x2 = x1
		mpz_set(x1, x); // x1 = x
		mpz_set(y2, y1); // y2 = y1
		mpz_set(y1, y); // y1 = y

#ifdef DEBUG
		fout << i << "\t" << q << "\t" << r << "\t" << x << "\t" << y << "\t" << a << "\t" << b << "\t" << x1 << "\t" << x2 << "\t" << y1 << "\t" << y2 << std::endl;
#endif

		i++;
	}

	mpz_set(d, a); // d = a
	mpz_set(x, x2); // x = x2
	mpz_set(y, y2); // y = y2

	t = clock() - t;
	t_1 += t;

#ifdef DEBUG
	std::cout << "x = " << x << std::endl;
	std::cout << "y = " << y << std::endl;
	std::cout << "d = " << d << std::endl;
	std::cout << "Time: " << t << " clicks (" << ((float)t) / CLOCKS_PER_SEC << " seconds)" << std::endl<<std::endl;

	fout << "--------------------------------------------------" << std::endl;
	fout << "x = " << x << std::endl;
	fout << "y = " << y << std::endl;
	fout << "d = " << d << std::endl;
	fout << "Time: " << t << " clicks (" << ((float)t) / CLOCKS_PER_SEC << " seconds)" << std::endl << std::endl;

	fout.close();
#endif
}

int check_right_bit(mpz_t * res, mpz_t * bit, mpz_t * a)
{
	mpz_and(*res, *a, *bit); // res = a & bit

	return (mpz_sgn(*res) == 0) ? 1 : 0;
}

void second(std::string file, const char * ch_a,const char * ch_b) // Расширенный бинарный алгоритм Евклида
{
	clock_t t;

#ifdef DEBUG
	std::ofstream fout(file);
#endif

	mpz_t for_check, for_check_res;
	mpz_init_set_str(for_check, "1", 10);
	mpz_init(for_check_res);

	mpz_t a, b, d, x, y;
	mpz_init_set_str(a, ch_a, 10);
	mpz_init_set_str(b, ch_b, 10);
	mpz_init(d);
	mpz_init(x);
	mpz_init(y);

#ifdef DEBUG
	std::cout << "Binary Evklid" << std::endl;
	std::cout << "a = " << a << std::endl;
	std::cout << "b = " << b << std::endl << std::endl;
		
	fout << "a = " << a << std::endl;
	fout << "b = " << b << std::endl;
#endif

	mpz_t g, u, v, A, B, C, D;

	mpz_init_set_str(g, "1", 10); // g = 1
	mpz_init(u);
	mpz_init(v);
	mpz_init(A);
	mpz_init(B);
	mpz_init(C);
	mpz_init(D);

	unsigned long long i = 0;

#ifdef DEBUG
	fout <<"i\t u\t v\t A\t B\t C\t D" << std::endl;
#endif

	t = clock();

	while (check_right_bit(&for_check_res, &for_check, &a) && check_right_bit(&for_check_res, &for_check, &b)) // while (a & 1 == 0 and b & 1 == 0)
	{
		mpz_tdiv_q_2exp(a, a, 1); // a = a / 2
		mpz_tdiv_q_2exp(b, b, 1); // b = b / 2
		mpz_mul_2exp(g, g, 1); // g = g * 2

		std::cout << "a = " << a << std::endl;
		std::cout << "g = " << g << std::endl;
	}

	mpz_set(u, a); // u = a
	mpz_set(v, b); // v = b
	mpz_set_str(A, "1", 10); // A = 1
	mpz_set_str(B, "0", 10); // B = 0
	mpz_set_str(C, "0", 10); // C = 0
	mpz_set_str(D, "1", 10); // D = 1

	while (mpz_sgn(u)) // while u != 0
	{
		while (check_right_bit(&for_check_res, &for_check, &u)) // while u % 2 == 0 :
		{
			mpz_tdiv_q_2exp(u, u, 1); // u = u / 2
			if (check_right_bit(&for_check_res, &for_check, &A) && check_right_bit(&for_check_res, &for_check, &B)) // if (A % 2 == 0 and B % 2 == 0)
			{
				mpz_tdiv_q_2exp(A, A, 1); // A = A / 2
				mpz_tdiv_q_2exp(B, B, 1); // B = B / 2
			}
			else // else
			{
				mpz_add(A, A, b); mpz_tdiv_q_2exp(A, A, 1);	// A = (A + b) / 2
				mpz_sub(B, B, a); mpz_tdiv_q_2exp(B, B, 1); // B = (B - a) / 2
			}
		}

		while (check_right_bit(&for_check_res, &for_check, &v)) // while v % 2 == 0
		{
			mpz_tdiv_q_2exp(v, v, 1); // v = v / 2
			if (check_right_bit(&for_check_res, &for_check, &C) && check_right_bit(&for_check_res, &for_check, &D)) // if (C % 2 == 0 and D % 2 == 0)
			{
				mpz_tdiv_q_2exp(C, C, 1); // C = C / 2
				mpz_tdiv_q_2exp(D, D, 1); // D = D / 2
			}
			else // else
			{
				mpz_add(C, C, b); mpz_tdiv_q_2exp(C, C, 1);	// C = (C + b) / 2
				mpz_sub(D, D, a); mpz_tdiv_q_2exp(D, D, 1); // D = (D - a) / 2
			}
		}

		if (mpz_cmp(u, v)>=0) // if u >= v
		{
			mpz_sub(u, u, v); // u = u - v
			mpz_sub(A, A, C); // A = A - C
			mpz_sub(B, B, D); // B = B - D
			
		}
		else // else
		{
			mpz_sub(v, v, u); // v = v - u
			mpz_sub(C, C, A); // C = C - A
			mpz_sub(D, D, B); // D = D - B
		}

#ifdef DEBUG
		fout << i << "\t" << u << "\t" << v << "\t" << A << "\t" << B << "\t" << C << "\t" << D << std::endl;
#endif

		mpz_mul(d, g, v); // d = g * v
		mpz_set(x, C); // x = C
		mpz_set(y, D); // y = D

		i++;
	}
	
	t = clock() - t;
	t_2 += t;

#ifdef DEBUG
	std::cout << "x = " << x << std::endl;
	std::cout << "y = " << y << std::endl;
	std::cout << "d = " << d << std::endl;
	std::cout << "Time: " << t << " clicks (" << ((float)t) / CLOCKS_PER_SEC << " seconds)" << std::endl << std::endl;

	fout << "--------------------------------------------------" << std::endl;
	fout << "x = " << x << std::endl;
	fout << "y = " << y << std::endl;
	fout << "d = " << d << std::endl;
	fout << "Time: " << t << " clicks (" << ((float)t) / CLOCKS_PER_SEC << " seconds)" << std::endl << std::endl;

	fout.close();
#endif
}

int check_thrid(mpz_t * r, mpz_t * b, mpz_t * add_var)
{
	mpz_tdiv_q_2exp(*add_var, *b, 1);
	mpz_sub(*add_var, *r, *add_var);
	// if((r - b / 2) == 0) // Correct  a * x + b * y = GCD
	// if((r - b / 2) > 0) // High speed (Верный алгоритм с "усеченными"), but  a * x + b * y <> GCD
	return (mpz_sgn(*add_var) > 0) ? 1 : 0;
}

void third(std::string file, const char * ch_a, const char * ch_b) // Расширенный алгоритм Евклида с "усечёнными" остатками
{
	clock_t t;

#ifdef DEBUG
	std::ofstream fout(file);
#endif

	mpz_t a, b, d, x, y;

	mpz_init_set_str(a, ch_a, 10);
	mpz_init_set_str(b, ch_b, 10);
	mpz_init(d);
	mpz_init(x);
	mpz_init(y);

#ifdef DEBUG
	std::cout << "Truncated Evklid" << std::endl;
	std::cout << "a = " << a << std::endl;
	std::cout << "b = " << b << std::endl << std::endl;

	fout << "a = " << a << std::endl;
	fout << "b = " << b << std::endl;
#endif

	mpz_t q, r, x1, x2, y1, y2;

	mpz_init(q);
	mpz_init(r);
	mpz_init(x1);
	mpz_init(x2);
	mpz_init(y1);
	mpz_init(y2);

	unsigned long long i = 0;

	// add var
	mpz_t var;
	mpz_init(var);

	t = clock();

#ifdef DEBUG
	fout << "i\t q\t r\t x\t y\t a\t b\t x1\t x2\t y1\t y2" << std::endl;
#endif

	if (mpz_sgn(b) == 0) // if (b == 0)
	{
		mpz_set(d, a); // d = a
		mpz_set_str(x, "1", 10); //	x = 1
		mpz_set_str(y, "0", 10); // y = 0
	}
	else
	{
		mpz_set_str(x2, "1", 10); // x2 = 1
		mpz_set_str(x1, "0", 10); // x1 = 0
		mpz_set_str(y2, "0", 10); // y2 = 0
		mpz_set_str(y1, "1", 10); // y1 = 1
	}

	while (mpz_sgn(b) > 0) // while (b>0)
	{
		mpz_tdiv_q(q, a, b); // q = a / b
		mpz_mul(var, q, b);	mpz_sub(r, a, var); // r = a - q * b
		mpz_mul(var, q, x1); mpz_sub(x, x2, var); // x = x2 - q * x1
		mpz_mul(var, q, y1); mpz_sub(y, y2, var); // y = y2 - q * y1

		if (check_thrid(&r, &b, &var)) // if((r - b / 2) > 0)
		{
			mpz_sub(r, b, r); // r = b - r
			mpz_sub(x, x1, x); // x = x1 - x
			mpz_sub(y, y1, y); // y = y1 - y
		}

		mpz_set(a, b); // a = b
		mpz_set(b, r); // b = r
		mpz_set(x2, x1); // x2 = x1
		mpz_set(x1, x); // x1 = x
		mpz_set(y2, y1); // y2 = y1
		mpz_set(y1, y); // y1 = y

#ifdef DEBUG
		fout << i << "\t" << q << "\t" << r << "\t" << x << "\t" << y << "\t" << a << "\t" << b << "\t" << x1 << "\t" << x2 << "\t" << y1 << "\t" << y2 << std::endl;
#endif

		i++;
	}

	mpz_set(d, a); // d = a
	mpz_set(x, x2); // x = x2
	mpz_set(y, y2); // y = y2

	t = clock() - t;
	t_3 += t;

#ifdef DEBUG
	mpz_t a_check, b_check, d_check;

	mpz_init_set_str(a_check, ch_a, 10);
	mpz_init_set_str(b_check, ch_b, 10);
	mpz_init(d_check);

	mpz_mul(a_check, a_check, x);
	mpz_mul(b_check, b_check, y);
	mpz_add(d_check, a_check, b_check);

	std::cout << "check: \n" << x << " * " << ch_a << " + " << y << " * " << ch_b << " =\n" << d_check << std::endl << std::endl;

	std::cout << "x = " << x << std::endl;
	std::cout << "y = " << y << std::endl;
	std::cout << "d = " << d << std::endl;
	std::cout << "Time: " << t << " clicks (" << ((float)t) / CLOCKS_PER_SEC << " seconds)" << std::endl << std::endl;

	fout << "--------------------------------------------------" << std::endl;
	fout << "x = " << x << std::endl;
	fout << "y = " << y << std::endl;
	fout << "d = " << d << std::endl;
	fout << "Time: " << t << " clicks (" << ((float)t) / CLOCKS_PER_SEC << " seconds)" << std::endl << std::endl;

	fout.close();
#endif
}

void for_task(char *file, const char *a, const char *b)
{
	std::string name_file;
	name_file = file;

#ifdef RELEASE
	for (int i = 0; i < 100000; i++)
#endif
	{

	first(name_file + "_1.txt", a, b);
	second(name_file + "_2.txt", a, b);
	third(name_file + "_3.txt", a, b);
	}
	std::cout << "Time t_1: " << t_1 << " clicks (" << ((float)t_1) / CLOCKS_PER_SEC << " seconds)" << std::endl << std::endl;
	std::cout << "Time t_2: " << t_2 << " clicks (" << ((float)t_2) / CLOCKS_PER_SEC << " seconds)" << std::endl << std::endl;
	std::cout << "Time t_3: " << t_3 << " clicks (" << ((float)t_3) / CLOCKS_PER_SEC << " seconds)" << std::endl << std::endl;

	mpz_t _a, _b, d;
	mpz_init_set_str(_a, a, 10);
	mpz_init_set_str(_b, b, 10);
	mpz_init(d);
	mpz_gcd(d, _a, _b);
	std::cout << "gcd: d = " << d << std::endl << std::endl;
}

int main(void)
{
	//for_task("4",	"17", "13");
	
	for_task("1",
		"9679032853442961649",
		"9679032648109410787");
	
	t_1 = t_2 = t_3 = 0;
	for_task("2",
		"177062133675607304176519220617962787321",
		"780496351258916797101051468810013169971");

	t_1 = t_2 = t_3 = 0;
	for_task("3", 
		"27698142912540626397833885903060828427690248156766927205949931751583144803164601",
		"7744912493415526406594555732735450703202025821016403443002080380116758393680983");

	getchar();
	return 0;
}