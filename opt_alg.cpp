#include"opt_alg.h"

solution MC(matrix(*ff)(matrix, matrix, matrix), int N, matrix lb, matrix ub, double epsilon, int Nmax, matrix ud1, matrix ud2)
{
	try
	{
		solution Xopt;
		while (true)
		{
			Xopt = rand_mat(N);
			for (int i = 0; i < N; ++i)
				Xopt.x(i) = (ub(i) - lb(i)) * Xopt.x(i) + lb(i);
			Xopt.fit_fun(ff, ud1, ud2);
			if (Xopt.y < epsilon)
			{
				Xopt.flag = 1;
				break;
			}
			if (solution::f_calls > Nmax)
			{
				Xopt.flag = 0;
				break;
			}
		}
		return Xopt;
	}
	catch (string ex_info)
	{
		throw ("solution MC(...):\n" + ex_info);
	}
}

double* expansion(matrix(*ff)(matrix, matrix, matrix), double x0, double d, double alpha, int Nmax, matrix ud1, matrix ud2)
{
	try
	{
		int i = 0;
		double ret = 0;
		double* p = new double[2]{ 0,0 };
		solution X0(x0),X1(x0+d);
		X0.fit_fun(ff, ud1, ud2);
		X1.fit_fun(ff, ud1, ud2);
		if(X0.y == X1.y) {
			p[0] = m2d(X0.x);
			p[1] = m2d(X1.x);
			return p;
		}
		if (X1.y > X0.y){
			d = -d;
			X1.x = X0.x + d;
			X1.fit_fun(ff, ud1, ud2);
			if (X1.y >= X0.y) {
				p[0] = m2d(X1.x);
				p[1] = m2d(X0.x)-d;
				return p;
			}
		}
		do {
			if (solution::f_calls > Nmax) {
				X0.flag = 0;
				return 0;
			}
			i = i + 1;
			ret = m2d(X0.x);
			X0 = X1;
			X1.x = x0 + d * pow(alpha,i);
			X1.fit_fun(ff, ud1, ud2);
		} while (X1.y <= X0.y);
		if (d > 0) {
			p[0] = ret;
			p[1] = m2d(X1.x);
		}
		else {
			p[1] = ret;
			p[0] = m2d(X1.x);
		}
		return p;
	}
	catch (string ex_info)
	{
		throw ("double* expansion(...):\n" + ex_info);
	}
}

solution fib(matrix(*ff)(matrix, matrix, matrix), double a, double b, double epsilon, matrix ud1, matrix ud2)
{
	try
	{
		solution Xopt,Xopt2;
		if (a > b) {
			double help = a;
			a = b;
			b = help;
		}
		double d,c;
		double k=0; 
		while (GetFib(k) < (b-a)/epsilon)
		{	
			k++;
		}
		cout << "a: " << a << " b: " << b << " k : " << k << endl;
		c = b - GetFib(k - 1) / GetFib(k) * (b - a);
		d = a + b - c;
		for (int i = 0; i < k - 3; i++) {
			Xopt.x = c;
			Xopt2.x = d;
			Xopt.fit_fun(ff, ud1, ud2);
			Xopt2.fit_fun(ff, ud1, ud2);
			if (Xopt.y < Xopt2.y) {
				b = d;
			}
			else {
				a = c;
			}
			c = b - GetFib(k - i - 2) / GetFib(k - i - 1) * (b - a);
			d = a + b - c;
		}
		Xopt.x = c;
		Xopt.fit_fun(ff, ud1, ud2);
		return Xopt;
	}
	catch (string ex_info)
	{
		throw ("solution fib(...):\n" + ex_info);
	}

}

solution lag(matrix(*ff)(matrix, matrix, matrix), double a, double b, double c, double epsilon, double gamma, int Nmax, matrix ud1, matrix ud2)
{
	try
	{
		solution Xopt,Xopt2,Xopt3;
		double d=0,dprev;
		int i = 0;
		Xopt2.flag = 1;
		do {
			Xopt.x = a;
			Xopt2.x = b;
			Xopt3.x = c;
			cout << "a: " << a << " b: " << b << " c: " << c << endl;
			Xopt.fit_fun(ff, ud1, ud2);
			Xopt2.fit_fun(ff, ud1, ud2);
			Xopt3.fit_fun(ff, ud1, ud2);
			double l = m2d(Xopt.y) * (pow(b, 2) - pow(c, 2)) + 
					   m2d(Xopt2.y) * (pow(c, 2) - pow(a, 2)) + 
					   m2d(Xopt3.y) * (pow(a, 2) - pow(b, 2));
			double m = m2d(Xopt.y) * (b - c) + m2d(Xopt2.y) * (c - a) + m2d(Xopt3.y) * (a - b);
			cout << "licznik: " << l << " mianownik: " << m << endl;
			if (m <= 0) {
				Xopt.flag = 0;
				return 0;
			}
			dprev = d;
			d = 0.5 * l / m;
			Xopt2.x = d;
			Xopt2.fit_fun(ff, ud1, ud2);
			if (a < d && d < c) {
				if (Xopt2.y < Xopt3.y) {
					c = d;
					b = c;
				}
				else {
					a = d;
				}
			}
			else if (c < d < b) {
				if (Xopt2.y < Xopt3.y) {
					c = d;
					b = c;
				}
				else {
					b = d;
				}
			}
			i++;
			if (solution::f_calls > Nmax) {
				Xopt.flag = 0;
				return 0;
			}
		} while (b - a < epsilon || fabs(m2d(Xopt2.x) - dprev) < gamma);
		return Xopt2;
	}
	catch (string ex_info)
	{
		throw ("solution lag(...):\n" + ex_info);
	}
}

solution HJ(matrix(*ff)(matrix, matrix, matrix), matrix x0, double s, double alpha, double epsilon, int Nmax, matrix ud1, matrix ud2)
{
	try
	{
		solution Xopt;
		//Tu wpisz kod funkcji

		return Xopt;
	}
	catch (string ex_info)
	{
		throw ("solution HJ(...):\n" + ex_info);
	}
}

solution HJ_trial(matrix(*ff)(matrix, matrix, matrix), solution XB, double s, matrix ud1, matrix ud2)
{
	try
	{
		//Tu wpisz kod funkcji

		return XB;
	}
	catch (string ex_info)
	{
		throw ("solution HJ_trial(...):\n" + ex_info);
	}
}

solution Rosen(matrix(*ff)(matrix, matrix, matrix), matrix x0, matrix s0, double alpha, double beta, double epsilon, int Nmax, matrix ud1, matrix ud2)
{
	try
	{
		solution Xopt;
		//Tu wpisz kod funkcji

		return Xopt;
	}
	catch (string ex_info)
	{
		throw ("solution Rosen(...):\n" + ex_info);
	}
}

solution pen(matrix(*ff)(matrix, matrix, matrix), matrix x0, double c, double dc, double epsilon, int Nmax, matrix ud1, matrix ud2)
{
	try {
		solution Xopt;
		//Tu wpisz kod funkcji

		return Xopt;
	}
	catch (string ex_info)
	{
		throw ("solution pen(...):\n" + ex_info);
	}
}

solution sym_NM(matrix(*ff)(matrix, matrix, matrix), matrix x0, double s, double alpha, double beta, double gamma, double delta, double epsilon, int Nmax, matrix ud1, matrix ud2)
{
	try
	{
		solution Xopt;
		//Tu wpisz kod funkcji

		return Xopt;
	}
	catch (string ex_info)
	{
		throw ("solution sym_NM(...):\n" + ex_info);
	}
}

solution SD(matrix(*ff)(matrix, matrix, matrix), matrix(*gf)(matrix, matrix, matrix), matrix x0, double h0, double epsilon, int Nmax, matrix ud1, matrix ud2)
{
	try
	{
		solution Xopt;
		//Tu wpisz kod funkcji

		return Xopt;
	}
	catch (string ex_info)
	{
		throw ("solution SD(...):\n" + ex_info);
	}
}

solution CG(matrix(*ff)(matrix, matrix, matrix), matrix(*gf)(matrix, matrix, matrix), matrix x0, double h0, double epsilon, int Nmax, matrix ud1, matrix ud2)
{
	try
	{
		solution Xopt;
		//Tu wpisz kod funkcji

		return Xopt;
	}
	catch (string ex_info)
	{
		throw ("solution CG(...):\n" + ex_info);
	}
}

solution Newton(matrix(*ff)(matrix, matrix, matrix), matrix(*gf)(matrix, matrix, matrix),
	matrix(*Hf)(matrix, matrix, matrix), matrix x0, double h0, double epsilon, int Nmax, matrix ud1, matrix ud2)
{
	try
	{
		solution Xopt;
		//Tu wpisz kod funkcji

		return Xopt;
	}
	catch (string ex_info)
	{
		throw ("solution Newton(...):\n" + ex_info);
	}
}

solution golden(matrix(*ff)(matrix, matrix, matrix), double a, double b, double epsilon, int Nmax, matrix ud1, matrix ud2)
{
	try
	{
		solution Xopt;
		//Tu wpisz kod funkcji

		return Xopt;
	}
	catch (string ex_info)
	{
		throw ("solution golden(...):\n" + ex_info);
	}
}

solution Powell(matrix(*ff)(matrix, matrix, matrix), matrix x0, double epsilon, int Nmax, matrix ud1, matrix ud2)
{
	try
	{
		solution Xopt;
		//Tu wpisz kod funkcji

		return Xopt;
	}
	catch (string ex_info)
	{
		throw ("solution Powell(...):\n" + ex_info);
	}
}

solution EA(matrix(*ff)(matrix, matrix, matrix), int N, matrix lb, matrix ub, int mi, int lambda, matrix sigma0, double epsilon, int Nmax, matrix ud1, matrix ud2)
{
	try
	{
		solution Xopt;
		//Tu wpisz kod funkcji

		return Xopt;
	}
	catch (string ex_info)
	{
		throw ("solution EA(...):\n" + ex_info);
	}
}
