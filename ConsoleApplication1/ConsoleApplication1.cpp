#include "math.h"
#include "stdafx.h"
#include "iostream"

#define PI 3.1415926535

using namespace std;

double* run(int n, double a, double b, double h ,double* y)
{
	double **arr = new double*[4];
	double *m = new double[n + 1];
	
	for (int i = 0; i < 4; i++)
	{
		arr[i] = new double[n + 1];

	}

	arr[0][0] = 1;
	arr[1][0] = 0;
	arr[3][0] = 0;

	for (int i = 1; i < n; i++)
	{
		arr[0][i] = 2 * h / 3;
		arr[1][i] = h / 6;
		arr[2][i] = h / 6;
		arr[3][i] = (y[i + 1] - y[i]) / h - (y[i] - y[i - 1]) / h;

	}

	arr[0][n] = 1;
	arr[2][n] = 0;
	arr[3][n] = 0 ;

	double *mu = new double[n + 1];
	double *l = new double[n + 1];
	
	mu[0] = arr[3][0] / arr[0][0];
	l[0] = -arr[1][0] / arr[0][0];
 
	for (int i = 1; i < n+1; i++)
	{
		mu[i] = (arr[3][i] - arr[2][i] * mu[i - 1]) / (arr[0][i] + arr[2][i] * l[i - 1]);
		l[i] = -arr[1][i] / (arr[0][i] + arr[2][i] * l[i - 1]);

	}

	m[n] = mu[n];

	for (int i = n - 1; i >= 0; i--)
	{
		m[i] = l[i] * m[i + 1] + mu[i];

	}

	delete[] mu;
	delete[] l;

	for (int i = 0; i < 4; i++)
	{
		delete[] arr[i];

	}

	return m;
}

double s3(int n,double* x,double* y ,double* m ,double h,double xp)
{
	for (int i = 1; i < n + 1; i++)
	{
		if (x[i - 1] <= xp && xp >= x[i])
			return (pow(x[i] - xp, 3) - pow(h, 2)*(x[i] - xp))*m[i - 1] / (6 * h) + 
				(pow(xp - x[i - 1], 3) - pow(h, 2)*(xp - x[i - 1]))*m[i] / (6 * h) + 
					(x[i] - xp)*y[i - 1] / h + (xp - x[i - 1])*y[i] / h;

	}
}

int main()
{
	double a = 0;
	double b = PI ;
	double deltaMaxP;
	
	for (int n = 5; n <= 10240; n *= 2)
	{
		double *x = new double[n + 1];
		double *y = new double[n + 1];
		double h = (b - a) / n;

		x[0] = a;
		x[n] = b;
		y[0] = y[n] = sin(a);

		for (int i = 1; i < n; i++)
		{
			x[i] = a + h*i;
			y[i] = sin(x[i]);

		}

		double deltaMax = 0;
		double *m = run(n, a, b, h, y);

		for (int i = 0; i < n; i++)
		{
			double delta = abs(s3(n, x, y, m, h, x[i] + h / 2) - sin(x[i] + h / 2));

			if (delta > deltaMax)
				deltaMax = delta;

		}

		if (n == 5)
		{
			cout << n << " " << deltaMax << endl;
		}
		else
		{
			cout << n << " " << deltaMax << " " << deltaMaxP / deltaMax <<endl;

		}

		deltaMaxP = deltaMax;

		delete[] m;
		delete[] x;
		delete[] y;

	}

    return 0;
}

