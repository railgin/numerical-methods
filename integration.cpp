#include <iostream>
#include <vector>
#include <math.h>
using namespace std;
#define PI 3.1415926535

double Function(double m)
{
	return sin(m);
}

double IntTr(int n, double a, double b, double e)
{
	double I = 0;
	double h = (b - a) / n;//шаг интегрирования
	for (int i = 0; i < n; i++)
	{
		double x1 = a + i * h;
		double x2 = a + (i + 1) * h;
		I += 0.5 * (x2 - x1) * (Function(x2) + Function(x1));
	}
	return I;
}


double IntSimps(int n, double a, double b, double e)
{
	double I = 0;
	double h = (b - a) / n;//шаг интегрирования
	double Inext = I;
	for (int i = 0; i < n; i++)
	{
		double x1 = a + i * h;
		double x2 = a + (i + 1) * h;
		I += (x2 - x1) / 6.0 * (Function(x1) + 4.0 * Function(0.5 * (x1 + x2)) + Function(x2));
	}

	return I;
}



int main()
{
	int N;
	cout << "Enter quantity of sections:";
	cin >> N;
	int copyN = N;
	cout << "Enter integration limits:";
	cout << endl;
	double a, b;
	cout << "a = ";
	cin >> a;
	cout << "b = ";
	cin >> b;
	cout << endl;
	double E;
	cout << "Enter epsilon:";
	cin >> E;
	cout << endl << endl;


	cout.precision(15);
	cout.setf(ios::fixed, ios::floatfield);

	cout << "Methods without double re-compute:" << endl << endl;

	cout << "Trapezoidal formula: " << endl;
	cout << "Int = " << IntTr(N, a, b, E) << endl;
	cout << endl;

	cout << "Simpson's formula: " << endl;
	cout << "Int = " << IntSimps(N, a, b, E) << endl;
	cout << endl;

	double T = IntTr(N, a, b, E);
	double S = IntSimps(N, a, b, E);
	cout << "|IntTr - IntSim| = " << abs(IntTr(N, a, b, E) - IntSimps(N, a, b, E)) << endl << endl << endl;

	cout << "Double recompute:" << endl << endl;

	int counter1 = 0;
	int counter2 = 0;

	while (abs(IntTr(N, a, b, E) - IntTr(2*N, a, b, E)) > E)
	{
		N = N * 2;
		counter1++;
	}
	cout << "Final N for trapezoidal formula = " << N << endl;
	cout << "Quantity of double-recompute iterations, trapezoidal formula: " << counter1 << endl;

	cout << "Trapezoidal integral after double recompute: " << endl;
	cout << "Int = " << IntTr(N, a, b, E) << endl;
	cout << endl;

	N = copyN;

	while (abs(IntSimps(N, a, b, E) - IntSimps(2 * N, a, b, E)) > E)
	{
		N = N * 2;
		counter2++;
	}
	cout << "Final N for Simpson's formula = " << N << endl;
	cout << "Quantity of double-recompute iterations, Simpsons formula: " << counter2 << endl;
	
	cout << "Simpsons integral after double recompute: " << endl;
	cout << "Int = " << IntSimps(N, a, b, E) << endl;
	cout << endl << endl;

	cout << "Difference between basic-N formulas and double-recompute method:" << endl;
	cout << "|IntTr* - IntTr| = " << abs(IntTr(N, a, b, E) - T) << endl;
	cout << "|IntSim* - IntSim| = " << abs(IntTr(N, a, b, E) - S) << endl;

	return 0;
}
