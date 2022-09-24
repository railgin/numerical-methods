#include <iostream>
#include <vector>
#include <math.h>
#define PI 3.1415926535
using namespace std;

double Function(double m)
{
	return exp(-m * m);
}

double Splain(const vector<double>& Vx, const vector<double>& Vy, double x, int n, double m[], double h)
{
	double S = 0; 
	double x1, x2;
	for (int i = 1; i <= n; i++)
	{
		x1 = Vx[i] - x;
		x2 = x - Vx[i - 1];
		if (x <= Vx[i] && x >= Vx[i - 1])
		{
			S = m[i - 1] * x1 * x1 * x1 / (6.0 * h) + m[i] * x2 * x2 * x2 / (6.0 * h) + (Vy[i] - m[i] * h * h / 6.0) * x2 / h +
				+(Vy[i - 1] - m[i - 1] * h * h / 6.0) * x1 / h;
		}
	}
	return S;
}

void methodProgonki(const vector<double>& Vx, const vector<double>& Vy, int n, double h, double m[])
{
	double* g = new double[n];
	g[0] = 0.0;
	g[n] = 0.0;
	for (int i = 1; i <= n - 1; i++)
	{
		g[i] = -(Vy[i] - Vy[i - 1]) / h + (Vy[i + 1] - Vy[i]) / h;
	}
	m[0] = 0.0;
	m[n] = 0.0;

	double* alpha = new double[n+1];
	double* beta = new double[n+1];
	
	alpha[0] = 0.0;
	beta[0] = g[0];

	double ai = h / 6.0;
	double bi = 2 * h / 3.0;
	double ci = h / 6.0;

	for (int i = 0; i <= n; i++)
	{
		alpha[i + 1] = -ci / (alpha[i]*ai + bi);
		beta[i + 1] = (g[i] - ai * beta[i]) / (ai * alpha[i] + bi);
	}
	for (int i = n - 1; i > 0; i--)
	{
		m[i] = alpha[i + 1] * Vx[i + 1] + beta[i + 1];
	}
}


int main()
{

	int n;
	cout << "Enter quantity of nodes: ";
	cin >> n;
	n -= 1;
	cout << "Enter interval boundaries: ";
	double a, b;
	cin >> a >> b;
	double h = (b - a) / n;
	double* m1 = new double[n];//коэффициенты сплайна узлов
	double* m2 = new double[n];


	double xi, xi_;//переменные для стандартных смещённых\несмещённых
	vector<double> v1;//вектора для переменных и значений функции в них
	vector<double> fi;
	vector<double> SmeshRavn;
	vector<double> FuncSmeshRavn;

	double Alpha;
	cout << "Enter Alpha from (0,1):";
	cin >> Alpha;
	
	for (int i = 0; i <= n; i++)
	{
		xi = a + i * h;//несмещённые узлы
		v1.push_back(xi);
		fi.push_back(Function(xi));
	}

	for (int i = 0; i <= n - 1; i++)//отдельный цикл для н-1 смещённых узлов
	{

		xi_ = v1[i] + Alpha * h;
		SmeshRavn.push_back(xi_);
		FuncSmeshRavn.push_back(Function(xi_));

	}


	cout.setf(ios::fixed, ios::floatfield);
	cout.precision(15);
	
	cout << endl << endl;
	methodProgonki(v1, fi, n, h, m1);
	methodProgonki(SmeshRavn, FuncSmeshRavn, n - 1, h, m2);
	
	

	cout << "Unbiased nodes:" << endl;
	cout << "xi\t\t\t\t Y(xi)\t\t\t\t YS(xi)\t\t\t\t |YSi - Yi|\n ";
	cout << "---------------------------------------------------------------------------------------------------------------" << endl;
	for (int i = 0; i <= n; i++)
	{
		cout << v1[i] << "\t\t" << fi[i] << "\t\t" << Splain(v1, fi, v1[i], n, m1, h) << "\t\t" << abs(fi[i] - Splain(v1, fi, v1[i], n, m1, h)) << endl;
	}

	cout << endl << endl;

	cout << "Displaced nodes:" << endl;
	cout << "xi~\t\t\t\t YS(xi~)\t\t\t YS(xi~)\t\t\t |YSi - Yi|\n ";
	cout << "---------------------------------------------------------------------------------------------------------------" << endl;
	for (int i = 0; i <= n - 1; i++)
	{
		cout << SmeshRavn[i] << "\t\t" << FuncSmeshRavn[i] << "\t\t" << Splain(v1, fi, SmeshRavn[i], n, m2, h) << "\t\t" << abs(FuncSmeshRavn[i] - Splain(v1, fi, SmeshRavn[i], n, m2, h)) << endl;
	}


	cout << endl << endl;
	
	return 0;
}
