#include <iostream>
#include <vector>
#include <math.h>
//8 вариант a = -2, b = 2
#define PI 3.1415926535
using namespace std;

double Function(double m)
{
	return m * exp(sin(m)) - 0.5 * m * m * m;
	//abs(m);
}

double PolyanomL(const vector<double>& Vx, const vector<double>& Vy, double x)
{
	double Lx = 0.0;
	for (int i = 0; i < Vx.size(); i++)
	{
		double l = 1.0;
		for (int j = 0; j < Vx.size(); j++)
		{
			if (j != i)
			{
				l *= (x - Vx[j]) / (Vx[i] - Vx[j]);
			}
		}
		Lx += Vy[i] * l;
	}
	return Lx;
}

double PolyanomNewtonForward(const vector<double> & Vx, const vector<double> & Vy, double x)
{

	double DeltaYi, y0, Resid_i, Resid_1 = 1;
	int i, j, k;
	for (i = 1, y0 = Vy[0]; i < Vy.size(); i++)
	{
		Resid_1 *= (x - Vx[i - 1]);
		for (j = 0, DeltaYi = 0; j <= i; j++)
		{
			for (k = 0, Resid_i = 1; k <= i; k++)
			{
				if (k != j)
				{
					Resid_i *= Vx[j] - Vx[k];
				}
			}
			DeltaYi += Vy[j] / Resid_i;
		}
		y0 += Resid_1 * DeltaYi;
	}
	return y0;
}


/*double PolyanomNewtonBack(vector<double> Vx, vector<double> Vy, double x)
{
	double DeltaYi, yn, Product_i, Product_1 = 1;
	int i, j, k;
	for (i = 1, yn = Vy[Vx.size() - i - 1]; i < Vx.size(); i++)
	{
		Product_1 *= (x - Vx[Vx.size() - i - 1]);
		for (j = 0, DeltaYi = 0; j <= i; j++)
		{
			for (k = 0, Product_i = 1; k <= i; k++)
			{
				if (k != j)
				{
					Product_i *= Vx[Vx.size() - i - j] - Vx[Vx.size() - i - k];
				}
			}
			DeltaYi += Vy[Vx.size() - i - j] / Product_i;
		}
		yn += Product_1 * DeltaYi;
	}
	return yn;
}*/


int main()
{
	int n;
	cout << "Enter quantity of nodes: ";
	cin >> n;
	n -= 1;//так надо..
	cout << "Enter interval boundaries: ";
	double a, b;
	cin >> a >> b;
	double h = (b - a) / n;

	double xi, xi_, Xi, Xi_;//переменные для стандартных смещённых\несмещённых, чебышевских смещённых\несмещённых
	vector<double> v1;//вектора для переменных и значений в них
	vector<double> v2;
	vector<double> v3;
	vector<double> v4;
	vector<double> v5;
	vector<double> v6;
	vector<double> v7;
	vector<double> v8;

	double Alpha;
	cout << "Enter Alpha from (0,1):";
	cin >> Alpha;

	for (int i = 0; i <= n; i++)
	{
		xi = a + i * h;//несмещённые узлы
		v1.push_back(xi);
		v2.push_back(Function(xi));

		Xi = (b + a + (b - a) * cos(PI * (2 * i + 1) / (2 * (n + 1)))) / 2;//несмещённые узлы Чебышева
		v5.push_back(Xi);
		v6.push_back(Function(Xi));
	}

	for (int i = 0; i <= n - 1; i++)//отдельный цикл для н-1 смещённых узлов
	{

		xi_ = v1[i] + Alpha * h;//смещённые узлы
		v3.push_back(xi_);
		v4.push_back(Function(xi_));

		Xi_ = v5[i] + Alpha * (v5[i + 1] - v5[i]);//смещённые узлы Чебышева, смещённые равноотстоящие 
		v7.push_back(Xi_);
		v8.push_back(Function(Xi_));
	}

	cout << endl;

	cout << "Lagrange:" << endl << endl;

	cout << "Unbiased nodes:" << endl;
	cout << "xi\t\t\t\t Y(xi)\t\t\t\t Y*(xi)\t\t\t\t |Y*i - Yi|\n ";
	cout << "---------------------------------------------------------------------------------------------------------------" << endl;
	for (int i = 0; i <= n; i++)
	{
		cout.precision(15);
		cout.setf(ios::fixed, ios::floatfield);
		cout << v1[i] << "\t\t" << v2[i] << "\t\t" << PolyanomL(v1, v2, v1[i]) << "\t\t" << abs(v2[i] - PolyanomL(v1, v2, v1[i])) << endl;
	}
	cout << endl << endl;

	cout << "Displaced nodes:" << endl;
	cout << "xi~\t\t\t\t Yi(xi~)\t\t\t Y*(xi~)\t\t\t |Y*i - Yi|\n ";
	cout << "---------------------------------------------------------------------------------------------------------------" << endl;
	for (int i = 0; i <= n - 1; i++)
	{
		cout.precision(15);
		cout.setf(ios::fixed, ios::floatfield);
		cout << v3[i] << "\t\t" << v4[i] << "\t\t" << PolyanomL(v1, v2, v3[i]) << "\t\t" << abs(v4[i] - PolyanomL(v1, v2, v3[i])) << endl;
	}

	cout << endl << endl;

	cout << "Chebyshev nodes:" << endl;
	cout << "Xi\t\t\t\t Y(Xi)\t\t\t\t Y*(Xi)\t\t\t\t |Y*i - Yi|\n ";
	cout << "---------------------------------------------------------------------------------------------------------------" << endl;
	for (int i = 0; i <= n; i++)
	{
		cout.precision(15);
		cout.setf(ios::fixed, ios::floatfield);
		cout << v5[i] << "\t\t" << v6[i] << "\t\t" << PolyanomL(v5, v6, v5[i]) << "\t\t" << abs(v6[i] - PolyanomL(v5, v6, v5[i])) << endl;
	}

	cout << endl << endl;

	cout << "Chebyshev displaced nodes:" << endl;
	cout << "Xi~\t\t\t\t Y(Xi~)\t\t\t\t Y*(Xi~)\t\t\t |Y*i - Yi|\n ";
	cout << "---------------------------------------------------------------------------------------------------------------" << endl;
	for (int i = 0; i <= n - 1; i++)
	{
		cout.precision(15);
		cout.setf(ios::fixed, ios::floatfield);
		cout << v7[i] << "\t\t" << v8[i] << "\t\t" << PolyanomL(v5, v6, v3[i]) << "\t\t" << abs(v4[i] - PolyanomL(v5, v6, v3[i])) << endl;
	}

	cout << endl << endl;
	cout << "===============================================================================================================" << endl;
	cout << endl << endl;

	cout << "NewtonForward:" << endl << endl;

	cout << "Unbiased nodes:" << endl;
	cout << "xi\t\t\t\t Y(xi)\t\t\t\t Y*(xi)\t\t\t\t |Y*i - Yi|\n ";
	cout << "---------------------------------------------------------------------------------------------------------------" << endl;
	for (int i = 0; i <= n; i++)
	{
		cout.precision(15);
		cout.setf(ios::fixed, ios::floatfield);
		cout << v1[i] << "\t\t" << v2[i] << "\t\t" << PolyanomNewtonForward(v1, v2, v1[i]) << "\t\t" << abs(v2[i] - PolyanomNewtonForward(v1, v2, v1[i])) << endl;
	}

	cout << endl << endl;

	cout << "Displaced nodes:" << endl;
	cout << "xi~\t\t\t\t Yi(xi~)\t\t\t Y*(xi~)\t\t\t |Y*i - Yi|\n ";
	cout << "---------------------------------------------------------------------------------------------------------------" << endl;
	for (int i = 0; i <= n - 1; i++)
	{
		cout.precision(15);
		cout.setf(ios::fixed, ios::floatfield);
		cout << v3[i] << "\t\t" << v4[i] << "\t\t" << PolyanomNewtonForward(v1, v2, v3[i]) << "\t\t" << abs(v4[i] - PolyanomNewtonForward(v1, v2, v3[i])) << endl;
	}

	cout << endl << endl;

	cout << "Chebyshev nodes:" << endl;
	cout << "Xi\t\t\t\t Y(Xi)\t\t\t\t Y(Xi)\t\t\t\t |Y*i - Yi|\n ";
	cout << "---------------------------------------------------------------------------------------------------------------" << endl;
	for (int i = 0; i <= n; i++)
	{
		cout.precision(15);
		cout.setf(ios::fixed, ios::floatfield);
		cout << v5[i] << "\t\t" << v6[i] << "\t\t" << PolyanomNewtonForward(v5, v6, v5[i]) << "\t\t" << abs(v6[i] - PolyanomNewtonForward(v5, v6, v5[i])) << endl;
	}

	cout << endl << endl;

	cout << "Chebyshev displaced nodes:" << endl;
	cout << "Xi~\t\t\t\t Yi(Xi~)\t\t\t Y*(Xi~)\t\t\t |Y*i - Yi|\n ";
	cout << "---------------------------------------------------------------------------------------------------------------" << endl;
	for (int i = 0; i <= n - 1; i++)
	{
		cout.precision(15);
		cout.setf(ios::fixed, ios::floatfield);
		cout << v7[i] << "\t\t" << v8[i] << "\t\t" << PolyanomNewtonForward(v5, v6, v3[i]) << "\t\t" << abs(v4[i] - PolyanomNewtonForward(v5, v6, v3[i])) << endl;
	}

	cout << endl << endl;
	cout << "================================================================================================================" << endl;
	cout << endl << endl;

	vector<double> inverse_v1;//развёрнутые векторы
	vector<double> inverse_v2;
	vector<double> inverse_v3;
	vector<double> inverse_v4;
	vector<double> inverse_v5;
	vector<double> inverse_v6;
	vector<double> inverse_v7;
	vector<double> inverse_v8;

	for (int i = n; i >= 0; i--)
	{
		inverse_v1.push_back(v1[i]);//несмещённые узлы в обратном порядке
		inverse_v2.push_back(v2[i]);


		inverse_v5.push_back(v5[i]);//несмещённые узлы чебышева в обратном порядке
		inverse_v6.push_back(v6[i]);
	}

	for (int i = n - 1; i >= 0; i--)// для н-1 смещённых узлов
	{
		inverse_v3.push_back(v3[i]);//смещённые узлы в обратном порядке
		inverse_v4.push_back(v4[i]);

		inverse_v7.push_back(v5[i]);//смещённые узлы чебышева в обратном порядке
		inverse_v8.push_back(v6[i]);
	}



	cout << "NewtonBack:" << endl << endl;
	cout << "Unbiased nodes:" << endl;
	cout << "xi\t\t\t\t Y(xi)\t\t\t\t Y*(xi)\t\t\t\t |Y*i - Yi|\n ";
	cout << "---------------------------------------------------------------------------------------------------------------" << endl;
	for (int i = 0; i <= n; i++)
	{
		cout.precision(15);
		cout.setf(ios::fixed, ios::floatfield);
		cout << inverse_v1[i] << "\t\t" << inverse_v2[i] << "\t\t" << PolyanomNewtonForward(inverse_v1, inverse_v2, inverse_v1[i]) << "\t\t" << abs(inverse_v2[i] - PolyanomNewtonForward(inverse_v1, inverse_v2, inverse_v1[i])) << endl;//векторы в обратном порядке
	}

	cout << endl << endl;

	cout << "Displaced nodes:" << endl;
	cout << "xi~\t\t\t\t Y(xi~)\t\t\t\t Y*(xi~)\t\t\t |Y*i - Yi|\n ";
	cout << "---------------------------------------------------------------------------------------------------------------" << endl;
	for (int i = 0; i <= n - 1; i++)
	{
		cout.precision(15);
		cout.setf(ios::fixed, ios::floatfield);
		cout << inverse_v3[i] << "\t\t" << inverse_v4[i] << "\t\t" << PolyanomNewtonForward(inverse_v1, inverse_v2, inverse_v3[i]) << "\t\t" << abs(inverse_v4[i] - PolyanomNewtonForward(inverse_v1, inverse_v2, inverse_v3[i])) << endl;
	}

	cout << endl << endl;

	cout << "Chebyshev nodes:" << endl;
	cout << "Xi\t\t\t\t Y(Xi)\t\t\t\t Y*(Xi)\t\t\t\t |Y*i - Yi|\n ";
	cout << "---------------------------------------------------------------------------------------------------------------" << endl;
	for (int i = 0; i <= n; i++)
	{
		cout.precision(15);
		cout.setf(ios::fixed, ios::floatfield);
		cout << inverse_v5[i] << "\t\t" << inverse_v6[i] << "\t\t" << PolyanomNewtonForward(inverse_v5, inverse_v6, inverse_v5[i]) << "\t\t" << abs(inverse_v6[i] - PolyanomNewtonForward(inverse_v5, inverse_v6, inverse_v5[i])) << endl;
	}

	cout << endl << endl;

	cout << "Chebyshev displaced nodes:" << endl;
	cout << "Xi~\t\t\t\t Yi(Xi~)\t\t\t Y*(Xi~)\t\t\t |Y*i - Yi|\n ";
	cout << "---------------------------------------------------------------------------------------------------------------" << endl;
	for (int i = 0; i <= n - 1; i++)
	{
		cout.precision(15);
		cout.setf(ios::fixed, ios::floatfield);
		cout << inverse_v7[i] << "\t\t" << inverse_v8[i] << "\t\t" << PolyanomNewtonForward(inverse_v5, inverse_v6, inverse_v3[i]) << "\t\t" << abs(inverse_v4[i] - PolyanomNewtonForward(inverse_v5, inverse_v6, inverse_v3[i])) << endl;
	}

	cout << endl << endl;
	cout << "================================================================================================================" << endl;
	cout << endl << endl;

	return 0;
}
