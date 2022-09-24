#include <iostream>
#include <vector>
using namespace std;


double* Iteration(int n, const vector<double> &B, const vector<vector<double>> &matrix, double E)
{
	
	vector<double> result1;
	vector<double> Approximation = B;
	double* result = new double[n];
	
	do
	{
		for (int i = 0; i < n; i++)
		{
			result[i]=(B[i] / matrix[i][i]);
			for (int j = 0; j < n; j++)
			{
				if (i == j)
				{
					continue;
				}
				else
				{
					result[i] -= (matrix[i][j] / matrix[i][i]) * Approximation[j];
				}
			}
		}

		bool pointer = true;
		for (int i = 0; i < n - 1; i++)
		{
			if (abs(result[i] - Approximation[i]) > E)
			{
				pointer = false;
				break;
			}
		}

		for (int i = 0; i < n; i++)
		{
			Approximation[i] = result[i];
		}

		if (pointer == 1) 
		{
			break;
		}

	} while (true); 
	

	return result;
}

int main()
{
	int n;
	cout << "Enter quantity of equations:";
	cin >> n;
	vector<vector<double>> matrix;
	vector<double> storage;
	vector<double> B;
	double a, b;
	
	cout << "Enter matrix coefficients:" << endl;
	for (int i = 0; i < n; i++)
	{
		for (int j = 0; j < n; j++)
		{
			cin >> a;
			storage.push_back(a);
		}
		matrix.push_back(storage);
		storage.clear();
	}

	cout << "Enter free member column :" << endl;
	for (int i = 0; i < n; i++)
	{
		cin >> b;
		B.push_back(b);
	}
	
	cout<<"Entered matrix with free member column:"<<endl;
	for (int i = 0; i < n; i++)
	{
		for (int j = 0; j < n; j++)
		{
			cout << matrix[i][j] << " ";
		}
		cout << "|" << B[i];
		cout << endl;
	}
	double eps;
	cout << "Enter epsilon:";
	cin >> eps;

	cout << endl;
	cout << "Solution:" << endl;
	Iteration(n, B, matrix, eps);
	for (int i = 0; i < n; i++)
	{
		cout << Iteration(n, B, matrix, eps)[i] << "|";
	}

	return 0;
}
