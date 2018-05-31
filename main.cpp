#include<iostream>
#include<cmath>
#include<stdlib.h>
#include <time.h>
#include <windows.h>
#include<fstream>
using namespace std;

double norma_residuum(int** A, int N, double * b, double *x)
{
	//r=A*x-b
	double *r = (double*)calloc(N, sizeof(double));
	for (int kolumna = 0; kolumna < N; kolumna++)//i
	{
		for (int wiersz = 0; wiersz < N; wiersz++) //j
		{
			r[kolumna] += A[wiersz][kolumna] * x[wiersz];
		}

		r[kolumna] -= b[kolumna];
	}
	//residuum
	double norma = 0;
	for (int i = 0; i < N; i++)
		norma += pow(r[i], 2);
	norma = sqrt(norma);
	return norma;

}

void gauss_seidel(int** A, int N, double * b, double *x)
{
	for (int kolumna = 0; kolumna < N; kolumna++)
	{
		double suma_poprzednich_el = 0;
		for (int wiersz = 0; wiersz < N; wiersz++)
		{
			if (wiersz == kolumna) continue;
			suma_poprzednich_el += A[kolumna][wiersz] * x[wiersz];
		}
		x[kolumna] = (b[kolumna] - suma_poprzednich_el) / A[kolumna][kolumna];
	}
}

void jacobi(int** A, int N, double * b, double *x)
{
	int ilosc = 0;
	double *x2 = (double*)calloc(N, sizeof(double));
	double suma_poprzednich_el = 0, suma_nastepnych_el = 0;

	for (int kolumna = 0; kolumna < N; kolumna++)
	{
		double suma_poprzednich_el = 0;
		for (int wiersz = 0; wiersz < N; wiersz++)
		{
			if (wiersz == kolumna) continue;
			suma_poprzednich_el += A[kolumna][wiersz] * x[wiersz];
		}
		x2[kolumna] = (b[kolumna] - suma_poprzednich_el) / A[kolumna][kolumna];
	}

	for (int i = 0; i < N; i++)
	{
		x[i] = x2[i];
	}
}
void faktoryzacjaLU(int** A, int N, double * b, double *x)
{
	double **L = (double**)malloc(sizeof(double)*N);
	double **U = (double**)malloc(sizeof(double)*N);
	double start = clock();
	for (int i = 0; i < N; i++)
	{
		L[i] = (double*)malloc(sizeof(double)*N);
		U[i] = (double*)malloc(sizeof(double)*N);
	}


	//przygotowanie L i U
	for (int i = 0; i < N; i++)
		for (int j = 0; j < N; j++)
		{
			U[j][i] = (double)A[j][i]; //A = U
			if (i == j) //L - macierz jednostkowa
				L[j][i] = 1.0;
			else L[j][i] = 0;
		}

	/*printf("U na poczatku\n");
	for (int i = 0; i < N; i++) {
		for (int j = 0; j < N; j++) {
			if (j % 10 == 0)
				printf("\n");
			printf("%f ", U[j][i]);
		}
	}

	printf("\nL na poczatku \n");
	for (int i = 0; i < N; i++) {
		for (int j = 0; j < N; j++) {
			if (j % 10 == 0)
				printf("\n");
			printf("%f ", L[j][i]);
		}
	}*/
	for (int k = 0; k < N - 1; k++) {
		for (int j = k + 1; j < N; j++) {
			L[k][j] = U[k][j] / U[k][k];
			for (int i = k; i < N; i++) {
				//printf("%f - %f * %f =", U[i][j],L[k][j],U[i][k]);
				U[i][j] -= (L[k][j] * U[i][k]);
				//	printf("%f \n", U[i][j]);
			}
		}
	}
	
	//y=Ux
	double *y = (double*)malloc(sizeof(double)*N);

	//podstawienie w przod
	for (int i = 0; i < N; i++)
	{
		double suma = 0;
		for (int j = 0; j < i; j++) {
			suma += L[j][i] * y[j];
		}

		y[i] = (b[i] - suma) / L[i][i];
	}

	// Podstawienie w tyl Ux = y
	for (int i = N - 1; i >= 0; i--) {
		double suma = 0;
		for (int j = i; j < N; j++) {
			suma += U[j][i] * x[j];
		}
		x[i] = (y[i] - suma) / U[i][i];

	}

	//printf("\nX\n");
	//for (int i = 0; i < N; i++) {
	//	printf("%f ",x[i]);
	//}
	double stop = clock();
	printf("\nCzas: %f\n", double(stop - start) / CLOCKS_PER_SEC);
	cout << "Norma: " << norma_residuum(A, N, b, x) << endl;
}
void init(int f, int N, double *x, double *b)
{
	for (int i = 0; i < N; i++)
	{
		b[i] = sin(i*(f + 1));
		x[i] = 0;
	}
}

int main() {
	int  index[] = { 1,6,5,6,0,9 };
	int a1 = index[3] + 5;
	int a2 = -1, a3 = -1;
	int N = 900 + index[4] * 10 + index[5];
	int f = index[2];
	//a1 = 3; a2 = a3 = -1; //zadanie C
	N = 3200;
	int **A = (int**)malloc(sizeof(int)*N);
	for (int i = 0; i < N; i++)
		A[i] = (int*)malloc(sizeof(int)*N);
	double *b = (double*)malloc(sizeof(double)*N);
	double *x = (double*)calloc(N, sizeof(double));

	//zapelnianie macierzy
	for (int wiersz = 0; wiersz < N; wiersz++)
	{
		for (int kolumna = 0; kolumna < N; kolumna++)
		{
			if (wiersz == kolumna)
				A[kolumna][wiersz] = a1;
			else if (abs(wiersz - kolumna) == 1)
				A[kolumna][wiersz] = a2;
			else if (abs(wiersz - kolumna) == 2)
				A[kolumna][wiersz] = a3;
			else
				A[kolumna][wiersz] = 0;

			//	if (A[kolumna][wiersz] == 0) printf << " ";
			//printf << A[kolumna][wiersz] << "  ";
		}
		//printf << endl;
	}

	//ZadanieC
	ofstream outJacobi("outJacobi.txt");
	ofstream outGauss("outGauss.txt");


	//Zadanie A-B

	init(f, N, x, b);

	double koniec = pow(10, -9);
	double norma = koniec + 1;
		int i = 0;
	double start = clock();
	printf( "Gauss-Seidel: \n");
	while (i < 100 && norma > koniec)
	{
	gauss_seidel(A, N, b, x);
	norma = norma_residuum(A, N, b, x);
	if (norma <= koniec)
	break;
	//outGauss << norma << endl;
	i++;
	}
	//printf("Liczba iteracji: %d \n", i);

	double stop = clock();
	cout << "Czas: " << double(stop - start) / CLOCKS_PER_SEC;

	//JACOBI
	init(f, N, x, b);
	i = 0; norma = koniec + 1; start = clock();
	printf("\nJacobi: \n" );

	while (i < 100 && norma > koniec)
	{
	jacobi(A, N, b, x);
	norma = norma_residuum(A, N, b, x);
		if (norma <= koniec)
			break;
		i++;
	//outJacobi << norma << endl;
		}
	//printf("Liczba iteracji: %d \n", i);
	stop = clock();
	cout << "Czas: " << double(stop - start) / CLOCKS_PER_SEC;


	//ZadanieD
	init(f, N, x, b);
	printf("\nFaktoryzacja ");
	faktoryzacjaLU(A, N, b, x);

	getchar();
	return 0;
}


//cout <<"Norma po"<< i <<"iteracjach: "<< norma << endl;	


/*printf("U\n");
for (int i = 0; i < N; i++) {
for (int j = 0; j < N; j++) {
if (j % 10 == 0)
printf("\n");
printf("%f ", U[j][i]);
}
}

printf("\nL\n");
for (int i = 0; i < N; i++) {
for (int j = 0; j < N; j++) {
if (j % 10 == 0)
printf("\n");
printf("%f ", L[j][i]);
}
}*/
