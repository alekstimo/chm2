#include <fstream>
#include <iostream>
#include <iomanip>
using namespace std;

void readMatrixFromFile(ifstream& fin, int n, int l, double** a) {
	for (int i = 0; i < n; i++)
		for (int j = 0; j < l; j++)
			fin >> a[i][j];
}

void readVectorFromFile(ifstream& fin, int n, double* f) {
	for (int i = 0; i < n; i++)
		fin >> f[i];
}

void readFromFile(ifstream& fin, int n, int l, double** a, double* f) {
	readMatrixFromFile(fin, n, l, a);
	readVectorFromFile(fin, n, f);
}

bool matrixFromFileIsCorrect(int n, int l, double** a) {
	bool flagOnExit = false;
	for (int i = n - 1, indexOfNull = 1; i > n - l && !flagOnExit; i--, indexOfNull++)
		for (int j = indexOfNull; j < l && !flagOnExit; j++)
			if (a[i][j])
				flagOnExit = true;
	return !flagOnExit;
}

void toMirrorMatrix(int n, int l, double** topA, double** bottomA) {
	int nullElementsInString = l - 1;
	for (int i = 0; i < l - 1; i++, nullElementsInString--)
		for (int j = 0; j < nullElementsInString; j++)
			bottomA[i][j] = 0;
	nullElementsInString = 0;
	for (int j = 0; j < l; j++, nullElementsInString++)
		for (int i = 0; i < n - nullElementsInString; i++)
			bottomA[i + nullElementsInString][l - 1 - j] = topA[i][j];
}

void printMatrix(int n, int l, double** topA, double** bottomA, double* f, int setPrecision) {
	int nullElementsInStringInBottom = l - 1;
	int nullElementsInStringInTop = 0;
	int tabCounterFromBegin = 0 - (l - 1);
	int tabCounterFromEnd = n - l;

	for (int i = 0; i < n; i++) {
		for (int k = 0; k < tabCounterFromBegin; k++)
			std::cout << '\t';
		tabCounterFromBegin++;

		if (nullElementsInStringInBottom > -1) {
			for (int j = nullElementsInStringInBottom; j < l - 1; j++)
				std::cout << setprecision(setPrecision) << bottomA[i][j] << '\t';
			nullElementsInStringInBottom--;
		}
		else {
			for (int j = 0; j < l - 1; j++)
				std::cout << setprecision(setPrecision) << bottomA[i][j] << '\t';
			nullElementsInStringInBottom--;
		}

		if (i > n - l)
			nullElementsInStringInTop++;

		for (int j = 0; j < l - nullElementsInStringInTop; j++)
			std::cout << setprecision(setPrecision) << topA[i][j] << '\t';

		for (int k = 0; k < tabCounterFromEnd; k++)
			std::cout << '\t';
		tabCounterFromEnd--;

		std::cout << "= " << f[i] << endl;
	}
	std::cout << endl;
}

void printVector(int n, double* vector, int setPrecision, bool isSolution) {
	if (isSolution) {
		for (int i = 0; i < n; i++)
			std::cout << "x" << i + 1 << " = " << setprecision(setPrecision) << vector[i] << endl;
		std::cout << endl;
	}
	else
		for (int i = 0; i < n; i++)
			std::cout << setprecision(setPrecision) << vector[i] << '\t';
}

void findSolution(int n, int l, double** topA, double** bottomA, double* f, double* x) {
	int t1 = l;
	double koef;
	for (int i = 0; i < n; i++) {
		koef = topA[i][0];
		topA[i][0] = bottomA[i][l - 1] = 1;
		if (n - i < l)
			t1--;
		for (int j = 1; j < t1; j++)
			topA[i][j] /= koef;
		f[i] /= koef;
		int t2 = t1;
		for (int j = l - 2, k = i + 1; j >= 0 && k < n; j--, k++) {
			koef = bottomA[k][j];
			bottomA[k][j] = 0;
			//int t3 = j;
			int t3 = 1;
			for (int p = j + 1; p < l - 1; p++) {
				bottomA[k][p] -= topA[i][p] * koef;
				//t3 = p + 1;
				t3++;
			}
			if (n - k < l)
				t2--;
			for (int p = 0; p < t2 && t3 < l; p++, t3++) {
				topA[k][p] -= topA[i][t3] * koef;
				if (p == 0)
					bottomA[k][l - 1] = topA[k][p];
			}
			f[k] -= f[i] * koef;
		}
		//printMatrix(n, l, topA, bottomA, f, 5);
	}
	koef = 0;
	for (int k = n - 1; k >= 0; k--) {
		for (int i = k - 1, j = 1; i > k - l && i >= 0; i--, j++) {
			koef = topA[i][j];
			topA[i][j] = 0;
			f[i] -= f[k] * koef;
		}
		//printMatrix(n, l, topA, bottomA, f, 5);
	}
}

void genVector(int n, double* vector, int e, int nullElementsCounter) {
	for (int i = 0; i < n - nullElementsCounter; i++) {
		double el = 1.0137 * (2 * (rand() % e) - e);
		while (el == 0) {
			el = 1.0137 * (2 * (rand() % e) - e);
		}
		vector[i] = el;
	}
	if (nullElementsCounter != 0)
		for (int i = n - nullElementsCounter; i < n; i++)
			vector[i] = 0;
}

void genMatrix(int n, int l, double** a, double* f, double* randomX, int e) {
	for (int i = 0; i < n - l; i++)
		genVector(l, a[i], e, 0);
	for (int i = n - l, nullElementsCounter = 0; i < n; i++, nullElementsCounter++)
		genVector(l, a[i], e, nullElementsCounter);
	genVector(n, f, e, 0);
	genVector(n, randomX, e, 0);
}

void countResultF(int n, int l, double** bottomA, double** topA, double* randomX, double* resultF) {
	/*int nullsInBottom = l - 1;
	int nullsInTop = 0;
	for (int i = 0; i < n; i++) {
		int k = l - 1;
		if (nullsInBottom > -1) {
			for (int j = nullsInBottom; j < l - 1; j++)
				if (i - k > -1) {
					resultF[i] += bottomA[i][j] * randomX[i - k];
					k++;
				}
			nullsInBottom--;
		}
		else {
			for (int j = 0; j < l - 1; j++)
				if (i - k > -1) {
					resultF[i] += bottomA[i][j] * randomX[i - k];
					k++;
				}
			nullsInBottom--;
		}

		if (i > n - l)
			nullsInTop++;

		for (int j = 0; j < l - nullsInTop; j++)
			if (i + k < n) {
				if (i - k > -1) {
					resultF[i] += topA[i][j] * randomX[i + k];
					k++;
				}
			}
	}*/
	for (int j = 0; j < l; j++) {
		for (int i = l - j - 1, k = 0; i < n; i++, k++) {
			resultF[i] += bottomA[i][j] * randomX[k];
		}
	}
	
	for (int i = 0; i < n-1; i++) {
		for (int j = 1,p=0; j < l-p; j++) {
			resultF[i] += topA[i][j] * randomX[i+1];
			if (n - i < l)
				p++;
		}
	}

}

void findingSolution2(int n, int l, double** topA, double** bottomA, double* f, double* x, double* resultF) {
	int t1 = l;
	double koef;
	for (int i = 0; i < n; i++) {
		koef = topA[i][0];
		topA[i][0] = bottomA[i][l - 1] = 1;
		if (n - i < l)
			t1--;
		for (int j = 1; j < t1; j++)
			topA[i][j] /= koef;
		f[i] /= koef;
		resultF[i] /= koef;
		int t2 = t1;
		for (int j = l - 2, k = i + 1; j >= 0 && k < n; j--, k++) {
			koef = bottomA[k][j];
			bottomA[k][j] = 0;
			int t3 =1;
			for (int p = j + 1; p < l - 1; p++) {
				bottomA[k][p] -= topA[i][t3] * koef;
				t3++;
			}
			if (n - k < l)
				t2--;
			for (int p = 0; p < t2 && t3 < l; p++, t3++) {
				topA[k][p] -= topA[i][t3] * koef;
				if (p == 0)
					bottomA[k][l - 1] = topA[k][p];
			}
			f[k] -= f[i] * koef;
			resultF[k] -= resultF[i] * koef;
		}
		//printMatrix(n, l, topA, bottomA, f, 5);
	}
	koef = 0;
	for (int k = n - 1; k >= 0; k--) {
		for (int i = k - 1, j = 1; i > k - l && i >= 0; i--, j++) {
			koef = topA[i][j];
			topA[i][j] = 0;
			f[i] -= f[k] * koef;
			resultF[i] -= resultF[k] * koef;
		}
		//printMatrix(n, l, topA, bottomA, f, 5);
	}
}

double searchingQ(int n, double* f) {
	double q = f[0];
	for (int i = 0; i < n; i++) {
		if (q < f[i])
			q = f[i];
	}
	return q;
}

double searchingDelta(int n, double q, double* resultF, double* randomX) {
	double d = 0;
	for (int i = 0; i < n; i++) {
		if (abs(resultF[i]) > q) {
			double pr = abs((resultF[i] - randomX[i]) / randomX[i]);
			if (pr > d)
				d = pr;
		}
		else {
			double pr = abs(resultF[i] - randomX[i]);
			if (pr > d)
				d = pr;
		}
	}
	return d;
}
void InitUnitVector(double* x, int n) {
	for (int i = 0; i<n; i++) {
		x[i] = 1;
	}
}

int main() {
	srand(time(NULL));
	setlocale(LC_ALL, "Russian");
	ifstream fin("test.txt");

	if (!fin.is_open())
		std::cout << "Ошибка открытия файла" << endl;
	else {
		int n; fin >> n;
		int l; fin >> l;
		double** topA = new double* [n], ** bottomA = new double* [n];
		for (int i = 0; i < n; i++) {
			topA[i] = new double[l];
			bottomA[i] = new double[l];
		}
		double* f = new double[n];
		double* x = new double[n];

		readFromFile(fin, n, l, topA, f);
		if (matrixFromFileIsCorrect(n, l, topA)) {
			toMirrorMatrix(n, l, topA, bottomA);
			printMatrix(n, l, topA, bottomA, f, 0);
			findSolution(n, l, topA, bottomA, f, x);
			printVector(n, f, 0, true);
		}
		else std::cout << "Матрица из файла не подходит по формату входных данных." << endl;
	}

	std::cout << "_______________________________Далее работаем с генерацией_______________________________" << endl << endl;
	int n = 100;
	int l = 10;
	int e = 10;
	double** topA = new double* [n], ** bottomA = new double* [n];
	for (int i = 0; i < n; i++) {
		topA[i] = new double[l];
		bottomA[i] = new double[l];
	}
	double* f = new double[n];
	double* x = new double[n];
	double* randomX = new double[n];
	double* unitX = new double[n];
	double* resultF = new double[n];
	for (int i = 0; i < n; resultF[i] = 0, i++);

	genMatrix(n, l, topA, f, randomX, e);
	toMirrorMatrix(n, l, topA, bottomA);
	InitUnitVector(unitX, n);
	//printMatrix(n, l, topA, bottomA, f, 5);
	//printVector(n, randomX, 4, false);
	//printVector(n, unitX, 4, false);
	//countResultF(n, l, bottomA, topA, randomX, resultF);
	countResultF(n, l, bottomA, topA, unitX, resultF);
	
	
	std::cout << endl;
	//printVector(n, resultF, 4, false);
	findingSolution2(n, l, topA, bottomA, f, x, resultF);
	std::cout << endl;
	//printVector(n, resultF, 4, false);

	//double q = searchingQ(n, randomX);
	double q = searchingQ(n, unitX);
	//double delta = searchingDelta(n, q, resultF, randomX);
	double delta = searchingDelta(n, q, resultF, unitX);
	std::cout << "Средняя относительная погрешность: " << endl;
	std::cout << " deltaX = " << setprecision(20) << delta << endl;

	std::system("Pause");
	return 0;
}