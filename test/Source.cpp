#include <iostream>
#include <cstdlib>
#include <omp.h>
#include <vector>
#include <cmath>
#include <iomanip>
using namespace std;


//Задача №1
//
//const double pi = 3.1415;
//vector <double> sinus;
//
//int main() {
//	double t = omp_get_wtime();
//	int n=100000000;
//	sinus.resize(n+1);
//#pragma omp parallel for schedule(dynamic, 1000)
//	for (int i = 0; i <= n; i++) {
//		sinus[i] = sin(i * pi / (2 * n));
//	}
//	cout << "Time: " << (omp_get_wtime() - t )*1000<<" ms"<< endl;
//	return 0;
//}

//Задача №2

//int main() {
//	double t = omp_get_wtime();
//	int n = 100000000;
//	double sum = 0;
//	double pi;
//#pragma omp parallel for schedule(guided, 1000) reduction(+:sum) 
//	for (int i = 1; i <= n; i++) {
//		sum += 1 / (1 + (2. * i - 1.) / (2. * n) * (2. * i - 1.) / (2. * n));
//	}
//
//	pi = 4*sum / n;
//	cout << pi << endl;
//	cout << "Time: " << (omp_get_wtime() - t) * 1000 << " ms" << endl;
//	return 0;
//}


//Задача №3


//int main() {
//	double t = omp_get_wtime();
//	int n = 100000000;
//	int k = 1;
//	bool f1;
//
//#pragma omp parallel for schedule(guided, 100) reduction(+:k)
//	for (int i = 3; i <= n ; i += 2) {
//		f1 = true;
//		for (int j = 2; j * j <= i; j++) {
//			if (i % j == 0) {
//				f1 = false;
//				break;
//			}
//		}
//		if (f1) {
//			k++;
//		}
//	}
//
//
//	cout << k << endl;
//	cout << "Time: " << (omp_get_wtime() - t) * 1000 << " ms" << endl;
//}



//Задача №4

struct Point {
	double x, y, z;
	Point() {
		x = 0; y = 0; z = 0;
	}
	Point(double a, double b, double c) {
		x = a; y = b; z = c;
	}
};
double dist(Point a, Point b) {
	return sqrt((a.x - b.x) * (a.x - b.x) + (a.y - b.y) * (a.y - b.y) + (a.z - b.z) * (a.z - b.z));
}
int main() {
	double t = omp_get_wtime();
	int n = 10000; 
	double max=0;
	vector <Point> setOfPoints(n);
	vector <vector <double>> matrix(n, vector <double>(n));
	for (int i = 0; i < n; i++) {
		setOfPoints[i].x = rand() % 100;
		setOfPoints[i].y = rand() % 100;
		setOfPoints[i].z = rand() % 100;
		//cout << setOfPoints[i].x << " " << setOfPoints[i].y << " " << setOfPoints[i].z<<endl;
	}

#pragma omp parallel for schedule(guided) collapse(2)
	for (int i = 0; i < n; i++) {
		for (int j = i; j < n; j++) {
			matrix[i][j] = dist(setOfPoints[i], setOfPoints[j]);
			if (matrix[i][j] > max) {
				max = matrix[i][j];
			}
			//cout << setw(10)<<matrix[i][j] << " ";
		}
		//cout << endl;
	}



	cout << max << endl;
	cout << "Time: " << (omp_get_wtime() - t) * 1000 << " ms" << endl;
}

//Задача №5

//int main() {
//	double t = omp_get_wtime();
//	int n, k, m;
//	n = 100;
//	k =100;
//	vector <vector <bool>> life(n+2, vector <bool>(k+2)); //таблица игрового поля
//	vector <vector <bool>> boof(n + 2, vector <bool>(k + 2)); //буферное игровое поле, чтобы правильно считать клетки вокруг
//
//	//Заполнение границы поля нулями
//	for (int i = 0; i < n + 2; i++) {
//		life[i][0] = 0;
//		life[i][k + 1] = 0;
//	}
//	for (int i = 0; i < k + 2; i++) {
//		life[0][i] = 0;
//		life[n + 1][i] = 0;
//	}
//
//	//Заполнение поля рандомными значениями
//	for (int i = 1; i < n + 1; i++) {
//		for (int j = 1; j < k + 1; j++) {
//			life[i][j] = rand() % 2;
//			boof[i][j] = life[i][j];
//			if (life[i][j] == true) {
//				//cout << setw(2)<<"1";
//			}
//			else {
//				//cout <<setw(2)<< "0";
//			}
//		}
//		//cout << endl;
//	}
//
//	//Расчёт следующего хода
//	for (int i = 1; i < n + 1; i++) {
//		for (int j = 1; j < k + 1; j++) {
//			m = boof[i - 1][j - 1] + boof[i - 1][j] + boof[i - 1][j + 1] + boof[i][j + 1] + boof[i + 1][j + 1] + boof[i + 1][j] + boof[i + 1][j - 1] + boof[i][j - 1];
//			if (life[i][j] == 1) {
//				if (m < 2 || m>3) {
//					life[i][j] = 0;
//				}
//			}
//			else {
//				if (m == 3) {
//					life[i][j] = 1;
//				}
//			}
//		}
//	}
//	//Параллельная на 4 секции
//
//	//#pragma omp parallel sections
//	//	{
//	//#pragma omp section
//	//		{
//	//			for (int i = 1; i < (n + 1)/2; i++) {
//	//				for (int j = 1; j < (k + 1)/2; j++) {
//	//					m = boof[i - 1][j - 1] + boof[i - 1][j] + boof[i - 1][j + 1] + boof[i][j + 1] + boof[i + 1][j + 1] + boof[i + 1][j] + boof[i + 1][j - 1] + boof[i][j - 1];
//	//					if (life[i][j] == 1) {
//	//						if (m < 2 || m>3) {
//	//							life[i][j] = 0;
//	//						}
//	//					}
//	//					else {
//	//						if (m == 3) {
//	//							life[i][j] = 1;
//	//						}
//	//					}
//	//				}
//	//			}
//	//		}
//	//#pragma omp section
//	//		{
//	//			for (int i = (n+1)/2; i < n + 1; i++) {
//	//				for (int j = 1; j < (k + 1)/2; j++) {
//	//					m = boof[i - 1][j - 1] + boof[i - 1][j] + boof[i - 1][j + 1] + boof[i][j + 1] + boof[i + 1][j + 1] + boof[i + 1][j] + boof[i + 1][j - 1] + boof[i][j - 1];
//	//					if (life[i][j] == 1) {
//	//						if (m < 2 || m>3) {
//	//							life[i][j] = 0;
//	//						}
//	//					}
//	//					else {
//	//						if (m == 3) {
//	//							life[i][j] = 1;
//	//						}
//	//					}
//	//				}
//	//			}
//	//		}
//	//#pragma omp section
//	//		{
//	//			for (int i = 1; i < (n + 1)/2; i++) {
//	//				for (int j = (k+1)/2; j < k + 1; j++) {
//	//					m = boof[i - 1][j - 1] + boof[i - 1][j] + boof[i - 1][j + 1] + boof[i][j + 1] + boof[i + 1][j + 1] + boof[i + 1][j] + boof[i + 1][j - 1] + boof[i][j - 1];
//	//					if (life[i][j] == 1) {
//	//						if (m < 2 || m>3) {
//	//							life[i][j] = 0;
//	//						}
//	//					}
//	//					else {
//	//						if (m == 3) {
//	//							life[i][j] = 1;
//	//						}
//	//					}
//	//				}
//	//			}
//	//		}
//	//#pragma omp section
//	//		{
//	//			for (int i = (n+1)/2; i < n + 1; i++) {
//	//				for (int j = (k+1)/2; j < k + 1; j++) {
//	//					m = boof[i - 1][j - 1] + boof[i - 1][j] + boof[i - 1][j + 1] + boof[i][j + 1] + boof[i + 1][j + 1] + boof[i + 1][j] + boof[i + 1][j - 1] + boof[i][j - 1];
//	//					if (life[i][j] == 1) {
//	//						if (m < 2 || m>3) {
//	//							life[i][j] = 0;
//	//						}
//	//					}
//	//					else {
//	//						if (m == 3) {
//	//							life[i][j] = 1;
//	//						}
//	//					}
//	//				}
//	//			}
//	//		}
//	//	}
//	
//
//	//cout << "-----------------" << endl;
//
//	//for (int i = 1; i < n + 1; i++) {
//	//	for (int j = 1; j < k + 1; j++) {
//	//		if (life[i][j] == true) {
//	//			//cout << setw(2) << "1";
//	//		}
//	//		else {
//	//			//cout << setw(2) << "0";
//	//		}
//	//	}
//	//	//cout << endl;
//	//}
//	cout << "Time: " << (omp_get_wtime() - t) * 1000 << " ms" << endl;
//}