#include <iostream>
#include <ctime>
#include <vector>
using namespace std;
#define LIM 10
class Matrix {
	
	int ca, cb, ra, rb;

public :
	vector <vector<int>> A, B, C;
	vector <vector<int>> matrix_iterative(vector <vector<int>>a, vector <vector<int>>b);
	int arrayinput();
	void initM();
	void display();
	Matrix();

	
};

Matrix::Matrix() {

}

int Matrix::arrayinput() {
	
	cout << "Enter the number of rows and columns for Matrix A"<<endl;
	cin >> ra >> ca;
	cout << "Enter the number of rows and columns for Matrix B" << endl;
	cin >> rb >> cb;
	if (ca != rb) {
		cout << "These Matrices cannot be multiplied" << endl;
		return -1;
	}
	else {
		return 1;
	}
}
void Matrix::initM() {
	srand(time(NULL));
	A.resize(ra);
	for (int i = 0; i < ra; ++i) {
		A[i].resize(ca);
	}
	B.resize(rb);
	for (int i = 0; i < rb; ++i)
		B[i].resize(cb);
	C.resize(ra);
	for (int i = 0; i < ra; ++i)
		C[i].resize(cb);
	for (int i = 0; i < ra; ++i) {
		for (int j = 0; j < ca; j++) {
			A[i][j] = rand() %LIM +1;
			//cout << a[i][j];
		}
	}
	for (int i = 0; i < rb; ++i) {
		for (int j = 0; j < cb; j++) {
			B[i][j] = rand() % LIM + 1;
			//cout << b[i][j];
		}
	}
}
vector <vector<int>> Matrix::matrix_iterative(vector <vector<int>>ia, vector <vector<int>>ib) {

	//cout << ra << cb << ca;
	for (int i = 0; i<ra; ++i)
		for (int j = 0; j < cb; ++j) {
			C[i][j] = 0;
			for (int k = 0; k < ca; ++k)
			{
				C[i][j] += ia[i][k] * ia[k][j];
				//cout << a[i][k] * b[k][j];
			}
		}
	return C;
}
void Matrix::display() {
	cout << endl << "Output Matrix: " << endl;
	for (int i = 0; i<ra; ++i)
		for (int j = 0; j<cb; ++j)
		{
			cout << "\t" << C[i][j];
			if (j == cb - 1)
				cout << endl;
		}
}
int main() {

	Matrix m;
	
	if (m.arrayinput() != 1) {

	}
	m.initM();
	m.matrix_iterative(m.A,m.B);
	m.display();
	
	
	system("pause");
	return 0;
}
