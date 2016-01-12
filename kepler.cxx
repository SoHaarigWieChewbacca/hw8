#include <iostream>
#include <fstream>
#include <cmath>

using namespace std;

void Euler(double* p1, double* p2, double* q1, double* q2, double* H, const double dt, const int N);
void Output(double* p1, double* p2, double* q1, double* q2, double* H, const double dt, const int N);

int main(){

	const double dt = 0.0005;
	const double period = 2*M_PI;
	const int N = (10*period/dt) + 1;
	
	double q1[N];
	double q2[N];
	double p1[N];
	double p2[N];
	double H[N];

	const double e = 0.6;
	q1(0) = 1.0 - e;
	q2(0) = 0.0;
	p1(0) = 0.0;
	p2(0) = sqrt((1.0+e)/(1.0-e));
	H(0) = -1/2.0;

	Euler(p1, p2, q1, q2, H, dt, N);
	Output(p1, p2, q1, q2, H, dt, N);
	
	return 0;
}


void Euler(double* p1, double* p2, double* q1, double* q2, double* H, const double dt, const int N){

	for(int i=1; i<N; i++){
			p1(i) = p1(i-1) - dt * (q1(i-1)/(q1(i-1)*q1(i-1) + q2(i-1)*q2(i-1))^(3/2));
			p2(i) = p2(i-1) - dt * (q2(i-1)/(q1(i-1)*q1(i-1) + q2(i-1)*q2(i-1))^(3/2));
			q1(i) = q1(i-1) + dt * p1(i);
			q2(i) = q2(i-1) + dt * p2(i);
			H(i) = 1/2.0 * (p1(i)*p1(i) + p2(i)*p2(i)) - 1.0/sqrt(q1(i)*q1(i) + q2(i)*q2(i));
		
	}

}

void Output(double* p1, double* p2, double* q1, double* q2, double* H, const double dt, const int N){

	ofstream out("dt_0.0005.dat");
	
	for(int i=0; i<N; i++)
		out << dt << "\t" << p1(i) << "\t" << p2(i) << "\t" << q1(i) << "\t" << q2(i) << "\t" << H(i) << endl;

	out.close();

}


