#include <stdio.h>
#include <math.h>
double fun(double x){
        double out;
        out = x*exp(-x)*sin(x);
        return out;
}

int main(){
        int i,j;
	double a = 1;
	double b = 3;
	int n = 500;
	double h = (b-a)/n;
	double integralf=0;	
        for(double x=a;x<3;x=x+h){
		integralf = integralf + fun(x)*h;
	}
        printf("integral=");
        printf("%f\n",integralf);
}
