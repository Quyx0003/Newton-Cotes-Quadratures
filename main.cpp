#include <iostream>
#include <iomanip>
#include <cmath>

// Midpoint method
void midpoint(double a, double b, int iterations, double (*func)(double)) {
    double prev_integral = 0.0; 
    double prev_prev_integral = 0.0; 

    std::cout << "MidPoint result:" << std::endl;
    std::cout   << std::setw(2)     << "i" 
                << std::setw(15)    << "A(h_i)" 
                << std::setw(20)    << "A(h_(i-1))-A(h_i)" 
                << std::setw(15)    << "alp^k" 
                << std::setw(15)    << "Rich-error" 
                << std::setw(15)    << "f-calc" 
                << std::endl;

    for (int i = 0; i <= iterations; i++) {
        double n = pow(2, i);
        double h = (b - a) / n;
        double sum = 0.0;
        for (int j = 0; j < n; j++) {
            double x_mid = a + (j + 0.5) * h; 
            sum += func(x_mid); 
        }
        double integral = h * sum; 
        std::cout << std::setw(2) << i+1 << std::setw(15) << integral; 

        if (i > 0) {
            double diff = prev_integral - integral; 
            std::cout << std::setw(20) << diff;
            if (i > 1) {
                double diff1 = prev_prev_integral - prev_integral; 
                double alpha_k = diff1 / diff;           
                double richardson = (diff) / (alpha_k - 1.0); 
                std::cout << std::setw(15) << alpha_k; 
                std::cout << std::setw(15) << richardson; 
                std::cout << std::setw(15) << n;  
            }
        }
        
        std::cout << std::endl;
        prev_prev_integral = prev_integral;
        prev_integral = integral;
    }
}

// Trapezoidal method
void trapezoidal(double a, double b, int iterations, double (*func)(double)) {
    double prev_integral = 0.0; 
    double prev_prev_integral = 0.0; 
    
    std::cout << "Trapezoidal result:" << std::endl;
    std::cout   << std::setw(2)     << "i" 
                << std::setw(15)    << "A(h_i)" 
                << std::setw(20)    << "A(h_(i-1))-A(h_i)" 
                << std::setw(15)    << "alp^k" 
                << std::setw(15)    << "Rich-error" 
                << std::setw(15)    << "f-calc" 
                << std::endl;
    
    for (int i = 0; i <= iterations; i++) {
        double n = pow(2, i);
        double h = (b - a) / n;
        double sum = 0.0;
        for (int j = 1; j < n; j++) {
            double x = a + j * h;
            sum += func(x);
        }
        double integral = h * (0.5 * (func(a) + func(b)) + sum); 
        std::cout << std::setw(2) << i+1 << std::setw(15) << integral;
        
        if (i > 0) {
            double diff = prev_integral - integral;
            std::cout << std::setw(20) << diff;
            if (i > 1) {
                double diff1 = prev_prev_integral - prev_integral; 
                double alpha_k = diff1 / diff;           
                double richardson = (diff) / (alpha_k - 1.0); 
                std::cout << std::setw(15) << alpha_k; 
                std::cout << std::setw(15) << richardson; 
                std::cout << std::setw(15) << n + 1;  
            }
        }
        
        std::cout << std::endl;
        prev_prev_integral = prev_integral;
        prev_integral = integral;
    }
}

// Simpson's method
void simpsons(double a, double b, int iterations, double (*func)(double)) {
    double prev_integral = 0.0; 
    double prev_prev_integral = 0.0; 
    
    std::cout << "Simpson's result:" << std::endl;
    std::cout   << std::setw(2)     << "i" 
                << std::setw(15)    << "A(h_i)" 
                << std::setw(20)    << "A(h_(i-1))-A(h_i)" 
                << std::setw(15)    << "alp^k" 
                << std::setw(15)    << "Rich-error" 
                << std::setw(15)    << "f-calc" 
                << std::endl;
    
    for (int i = 0; i <= iterations; i++) {
        double n = pow(2, i);
        double h = (b - a) / n;
        double sum_even = 0.0;
        double sum_odd = 0.0;
        for (int j = 1; j < n; j++) {
            double x = a + j * h;
            if (j % 2 == 0) {
                sum_even += func(x);
            } else {
                sum_odd += func(x);
            }
        }
        double integral = (h / 3.0) * (func(a) + func(b) + 4.0 * sum_odd + 2.0 * sum_even);
        std::cout << std::setw(2) << i+1 << std::setw(15) << integral;
        
        if (i > 0) {
            double diff = prev_integral - integral;
            std::cout << std::setw(20) << diff;
            if (i > 1) {
                double diff1 = prev_prev_integral - prev_integral; 
                double alpha_k = diff1 / diff;           
                double richardson = (diff) / (alpha_k - 1.0); 
                std::cout << std::setw(15) << alpha_k; 
                std::cout << std::setw(15) << richardson; 
                std::cout << std::setw(15) << n + 1;  
            }
        }
        
        std::cout << std::endl;
        prev_prev_integral = prev_integral;
        prev_integral = integral;
    }
}

double f1(double x) {
    return cos(x*x) * exp(-x);
}

double f2(double x) {
    return sqrt(x) * cos(x*x) * exp(-x);
}

double f3(double x) {
    return 1/sqrt(x) * cos(x*x) * exp(-x);
}

double f4(double x) {
    return 1000 * exp(-1/x) * exp(-1/(1-x));
}

int main(int, char**){
    // Exercise 1
    std::cout << "Exercise 1 - cos(x*x)*exp(-x)" << std::endl;
    midpoint(0.0, 1.0, 15, f1);
    std::cout << std::endl;

    std::cout << "Exercise 1 - cos(x*x)*exp(-x)" << std::endl;
    trapezoidal(0.0, 1.0, 15, f1);
    std::cout << std::endl;

    std::cout << "Exercise 1 - cos(x*x)*exp(-x)" << std::endl;
    simpsons(0.0, 1.0, 9, f1);
    std::cout << std::endl;

    // Exercise 2
    std::cout << "Exercise 2 - sqrt(x)*cos(x*x)*exp(-x)" << std::endl;
    simpsons(0.0, 1.0, 22, f2);
    std::cout << std::endl;

    // Exercise 3
    std::cout << "Exercise 3 - 1/sqrt(x)*cos(x*x)*exp(-x)" << std::endl;
    midpoint(0.0, 1.0, 23, f3);
    std::cout << std::endl;

    // Exercise 4
    std::cout << "Exercise 4 - 1000*exp(-1.0/x)*exp(-1.0/(1.0-x))" << std::endl;
    trapezoidal(0.0, 1.0, 5, f4);
    std::cout << std::endl;
}


