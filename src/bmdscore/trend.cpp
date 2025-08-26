include "trend.h";

double approximate_pi(int terms) {
    double pi = 0.0;
    double sign = 1.0;

    for (int i = 0; i < terms; ++i) {
        pi += sign / (2.0 * i + 1.0);
        sign *= -1;
    }
    return 4.0 * pi;

}
