#include "convex_hull.h"
#include "libBMP.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <math.h>

#define VERBOSE 0

const double randMin = -10.0;
const double randMax = 10.0;

bool set_points(int* n) {
    printf("Enter number of points: ");

    char term;
    if (scanf("%i%c", n, &term) != 2 || term != '\n') {
        printf("Error setting points.\n");
        return false;
    }

    return true;
}

bool check_valid(int n) {
    if (n <= 0) {
        printf("Number of points has not been properly set yet.\n");
        return false;
    }

    return true;
}

void print_points(int n, double x[], double y[]) {
    printf("Points:\n");

    for (int i = 0; i < n; i++) {
        printf("\t%i: (%lf, %lf)\n", i, x[i], y[i]);
    }
}

// 1
bool read_points(int n, double x[], double y[]) {
    if (!check_valid(n)) {
        return false;
    }

    printf("Enter points in following format: (x0,y0) (x1,y1) (x2,y2) ... : ");
    char str[256];
    scanf("%255[^\n\r]", str);

    char* token = strtok(str, " ");
    for (int i = 0; i < n; i++) {
        if (token == NULL) {
            printf("Error reading points. Not enough points given to match number of points set.\n");
            return false;
        }

        if (sscanf(token, "(%lf,%lf)", (x + i), (y + i)) != 2) {
            printf("Error reading points. Invalid format.\n");
            return false;
        }

        token = strtok(NULL, " ");
    }

    return true;
}

// 2
bool rand_points(int n, double x[], double y[]) {
    if (!check_valid(n)) {
        return false;
    }

    double div = RAND_MAX / (randMax - randMin);

    for (int i = 0; i < n; i++) {
        x[i] = randMin + (rand() / div);
        y[i] = randMin + (rand() / div);
    }

    return true;
}

// 3
bool display_corners(int m, double x[], double y[], int c[]) {
    printf("Corners:\n");

    for (int i = 0; i < m; i++) {
        printf("\t%i: (%lf, %lf)\n", c[i], x[c[i]], y[c[i]]);
    }

    return true;
}

// 4
bool switch_point(int n, double x[], double y[], int* i_start, int* i_switch) {
    if (!check_valid(n)) {
        return false;
    }

    double minY = y[0];
    double maxX = x[0];
    *i_start = 0;

    double maxY = y[0];
    double minX = x[0];
    *i_switch = 0;

    for (int i = 1; i < n; i++) {
        if (y[i] < minY) {
            minY = y[i];
            maxX = x[i];
            *i_start = i;
        }
        else if (y[i] == minY && x[i] > maxX) {
            maxX = x[i];
            *i_start = i;
        }
        
        if (y[i] > maxY) {
            maxY = y[i];
            minX = x[i];
            *i_switch = i;
        }
        else if (y[i] == maxY && x[i] < minX) {
            minX = x[i];
            *i_switch = i;
        }
    }

    return true;
}

void printOrder(int n, int c[]) {
    printf("Corner order:\n");
    for (int i = 0; i < n; i++) {
        printf("\t%i\n", c[i]);
    }
}

// 5
int hull(int n, double x[], double y[], int c[]) {
    if (!check_valid(n)) {
        return -1;
    }

    // (1)
    for (int i = 0; i < n; i++) {
        c[i] = i;
    }

    // (2)
    if (n == 1) {
        return 1;
    }

    // (3)
    int i_start, i_switch;
    switch_point(n, x, y, &i_start, &i_switch);

    // (4)
    if (i_start != 0) {
        c[0] = c[i_start];
        c[i_start] = 0;
    }

    if (VERBOSE) printOrder(n, c);

    // (5)
    bool up = true;
    for (int m = 1; m <= n; m++) {
        // (9)
        double phiMin = 4.0;
        double lMax = 0.0;
        int jMin = m;

        // (6)
        for (int jIter = m; jIter <= n; jIter++) {
            int j = jIter % n;

            if (j == 0 && m == 1) {
                continue;
            }

            // (7)
            double dx = x[c[j]] - x[c[m-1]];
            double dy = y[c[j]] - y[c[m-1]];

            // (8)
            double phi = atan2(dy, dx);
            double l = hypot(dx, dy);

            // (9)
            if (up && (phi >= 0.0) && phi < phiMin) {
                phiMin = phi;
                lMax = l;
                jMin = j;
            }
            else if (up && phi == phiMin && l > lMax) {
                lMax = l;
                jMin = j;
            }
            else if (!up && (phi <= 0.0) && phi < phiMin) {
                phiMin = phi;
                lMax = l;
                jMin = j;
            }
            else if (!up && phi == phiMin && l > lMax) {
                lMax = l;
                jMin = j;
            }
        }

        if (VERBOSE) printf("m: %i, jMin: %i\n", m, jMin);

        // (10)
        if (x[c[jMin]] == x[c[0]] && y[c[jMin]] == y[c[0]]) {
            return m;
        }

        // (11)
        if (jMin != m) {
            int tmp = c[m];
            c[m] = c[jMin];
            c[jMin] = tmp;
        }

        if (VERBOSE) printOrder(n, c);

        // (12)
        if (c[m] == i_switch) {
            up = false;
        }
    }

    return -1;
}

// 6
const int W = 1000;
const int H = 1000;

void toMath(int X, int Y, double* x, double* y, double xMin, double xMax, double yMin, double yMax) {
    *x = xMin + X * (xMax - xMin) / W;
    *y = yMax - Y * (yMax - yMin) / H;
    return;
}

void toBMP(double x, double y, int* X, int* Y, double xMin, double xMax, double yMin, double yMax) {
    *X = (x - xMin) * W / (xMax - xMin);
    *Y = (yMax - y) * H / (yMax - yMin);
    return;
}

void white(uint32_t* data) {
    for (int i = 0; i < W * H; i++) {
        data[i] = COLOR_WHITE;
    }
    return;
}

void addAxes(uint32_t* data, double xMin, double xMax, double yMin, double yMax) {
    int XOrigin;
    int YOrigin;
    toBMP(0.0, 0.0, &XOrigin, &YOrigin, xMin, xMax, yMin, yMax);

    for (int Y = 0; Y < H; Y++) {
        data[Y * W + XOrigin] = COLOR_BLUE;
    }
    for (int X = 0; X < W; X++) {
        data[YOrigin * W + X] = COLOR_BLUE;
    }
}

void plot_hull(int m, int n, double x[], double y[], int c[]) {
    if (m < 0 || n < 0) {
        printf("Error creating plot.\n");
        return;
    }

    double xMin = x[0];
    double xMax = x[0];
    double yMin = y[0];
    double yMax = y[0];
    for (int i = 0; i < n; i++) {
        if (x[i] < xMin) {
            xMin = x[i];
        }
        else if (x[i] > xMax) {
            xMax = x[i];
        }
        
        if (y[i] < yMin) {
            yMin = y[i];
        }
        else if (y[i] > yMax) {
            yMax = y[i];
        }
    }
    const int borderSteps = 10;
    const int crossSteps = 5;
    double xStepSize = (xMax - xMin) / W;
    double yStepSize = (yMax - yMin) / H;
    xMin -= xStepSize * borderSteps;
    xMax += xStepSize * borderSteps;
    yMin -= yStepSize * borderSteps;
    yMax += yStepSize * borderSteps;

    uint32_t* data = (uint32_t*) malloc(W * H * sizeof(uint32_t));

    white(data);

    int X, Y;
    for (int i = 0; i < n; i++) {
        double xP = x[i];
        double yP = y[i];
        toBMP(xP, yP, &X, &Y, xMin, xMax, yMin, yMax);
        for (int j = -crossSteps; j <= crossSteps; j++) {
            int XCross = X + j;
            int YCross = Y + j;

            if (XCross < 0 || XCross > W || YCross < 0 || YCross > H) {
                continue;
            }

            data[Y * W + XCross] = COLOR_BLUE;
            data[YCross * W + X] = COLOR_BLUE;
        }
    }

    for (int i = 0; i < m; i++) {
        double xP = x[c[i]];
        double yP = y[c[i]];
        toBMP(xP, yP, &X, &Y, xMin, xMax, yMin, yMax);
        for (int j = -crossSteps; j <= crossSteps; j++) {
            int XCross = X + j;
            int YCross = Y + j;

            if (XCross < 0 || XCross > W || YCross < 0 || YCross > H) {
                continue;
            }

            data[Y * W + XCross] = COLOR_RED;
            data[YCross * W + X] = COLOR_RED;
        }

        int next = (i + 1) % m;
        double xPNext = x[c[next]];
        double yPNext = y[c[next]];
        double dx = xP - xPNext;
        double dy = yP - yPNext;

        double l = hypot(dx, dy);
        double diag = hypot(xMax - xMin, yMax - yMin);
        int lineSteps = 1000 * l / diag;

        for (int j = 0; j < lineSteps; j++) {
            toBMP(xPNext + j * dx / lineSteps, yPNext + j * dy / lineSteps, &X, &Y, xMin, xMax, yMin, yMax);
            data[Y * W + X] = COLOR_RED;
        }
    }

    addAxes(data, xMin, xMax, yMin, yMax);

    bmp_create("hull.bmp", data, W, H);
}

int main() {
    srand(time(NULL));

    int n = 0;
    double* x;
    double* y;

    if (set_points(&n)) {
        x = (double*) malloc(n * sizeof(double));
        y = (double*) malloc(n * sizeof(double));

        //read_points(n, x, y);

        //print_points(n, x, y);

        rand_points(n, x, y);

        print_points(n, x, y);

        /*
        int i_start, i_switch;
        switch_point(n, x, y, &i_start, &i_switch);
        printf("%i, %i", i_start, i_switch);
        */

       int* c = (int*) malloc(n * sizeof(int));
       int m = hull(n, x, y, c);

       display_corners(m, x, y, c);

       plot_hull(m, n, x, y, c);
    }
}