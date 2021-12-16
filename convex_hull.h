#include <stdbool.h>

// 'bool' statt 'void' als Rückgabewert, um Information über mögliche Errors zu liefern

bool read_points(int n, double x[], double y[]);

bool rand_points(int n, double x[], double y[]);

bool display_corners(int m, double x[], double y[], int c[]);

bool switch_point(int n, double x[], double y[], int* i_start, int* i_switch);

int hull(int n, double x[], double y[], int c[]);

void plot_hull(int m, int n, double x[], double y[], int c[]);