#pragma once
#include <gtk/gtk.h>
#include <math.h>

bool loadData(const char* filename, double A[1000][10], int* nRows, int* nCols);
void showData(GtkWindow* parent, const char* title, double A[1000][10], int nRows, int nCols);