#include "GraphFromFile.h"

bool saveData(const char* filename, double A[1000][10], int nRows, int nCols){
    FILE* dat = fopen(filename, "w");
    if (!dat) {
        printf("Greska prilikom otvaranja/kreiranja datoteke!");
    }
    else {
        fprintf(dat, "43 1\n");
        for (int i = 0; i < nRows; i++) {
            for (int j = 0; j < nCols; j++) {
                fprintf(dat, "%f\n", A[i][j]);
            }
        }
        fclose(dat);
    }
    return true;
}

bool loadData(const char* filename, double A[1000][10], int* nRows, int* nCols)
{
    FILE* fp;
    char c;
    int redovi = 1;

    *nRows = 0;
    *nCols = 0;
    fp = fopen(filename, "r");
    if (fp == NULL)
    {
        printf("Could not open file %s", filename);
        return false;
    }
    int I, J;
    fscanf_s(fp, "%d %d", &I, &J);
    *nRows = I;
    *nCols = J;

    for (int i = 0; i < I; i++)
    {
        for (int j = 0; j < J; j++)
        {
            double val;
            fscanf_s(fp, "%lf", &val);
            A[i][j] = val;
        }
    }
    fclose(fp);
    return true;
}

void showData(GtkWindow* parent, const char* title, double A[1000][10], int nRows, int nCols)
{
    GtkWidget* dialog, * label, * content_area;
    GtkDialogFlags flags;

    // Create the widgets
    flags = GTK_DIALOG_DESTROY_WITH_PARENT;
    dialog = gtk_dialog_new_with_buttons(title, parent, flags, ("_OK"), GTK_RESPONSE_NONE, NULL);
    content_area = gtk_dialog_get_content_area(GTK_DIALOG(dialog));


    // Ensure that the dialog box is destroyed when the user responds

    g_signal_connect_swapped(dialog,
        "response",
        G_CALLBACK(gtk_window_destroy),
        dialog);

    // Add the label, and show everything we’ve added

    char tmp[1024];
    char* pTmp = tmp;
    for (int i = 0; i < nRows; i++)
    {
        if (i == 0) {
            size_t tmpara = sprintf(pTmp, "\n\nGenerator podaci: \n   ");
            pTmp += tmpara; sprintf(pTmp, "\r\n");
        }
        else if (i == 4) {
            size_t tmpara = sprintf(pTmp, "\n\nLinija podaci: \n   ");
            pTmp += tmpara; sprintf(pTmp, "\r\n");
        }
        else if (i == 11) {
            size_t tmpara = sprintf(pTmp, "\n\nTransformator podaci: \n   ");
            pTmp += tmpara; sprintf(pTmp, "\r\n");
        }
        else if (i == 28) {
            size_t tmpara = sprintf(pTmp, "\n\nPotrosac podaci: \n   ");
            pTmp += tmpara; sprintf(pTmp, "\r\n");
        }
        else if (i == 37) {
            size_t tmpara = sprintf(pTmp, "\n\nKond. baterija podaci: \n   ");
            pTmp += tmpara; sprintf(pTmp, "\r\n");
        }
        else if (i == 39) {
            size_t tmpara = sprintf(pTmp, "\n\nOstalo podaci: \n   ");
            pTmp += tmpara; sprintf(pTmp, "\r\n");
        }

        

        for (int j = 0; j < nCols; j++)
        {
            size_t nCh = sprintf(pTmp, " | %.2f |", A[i][j]);
            pTmp += nCh;
        }
        sprintf(pTmp, "\r\n");
    }

    label = gtk_label_new(tmp);

    gtk_box_append(GTK_BOX(content_area), label);
    gtk_widget_show(dialog);
}
