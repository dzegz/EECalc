#include <gtk/gtk.h>
#include <math.h>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <Eigen/Dense>
#include "GraphFromFile.h"
#include <vector>

using namespace Eigen;
using namespace std::complex_literals;
int izborKvar = 6, izborTranSpoj = 0, izborPotSpoj = 1, izborKondSpoj = 1;
// konstante

const double PI = std::acos(-1);
const double w = 2 * 50 * PI;
std::complex<double> a = std::exp(1i * (2. / 3. * PI));

//combo boxovi - izbornici

GtkWidget* tranMjerenje, * potSpoj, * tranSpoj, *tipKvar, *kondSpoj;

// prozori

GtkCssProvider* cssProvider, *rezCssProvider;
GtkStyleContext* cssStyleContext;
GdkDisplay* display, *rezDisplay;
GtkWidget* grid, * rezGrid;

// regularni unosi

GtkWidget* unosLinR0, * unosLinR1, * unosLinX0, * unosLinX1, * unosLinC0, * unosLinC1, * unosLinL;
GtkWidget* unosGenLinNap, * unosGenFi, * unosGenR, * unosGenX, * unosf_x, * unosf_y, * unosfi_x, * unosfi_y;
GtkWidget* unosPotZ1, * unosPotZ2, * unosPotZ3, * unosPotZg, * unosTranNomSn, *unosTranNomPr, *unosTranNomSek;
GtkWidget* unosTranUk0, * unosTranUk1, * unosTranPcu0, * unosTranPcu1, * unosTranPfe0, * unosTranPfe1, * unosTranI00, * unosTranI01, * unosTranC, * unosTranZgi, * unosTranZgj;
GtkWidget* unosKvarZ, *unosKvarL;
GtkWidget* unosKondC;

// ostalo

GtkWidget* colorBtnX, *colorBtnY;
GtkWidget* area, * btnPng;
bool exportToPng = FALSE;
GtkWidget* win;
char pngPathName[512];
char filePathName[512];
double TextData[1000][10];
bool fault = true;

VectorXcd naponiCvorova(12);
//naponiCvorova.resize(9,1);

MatrixXcd xfmr_model_phi_equiv(Matrix2cd y0, Matrix2cd y1, Matrix2cd y2) {
    MatrixXcd y;
    y.resize(6, 6);

    Matrix3cd T, Ti;
    T << 1, 1, 1,
        1, a, a* a,
        1, a* a, a;
    T *= (1 / 3.);
    Ti = T.inverse();

    y.block(0, 0, 3, 3) = Ti * Vector3cd(y0(0, 0), y1(0, 0), y2(0, 0)).asDiagonal() * T;
    y.block(0, 3, 3, 3) = Ti * Vector3cd(y0(0, 1), y1(0, 1), y2(0, 1)).asDiagonal() * T;
    y.block(3, 0, 3, 3) = Ti * Vector3cd(y0(1, 0), y1(1, 0), y2(1, 0)).asDiagonal() * T;
    y.block(3, 3, 3, 3) = Ti * Vector3cd(y0(1, 1), y1(1, 1), y2(1, 1)).asDiagonal() * T;

    return y;


}
Matrix2cd xfmr_model_symm_12(int  c, std::complex<double> zs1, std::complex<double> ysh1, int x, std::string side) {
    std::complex<double> nc = std::exp(1i * (c * PI / 6.));

    Matrix2cd T, T2, yx;

    if (side.compare("i") == 0) {
        T << 1, 0,
            0, nc;
        T2 << 1, 0,
            0, std::conj(nc);

    }
    else if (side.compare("j") == 0) {
        T << std::pow(nc, -1), 0,
            0, 1;
        T2 << std::pow(std::conj(nc), -1), 0,
            0, 1;
    }

    yx << std::pow(zs1, -1) + ysh1 * 0.5, -std::pow(zs1, -1),
        -std::pow(zs1, -1), std::pow(zs1, -1) + ysh1 * 0.5;

    if (x == 1) return T.conjugate() * yx * T;
    return T2.conjugate() * yx * T2;

}
Matrix2cd xfmr_model_symm_0(int  c, std::complex<double> zs0, std::complex<double> ysh0, std::complex<double> zgi, std::complex<double> zgj, std::string side, double zb1, double  zb2) {
    double zb;
    std::complex<double> nc = std::exp(1i * (c * PI / 6.)), zi;
    Matrix2cd T, yx;
    if (side.compare("i") == 0) {
        zb = zb1;
        T << 1, 0,
            0, nc;

    }
    else if (side.compare("j") == 0) {
        zb = zb2;
        T << std::pow(nc, -1), 0,
            0, 1;
    }

    if (zgj == std::complex<double>(-1, 0)) {
        zgi /= zb;
        zi = zgi * std::complex < double>(3, 0);
        yx << std::pow(zs0 + zi, -1), 0,
            0, 0;
    }
    else if (zgi == std::complex<double>(-1, 0)) {
        zgj /= zb;
        zi = zgj * std::complex < double>(3, 0);
        yx << 0, 0,
            0, std::pow(zs0 + zi, -1);
    }
    else {
        zi = std::complex < double>(3, 0) * zgi / zb + std::complex < double>(3, 0) * zgj / zb;
        yx << std::pow(zs0 + zi, -1) + ysh0 * 0.5, -std::pow(zs0 + zi, -1),
            -std::pow(zs0 + zi, -1), std::pow(zs0 + zi, -1) + ysh0 * 0.5;
    }

    return T.conjugate() * yx * T;

}
Vector2cd xfmr_sequence_params(double Uni, double Unj, double Sn, double uk, double Pcu, double i0, double Pfe, std::string side, double zb1, double zb2) {
    double Un, zb, zsm, rs, xs, yshm, gsh, bsh;

    std::complex<double> zs, ysh;

    if (side.compare("i") == 0) {
        Un = Uni; zb = zb1;
    }
    else if (side.compare("j") == 0) {
        Un = Unj; zb = zb2;
    }

    zsm = uk * Un * Un / (100. * Sn);
    rs = Un * Un * Pcu / (Sn * Sn);
    xs = std::sqrt(zsm * zsm - rs * rs);

    yshm = i0 * Sn / (100 * Un * Un);
    gsh = Pfe / (Un * Un);
    bsh = std::sqrt(yshm * yshm - gsh * gsh);

    return Vector2cd(std::complex<double>(rs, xs) / zb, std::complex<double>(gsh, -bsh) * zb);

}

Matrix3cd breakageModel(std::vector<std::string> type, std::complex<double> z, double zb) {
    //ovdje kontam da je uredu da svaki kvar bude ista impedansa
    //isto kontam da ce biti viska ova obrnuta provjera jer cemo tamo iz ček boksa uzimati kkvarove al nejse
    Matrix3cd f = Matrix3cd::Constant(0);
    for (int i = 0; i < type.size(); i++) {
        if (type.at(i).compare("ab") == 0 || type.at(i).compare("ba") == 0) {
            f(0, 1) -= std::pow(z, -1);
            f(1, 0) -= std::pow(z, -1);
            f(0, 0) += std::pow(z, -1);
            f(1, 1) += std::pow(z, -1);
        }
        else if (type.at(i).compare("ac") == 0 || type.at(i).compare("ca") == 0) {
            f(0, 2) -= std::pow(z, -1);
            f(2, 0) -= std::pow(z, -1);
            f(0, 0) += std::pow(z, -1);
            f(2, 2) += std::pow(z, -1);
        }
        else if (type.at(i).compare("bc") == 0 || type.at(i).compare("cb") == 0) {
            f(2, 1) -= std::pow(z, -1);
            f(1, 2) -= std::pow(z, -1);
            f(1, 1) += std::pow(z, -1);
            f(2, 2) += std::pow(z, -1);
        }
        else if (type.at(i).compare("a0") == 0) {
            f(0, 0) += std::pow(z, -1);
        }
        else if (type.at(i).compare("b0") == 0) {
            f(1, 1) += std::pow(z, -1);
        }
        else if (type.at(i).compare("c0") == 0) {
            f(2, 2) += std::pow(z, -1);
        }
    }
    f *= zb;
    return f;

}
Matrix3cd loadModel(std::string type, std::complex<double> za, std::complex<double> zb, std::complex<double> zc, std::complex<double> zg, double  zbase) {
    //ovdje kontam da treba napraviti preklapanje funckcija jer ne treba parametar zg ako je delta al haj opet kontam svakako se cita tamo iz gui-a pa
    //samo mozemo ignorisati vrijednost, također ova fja ide i za kondenzatorsku bateriju, samosto je zg vazda nula ako je "star"

    za /= zbase;
    zb /= zbase;
    zc /= zbase;
    zg /= zbase;
    Matrix3cd r;
    if (type.compare("delta") == 0) {
        r << std::pow(za, -1) + std::pow(zc, -1), -std::pow(za, -1), -std::pow(zc, -1),
            -std::pow(za, -1), std::pow(za, -1) + std::pow(zb, -1), -std::pow(zb, -1),
            -std::pow(zc, -1), -std::pow(zb, -1), std::pow(zc, -1) + std::pow(zb, -1);
        return r;
    }
    else{
        r << za + zg, zg, zg,
            zg, zb + zg, zg,
            zg, zg, zc + zg;
        return r.inverse();
    }

}
MatrixXcd lineModel(double r0, double x0, double r1, double x1, double c0, double c1, double l, double zb) {
    Matrix3cd T, Ti;
    T << 1, 1, 1,
        1, a, a* a,
        1, a* a, a;
    T *= (1 / 3.);
    Ti = T.inverse();

    std::complex<double> z0(r0, x0);
    std::complex<double> z1(r1, x1);

    auto zs0 = z0 * l; zs0 /= zb;
    auto zs1 = z1 * l; zs1 /= zb;

    auto ysh0 = std::complex<double>(0, w * c0 * l); ysh0 *= zb;
    auto ysh1 = std::complex<double>(0, w * c1 * l); ysh1 *= zb;

    auto yii = Ti * Vector3cd(std::pow(zs0, -1) + 0.5 * ysh0, std::pow(zs1, -1) + 0.5 * ysh1, std::pow(zs1, -1) + 0.5 * ysh1).asDiagonal() * T;
    auto yij = Ti * Vector3cd(-std::pow(zs0, -1), -std::pow(zs1, -1), -std::pow(zs1, -1)).asDiagonal() * T;
    MatrixXcd r; r.resize(6, 6);
    r.block(0, 0, 3, 3) = yii;
    r.block(3, 3, 3, 3) = yii;
    r.block(0, 3, 3, 3) = yij;
    r.block(3, 0, 3, 3) = yij;

    return r;

}
Vector3cd generatorVoltage(double egll, double phase, double vb) {
    return egll / sqrt(3) * Vector3cd(1, a * a, a) * std::exp(1i * (phase * PI / 180.)) / vb;
}
Matrix3cd generatorAdmittance(double rg, double xg, double zb) {
    std::complex<double> zg(rg, xg);
    zg /= zb;
    return Vector3cd(1. / zg, 1. / zg, 1. / zg).asDiagonal();

}
Vector3d get_base_values(double sb, double vb) {

    double ib = sb / vb;
    double zb = vb * vb / sb;

    return Vector3d(vb, ib, zb);
}

const char* getEntryText(GtkWidget* pEntry)
{
    GtkEntryBuffer* buffer = gtk_entry_get_buffer(GTK_ENTRY(pEntry));
    const char* text = gtk_entry_buffer_get_text(buffer);
    return text;
}

void setEntryText(GtkWidget* pEntry, const char* text, int len)
{
    GtkEntryBuffer* buffer = gtk_entry_get_buffer(GTK_ENTRY(pEntry));

    size_t nch = len;
    if (nch <= 0)
        nch = strlen(text);

    gtk_entry_buffer_set_text(buffer, text, (int)nch);
}

void dajIzborKvar(GtkWidget* Izbornik) {
    izborKvar = gtk_combo_box_get_active(GTK_COMBO_BOX(Izbornik));
}

void dajIzborTran(GtkWidget* Izbornik) {
    izborTranSpoj = gtk_combo_box_get_active(GTK_COMBO_BOX(Izbornik));
}

void dajIzborPot(GtkWidget* Izbornik) {
    izborPotSpoj = gtk_combo_box_get_active(GTK_COMBO_BOX(Izbornik));
}

void dajIzborKond(GtkWidget* Izbornik) {
    izborKondSpoj = gtk_combo_box_get_active(GTK_COMBO_BOX(Izbornik));
}

std::complex<double> dajEntryCmplx(GtkWidget* pEntry) {
    const char* val1 = getEntryText(pEntry);
    std::string s = val1;
    std::istringstream is('(' + s + ')');
    std::complex<double> c;
    is >> c;
    return c;
}

double dajEntryDbl(GtkWidget* pEntry)
{
    const char* val1 = getEntryText(pEntry);
    double br = atof(val1);
    return br;
}

int dajEntryInt(GtkWidget* pEntry)
{
    const char* val1 = getEntryText(pEntry);
    int br = atoi(val1);
    return br;
}

void setEntryDbl(GtkWidget* pEntry, double val)
{
    char tmp[64];
    int n = sprintf(tmp, "%.4f", val);
    setEntryText(pEntry, tmp, n);
}

void setEntryCmplx(GtkWidget* pEntry, double real, double imag) {
    char tmp[64];
    int n = sprintf(tmp, "%.4f,%.4f", real, imag);
    setEntryText(pEntry, tmp, n);
}

static void brisi(int width, int height)
{
    //empty = true;
    gtk_widget_queue_draw(area);
    return;
}


static void
on_save_response(GtkNativeDialog* dialog,
    int        response)
{
    if (response == GTK_RESPONSE_ACCEPT)
    {
        GtkFileChooser* chooser = GTK_FILE_CHOOSER(dialog);
        GFile* file = gtk_file_chooser_get_file(chooser);
        //      pngPathName = gtk_file_chooser_get_current_name(GTK_FILE_CHOOSER(dialog));
        const char* path = g_file_get_path(file);
        strcpy(pngPathName, path);
        g_object_unref(file);
        exportToPng = TRUE;
        gtk_widget_queue_draw(area);
    }
    g_object_unref(dialog);
}



//static void savePng(int width, int height)
//{
//    //cairo_surface_t *surface;
//    //surface = cairo_image_surface_create (CAIRO_FORMAT_RGB24, width, height);
//    //cairo_surface_write_to_png (surface, "/Users/dzeni/Desktop/sine.png");
//    GtkWidget *dialog;
//
//      dialog = gtk_file_chooser_dialog_new ("Select file",
//                                            GTK_WINDOW (win),
//                                            GTK_FILE_CHOOSER_ACTION_SAVE,
//                                            "_Cancel", GTK_RESPONSE_CANCEL,
//                                            "_Save", GTK_RESPONSE_OK,
//                                            NULL);
//      gtk_dialog_set_default_response (GTK_DIALOG (dialog), GTK_RESPONSE_OK);
//      gtk_window_set_modal (GTK_WINDOW (dialog), TRUE);
//      gtk_widget_show (dialog);
//
//      g_signal_connect (dialog, "response",
//                        G_CALLBACK (on_save_response),
//                        NULL);
//
//
//}

static void savePng()
{

    //GtkFileChooserNative* dlg;
    //GtkFileChooserAction action = GTK_FILE_CHOOSER_ACTION_SAVE;
    //dlg = gtk_file_chooser_native_new("Snimi graf", (GtkWindow*)win, action, "_Save", "_Cancel");
    //g_signal_connect(dlg, "response",
    //    G_CALLBACK(on_save_response),
    //    NULL);
    //gtk_native_dialog_show((GtkNativeDialog*)dlg);
    GtkFileChooserAction action = GTK_FILE_CHOOSER_ACTION_SAVE;
    GtkFileChooserNative* dlg = gtk_file_chooser_native_new("Snimi graf", GTK_WINDOW(win), action, NULL, NULL);
    //dlg = gtk_file_chooser_new("Snimi graf", (GtkWindow*)win, action, "_Save", "_Cancel");
    g_signal_connect(dlg, "response",
        G_CALLBACK(on_save_response),
        NULL);
    gtk_native_dialog_show((GtkNativeDialog*)dlg);

}


//-------------------------------------------------------------------------------------------------------------------------------------------
// READING FROM THE FILE 
//-------------------------------------------------------------------------------------------------------------------------------------------

static void on_open_response(GtkDialog* dialog, int response)
{
    int nRows, nCols;
    if (response == GTK_RESPONSE_ACCEPT)
    {
        GtkFileChooser* chooser = GTK_FILE_CHOOSER(dialog);
        GFile* file = gtk_file_chooser_get_file(chooser);
        const char* path = g_file_get_path(file);
        if (loadData(path, TextData, &nRows, &nCols))
        {
            g_object_unref(file);
            g_object_unref(dialog);
            showData((GtkWindow*)win, "Podaci iz file-a", TextData, nRows, nCols);
            //
            //varijable

            // prenosna linija
            setEntryDbl(unosLinR0, TextData[4][0]);
            setEntryDbl(unosLinR1, TextData[6][0]);
            setEntryDbl(unosLinX0, TextData[5][0]);
            setEntryDbl(unosLinX1, TextData[7][0]);
            setEntryDbl(unosLinC0, TextData[8][0]);
            setEntryDbl(unosLinC1, TextData[9][0]);
            setEntryDbl(unosLinL, TextData[10][0]);

            //generator
            setEntryDbl(unosGenLinNap, TextData[0][0]);
            setEntryDbl(unosGenFi, TextData[1][0]);
            setEntryDbl(unosGenR, TextData[2][0]);
            setEntryDbl(unosGenX, TextData[3][0]);

            izborTranSpoj = TextData[11][0];
            gtk_combo_box_set_active(GTK_COMBO_BOX(tranSpoj), izborTranSpoj);
            //transformator
            setEntryDbl(unosTranNomSn, TextData[12][0]);
            setEntryDbl(unosTranNomPr, TextData[13][0]);
            setEntryDbl(unosTranNomSek, TextData[14][0]);
            setEntryDbl(unosTranUk0, TextData[15][0]);
            setEntryDbl(unosTranUk1, TextData[19][0]);
            setEntryDbl(unosTranPcu0, TextData[16][0]);
            setEntryDbl(unosTranPcu1, TextData[20][0]);
            setEntryDbl(unosTranPfe0, TextData[18][0]);
            setEntryDbl(unosTranPfe1, TextData[22][0]);
            setEntryDbl(unosTranI00, TextData[17][0]);
            setEntryDbl(unosTranI01, TextData[21][0]);
            setEntryDbl(unosTranC, TextData[23][0]);
            setEntryCmplx(unosTranZgi, TextData[24][0], TextData[25][0]);
            setEntryCmplx(unosTranZgj, TextData[26][0], TextData[27][0]);

            izborPotSpoj = TextData[28][0];
            izborKondSpoj = TextData[37][0];
            gtk_combo_box_set_active(GTK_COMBO_BOX(potSpoj), izborPotSpoj);
            gtk_combo_box_set_active(GTK_COMBO_BOX(kondSpoj), izborKondSpoj);

            //potrosac
            setEntryCmplx(unosPotZ1, TextData[29][0], TextData[30][0]);
            setEntryCmplx(unosPotZ2, TextData[31][0], TextData[32][0]);
            setEntryCmplx(unosPotZ3, TextData[33][0], TextData[34][0]);
            setEntryCmplx(unosPotZg, TextData[35][0], TextData[36][0]);

            izborKvar = TextData[39][0];
            gtk_combo_box_set_active(GTK_COMBO_BOX(tipKvar), izborKvar);
            //ostalo
            setEntryCmplx(unosKvarZ, TextData[40][0], TextData[41][0]);
            setEntryDbl(unosKvarL, TextData[42][0]);
            setEntryDbl(unosKondC, TextData[38][0]);
            //
            return;
        }
        g_object_unref(file);
    }
    g_object_unref(dialog);
}

static void activate_open()
{
    GtkFileChooserNative* dlg;
    GtkFileChooserAction action = GTK_FILE_CHOOSER_ACTION_OPEN;
    dlg = gtk_file_chooser_native_new("Open file", (GtkWindow*) win, action, "_Open", "_Cancel");
    g_signal_connect(dlg, "response",
        G_CALLBACK(on_open_response),
        NULL);
    gtk_native_dialog_show((GtkNativeDialog*)dlg);
}


//-------------------------------------------------------------------------------------------------------------------------------------------
// CALCULATIONS 
//-------------------------------------------------------------------------------------------------------------------------------------------

void calculate_results() {
    //varijable

    // prenosna linija
    double R0 = dajEntryDbl(unosLinR0);
    double R1 = dajEntryDbl(unosLinR1);
    double X0 = dajEntryDbl(unosLinX0);
    double X1 = dajEntryDbl(unosLinX1);
    double C0 = dajEntryDbl(unosLinC0)/1000000;
    double C1 = dajEntryDbl(unosLinC1)/1000000;
    double l = dajEntryDbl(unosLinL)*1000;

    //generator
    double eg_ll = dajEntryDbl(unosGenLinNap)*1000;
    double fi = dajEntryDbl(unosGenFi);
    double rg = dajEntryDbl(unosGenR);
    double xg = dajEntryDbl(unosGenX);

    //transformator
    double Sn = dajEntryDbl(unosTranNomSn)*1000000;
    double Uni = dajEntryDbl(unosTranNomPr)*1000;
    double Unj = dajEntryDbl(unosTranNomSek)*1000;
    double uk0 = dajEntryDbl(unosTranUk0);
    double uk1 = dajEntryDbl(unosTranUk1);
    double pcu0 = dajEntryDbl(unosTranPcu0)*1000;
    double pcu1 = dajEntryDbl(unosTranPcu1)*1000;
    double pfe0 = dajEntryDbl(unosTranPfe0)*1000;
    double pfe1 = dajEntryDbl(unosTranPfe1)*1000;
    double i00 = dajEntryDbl(unosTranI00);
    double i01 = dajEntryDbl(unosTranI01);
    int satni_broj = dajEntryInt(unosTranC);
    std::complex<double> zgi = dajEntryCmplx(unosTranZgi);
    std::complex<double> zgj = dajEntryCmplx(unosTranZgj);

    bool pokus_primar = true;
    if (izborTranSpoj == 1) { //YD
        zgj = std::complex<double>(-1, 0);
    }
    else if (izborTranSpoj == 2) { //DY
        zgi = std::complex<double>(-1, 0);
        pokus_primar = false;
    }

    //potrosac
    std::complex<double> Z1 = dajEntryCmplx(unosPotZ1);
    std::complex<double> Z2 = dajEntryCmplx(unosPotZ2);
    std::complex<double> Z3 = dajEntryCmplx(unosPotZ3);
    std::complex<double> Zg = dajEntryCmplx(unosPotZg);

    //ostalo
    std::complex<double> Zf = dajEntryCmplx(unosKvarZ);
    double udaljenostKvar = dajEntryDbl(unosKvarL)*1000;
    double kondC = dajEntryDbl(unosKondC)/1000000000;

    //-------------------------------------------------------------------------------
    // racun    
    Matrix3cd T, Ti;
    T << 1, 1, 1,
        1, a, a* a,
        1, a* a, a;
    T *= (1 / 3.);
    Ti = T.inverse();
    ////////////////////////////
    double sb = 1e6, v1b = Uni / sqrt(3), v2b = Unj / sqrt(3);
    //base values
    auto bv1 = get_base_values(sb, v1b);
    auto bv2 = get_base_values(sb, v2b);

    ////model generatora
    auto yg = generatorAdmittance(rg, xg, bv1(2));
    auto eg = generatorVoltage(eg_ll, fi, bv1(0));
    auto ig = yg * eg;

    //// linija ako kvar nije prisutan
    auto yl = lineModel(R0,X0,R1,X1,C0,C1,l,bv2(2));

    //linija ako je prisutan kvar
    auto yl1 = lineModel(R0, X0, R1, X1, C0, C1, udaljenostKvar, bv2(2));
    auto yl2 = lineModel(R0, X0, R1, X1, C0, C1, l-udaljenostKvar, bv2(2));

    //potrosac
    Matrix3cd yp, yc;
    if (izborPotSpoj == 1) { //trokut
        yp = loadModel("delta", Z1, Z2, Z3, Zg, bv2(2));
    }
    else yp = loadModel("star", Z1, Z2, Z3, Zg, bv2(2));

    //kondenzator
    if (izborKondSpoj == 1) { //trokut
        yc = loadModel("delta", std::complex<double>(0, std::pow(-w*kondC, -1)), std::complex<double>(0, std::pow(-w * kondC, -1)), std::complex<double>(0, std::pow(-w * kondC, -1)), 0, bv2(2));
    }
    else yc = loadModel("star", std::complex<double>(0, std::pow(-w * kondC, -1)), std::complex<double>(0, std::pow(-w * kondC, -1)), std::complex<double>(0, std::pow(-w * kondC, -1)), 0, bv2(2));

    // test parametara zs i ysh
    Vector2cd p, p2;
    Matrix2cd y, y2, yy, yyy;
    if (pokus_primar) {
        p = xfmr_sequence_params(Uni, Unj, Sn, uk0, pcu0, i00, pfe0, "i", bv1(2), bv2(2)); // zs_x, ysh_x
        //test Yij_0
        y = xfmr_model_symm_0(satni_broj, p(0), p(1), zgi, zgj, "i", bv1(2), bv2(2)); // Ytransf_0
        //test yij_12
        yy = xfmr_model_symm_12(satni_broj, p(0), p(1), 1, "i"); // Ytransf_1
        yyy = xfmr_model_symm_12(satni_broj, p(0), p(1), 2, "i"); // Ytransf_2
    }
    else {
        p2 = xfmr_sequence_params(Uni, Unj, Sn, uk1, pcu1, i01, pfe1, "j", bv1(2), bv2(2));
        //test Yij_0
        y2 = xfmr_model_symm_0(satni_broj, p2(0), p2(1), zgi, zgj, "j", bv1(2), bv2(2));
        //test yij_12
        yy = xfmr_model_symm_12(satni_broj, p2(0), p2(1), 1, "j");
        yyy = xfmr_model_symm_12(satni_broj, p2(0), p2(1), 2, "j");
    }


    //pi ekv trafoa
    auto Yt = xfmr_model_phi_equiv(y, yy, yyy); //Ytransf11, Ytransf12, Ytransf21, Ytransf22

    //funkcija breakageModel za model yf
    //breakageModel(std::vector<std::string> type, std::complex<double> z, double zb) {
    std::vector<std::string> kvarovi;
    fault = true;
    if (izborKvar == 0) {
        kvarovi.push_back(std::string("ab"));
    }
    else if (izborKvar == 1) {
        kvarovi.push_back(std::string("bc"));
    }
    else if (izborKvar == 2) {
        kvarovi.push_back(std::string("ac"));
    }
    else if (izborKvar == 3) {
        kvarovi.push_back(std::string("a0"));
    }
    else if (izborKvar == 4) {
        kvarovi.push_back(std::string("b0"));
    }
    else if (izborKvar == 5) {
        kvarovi.push_back(std::string("c0"));
    }
    else {
        fault = false;
    }

    Matrix3cd yf;
    if(fault) yf = breakageModel(kvarovi, Zf, bv2(2));
    
    Matrix3cd zm = Matrix3cd::Constant(0);
    Vector3cd zv = Vector3cd::Constant(0);
    
    if (!fault) {
        MatrixXcd Y; Y.resize(9, 9);
        VectorXcd I; I.resize(9);
        I << ig, zv, zv;
        Y << yg + Yt.block(0, 0, 3, 3), Yt.block(0, 3, 3, 3), zm,
            Yt.block(3, 0, 3, 3), Yt.block(3, 3, 3, 3) + yl.block(0, 0, 3, 3), yl.block(0, 3, 3, 3),
            zm, yl.block(3, 0, 3, 3), yl.block(3, 3, 3, 3) + yc + yp;
        naponiCvorova = Y.inverse() * I;
        for (int i = 0; i < 3; i++) naponiCvorova(i) *= v1b;
        for (int i = 3; i < 9; i++) naponiCvorova(i) *= v2b;
    }
    else {// if ima kvara
        MatrixXcd Y; Y.resize(12, 12);
        VectorXcd I; I.resize(12);
        I << ig, zv, zv, zv; //ovjde nema struja kod cvora potrošača jer zadajemo preko impedanse beli
        Y << yg + Yt.block(0, 0, 3, 3), Yt.block(0, 3, 3, 3), zm, zm,
            Yt.block(3, 0, 3, 3), Yt.block(3, 3, 3, 3) + yl1.block(0, 0, 3, 3), yl1.block(0, 3, 3, 3), zm,
            zm, yl1.block(3, 0, 3, 3), yl1.block(3, 3, 3, 3) + yl2.block(0, 0, 3, 3) + yf, yl2.block(0, 3, 3, 3),
            zm, zm, yl2.block(3, 0, 3, 3), yl2.block(3, 3, 3, 3) + yc + yp;
        naponiCvorova = Y.inverse() * I;
        for (int i = 0; i < 3; i++) naponiCvorova(i) *= v1b;
        for (int i = 3; i < 12; i++) naponiCvorova(i) *= v2b;
    }
}

//-------------------------------------------------------------------------------------------------------------------------------------------
// RESULTS 
//-------------------------------------------------------------------------------------------------------------------------------------------


static void results_close(GtkButton* btn, gpointer user_data) {
    GtkWindow* tempWin = GTK_WINDOW(user_data);
    gtk_window_destroy(tempWin);
}

void results_show(GtkWidget* p_widget, gpointer user_data) { 
     
    
    //naponiCvorova = VectorXcd::Constant(12,0);
    calculate_results();
    //gtk_widget_queue_draw(area);
    GApplication* app = G_APPLICATION(user_data);
    
    //gtk_button_released(p_widget);
    GtkWidget* rezWin;
    GtkWidget* rezMainVertBox, * rezBoxToolBar;
    rezMainVertBox = gtk_box_new(GTK_ORIENTATION_VERTICAL, 5); //5px izmedju boxova
    rezBoxToolBar = gtk_box_new(GTK_ORIENTATION_HORIZONTAL, 2);

    rezWin = gtk_window_new();
    gtk_window_set_title(GTK_WINDOW(rezWin), "EECalc | Rezultati proracuna");
    gtk_window_set_default_size(GTK_WINDOW(rezWin), 320, 200);


    rezDisplay = gtk_widget_get_display(GTK_WIDGET(rezWin));
    rezCssProvider = gtk_css_provider_new();
    gtk_css_provider_load_from_data(rezCssProvider, "entry { min-height: 0pt; }", -1);
    gtk_style_context_add_provider_for_display(rezDisplay, GTK_STYLE_PROVIDER(rezCssProvider), GTK_STYLE_PROVIDER_PRIORITY_APPLICATION);

    
    gtk_window_set_transient_for(GTK_WINDOW(rezWin), GTK_WINDOW(win));
    //gtk_window_set_position(GTK_WINDOW(rezWin), GTK_WIN_POS_CENTER_ON_PARENT);
    
    //createToolBar(rezBoxToolBar);

    gtk_box_append((GtkBox*)rezMainVertBox, rezBoxToolBar);

    rezGrid = gtk_grid_new();
    gtk_grid_set_column_homogeneous(GTK_GRID(rezGrid), TRUE);
    gtk_window_set_child(GTK_WINDOW(rezWin), rezMainVertBox);
    gtk_box_append((GtkBox*)rezMainVertBox, rezGrid);
    
    // SASTAV
    int i = 0;

    GtkWidget* textRezultati, *btnZatvoriRez;
    textRezultati = gtk_label_new("\nRezultati proracuna\n");
    gtk_grid_attach(GTK_GRID(rezGrid), textRezultati, 0, i, 2, 2);
    i += 2;

    GtkWidget* nc1, * nc2, *nc3, *nc4;
    gchar* display;

    nc1 = gtk_label_new("Napon V1 | | | | ");
    display = g_strdup_printf("Naponi cvora 1 \n  %.2f < %.2f V \n  %.2f < %.2f V \n  %.2f < %.2f V\n",
        std::abs(naponiCvorova(0)), std::arg(naponiCvorova(0))*180/PI, std::abs(naponiCvorova(1)), std::arg(naponiCvorova(1))*180/PI, std::abs(naponiCvorova(2)), std::arg(naponiCvorova(2))*180/PI);
    gtk_label_set_text(GTK_LABEL(nc1), display);
    g_free(display);

    nc2 = gtk_label_new("Napon V2 | | | | ");
    display = g_strdup_printf("Naponi cvora 2 \n  %.2f < %.2f V \n  %.2f < %.2f V \n  %.2f < %.2f V\n",
        std::abs(naponiCvorova(3)), std::arg(naponiCvorova(3)) * 180 / PI, std::abs(naponiCvorova(4)), std::arg(naponiCvorova(4)) * 180 / PI, std::abs(naponiCvorova(5)), std::arg(naponiCvorova(5)) * 180 / PI);
    gtk_label_set_text(GTK_LABEL(nc2), display);
    g_free(display);

    nc3 = gtk_label_new("Napon V3 | | | | ");
    display = g_strdup_printf("Naponi cvora 3 \n  %.2f < %.2f V \n  %.2f < %.2f V \n  %.2f < %.2f V\n",
        std::abs(naponiCvorova(6)), std::arg(naponiCvorova(6)) * 180 / PI, std::abs(naponiCvorova(7)), std::arg(naponiCvorova(7)) * 180 / PI, std::abs(naponiCvorova(8)), std::arg(naponiCvorova(8)) * 180 / PI);
    gtk_label_set_text(GTK_LABEL(nc3), display);
    g_free(display);

    nc4 = gtk_label_new("Napon V4 | | | | ");
    if (fault) {
        display = g_strdup_printf("Naponi cvora 4 \n  %.2f < %.2f V \n  %.2f < %.2f V \n  %.2f < %.2f V\n",
            std::abs(naponiCvorova(9)), std::arg(naponiCvorova(9)) * 180 / PI, std::abs(naponiCvorova(10)), std::arg(naponiCvorova(10)) * 180 / PI, std::abs(naponiCvorova(11)), std::arg(naponiCvorova(11)) * 180 / PI);
        gtk_label_set_text(GTK_LABEL(nc4), display);
        g_free(display);
    }

   
    gtk_grid_attach(GTK_GRID(rezGrid), nc1, 0, i++, 1, 1); i++;
    gtk_grid_attach(GTK_GRID(rezGrid), nc2, 0, i++, 1, 1); i++;
    gtk_grid_attach(GTK_GRID(rezGrid), nc3, 0, i++, 1, 1); i++;
    gtk_grid_attach(GTK_GRID(rezGrid), nc4, 0, i++, 1, 1); i++;
    i++;
    btnZatvoriRez = gtk_button_new_with_label("Zatvori");
    //gtk_window_set_child(GTK_WINDOW(rezWin), btnZatvoriRez); 
    gtk_grid_attach(GTK_GRID(rezGrid), btnZatvoriRez, 0, i++, 2, 2);
    g_signal_connect(btnZatvoriRez, "clicked", G_CALLBACK(results_close),  rezWin);
    
    //------------------
    gtk_window_present(GTK_WINDOW(rezWin));
    gtk_widget_show(rezWin);
    
    if (!fault) {
        gtk_widget_hide(nc4);
    }
    

}

static void
quit_activated(GSimpleAction* action, GVariant* parameter, gpointer user_data)
{
    GApplication* app = G_APPLICATION(user_data);
    g_application_quit(app);
}

static void
test_activated(GSimpleAction* action, GVariant* parameter, gpointer user_data)
{
    printf("Test activated\n");
}
static void
edit_activated(GSimpleAction* action, GVariant* parameter, gpointer user_data)
{
    printf("Edit activated\n");
}
static void
view_activated(GSimpleAction* action, GVariant* parameter, gpointer user_data)
{
    printf("View activated\n");
}


static void createMenu(GApplication* app)
{
    GMenu* menubar = g_menu_new();
    //prvi submenu
    GSimpleAction* act_quit = g_simple_action_new("quit", NULL);
    g_action_map_add_action(G_ACTION_MAP(app), G_ACTION(act_quit));
    g_signal_connect(act_quit, "activate", G_CALLBACK(quit_activated), app);

    GSimpleAction* act_test = g_simple_action_new("test", NULL);
    g_action_map_add_action(G_ACTION_MAP(app), G_ACTION(act_test));
    g_signal_connect(act_test, "activate", G_CALLBACK(test_activated), app);

    GMenuItem* menu_item_menu = g_menu_item_new("Menu", NULL);
    GMenu* menu = g_menu_new();

    GMenuItem* menu_item_test = g_menu_item_new("Test", "app.test");
    g_menu_append_item(menu, menu_item_test);
    g_object_unref(menu_item_test);

    GMenuItem* menu_item_quit = g_menu_item_new("Quit", "app.quit");
    g_menu_append_item(menu, menu_item_quit);
    g_object_unref(menu_item_quit);
    g_menu_item_set_submenu(menu_item_menu, G_MENU_MODEL(menu));
    g_menu_append_item(menubar, menu_item_menu);
    g_object_unref(menu_item_menu);

    //drugi submenu
    GSimpleAction* act_edit = g_simple_action_new("edit", NULL);
    g_action_map_add_action(G_ACTION_MAP(app), G_ACTION(act_edit));
    g_signal_connect(act_edit, "activate", G_CALLBACK(edit_activated), app);

    GSimpleAction* act_view = g_simple_action_new("view", NULL);
    g_action_map_add_action(G_ACTION_MAP(app), G_ACTION(act_view));
    g_signal_connect(act_view, "activate", G_CALLBACK(view_activated), app);

    GMenuItem* menu_item_editSB = g_menu_item_new("Edit SB", NULL);
    GMenu* menu_editSB = g_menu_new();

    GMenuItem* menu_item_edit = g_menu_item_new("Edit", "app.edit");
    g_menu_append_item(menu_editSB, menu_item_edit);
    g_object_unref(menu_item_edit);

    GMenuItem* menu_item_view = g_menu_item_new("Test", "app.view");
    g_menu_append_item(menu_editSB, menu_item_view);
    g_object_unref(menu_item_view);


    g_menu_item_set_submenu(menu_item_editSB, G_MENU_MODEL(menu_editSB));
    g_menu_append_item(menubar, menu_item_editSB);
    g_object_unref(menu_item_editSB);

    gtk_application_set_menubar(GTK_APPLICATION(app), G_MENU_MODEL(menubar));
    gtk_application_window_set_show_menubar(GTK_APPLICATION_WINDOW(win), TRUE);
}

static void comboBoxSelChange()
{
    printf("Selection changed");
}

void scaleChanged(GtkRange* self, gpointer user_data)
{
    double value = gtk_range_get_value(self);
    if (user_data == 0)
        setEntryDbl(unosGenR, value);
    else
        setEntryDbl(unosLinR1, value);
}

static void createToolBar(GtkWidget* boxToolBar)
{
    //Ovdje izmijenite put do slika u Icons folderu!
    GtkWidget* imageExport;
    imageExport = gtk_image_new_from_file("C:/Users/Aldin/Posao/CppProjects/EECalc/Icons/export.png");
    gtk_image_set_icon_size((GtkImage*)imageExport, GTK_ICON_SIZE_LARGE);

    GtkWidget* imageDraw;
    imageDraw = gtk_image_new_from_file("C:/Users/Aldin/Posao/CppProjects/EECalc/Icons/chart.png");
    gtk_image_set_icon_size((GtkImage*)imageDraw, GTK_ICON_SIZE_LARGE);

    GtkWidget* imageDelete;
    imageDelete = gtk_image_new_from_file("C:/Users/Aldin/Posao/CppProjects/EECalc/Icons/delete.png");
    gtk_image_set_icon_size((GtkImage*)imageDelete, GTK_ICON_SIZE_LARGE);

    GtkWidget* imageExit;
    imageExit = gtk_image_new_from_file("C:/Users/Aldin/Posao/CppProjects/EECalc/Icons/openFile.png");
    gtk_image_set_icon_size((GtkImage*)imageExit, GTK_ICON_SIZE_LARGE);


    GtkWidget* btnCrtaj, * btnBrisi, * btnExit, * lblEmpty;
    lblEmpty = gtk_label_new(" ");
    gtk_widget_set_hexpand(lblEmpty, TRUE);

    btnCrtaj = gtk_button_new_with_label("Izvrsi proracun");
    btnBrisi = gtk_button_new_with_label("Brisi");
    btnPng  = gtk_button_new_with_label("Save as .png");
    btnExit = gtk_button_new_with_label("Exit");

    GtkWidget* boxExport, * lblExport, * lblDraw, * lblDelete, * lblExit;
    GtkWidget* boxDraw, * boxDelete, * boxExit;
    boxExport = gtk_box_new(GTK_ORIENTATION_VERTICAL, 3);
    boxDraw = gtk_box_new(GTK_ORIENTATION_VERTICAL, 3);
    boxDelete = gtk_box_new(GTK_ORIENTATION_VERTICAL, 3);
    boxExit = gtk_box_new(GTK_ORIENTATION_VERTICAL, 3);

    lblExport = gtk_label_new("Export to .png");
    lblDraw = gtk_label_new("Izvrsi proracun");
    lblDelete = gtk_label_new("Brisi");
    lblExit = gtk_label_new("Open File");
    gtk_box_append((GtkBox*)boxDraw, imageDraw);
    gtk_box_append((GtkBox*)boxDraw, lblDraw);
    gtk_box_append((GtkBox*)boxDelete, imageDelete);
    gtk_box_append((GtkBox*)boxDelete, lblDelete);
    gtk_box_append((GtkBox*)boxExport, imageExport);
    gtk_box_append((GtkBox*)boxExport, lblExport);
    gtk_box_append((GtkBox*)boxExit, imageExit);
    gtk_box_append((GtkBox*)boxExit, lblExit);

    //merganje
    gtk_button_set_child((GtkButton*)btnCrtaj, boxDraw);
    gtk_button_set_child((GtkButton*)btnBrisi, boxDelete);
    gtk_button_set_child((GtkButton*)btnPng, boxExport);
    gtk_button_set_child((GtkButton*)btnExit, boxExit);

    gtk_box_append((GtkBox*)boxToolBar, btnCrtaj);
    gtk_box_append((GtkBox*)boxToolBar, btnBrisi);
    gtk_box_append((GtkBox*)boxToolBar, btnPng);
    gtk_box_append((GtkBox*)boxToolBar, lblEmpty);
    gtk_box_append((GtkBox*)boxToolBar, btnExit);

    g_signal_connect(btnCrtaj, "clicked", G_CALLBACK(results_show), NULL);
    g_signal_connect(btnBrisi, "clicked", G_CALLBACK(brisi), NULL);
    g_signal_connect(btnExit, "clicked", G_CALLBACK(activate_open), NULL);
    g_signal_connect(btnPng, "clicked", G_CALLBACK(savePng), NULL);
}


//-------------------------------------------------------------------------------------------------------------------------------------------
// MAIN 
//-------------------------------------------------------------------------------------------------------------------------------------------

static void app_activate(GApplication* app, gpointer user_data)
{
    GtkWidget* mainVertBox, * boxToolBar;
    mainVertBox = gtk_box_new(GTK_ORIENTATION_VERTICAL, 5); //5px izmedju boxova
    boxToolBar = gtk_box_new(GTK_ORIENTATION_HORIZONTAL, 50);

    int i = 0;
    g_object_set(gtk_settings_get_default(), "gtk-application-prefer-dark-theme", TRUE, NULL);
    win = gtk_application_window_new(GTK_APPLICATION(app));
    gtk_window_set_title(GTK_WINDOW(win), "EECalc");
    gtk_window_set_default_size(GTK_WINDOW(win), 800, 800);

    display = gtk_widget_get_display(GTK_WIDGET(win));
    cssProvider = gtk_css_provider_new();
    gtk_css_provider_load_from_data(cssProvider, "entry { min-height: 0pt; }", -1);
    gtk_style_context_add_provider_for_display(display, GTK_STYLE_PROVIDER(cssProvider), GTK_STYLE_PROVIDER_PRIORITY_APPLICATION);

    createMenu(app); // izbornik glavni

    createToolBar(boxToolBar);

    gtk_box_append((GtkBox*)mainVertBox, boxToolBar);

    grid = gtk_grid_new();
    gtk_grid_set_column_homogeneous(GTK_GRID(grid), TRUE);
    gtk_window_set_child(GTK_WINDOW(win), mainVertBox);
    gtk_box_append((GtkBox*)mainVertBox, grid);

    //------------------------------------------------------------------------
    // GENERATOR & LINIJA
    //------------------------------------------------------------------------


    GtkWidget* textGenerator, * textLinija;
    textGenerator = gtk_label_new("Podaci o generatoru\n");
    textLinija = gtk_label_new("Podaci o liniji\n");
    gtk_grid_attach(GTK_GRID(grid), textGenerator, 0, i, 2, 2);
    gtk_grid_attach(GTK_GRID(grid), textLinija, 2, i, 2, 2);
    i += 2;  

    GtkWidget* genLinNap, * linR0;
    genLinNap = gtk_label_new("Linijski napon generatora (kV) ");
    linR0 = gtk_label_new("R0 ");
    unosGenLinNap = gtk_entry_new();
    unosLinR0 = gtk_entry_new();
    setEntryText(unosGenLinNap, "22.5", 4);
    setEntryText(unosLinR0, "0.158", 5);
    gtk_grid_attach(GTK_GRID(grid), genLinNap, 0, i, 1, 1);
    gtk_grid_attach(GTK_GRID(grid), unosGenLinNap, 1, i, 1, 1);
    gtk_grid_attach(GTK_GRID(grid), linR0, 2, i, 1, 1);
    gtk_grid_attach(GTK_GRID(grid), unosLinR0, 3, i++, 1, 1);
    
    GtkWidget* genFi, * linX0;
    genFi = gtk_label_new("Fazni pomak ");
    linX0 = gtk_label_new("X0 ");
    unosGenFi = gtk_entry_new();
    unosLinX0 = gtk_entry_new();
    setEntryText(unosGenFi, "0", 4);
    setEntryText(unosLinX0, "0.489", 5);
    gtk_grid_attach(GTK_GRID(grid), genFi, 0, i, 1, 1);
    gtk_grid_attach(GTK_GRID(grid), unosGenFi, 1, i, 1, 1);
    gtk_grid_attach(GTK_GRID(grid), linX0, 2, i, 1, 1);
    gtk_grid_attach(GTK_GRID(grid), unosLinX0, 3, i++, 1, 1);

    GtkWidget* genR, * linR1;
    genR = gtk_label_new("Realni dio impedanse generatora ");
    linR1 = gtk_label_new("R1 ");
    unosGenR = gtk_entry_new();
    unosLinR1 = gtk_entry_new();
    setEntryText(unosGenR, "0.0001", 6);
    setEntryText(unosLinR1, "0.026", 5);
    gtk_grid_attach(GTK_GRID(grid), genR, 0, i, 1, 1);
    gtk_grid_attach(GTK_GRID(grid), unosGenR, 1, i, 1, 1);
    gtk_grid_attach(GTK_GRID(grid), linR1, 2, i, 1, 1);
    gtk_grid_attach(GTK_GRID(grid), unosLinR1, 3, i++, 1, 1);
       
    GtkWidget* genX, * linX1;
    genX = gtk_label_new("Imaginarni dio impedanse generatora ");
    linX1 = gtk_label_new("X1 ");
    unosGenX = gtk_entry_new();
    unosLinX1 = gtk_entry_new();
    setEntryText(unosGenX, "0.0005", 6);
    setEntryText(unosLinX1, "0.188", 5);
    gtk_grid_attach(GTK_GRID(grid), genX, 0, i, 1, 1);
    gtk_grid_attach(GTK_GRID(grid), unosGenX, 1, i, 1, 1);
    gtk_grid_attach(GTK_GRID(grid), linX1, 2, i, 1, 1);
    gtk_grid_attach(GTK_GRID(grid), unosLinX1, 3, i++, 1, 1);

    GtkWidget* linC0, * linC1, *linL;
    linC0 = gtk_label_new("C0 (uF) ");
    linC1 = gtk_label_new("C1 (uF) ");
    linL = gtk_label_new("Duzina linije (km) ");
    unosLinC0 = gtk_entry_new();
    unosLinC1 = gtk_entry_new();
    unosLinL = gtk_entry_new();
    setEntryText(unosLinC0, "0.003", 5);
    setEntryText(unosLinC1, "0.006", 5);
    setEntryText(unosLinL, "50", 2);
    gtk_grid_attach(GTK_GRID(grid), linC0, 2, i, 1, 1);
    gtk_grid_attach(GTK_GRID(grid), unosLinC0, 3, i++, 1, 1);
    gtk_grid_attach(GTK_GRID(grid), linC1, 2, i, 1, 1);
    gtk_grid_attach(GTK_GRID(grid), unosLinC1, 3, i++, 1, 1);
    gtk_grid_attach(GTK_GRID(grid), linL, 2, i, 1, 1);
    gtk_grid_attach(GTK_GRID(grid), unosLinL, 3, i++, 1, 1);

    //------------------------------------------------------------------------
    // TRANSFORMATOR & POTROSAC/KOND BATERIJA
    //------------------------------------------------------------------------

    GtkWidget* textTransformator, * textPotrosac;
    textTransformator = gtk_label_new("\n\nPodaci o transformatoru\n");
    textPotrosac = gtk_label_new("\n\nPodaci o potrosacu\n");
    gtk_grid_attach(GTK_GRID(grid), textTransformator, 0, i, 2, 2);
    gtk_grid_attach(GTK_GRID(grid), textPotrosac, 2, i, 2, 2);
    i += 2;

   
    tranSpoj = gtk_combo_box_text_new();
    gtk_combo_box_text_append(GTK_COMBO_BOX_TEXT(tranSpoj), "1", "YY");
    gtk_combo_box_text_append(GTK_COMBO_BOX_TEXT(tranSpoj), "2", "YD");
    gtk_combo_box_text_append(GTK_COMBO_BOX_TEXT(tranSpoj), "3", "DY");
    gtk_combo_box_set_active(GTK_COMBO_BOX(tranSpoj), 0);
    g_signal_connect(tranSpoj, "changed", G_CALLBACK(dajIzborTran), NULL);
    
    potSpoj = gtk_combo_box_text_new();
    gtk_combo_box_text_append(GTK_COMBO_BOX_TEXT(potSpoj), "1", "Zvijezda");
    gtk_combo_box_text_append(GTK_COMBO_BOX_TEXT(potSpoj), "2", "Trokut");
    gtk_combo_box_set_active(GTK_COMBO_BOX(potSpoj), 1);
    g_signal_connect(potSpoj, "changed", G_CALLBACK(dajIzborPot), NULL);

    GtkWidget* textTranSpoj, * textPotSpoj;
    textTranSpoj = gtk_label_new("Odaberite spoj ");
    textPotSpoj = gtk_label_new("Odaberite spoj ");
    gtk_grid_attach(GTK_GRID(grid), textTranSpoj, 0, i, 1, 1);
    gtk_grid_attach(GTK_GRID(grid), tranSpoj, 1, i, 1, 1);
    gtk_grid_attach(GTK_GRID(grid), textPotSpoj, 2, i, 1, 1);
    gtk_grid_attach(GTK_GRID(grid), potSpoj, 3, i++, 1, 1);

    GtkWidget* tranNomSn, * potZ1;
    tranNomSn = gtk_label_new("Nominalna snaga Sn (MVA) ");
    potZ1 = gtk_label_new("Z1 ");
    unosTranNomSn = gtk_entry_new();
    unosPotZ1 = gtk_entry_new();
    setEntryText(unosTranNomSn, "600", 3);
    setEntryText(unosPotZ1, "300,50", 6);
    gtk_grid_attach(GTK_GRID(grid), tranNomSn, 0, i, 1, 1);
    gtk_grid_attach(GTK_GRID(grid), unosTranNomSn, 1, i, 1, 1);
    gtk_grid_attach(GTK_GRID(grid), potZ1, 2, i, 1, 1);
    gtk_grid_attach(GTK_GRID(grid), unosPotZ1, 3, i++, 1, 1);

    GtkWidget* tranNomPr, * potZ2;
    tranNomPr = gtk_label_new("Nominalni napon primara Uni (kV) ");
    potZ2 = gtk_label_new("Z2 ");
    unosTranNomPr = gtk_entry_new();
    unosPotZ2 = gtk_entry_new();
    setEntryText(unosTranNomPr, "22", 2);
    setEntryText(unosPotZ2, "350,80", 6);
    gtk_grid_attach(GTK_GRID(grid), tranNomPr, 0, i, 1, 1);
    gtk_grid_attach(GTK_GRID(grid), unosTranNomPr, 1, i, 1, 1);
    gtk_grid_attach(GTK_GRID(grid), potZ2, 2, i, 1, 1);
    gtk_grid_attach(GTK_GRID(grid), unosPotZ2, 3, i++, 1, 1);

    GtkWidget* tranNomSek, * potZ3;
    tranNomSek = gtk_label_new("Nominalni napon sekundara Unj (kV) ");
    potZ3 = gtk_label_new("Z3 ");
    unosTranNomSek = gtk_entry_new();
    unosPotZ3 = gtk_entry_new();
    setEntryText(unosTranNomSek, "220", 3);
    setEntryText(unosPotZ3, "350,30", 6);
    gtk_grid_attach(GTK_GRID(grid), tranNomSek, 0, i, 1, 1);
    gtk_grid_attach(GTK_GRID(grid), unosTranNomSek, 1, i, 1, 1);
    gtk_grid_attach(GTK_GRID(grid), potZ3, 2, i, 1, 1);
    gtk_grid_attach(GTK_GRID(grid), unosPotZ3, 3, i++, 1, 1);
    /////////////////////////////////

    GtkWidget* tranUk0, * potZg;
    tranUk0 = gtk_label_new("uk0 (%) ");
    potZg = gtk_label_new("Zg (opcionalno) ");
    unosTranUk0 = gtk_entry_new();
    unosPotZg = gtk_entry_new();
    setEntryText(unosTranUk0, "11", 2);
    setEntryText(unosPotZg, "0,0", 3);
    gtk_grid_attach(GTK_GRID(grid), tranUk0, 0, i, 1, 1);
    gtk_grid_attach(GTK_GRID(grid), unosTranUk0, 1, i, 1, 1);
    gtk_grid_attach(GTK_GRID(grid), potZg, 2, i, 1, 1);
    gtk_grid_attach(GTK_GRID(grid), unosPotZg, 3, i++, 1, 1);
   // i += 2;

    GtkWidget* textKondBat;
    textKondBat = gtk_label_new("\n\nPodaci o kondenzatorskoj bateriji");
    gtk_grid_attach(GTK_GRID(grid), textKondBat, 2, i, 2, 2);
    
    //
    /*kondSpoj = gtk_combo_box_text_new();
    gtk_combo_box_text_append(GTK_COMBO_BOX_TEXT(kondSpoj), "1", "Zvijezda");
    gtk_combo_box_text_append(GTK_COMBO_BOX_TEXT(kondSpoj), "2", "Trokut");
    gtk_combo_box_set_active(GTK_COMBO_BOX(kondSpoj), 0);

    GtkWidget* textKondSpoj;
    textKondSpoj = gtk_label_new("Odaberite spoj ");
    gtk_grid_attach(GTK_GRID(grid), textKondSpoj, 2, i, 1, 1);
    gtk_grid_attach(GTK_GRID(grid), kondSpoj, 3, i++, 1, 1);*/

    //
    GtkWidget* tranPcu0;
    tranPcu0 = gtk_label_new("Pcu0 (kW) ");
    unosTranPcu0 = gtk_entry_new();
    setEntryText(unosTranPcu0, "7", 1);
    gtk_grid_attach(GTK_GRID(grid), tranPcu0, 0, i, 1, 1);
    gtk_grid_attach(GTK_GRID(grid), unosTranPcu0, 1, i++, 1, 1);

    GtkWidget* tranI00;
    tranI00 = gtk_label_new("I00 (%) ");
    unosTranI00 = gtk_entry_new();
    setEntryText(unosTranI00, "0.04", 4);
    gtk_grid_attach(GTK_GRID(grid), tranI00, 0, i, 1, 1);
    gtk_grid_attach(GTK_GRID(grid), unosTranI00, 1, i++, 1, 1);

    GtkWidget* tranPfe0;
    tranPfe0 = gtk_label_new("Pfe0 (kW) ");
    unosTranPfe0 = gtk_entry_new();
    setEntryText(unosTranPfe0, "3", 1);
    gtk_grid_attach(GTK_GRID(grid), tranPfe0, 0, i, 1, 1);
    gtk_grid_attach(GTK_GRID(grid), unosTranPfe0, 1, i++, 1, 1);
    
    kondSpoj = gtk_combo_box_text_new();
    gtk_combo_box_text_append(GTK_COMBO_BOX_TEXT(kondSpoj), "1", "Zvijezda");
    gtk_combo_box_text_append(GTK_COMBO_BOX_TEXT(kondSpoj), "2", "Trokut");
    gtk_combo_box_set_active(GTK_COMBO_BOX(kondSpoj), 1);
    g_signal_connect(kondSpoj, "changed", G_CALLBACK(dajIzborKond), NULL);

    GtkWidget* tranUk1, *textKondSpoj;
    tranUk1 = gtk_label_new("uk1 (%) ");
    textKondSpoj = gtk_label_new("Odaberite spoj ");
    unosTranUk1 = gtk_entry_new();
    setEntryText(unosTranUk1, "10", 2);
    gtk_grid_attach(GTK_GRID(grid), tranUk1, 0, i, 1, 1);
    gtk_grid_attach(GTK_GRID(grid), unosTranUk1, 1, i, 1, 1);
    gtk_grid_attach(GTK_GRID(grid), textKondSpoj, 2, i, 1, 1);
    gtk_grid_attach(GTK_GRID(grid), kondSpoj, 3, i++, 1, 1);

    GtkWidget* tranPcu1, *kondC;
    tranPcu1 = gtk_label_new("Pcu1 (kW) ");
    kondC = gtk_label_new("C (nF) ");
    unosTranPcu1 = gtk_entry_new();
    unosKondC = gtk_entry_new();
    setEntryText(unosTranPcu1, "7", 1);
    setEntryText(unosKondC, "1", 1);
    gtk_grid_attach(GTK_GRID(grid), tranPcu1, 0, i, 1, 1);
    gtk_grid_attach(GTK_GRID(grid), unosTranPcu1, 1, i, 1, 1);
    gtk_grid_attach(GTK_GRID(grid), kondC, 2, i, 1, 1);
    gtk_grid_attach(GTK_GRID(grid), unosKondC, 3, i++, 1, 1);

    GtkWidget* tranI01;
    tranI01 = gtk_label_new("I01 (%) ");
    unosTranI01 = gtk_entry_new();
    setEntryText(unosTranI01, "0.05", 4);
    gtk_grid_attach(GTK_GRID(grid), tranI01, 0, i, 1, 1);
    gtk_grid_attach(GTK_GRID(grid), unosTranI01, 1, i++, 1, 1);

    GtkWidget* tranPfe1;
    tranPfe1 = gtk_label_new("Pfe1 (kW) ");
    unosTranPfe1 = gtk_entry_new();
    setEntryText(unosTranPfe1, "3", 1);
    gtk_grid_attach(GTK_GRID(grid), tranPfe1, 0, i, 1, 1);
    gtk_grid_attach(GTK_GRID(grid), unosTranPfe1, 1, i++, 1, 1);

    /*tranMjerenje = gtk_combo_box_text_new();
    gtk_combo_box_text_append(GTK_COMBO_BOX_TEXT(tranMjerenje), "1", "Primar");
    gtk_combo_box_text_append(GTK_COMBO_BOX_TEXT(tranMjerenje), "2", "Sekundar");
    gtk_combo_box_set_active(GTK_COMBO_BOX(tranMjerenje), 0);*/
    
    /*GtkWidget* textTranMjerenje;
    textTranMjerenje = gtk_label_new("Odaberite stranu mjerenja ");
    gtk_grid_attach(GTK_GRID(grid), textTranMjerenje, 0, i, 1, 1);
    gtk_grid_attach(GTK_GRID(grid), tranMjerenje, 1, i++, 1, 1);*/


    GtkWidget* tranC;
    tranC = gtk_label_new("Satni broj c ");
    unosTranC = gtk_entry_new();
    setEntryText(unosTranC, "0", 1);
    gtk_grid_attach(GTK_GRID(grid), tranC, 0, i, 1, 1);
    gtk_grid_attach(GTK_GRID(grid), unosTranC, 1, i++, 1, 1);

    GtkWidget* tranZgi;
    tranZgi = gtk_label_new("Zgi ");
    unosTranZgi = gtk_entry_new();
    setEntryText(unosTranZgi, "1,0.8", 5);
    gtk_grid_attach(GTK_GRID(grid), tranZgi, 0, i, 1, 1);
    gtk_grid_attach(GTK_GRID(grid), unosTranZgi, 1, i++, 1, 1 );

    GtkWidget* tranZgj;
    tranZgj = gtk_label_new("Zgj ");
    unosTranZgj = gtk_entry_new();
    setEntryText(unosTranZgj, "100,80", 6);
    gtk_grid_attach(GTK_GRID(grid), tranZgj, 0, i, 1, 1);
    gtk_grid_attach(GTK_GRID(grid), unosTranZgj, 1, i++, 1, 1);


    //------------------------------------------------------------------------
    // OSTALO
    //------------------------------------------------------------------------

        GtkWidget* textOstalo;
        textOstalo = gtk_label_new("\n\nOSTALO\n");
        gtk_grid_attach(GTK_GRID(grid), textOstalo, 1, i, 2, 2);
        i += 2;

        tipKvar = gtk_combo_box_text_new();
        gtk_combo_box_text_append(GTK_COMBO_BOX_TEXT(tipKvar), "1", "Faze 1-2");
        gtk_combo_box_text_append(GTK_COMBO_BOX_TEXT(tipKvar), "2", "Faze 2-3");
        gtk_combo_box_text_append(GTK_COMBO_BOX_TEXT(tipKvar), "3", "Faze 1-3");
        gtk_combo_box_text_append(GTK_COMBO_BOX_TEXT(tipKvar), "4", "Faza 1-masa");
        gtk_combo_box_text_append(GTK_COMBO_BOX_TEXT(tipKvar), "5", "Faza 2-masa");
        gtk_combo_box_text_append(GTK_COMBO_BOX_TEXT(tipKvar), "6", "Faza 3-masa");
        gtk_combo_box_text_append(GTK_COMBO_BOX_TEXT(tipKvar), "7", "NEMA KVARA");
        gtk_combo_box_set_active(GTK_COMBO_BOX(tipKvar), 6);

        g_signal_connect(tipKvar, "changed", G_CALLBACK(dajIzborKvar), NULL);
        

        GtkWidget* textKvar;
        textKvar = gtk_label_new("Odaberite tip kvara ");
        gtk_grid_attach(GTK_GRID(grid), textKvar, 0, i, 1, 1);
        gtk_grid_attach(GTK_GRID(grid), tipKvar, 1, i++, 1, 1);

        GtkWidget* kvarZ;
        kvarZ = gtk_label_new("Impedansa kvara ");
        unosKvarZ = gtk_entry_new();
        setEntryText(unosKvarZ, "1,0", 3);
        gtk_grid_attach(GTK_GRID(grid), kvarZ, 0, i, 1, 1);
        gtk_grid_attach(GTK_GRID(grid), unosKvarZ, 1, i++, 1, 1);

        GtkWidget* kvarL;
        kvarL = gtk_label_new("Udaljenost kvara (km) ");
        unosKvarL = gtk_entry_new();
        setEntryText(unosKvarL, "2", 1);
        gtk_grid_attach(GTK_GRID(grid), kvarL, 0, i, 1, 1);
        gtk_grid_attach(GTK_GRID(grid), unosKvarL, 1, i++, 1, 1);

        
    /*

    //klizac
    GtkWidget* scaleX, * scaleY;
    GtkAdjustment* adjustment1, * adjustment2;
    adjustment1 = gtk_adjustment_new(0.5, 0.01, 100.0, 0.01, 15.0, 0.0);
    scaleX = gtk_scale_new(GTK_ORIENTATION_HORIZONTAL, adjustment1);
    gtk_scale_set_digits(GTK_SCALE(scaleX), 2);
    gtk_scale_set_value_pos(GTK_SCALE(scaleX), GTK_POS_TOP);
    gtk_scale_set_draw_value(GTK_SCALE(scaleX), TRUE);
    gtk_grid_attach(GTK_GRID(grid), scaleX, 0, i, 2, 1);
    void* s1 = 0;
    g_signal_connect(scaleX, "value_changed", G_CALLBACK(scaleChanged), s1);


    adjustment2 = gtk_adjustment_new(0.5, 0.01, 100.0, 0.01, 15.0, 0.0);
    scaleY = gtk_scale_new(GTK_ORIENTATION_HORIZONTAL, adjustment2);
    gtk_scale_set_digits(GTK_SCALE(scaleY), 2);
    gtk_scale_set_value_pos(GTK_SCALE(scaleY), GTK_POS_TOP);
    gtk_scale_set_draw_value(GTK_SCALE(scaleY), TRUE);
    gtk_grid_attach(GTK_GRID(grid), scaleY, 2, i++, 2, 1);
    s1 = (void*)1;
    g_signal_connect(scaleY, "value_changed", G_CALLBACK(scaleChanged), s1);

    GtkWidget* f_x, * f_y;
    f_x = gtk_label_new("f_x = ");
    f_y = gtk_label_new("f_y = ");
    unosf_x = gtk_entry_new();
    unosf_y = gtk_entry_new();
    setEntryText(unosf_x, "50", 2);
    setEntryText(unosf_y, "150", 3);
    gtk_grid_attach(GTK_GRID(grid), f_x, 0, i, 1, 1);
    gtk_grid_attach(GTK_GRID(grid), unosf_x, 1, i, 1, 1);
    gtk_grid_attach(GTK_GRID(grid), f_y, 2, i, 1, 1);
    gtk_grid_attach(GTK_GRID(grid), unosf_y, 3, i++, 1, 1);

    GtkWidget* fi_x, * fi_y;
    fi_x = gtk_label_new("fi_x = ");
    fi_y = gtk_label_new("fi_y = ");
    unosfi_x = gtk_entry_new();
    unosfi_y = gtk_entry_new();
    gtk_grid_attach(GTK_GRID(grid), fi_x, 0, i, 1, 1);
    gtk_grid_attach(GTK_GRID(grid), unosfi_x, 1, i, 1, 1);
    gtk_grid_attach(GTK_GRID(grid), fi_y, 2, i, 1, 1);
    gtk_grid_attach(GTK_GRID(grid), unosfi_y, 3, i++, 1, 1);

    GtkWidget* lblClrX, * lblClrY;
    colorBtnX = gtk_color_button_new();
    GdkRGBA clr;
    clr.alpha = 1;
    clr.red = 1;
    clr.green = 0;
    clr.blue = 0;
    gtk_color_chooser_set_rgba(GTK_COLOR_CHOOSER(colorBtnX), &clr);
    colorBtnY = gtk_color_button_new();
    clr.red = 0;
    clr.green = 1;
    gtk_color_chooser_set_rgba(GTK_COLOR_CHOOSER(colorBtnY), &clr);
    lblClrX = gtk_label_new("Boja: ");
    lblClrY = gtk_label_new("Boja: ");
    gtk_grid_attach(GTK_GRID(grid), lblClrX, 0, i, 1, 1);
    gtk_grid_attach(GTK_GRID(grid), colorBtnX, 1, i, 1, 1);
    gtk_grid_attach(GTK_GRID(grid), lblClrY, 2, i, 1, 1);
    gtk_grid_attach(GTK_GRID(grid), colorBtnY, 3, i++, 1, 1);

    area = gtk_drawing_area_new();
    gtk_widget_set_size_request(area, 300, 300);
    gtk_widget_set_vexpand(area, TRUE);
    gtk_drawing_area_set_draw_func(GTK_DRAWING_AREA(area), draw_function, NULL, NULL);

    gtk_grid_attach(GTK_GRID(grid), area, 0, i++, 4, 1);*/
    
    GtkWidget* label;
    label = gtk_label_new("\n\n\n\n\n\n\nCopyright (c) 2022 ETF ");
    gtk_grid_attach(GTK_GRID(grid), label, 1, i++, 2, 1);

    gtk_widget_show(win);

}

int main(int argc, char** argv) {
    GtkApplication* app;
    int stat;
   
    app = gtk_application_new("com.example.GtkApplication", G_APPLICATION_FLAGS_NONE);
    g_signal_connect(app, "activate", G_CALLBACK(app_activate), NULL);
    stat = g_application_run(G_APPLICATION(app), argc, argv);
    g_object_unref(app);
    return stat;
}

