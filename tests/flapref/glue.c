/* Test glue for the upstream FLAP reference (xfoil-c tests only): the F77 externals it calls
 * (SEVAL, DEVAL, D2VAL, SINVRT, SCALC, SEGSPL, GEOPAR) are routed to the module procedures of
 * the reference library libxfoil_fortran.so (the same code the C port was verified against),
 * so the test isolates FLAP/SSS/SCHECK/INSIDE/GETXYF; plotting and prompts are stubs. */
#include <stddef.h>

typedef double R;
extern R __m_spline_MOD_seval(R *, R *, R *, R *, int *);
extern R __m_spline_MOD_deval(R *, R *, R *, R *, int *);
extern R __m_spline_MOD_d2val(R *, R *, R *, R *, int *);
extern void __m_spline_MOD_sinvrt(R *, R *, R *, R *, R *, int *);
extern void __m_spline_MOD_scalc(R *, R *, R *, int *);
extern void __m_spline_MOD_segspl(R *, R *, R *, int *);
extern void __m_xgeom_MOD_geopar(R *, R *, R *, R *, R *, int *, R *, R *, R *, R *, R *, R *, R *, R *, R *,
                                 R *, R *, R *, R *, R *, R *, R *);

R seval_(R *ss, R *x, R *xs, R *s, int *n) { return __m_spline_MOD_seval(ss, x, xs, s, n); }
R deval_(R *ss, R *x, R *xs, R *s, int *n) { return __m_spline_MOD_deval(ss, x, xs, s, n); }
R d2val_(R *ss, R *x, R *xs, R *s, int *n) { return __m_spline_MOD_d2val(ss, x, xs, s, n); }
void sinvrt_(R *si, R *xi, R *x, R *xs, R *s, int *n) { __m_spline_MOD_sinvrt(si, xi, x, xs, s, n); }
void scalc_(R *x, R *y, R *s, int *n) { __m_spline_MOD_scalc(x, y, s, n); }
void segspl_(R *x, R *xs, R *s, int *n) { __m_spline_MOD_segspl(x, xs, s, n); }
void geopar_(R *x, R *xp, R *y, R *yp, R *s, int *n, R *t, R *sle, R *chord, R *area, R *radle, R *angte,
             R *ei11a, R *ei22a, R *apx1a, R *apx2a, R *ei11t, R *ei22t, R *apx1t, R *apx2t, R *thick, R *cambr)
{
    __m_xgeom_MOD_geopar(x, xp, y, yp, s, n, t, sle, chord, area, radle, angte, ei11a, ei22a, apx1a, apx2a,
                         ei11t, ei22t, apx1t, apx2t, thick, cambr);
}

/* plotting and prompts: not used by the tests (LPLOT = .TRUE., all inputs given) */
void pltini_(void) {}
void getcolor_(int *icol) { *icol = 0; }
void newcolorname_(const char *name, size_t len) { (void)name; (void)len; }
void newcolor_(int *icol) { (void)icol; }
void plot_(R *x, R *y, int *ipen) { (void)x; (void)y; (void)ipen; }
void plsymb_(R *x, R *y, R *h, int *isym, R *ang, int *ninc) { (void)x; (void)y; (void)h; (void)isym; (void)ang; (void)ninc; }
void pltair_(void) {}
void plnewp_(const char *name, size_t len) { (void)name; (void)len; }
void askr_(const char *prompt, R *v, size_t len) { (void)prompt; (void)v; (void)len; }
