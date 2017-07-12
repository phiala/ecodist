#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* .C calls */
extern void bcdistc(void *, void *, void *, void *);
extern void bootstrap(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void jabs(void *, void *, void *, void *);
extern void jfirst(void *, void *, void *, void *);
extern void jpres(void *, void *, void *, void *);
extern void jsec(void *, void *, void *, void *);
extern void mrmperm(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void newpermone(void *, void *, void *, void *, void *, void *, void *, void *);
extern void newpermtwo(void *, void *, void *, void *, void *, void *, void *, void *);
extern void pdiff(void *, void *, void *, void *);
extern void permpart(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void permute(void *, void *, void *, void *, void *, void *, void *, void *);
extern void psum(void *, void *, void *, void *);
extern void xpermute(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void xpermpart(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);

static const R_CMethodDef CEntries[] = {
    {"bcdistc",    (DL_FUNC) &bcdistc,     4},
    {"bootstrap",  (DL_FUNC) &bootstrap,  11},
    {"jabs",       (DL_FUNC) &jabs,        4},
    {"jfirst",     (DL_FUNC) &jfirst,      4},
    {"jpres",      (DL_FUNC) &jpres,       4},
    {"jsec",       (DL_FUNC) &jsec,        4},
    {"mrmperm",    (DL_FUNC) &mrmperm,    15},
    {"newpermone", (DL_FUNC) &newpermone,  8},
    {"newpermtwo", (DL_FUNC) &newpermtwo,  8},
    {"pdiff",      (DL_FUNC) &pdiff,       4},
    {"permpart",   (DL_FUNC) &permpart,   13},
    {"permute",    (DL_FUNC) &permute,     8},
    {"psum",       (DL_FUNC) &psum,        4},
    {"xpermute",   (DL_FUNC) &xpermute,   10},
    {"xpermpart",  (DL_FUNC) &xpermpart,  12},
    {NULL, NULL, 0}
};

void R_init_ecodist(DllInfo *dll)
{
    R_registerRoutines(dll, CEntries, NULL, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
