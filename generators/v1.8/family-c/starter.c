/* Round 13, C port (Theo 2026-11-04, 0 tokens): focused local search for the starter families C3/C4/C5.
   usage: starter K F SECONDS SEED
   Two rich starters u, v in PG(1,5)^k and one binary starter b in {0,1}^k; the rows are the 6 constants plus
   g(rho^r w), g in PGL(2,5), r in Z_m (m = k - f, last f columns fixed). Cost = (column class, orbit) pairs with
   orbit two-symbol or rich that no starter meets. Prints the starters as JSON when the cost reaches 0. */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <stdint.h>

#define NORB 17
static int ORB[1296], NSYM[NORB];
static int K, F, M, NS = 3, ALPH[3] = {6, 6, 2};
static int NE;              /* entries = classes * m */
static int (*POS)[4];       /* starter positions read by entry e */
static int *ECLS;           /* class of entry e */
static int **TOUCH, *NTOUCH;
static int NC;
static int *CNT, *NEED, *EO[3];
static int W[3][64];
static long HOLES, WC;
static int *WT;
static uint64_t rs;
static inline uint64_t rnd(void) { rs ^= rs << 13; rs ^= rs >> 7; rs ^= rs << 17; return rs; }
static inline int ri(int n) { return (int)(rnd() % (uint64_t)n); }
static inline double rf(void) { return (rnd() >> 11) * (1.0 / 9007199254740992.0); }

static int inv5(int d) { int r = 1; for (int i = 0; i < 3; i++) r = r * d % 5; return r; }
static void build_orbits(void) {
    int G[120][6], ng = 0;
    for (int a = 0; a < 5; a++) for (int b = 0; b < 5; b++) for (int c = 0; c < 5; c++) for (int d = 0; d < 5; d++) {
        if (((a * d - b * c) % 5 + 5) % 5 == 0) continue;
        int img[6];
        for (int x = 0; x < 6; x++) {
            int num, den;
            if (x == 5) { num = a; den = c; } else { num = (a * x + b) % 5; den = (c * x + d) % 5; }
            img[x] = den == 0 ? 5 : num * inv5(den) % 5;
        }
        int dup = 0;
        for (int h = 0; h < ng && !dup; h++) if (!memcmp(G[h], img, sizeof img)) dup = 1;
        if (!dup) memcpy(G[ng++], img, sizeof img);
    }
    if (ng != 120) { fprintf(stderr, "PGL size %d\n", ng); exit(2); }
    for (int i = 0; i < 1296; i++) ORB[i] = -1;
    int no = 0;
    for (int code = 0; code < 1296; code++) {
        if (ORB[code] >= 0) continue;
        int t[4] = {code / 216, code / 36 % 6, code / 6 % 6, code % 6};
        int s = 0, seen[6] = {0};
        for (int j = 0; j < 4; j++) if (!seen[t[j]]) { seen[t[j]] = 1; s++; }
        for (int g = 0; g < 120; g++) ORB[G[g][t[0]] * 216 + G[g][t[1]] * 36 + G[g][t[2]] * 6 + G[g][t[3]]] = no;
        NSYM[no++] = s;
    }
    if (no != NORB) { fprintf(stderr, "orbits %d\n", no); exit(2); }
}

static int img(int c, int r) { return c < M ? (c + r) % M : c; }
static void build_classes(void) {
    int nsub = 0, cap = 1 << 20;
    char *seen = calloc((size_t)1 << 24, 1);
    POS = malloc(sizeof(int[4]) * cap); ECLS = malloc(sizeof(int) * cap);
    NC = 0; NE = 0;
    for (int a = 0; a < K; a++) for (int b = a + 1; b < K; b++) for (int c = b + 1; c < K; c++) for (int d = c + 1; d < K; d++) {
        int key = ((a * 64 + b) * 64 + c) * 64 + d;
        if (seen[key]) continue;
        for (int r = 0; r < M; r++) {
            int s[4] = {img(a, r), img(b, r), img(c, r), img(d, r)};
            for (int i = 1; i < 4; i++) for (int j = i; j > 0 && s[j - 1] > s[j]; j--) { int t = s[j]; s[j] = s[j - 1]; s[j - 1] = t; }
            seen[((s[0] * 64 + s[1]) * 64 + s[2]) * 64 + s[3]] = 1;
        }
        int S[4] = {a, b, c, d};
        for (int r = 0; r < M; r++) {         /* row rho^r w read at column S[j] is w[(S[j]-r) mod m] */
            for (int j = 0; j < 4; j++) POS[NE][j] = S[j] < M ? ((S[j] - r) % M + M) % M : S[j];
            ECLS[NE++] = NC;
        }
        NC++; nsub++;
    }
    free(seen);
    NTOUCH = calloc(K, sizeof(int)); TOUCH = malloc(K * sizeof(int *));
    for (int e = 0; e < NE; e++) for (int j = 0; j < 4; j++) NTOUCH[POS[e][j]]++;
    for (int c = 0; c < K; c++) { TOUCH[c] = malloc(NTOUCH[c] * sizeof(int)); NTOUCH[c] = 0; }
    for (int e = 0; e < NE; e++) for (int j = 0; j < 4; j++) TOUCH[POS[e][j]][NTOUCH[POS[e][j]]++] = e;
    CNT = calloc(NC * NORB, sizeof(int)); NEED = calloc(NC * NORB, sizeof(int));
    for (int c = 0; c < NC; c++) for (int o = 0; o < NORB; o++) NEED[c * NORB + o] = NSYM[o] >= 2;
    for (int i = 0; i < NS; i++) EO[i] = malloc(NE * sizeof(int));
    WT = malloc(NC * NORB * sizeof(int)); for (int x = 0; x < NC * NORB; x++) WT[x] = 1;
}

static inline int eorb(int i, int e) { int *p = POS[e]; return ORB[W[i][p[0]] * 216 + W[i][p[1]] * 36 + W[i][p[2]] * 6 + W[i][p[3]]]; }
static void load(void) {
    memset(CNT, 0, NC * NORB * sizeof(int));
    for (int i = 0; i < NS; i++) for (int e = 0; e < NE; e++) { EO[i][e] = eorb(i, e); CNT[ECLS[e] * NORB + EO[i][e]]++; }
    HOLES = 0; WC = 0;
    for (int x = 0; x < NC * NORB; x++) if (NEED[x] && !CNT[x]) { HOLES++; WC += WT[x]; }
}
static long set_entry(int i, int c, int a) {          /* returns change in weighted cost */
    long h0 = WC;
    W[i][c] = a;
    int *t = TOUCH[c], n = NTOUCH[c];
    for (int q = 0; q < n; q++) {
        int e = t[q];
        int on = eorb(i, e), oo = EO[i][e];
        if (on == oo) continue;
        int b = ECLS[e] * NORB;
        if (--CNT[b + oo] == 0 && NEED[b + oo]) { HOLES++; WC += WT[b + oo]; }
        if (CNT[b + on]++ == 0 && NEED[b + on]) { HOLES--; WC -= WT[b + on]; }
        EO[i][e] = on;
    }
    return WC - h0;
}
static void randomise(void) { for (int i = 0; i < NS; i++) for (int c = 0; c < K; c++) W[i][c] = ri(ALPH[i]); W[2][0] = 0; }

static void dump(const char *tag, long holes, double secs, long long it) {
    printf("{\"tag\": \"%s\", \"k\": %d, \"f\": %d, \"holes\": %ld, \"seconds\": %.1f, \"moves\": %lld, \"starters\": [", tag, K, F, holes, secs, it);
    for (int i = 0; i < NS; i++) { printf(i ? ", [" : "["); for (int c = 0; c < K; c++) printf(c ? ", %d" : "%d", W[i][c]); printf("]"); }
    printf("]}\n"); fflush(stdout);
}

int main(int argc, char **argv) {
    if (argc < 5) { fprintf(stderr, "usage: starter K F SECONDS SEED\n"); return 1; }
    K = atoi(argv[1]); F = atoi(argv[2]); M = K - F;
    double budget = atof(argv[3]); rs = strtoull(argv[4], 0, 10) * 2654435761ULL + 88172645463325252ULL;
    if (K > 60 || F < 0 || M < 2) return 1;
    build_orbits(); build_classes();
    int *holes = malloc(NC * NORB * sizeof(int));
    time_t t0 = time(0);
    long best = 1 << 30; int BW[3][64];
    long long it = 0, restarts = 0;
    while (difftime(time(0), t0) < budget) {
        for (int x = 0; x < NC * NORB; x++) WT[x] = 1;
        randomise(); load(); restarts++;
        long rbest = HOLES; int RW[3][64]; memcpy(RW, W, sizeof W);
        long long since = 0;
        while (HOLES > 0 && since < 40000000LL) {
            it++; since++;
            if ((it & 0xFFFF) == 0 && difftime(time(0), t0) >= budget) break;
            int i, c, a;
            if (rf() < 0.7) {                              /* focused: pick a hole, set one entry so that tuple lands in it */
                int nh = 0;
                for (int x = 0; x < NC * NORB && nh < 64; x++) if (NEED[x] && !CNT[x]) holes[nh++] = x;
                int x = holes[ri(nh)], cls = x / NORB, o = x % NORB;
                i = NSYM[o] >= 3 ? ri(2) : ri(3);
                int e = cls * M + ri(M);                   /* entries of a class are contiguous, M per class */
                int cand[6], nc = 0;
                int j = ri(4); c = POS[e][j];
                int old = W[i][c];
                for (int v = 0; v < ALPH[i]; v++) { W[i][c] = v; if (eorb(i, e) == o) cand[nc++] = v; }
                W[i][c] = old;
                if (!nc) continue;
                a = cand[ri(nc)];
                if (a == old) continue;
            } else {
                i = ri(3); c = ri(K);
                a = ri(ALPH[i] - 1); if (a >= W[i][c]) a++;
            }
            if (i == 2 && c == 0) continue;
            int old = W[i][c];
            long d = set_entry(i, c, a);
            if (d > 0 && rf() >= 0.01) set_entry(i, c, old);
            if (since % 20000 == 0) { for (int x = 0; x < NC * NORB; x++) if (NEED[x] && !CNT[x]) { WT[x]++; WC++; } }
            else if (HOLES < rbest) { rbest = HOLES; memcpy(RW, W, sizeof W); since = 0; }
        }
        if (HOLES == 0) { rbest = 0; memcpy(RW, W, sizeof W); }
        if (rbest < best) { best = rbest; memcpy(BW, RW, sizeof W); memcpy(W, BW, sizeof W); dump(best ? "best" : "witness", best, difftime(time(0), t0), it); }
        if (best == 0) break;
    }
    memcpy(W, BW, sizeof W); load();
    for (int x = 0; x < NC * NORB; x++) if (NEED[x] && !CNT[x]) {
        int e = (x / NORB) * M; fprintf(stderr, "hole class=%d cols=%d,%d,%d,%d orbit=%d nsym=%d\n", x / NORB, POS[e][0], POS[e][1], POS[e][2], POS[e][3], x % NORB, NSYM[x % NORB]); }
    fprintf(stderr, "k=%d f=%d classes=%d best=%ld moves=%lld restarts=%lld\n", K, F, NC, best, it, restarts);
    return 0;
}
