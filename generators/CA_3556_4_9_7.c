/* Generaliser harness for the coset-union construction of covering arrays CA(N; 4, k, q), q prime.
 *
 * Fixed part (not evolved): size loop, 4-set enumeration, expansion of the construction into rows,
 * greedy completion of every uncovered 4-tuple, greedy pruning of redundant rows, output.
 * Evolved part: the block between EVOLVE-BLOCK-START and EVOLVE-BLOCK-END, which must define
 *     void construct(int q, int k, int table_N, double seconds, Construction *out);
 *
 * GRID="q,k,N;q,k,N;..." (N = best known table size), SIZE_SECONDS = budget of construct() per size.
 * Each size runs in a forked child under a hard alarm, so a hung construct() costs only that size.
 * Output per size: "# size q k", then the rows, then "# end q k" (written by the parent).
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <unistd.h>
#include <signal.h>
#include <sys/time.h>
#include <sys/wait.h>

#define MAXK 16
#define MAXD 3
#define MAXM 2000
#define MAXNT 1820          /* C(16,4) */
#define MAXROWS 40000

typedef struct {
  int d;                       /* dimension of the linear code C (0..3); d = 0 means plain rows   */
  int G[MAXD][MAXK];           /* generator matrix of C over GF(q), entries 0..q-1                */
  int m;                       /* number of translates (cosets s_i + C)                          */
  unsigned char S[MAXM][MAXK]; /* the translates s_i, entries 0..q-1                             */
} Construction;               /* rows = { s_i + u G : i < m, u in GF(q)^d }, N = m q^d before fix-up */

/* ---- utilities available to the evolved block ---- */
static int Q, K, NT;                 /* current size, NT = C(K,4)                    */
static int cols[MAXNT][4];           /* the 4-subsets of columns, lexicographic      */
static int smask[MAXNT];             /* bitmask of the columns of each 4-subset      */
static unsigned long long rs = 0x9E3779B97F4A7C15ULL;
static inline unsigned long long rnd(void) { rs ^= rs << 13; rs ^= rs >> 7; rs ^= rs << 17; return rs; }
static double now(void) { struct timeval tv; gettimeofday(&tv, 0); return tv.tv_sec + tv.tv_usec * 1e-6; }
static int ipow(int b, int e) { int r = 1; while (e--) r *= b; return r; }

// EVOLVE-BLOCK-START
/* Seed strategy (port of evolve/coset/genq.c, the generator behind CA(825;4,7,5) and CA(1050;4,8,5)).
 * Code: a random [k,2] code over GF(q) of rank 2 on every 4 columns. On a 4-set S the code is a
 * 2-dim subspace of GF(q)^4, so a translate s is known on S by its label = (b1.s_S, b2.s_S) for a
 * basis b1, b2 of the null space: q^2 labels per 4-set. The union of the m cosets is a CA iff every
 * 4-set sees all q^2 labels. Weighted tabu search over the translates (moves change 1 or 2 entries
 * of one translate), m started just under table_N / q^2 and lowered after every success. */
static int inv_mod(int a, int q) { int r = 1; for (int e = q - 2; e > 0; e--) r = r * a % q; return r; }

/* null space of the d x 4 matrix P over GF(q): basis rows in B, returns the rank of P */
static int null_space(int q, int d, int P[MAXD][4], int B[4][4]) {
  int M[MAXD][4], r = 0, pc[MAXD];
  memcpy(M, P, sizeof M);
  for (int c = 0; c < 4 && r < d; c++) {
    int p = -1;
    for (int i = r; i < d; i++) if (M[i][c]) { p = i; break; }
    if (p < 0) continue;
    for (int j = 0; j < 4; j++) { int t = M[p][j]; M[p][j] = M[r][j]; M[r][j] = t; }
    int iv = inv_mod(M[r][c], q);
    for (int j = 0; j < 4; j++) M[r][j] = M[r][j] * iv % q;
    for (int i = 0; i < d; i++) if (i != r && M[i][c]) { int f = M[i][c]; for (int j = 0; j < 4; j++) M[i][j] = ((M[i][j] - f * M[r][j]) % q + q) % q; }
    pc[r++] = c;
  }
  int isp[4] = { 0, 0, 0, 0 }, nb = 0;
  for (int i = 0; i < r; i++) isp[pc[i]] = 1;
  for (int f = 0; f < 4; f++) if (!isp[f]) {
    int v[4] = { 0, 0, 0, 0 }; v[f] = 1;
    for (int i = 0; i < r; i++) v[pc[i]] = (q - M[i][f]) % q;
    memcpy(B[nb++], v, sizeof v);
  }
  return r;
}

static int GQ, GD, GE, GLS, GM_;                  /* q, d, e = 4 - d, labels per 4-set, m */
static int (*mulc)[4][2][7];                      /* mulc[t][j][comp][v] = b_comp[j] * v mod q (q <= 7) */
static int *gcnt, *gul, *gup, gnu; static unsigned *gW;
static unsigned char (*SH)[MAXK];
static int (*glab)[MAXNT];
static long long (*gtabu)[MAXK];

static int col_ok(int q, int d, Construction *out, int c) {
  for (int a = 0; a < c; a++) for (int b = a + 1; b < c; b++) for (int e = b + 1; e < c; e++) {
    int P[MAXD][4] = { { 0 } }, B[4][4] = { { 0 } }, cc[4] = { a, b, e, c };
    for (int i = 0; i < d; i++) for (int j = 0; j < 4; j++) P[i][j] = out->G[i][cc[j]];
    if (null_space(q, d, P, B) < d) return 0;
  }
  return 1;
}

static int setup_code(int q, int k, int d, Construction *out, double tend) {
  for (long tries = 0; now() < tend; tries++) {
    int okc = 1;
    for (int c = 0; c < k && okc; c++) {
      int got = 0;
      for (int r = 0; r < 2000 && !got; r++) {
        for (int i = 0; i < d; i++) out->G[i][c] = (int)(rnd() % q);
        got = col_ok(q, d, out, c);
      }
      if (!got) okc = 0;
    }
    if (!okc) continue;
    int ok = 1;
    for (int t = 0; t < NT && ok; t++) {
      int P[MAXD][4] = { { 0 } }, B[4][4] = { { 0 } };
      for (int i = 0; i < d; i++) for (int j = 0; j < 4; j++) P[i][j] = out->G[i][cols[t][j]];
      if (null_space(q, d, P, B) < d) { ok = 0; break; }
      for (int j = 0; j < 4; j++) for (int comp = 0; comp < 2; comp++) for (int v = 0; v < q; v++)
        mulc[t][j][comp][v] = comp < 4 - d ? B[comp][j] * v % q : 0;
    }
    if (ok) return 1;
  }
  return 0;
}

static inline int lab_of(const unsigned char *s, int t) {
  int l0 = 0, l1 = 0;
  for (int j = 0; j < 4; j++) { int v = s[cols[t][j]]; l0 += mulc[t][j][0][v]; l1 += mulc[t][j][1][v]; }
  return GE == 2 ? (l0 % GQ) * GQ + l1 % GQ : l0 % GQ;
}
static inline void ginc(int k) { if (gcnt[k]++ == 0) { int p = gup[k]; int l = gul[--gnu]; gul[p] = l; gup[l] = p; gup[k] = -1; } }
static inline void gdec(int k) { if (--gcnt[k] == 0) { gup[k] = gnu; gul[gnu++] = k; } }

static void gbuild(void) {
  int NTL = NT * GLS; gnu = 0;
  memset(gcnt, 0, sizeof(int) * NTL);
  for (int i = 0; i < GM_; i++) for (int t = 0; t < NT; t++) { glab[i][t] = lab_of(SH[i], t); gcnt[t * GLS + glab[i][t]]++; }
  for (int k = 0; k < NTL; k++) { if (!gcnt[k]) { gup[k] = gnu; gul[gnu++] = k; } else gup[k] = -1; }
}

/* tabu search on the translates until every (4-set, label) is hit or tend; returns the best count of misses */
static int gsearch(double tend) {
  long long it = 0; int bestnu = gnu, NTL = NT * GLS, q = GQ;
  for (int k = 0; k < NTL; k++) gW[k] = 1;
  memset(gtabu, 0, sizeof(gtabu[0]) * GM_);
  while (gnu > 0) {
    if ((it & 255) == 0 && now() > tend) break;
    it++;
    if (it % 20000 == 0) for (int k = 0; k < NTL; k++) gW[k] = 1 + (gW[k] - 1) * 3 / 4;
    int k = gul[rnd() % gnu], t = k / GLS, L = k % GLS;
    int L0 = GE == 2 ? L / q : L, L1 = GE == 2 ? L % q : 0;
    long long bw = 1LL << 62; int bq = -1, bc = 0, bcm = 0; unsigned char brow[MAXK];
    for (int i = 0; i < GM_; i++) {
      int a = glab[i][t], a0 = GE == 2 ? a / q : a, a1 = GE == 2 ? a % q : 0;
      int D0 = (L0 - a0 + q) % q, D1 = (L1 - a1 + q) % q;
      for (int j1 = 0; j1 < 4; j1++) for (int j2 = j1; j2 < 4; j2++) {
        int single = j2 == j1;
        for (int x1 = 1; x1 < q; x1++) for (int x2 = single ? 0 : 1; x2 < (single ? 1 : q); x2++) {
          int s0 = mulc[t][j1][0][x1] + (single ? 0 : mulc[t][j2][0][x2]);
          int s1 = mulc[t][j1][1][x1] + (single ? 0 : mulc[t][j2][1][x2]);
          if (s0 % q != D0 || s1 % q != D1) continue;
          unsigned char w[MAXK]; memcpy(w, SH[i], K);
          int c1 = cols[t][j1], c2 = cols[t][j2], cm = 1 << c1;
          w[c1] = (unsigned char)((w[c1] + x1) % q);
          if (!single) { w[c2] = (unsigned char)((w[c2] + x2) % q); cm |= 1 << c2; }
          long long dw = 0; int du = 0;
          for (int t2 = 0; t2 < NT; t2++) if (smask[t2] & cm) {
            int aa = t2 * GLS + glab[i][t2], bb = t2 * GLS + lab_of(w, t2);
            if (aa == bb) continue;
            if (gcnt[aa] == 1) { dw += gW[aa]; du++; }
            if (gcnt[bb] == 0) { dw -= gW[bb]; du--; }
          }
          int ist = 0;
          for (int c = 0; c < K; c++) if (((cm >> c) & 1) && gtabu[i][c] > it) { ist = 1; break; }
          if (ist && gnu + du >= bestnu) continue;
          if (dw < bw) { bw = dw; bq = i; bcm = cm; bc = 1; memcpy(brow, w, K); }
          else if (dw == bw && rnd() % (++bc) == 0) { bq = i; bcm = cm; memcpy(brow, w, K); }
        }
      }
    }
    if (bq < 0) continue;
    if (bw >= 0) for (int u = 0; u < gnu; u++) gW[gul[u]]++;
    for (int t2 = 0; t2 < NT; t2++) if (smask[t2] & bcm) {
      int nl = lab_of(brow, t2);
      if (nl != glab[bq][t2]) { gdec(t2 * GLS + glab[bq][t2]); ginc(t2 * GLS + nl); glab[bq][t2] = nl; }
    }
    memcpy(SH[bq], brow, K);
    for (int c = 0; c < K; c++) if ((bcm >> c) & 1) gtabu[bq][c] = it + 3 + (long long)(rnd() % 6);
    if (gnu < bestnu) bestnu = gnu;
  }
  return gnu < bestnu ? gnu : bestnu;
}

/* drop translates one by one (fewest uniquely-covered labels first) until GM_ == target; gcnt/glab must be valid */
static void drop_rows(int target) {
  while (GM_ > target) {
    int bi = 0, bl = 1 << 30;
    for (int i = 0; i < GM_; i++) {
      int l = 0;
      for (int t = 0; t < NT; t++) if (gcnt[t * GLS + glab[i][t]] == 1) l++;
      if (l < bl || (l == bl && rnd() % 2)) { bl = l; bi = i; }
    }
    for (int t = 0; t < NT; t++) gcnt[t * GLS + glab[bi][t]]--;
    memcpy(SH[bi], SH[GM_ - 1], MAXK); memcpy(glab[bi], glab[GM_ - 1], sizeof(glab[0]));
    GM_--;
  }
}

static int fill_mulc(int q, int d, Construction *out) {
  for (int t = 0; t < NT; t++) {
    int P[MAXD][4] = { { 0 } }, B[4][4] = { { 0 } };
    for (int i = 0; i < d; i++) for (int j = 0; j < 4; j++) P[i][j] = out->G[i][cols[t][j]];
    if (null_space(q, d, P, B) < d) return 0;
    for (int j = 0; j < 4; j++) for (int comp = 0; comp < 2; comp++) for (int v = 0; v < q; v++)
      mulc[t][j][comp][v] = comp < 4 - d ? B[comp][j] * v % q : 0;
  }
  return 1;
}
static int set_and_check(int q, int d, Construction *out, int m) {
  if (!fill_mulc(q, d, out)) return 0;
  GM_ = m; gbuild(); return gnu == 0;
}
static int sB[MAXNT][2][4], sph[MAXNT][4][2], sidx[MAXNT];
static int tl[MAXK][MAXNT], tj[MAXK][MAXNT], tn[MAXK];
static int sW[4][MAXK], SX[MAXM][4];
static int nB[MAXNT][2][4], nphG[MAXNT][4][2], nph[MAXNT][2], nidx[MAXNT];

static int sa_init(int q, int k, int d, int n, Construction *out, unsigned char *cov, int *pwn, double tend) {
  int e = 4 - d;
  if (!setup_code(q, k, d, out, tend)) return -1;
  for (int a = 0; a < n; a++) for (int c = 0; c < k; c++) sW[a][c] = (int)(rnd() % q);
  int unc = 0;
  for (int t = 0; t < NT; t++) {
    int P[MAXD][4] = { { 0 } }, B[4][4] = { { 0 } };
    for (int i = 0; i < d; i++) for (int j = 0; j < 4; j++) P[i][j] = out->G[i][cols[t][j]];
    null_space(q, d, P, B);
    for (int comp = 0; comp < e; comp++) for (int j = 0; j < 4; j++) sB[t][comp][j] = B[comp][j];
    sidx[t] = 0;
    for (int a = 0; a < n; a++) for (int comp = 0; comp < e; comp++) {
      int s = 0; for (int j = 0; j < 4; j++) s += sB[t][comp][j] * sW[a][cols[t][j]];
      sph[t][a][comp] = s % q; sidx[t] += sph[t][a][comp] * pwn[a * e + comp];
    }
    if (!cov[sidx[t]]) unc++;
  }
  return unc;
}
static inline int sa_acc(int delta) {
  if (delta <= 0) return 1;
  if (delta == 1) return rnd() % 1000 < 30;
  if (delta == 2) return rnd() % 1000 < 2;
  return 0;
}

/* structured translates: X (template in GF(q)^n) times W (n x k), code G searched jointly */
static int structured(int q, int k, int d, Construction *out, double pend, long best_est) {
  int e = 4 - d, P = ipow(q, d), GLSl = ipow(q, e);
  for (int c = 0; c < k; c++) tn[c] = 0;
  for (int t = 0; t < NT; t++) for (int j = 0; j < 4; j++) { int c = cols[t][j]; tl[c][tn[c]] = t; tj[c][tn[c]++] = j; }
  if (d == 3 && k <= q + 1 && (long)q * P < best_est) {
    int w[MAXK];
    for (int c = 0; c < k; c++) {
      if (c < q) { int x = c; out->G[0][c] = 1; out->G[1][c] = x * x % q; out->G[2][c] = x * x * x % q; w[c] = x; }
      else { out->G[0][c] = 0; out->G[1][c] = 0; out->G[2][c] = 1; w[c] = 0; }
    }
    for (int i = 0; i < q; i++) for (int c = 0; c < k; c++) SH[i][c] = (unsigned char)(i * w[c] % q);
    if (set_and_check(q, d, out, q)) return q;
  }
  int ntp = d == 3 ? 3 : 1;
  for (int tp = 0; tp < ntp; tp++) {
    int n, mX = 0;
    if (d == 3) {
      n = tp + 2;
      memset(SX[0], 0, sizeof SX[0]); mX = 1;
      for (int a = 0; a < n; a++) for (int t = 1; t < q; t++) { memset(SX[mX], 0, sizeof SX[0]); SX[mX][a] = t; mX++; }
    } else {
      n = 3;
      for (int x = 0; x < q * q * q; x++) {
        int v0 = x % q, v1 = x / q % q, v2 = x / q / q;
        if (v0 && v1 && v2) continue;
        SX[mX][0] = v0; SX[mX][1] = v1; SX[mX][2] = v2; SX[mX][3] = 0; mX++;
      }
    }
    if ((long)mX * P >= best_est || mX > MAXM) continue;
    if (now() >= pend) break;
    double tend = now() + (pend - now()) * (tp == 0 ? 0.6 : 0.4);
    int pwn[8], sz = ipow(q, n * e);
    for (int i = 0; i < n * e; i++) pwn[i] = ipow(q, i);
    unsigned char *cov = malloc(sz);
    for (int idx = 0; idx < sz; idx++) {
      int ph[4][2] = { { 0 } }, mark[64] = { 0 }, cm = 0;
      for (int a = 0; a < n; a++) for (int comp = 0; comp < e; comp++) ph[a][comp] = idx / pwn[a * e + comp] % q;
      for (int i = 0; i < mX && cm < GLSl; i++) {
        int l0 = 0, l1 = 0;
        for (int a = 0; a < n; a++) { l0 += ph[a][0] * SX[i][a]; l1 += ph[a][1] * SX[i][a]; }
        int L = e == 2 ? (l0 % q) * q + l1 % q : l0 % q;
        if (!mark[L]) { mark[L] = 1; cm++; }
      }
      cov[idx] = cm == GLSl;
    }
    long long it = 0, lastimp = 0;
    int unc = sa_init(q, k, d, n, out, cov, pwn, tend), bestu = unc;
    while (unc > 0) {
      if ((it & 63) == 0 && now() > tend) break;
      it++;
      if (it - lastimp > 300000) { unc = sa_init(q, k, d, n, out, cov, pwn, tend); if (unc < 0) break; bestu = unc; lastimp = it; continue; }
      int t = (int)(rnd() % NT); while (cov[sidx[t]]) t = (t + 1) % NT;
      int j = (int)(rnd() % 4), c = cols[t][j], delta = 0, L = tn[c];
      if (rnd() % 5 == 0) {
        int g[MAXD], nz = 0;
        for (int i = 0; i < d; i++) { g[i] = (int)(rnd() % q); nz |= g[i]; }
        if (!nz) continue;
        int bad = 0;
        for (int ii = 0; ii < L && !bad; ii++) {
          int t2 = tl[c][ii], j2 = tj[c][ii];
          int Pm[MAXD][4] = { { 0 } }, B[4][4] = { { 0 } };
          for (int i = 0; i < d; i++) for (int jj = 0; jj < 4; jj++) Pm[i][jj] = jj == j2 ? g[i] : out->G[i][cols[t2][jj]];
          if (null_space(q, d, Pm, B) < d) { bad = 1; break; }
          int ni = 0;
          for (int comp = 0; comp < e; comp++) for (int jj = 0; jj < 4; jj++) nB[ii][comp][jj] = B[comp][jj];
          for (int a = 0; a < n; a++) for (int comp = 0; comp < e; comp++) {
            int s = 0; for (int jj = 0; jj < 4; jj++) s += B[comp][jj] * sW[a][cols[t2][jj]];
            nphG[ii][a][comp] = s % q; ni += nphG[ii][a][comp] * pwn[a * e + comp];
          }
          nidx[ii] = ni; delta += (!cov[ni]) - (!cov[sidx[t2]]);
        }
        if (bad || !sa_acc(delta)) continue;
        for (int i = 0; i < d; i++) out->G[i][c] = g[i];
        for (int ii = 0; ii < L; ii++) {
          int t2 = tl[c][ii];
          memcpy(sB[t2], nB[ii], sizeof sB[0]); memcpy(sph[t2], nphG[ii], sizeof sph[0]); sidx[t2] = nidx[ii];
        }
      } else {
        int a = (int)(rnd() % n), old = sW[a][c], v = (old + 1 + (int)(rnd() % (q - 1))) % q, dv = (v - old + q) % q;
        for (int ii = 0; ii < L; ii++) {
          int t2 = tl[c][ii], j2 = tj[c][ii], ni = sidx[t2];
          for (int comp = 0; comp < e; comp++) {
            int np = (sph[t2][a][comp] + sB[t2][comp][j2] * dv) % q;
            nph[ii][comp] = np; ni += (np - sph[t2][a][comp]) * pwn[a * e + comp];
          }
          nidx[ii] = ni; delta += (!cov[ni]) - (!cov[sidx[t2]]);
        }
        if (!sa_acc(delta)) continue;
        sW[a][c] = v;
        for (int ii = 0; ii < L; ii++) {
          int t2 = tl[c][ii];
          for (int comp = 0; comp < e; comp++) sph[t2][a][comp] = nph[ii][comp];
          sidx[t2] = nidx[ii];
        }
      }
      unc += delta; if (unc < bestu) { bestu = unc; lastimp = it; }
    }
    free(cov);
    fprintf(stderr, "  structured d=%d n=%d m=%d: unc=%d\n", d, n, mX, unc);
    if (unc == 0) {
      for (int i = 0; i < mX; i++) for (int c = 0; c < k; c++) {
        int s = 0; for (int a = 0; a < n; a++) s += SX[i][a] * sW[a][c];
        SH[i][c] = (unsigned char)(s % q);
      }
      if (set_and_check(q, d, out, mX)) return mX;
    }
  }
  return 0;
}

static unsigned char BS[MAXM][MAXK], LS[MAXM][MAXK];
static int BG[MAXD][MAXK];

/* ---- fast label search for any d in 0..3 ---- */
#define FLMAX 625
static int FLS, Fm, Fnu, *Fcnt, *Ful, *Fup, *FW;
static unsigned short *FADD, Fneg[FLMAX], (*Fdl)[4][7], (*Flab)[MAXNT];
static unsigned char *Ffix, EXS[MAXM][MAXK];
static long Fbest = 1L << 40; static int Fbd = 2, Fbm = 0;
static int LG[MAXD][MAXK], Ld, Lm;

static int f_setup(int q, int d, int G[MAXD][MAXK]) {
  int e = 4 - d; FLS = ipow(q, e);
  for (int a = 0; a < FLS; a++) {
    int x = a, r = 0, pw = 1;
    for (int i = 0; i < e; i++) { r += ((q - x % q) % q) * pw; x /= q; pw *= q; }
    Fneg[a] = (unsigned short)r;
    for (int b = 0; b < FLS; b++) {
      int u = a, w = b; r = 0; pw = 1;
      for (int i = 0; i < e; i++) { r += ((u % q + w % q) % q) * pw; u /= q; w /= q; pw *= q; }
      FADD[a * FLS + b] = (unsigned short)r;
    }
  }
  memset(Ffix, 0, (size_t)NT * FLS * 4);
  for (int c = 0; c < MAXK; c++) tn[c] = 0;
  for (int t = 0; t < NT; t++) {
    int P[MAXD][4] = { { 0 } }, B[4][4] = { { 0 } };
    for (int i = 0; i < d; i++) for (int j = 0; j < 4; j++) P[i][j] = G[i][cols[t][j]];
    if (null_space(q, d, P, B) < d) return 0;
    for (int j = 0; j < 4; j++) {
      int c = cols[t][j]; tl[c][tn[c]] = t; tj[c][tn[c]++] = j;
      for (int dv = 0; dv < q; dv++) {
        int l = 0, pw = 1;
        for (int i = 0; i < e; i++) { l += (B[i][j] * dv % q) * pw; pw *= q; }
        Fdl[t][j][dv] = (unsigned short)l;
        if (dv && l) Ffix[(t * FLS + l) * 4 + j] = (unsigned char)dv;
      }
    }
  }
  return 1;
}
static inline void f_inc(int x) { if (Fcnt[x]++ == 0) { int p = Fup[x], l = Ful[--Fnu]; Ful[p] = l; Fup[l] = p; Fup[x] = -1; } }
static inline void f_dec(int x) { if (--Fcnt[x] == 0) { Fup[x] = Fnu; Ful[Fnu++] = x; } }
static void f_build(void) {
  int NTL = NT * FLS; Fnu = 0;
  memset(Fcnt, 0, sizeof(int) * NTL);
  for (int i = 0; i < Fm; i++) for (int t = 0; t < NT; t++) {
    int l = 0;
    for (int j = 0; j < 4; j++) l = FADD[l * FLS + Fdl[t][j][SH[i][cols[t][j]]]];
    Flab[i][t] = (unsigned short)l; Fcnt[t * FLS + l]++;
  }
  for (int x = 0; x < NTL; x++) { if (!Fcnt[x]) { Fup[x] = Fnu; Ful[Fnu++] = x; } else Fup[x] = -1; }
}
/* focused weighted repair: pick a hole, sample translates one entry away from it, take the best */
static int f_solve(int q, double te) {
  long long it = 0; int NTL = NT * FLS;
  for (int x = 0; x < NTL; x++) FW[x] = 1;
  while (Fnu > 0) {
    if ((it & 255) == 0 && now() > te) return 0;
    it++;
    if (it % 200000 == 0) for (int x = 0; x < NTL; x++) FW[x] = 1 + (FW[x] - 1) * 3 / 4;
    int u = Ful[rnd() % Fnu], t = u / FLS, L = u % FLS;
    int st = (int)(rnd() % Fm), found = 0, bi = -1, bc = -1, bdv = 0; long bw = 1L << 60;
    for (int n = 0; n < Fm && found < 6; n++) {
      int i = st + n; if (i >= Fm) i -= Fm;
      int D = FADD[L * FLS + Fneg[Flab[i][t]]];
      for (int jj = 0; jj < 4; jj++) {
        int j = (jj + n) & 3, dv = Ffix[(t * FLS + D) * 4 + j];
        if (!dv) continue;
        int c = cols[t][j]; long w = 0;
        for (int ii = 0; ii < tn[c]; ii++) {
          int t2 = tl[c][ii], dl = Fdl[t2][tj[c][ii]][dv];
          if (!dl) continue;
          int a = Flab[i][t2], x = t2 * FLS + a, y = t2 * FLS + FADD[a * FLS + dl];
          if (Fcnt[x] == 1) w += FW[x];
          if (Fcnt[y] == 0) w -= FW[y];
        }
        found++;
        if (w < bw || (w == bw && rnd() % found == 0)) { bw = w; bi = i; bc = c; bdv = dv; }
        break;
      }
    }
    if (bi < 0) { FW[u]++; continue; }
    if (bw > 0) { FW[u]++; if (rnd() % 32) continue; }
    for (int ii = 0; ii < tn[bc]; ii++) {
      int t2 = tl[bc][ii], dl = Fdl[t2][tj[bc][ii]][bdv];
      if (!dl) continue;
      int a = Flab[bi][t2], b = FADD[a * FLS + dl];
      f_inc(t2 * FLS + b); f_dec(t2 * FLS + a);
      Flab[bi][t2] = (unsigned short)b;
    }
    SH[bi][bc] = (unsigned char)((SH[bi][bc] + bdv) % q);
  }
  return 1;
}
static void f_drop(void) {
  int bi = 0, bl = 1 << 30, ties = 0;
  for (int i = 0; i < Fm; i++) {
    int l = 0;
    for (int t = 0; t < NT; t++) if (Fcnt[t * FLS + Flab[i][t]] == 1) l++;
    if (l < bl) { bl = l; bi = i; ties = 1; }
    else if (l == bl && rnd() % (++ties) == 0) bi = i;
  }
  for (int t = 0; t < NT; t++) f_dec(t * FLS + Flab[bi][t]);
  memcpy(SH[bi], SH[Fm - 1], MAXK); memcpy(Flab[bi], Flab[Fm - 1], sizeof(Flab[0]));
  Fm--;
}
static void f_record(int q, int k, int d, int G[MAXD][MAXK]) {
  long N = (long)Fm * ipow(q, d);
  if (N < Fbest || (N == Fbest && d < Fbd)) {
    Fbest = N; Fbd = d; Fbm = Fm; memcpy(BS, SH, sizeof(SH[0]) * Fm);
    for (int i = 0; i < d; i++) for (int c = 0; c < k; c++) BG[i][c] = G[i][c];
    fprintf(stderr, "  d=%d m=%d N=%ld\n", d, Fm, N);
  }
}
/* SH solves at Fm: drop the weakest translate and repair until te; LS/Lm always hold the last solution */
static void f_descend(int q, int k, int d, int G[MAXD][MAXK], double te) {
  memcpy(LS, SH, sizeof(SH[0]) * Fm); Lm = Fm;
  while (Fm > FLS && now() < te) {
    f_drop();
    if (!f_solve(q, te)) return;
    memcpy(LS, SH, sizeof(SH[0]) * Fm); Lm = Fm;
    f_record(q, k, d, G);
  }
}
static void f_scratch(int q, int k, int ds, int table_N, Construction *out, double te) {
  int G[MAXD][MAXK];
  memset(G, 0, sizeof G);
  if (ds == 1) { for (int c = 0; c < k; c++) G[0][c] = 1; }
  else if (ds > 1) {
    if (!setup_code(q, k, ds, out, now() + (te - now()) * 0.2)) return;
    for (int i = 0; i < ds; i++) for (int c = 0; c < k; c++) G[i][c] = out->G[i][c];
  }
  if (!f_setup(q, ds, G)) return;
  int P = ipow(q, ds), mmax = 40000 / P; if (mmax > MAXM) mmax = MAXM;
  int m = (int)(1.1 * table_N / P) + 1; if (m < FLS) m = FLS; if (m > mmax) m = mmax;
  for (;;) {
    Fm = m;
    for (int i = 0; i < m; i++) for (int c = 0; c < k; c++) SH[i][c] = (unsigned char)(rnd() % q);
    f_build();
    if (f_solve(q, now() + (te - now()) * 0.35)) { f_record(q, k, ds, G); f_descend(q, k, ds, G, te); return; }
    if (Fbest == 1L << 40) {            /* nothing valid yet: keep the near miss for the harness to complete */
      Fbd = ds; Fbm = Fm; memcpy(BS, SH, sizeof(SH[0]) * Fm);
      for (int i = 0; i < ds; i++) for (int c = 0; c < k; c++) BG[i][c] = G[i][c];
    }
    if (now() >= te || m >= mmax) return;
    m = m * 5 / 4 + 1; if (m > mmax) m = mmax;
  }
}
/* RS code over GF(q+1), q+1 = 2^r: projective codewords, logs of the nonzero entries, zeros are free cells.
 * Scalars act as the all-ones vector on the logs, so this is a coset union with d = 1 over Z_q. */
static int fusion_seed(int q, int k) {
  int Q = q + 1, poly = Q == 4 ? 7 : Q == 8 ? 11 : 0;
  if (!poly || k > Q + 1) return 0;
  int ex[16], lg[16] = { 0 }, x = 1, n = 0;
  for (int i = 0; i < Q - 1; i++) { ex[i] = x; lg[x] = i; x <<= 1; if (x & Q) x ^= poly; }
  #define GFM(a, b) (((a) && (b)) ? ex[(lg[a] + lg[b]) % (Q - 1)] : 0)
  for (int a0 = 0; a0 < Q; a0++) for (int a1 = 0; a1 < Q; a1++) for (int a2 = 0; a2 < Q; a2++) for (int a3 = 0; a3 < Q; a3++) {
    int f = a0 ? a0 : a1 ? a1 : a2 ? a2 : a3;
    if (f != 1) continue;
    if (n >= MAXM) return 0;
    for (int c = 0; c < k; c++) {
      int v;
      if (c == Q) v = a3;
      else { int p = c; v = a0; v ^= GFM(a1, p); p = GFM(p, c); v ^= GFM(a2, p); p = GFM(p, c); v ^= GFM(a3, p); }
      SH[n][c] = (unsigned char)(v ? lg[v] : (int)(rnd() % q));
    }
    n++;
  }
  #undef GFM
  return n;
}
static int can_split(int q) { return Ld > 0 && (long)Lm * q <= MAXM && ipow(q, 5 - Ld) <= FLMAX; }

void construct(int q, int k, int table_N, double seconds, Construction *out) {
  double t0 = now(), tend = t0 + seconds - 0.5;
  mulc = calloc(MAXNT, sizeof *mulc);
  gcnt = malloc(sizeof(int) * NT * q * q); gul = malloc(sizeof(int) * NT * q * q); gup = malloc(sizeof(int) * NT * q * q);
  gW = malloc(sizeof(unsigned) * NT * q * q);
  SH = calloc(MAXM, sizeof *SH); glab = malloc(sizeof *glab * MAXM); gtabu = malloc(sizeof *gtabu * MAXM);
  Fcnt = malloc(sizeof(int) * NT * FLMAX); Ful = malloc(sizeof(int) * NT * FLMAX); Fup = malloc(sizeof(int) * NT * FLMAX);
  FW = malloc(sizeof(int) * NT * FLMAX);
  FADD = malloc(sizeof(unsigned short) * FLMAX * FLMAX);
  Fdl = calloc(MAXNT, sizeof *Fdl); Flab = malloc(sizeof *Flab * MAXM);
  Ffix = malloc((size_t)NT * FLMAX * 4);
  Fbest = 1L << 40; Fbd = 2; Fbm = 0;
  static int G[MAXD][MAXK];
  int q4 = q * q * q * q, seeded = 0;

  /* exact seeds */
  GQ = q; GD = 3; GE = 1; GLS = q;
  int m0 = structured(q, k, 3, out, t0 + seconds * 0.2, 1L << 40);
  if (m0 > 0) {
    Fm = m0;
    for (int i = 0; i < 3; i++) for (int c = 0; c < k; c++) G[i][c] = out->G[i][c];
    f_record(q, k, 3, G); seeded = 1;
  }
  if (Fbest > q4) {
    int fm = fusion_seed(q, k);
    if (fm > 0 && (long)fm * q < Fbest) {
      memset(G, 0, sizeof G);
      for (int c = 0; c < k; c++) G[0][c] = 1;
      if (f_setup(q, 1, G)) {
        Fm = fm; f_build();
        fprintf(stderr, "  fusion m=%d holes=%d\n", fm, Fnu);
        if (Fnu == 0) { f_record(q, k, 1, G); seeded = 1; }
      }
    }
  }
  /* from-scratch search on a fine grid */
  if (Fbest > q4) {
    double frac = !seeded ? 0.6 : q >= 7 ? 0.0 : Fbest > table_N ? 0.5 : 0.2;
    if (frac > 0) f_scratch(q, k, q == 3 ? 1 : 2, table_N, out, now() + (tend - now()) * frac);
  }
  /* drop-and-repair on the best array, splitting each coset into q finer cosets between levels */
  if (Fbest > q4 && Fbest < 1L << 40) {
    Ld = Fbd; Lm = Fbm; memcpy(LS, BS, sizeof(BS[0]) * Fbm);
    for (int i = 0; i < Ld; i++) for (int c = 0; c < k; c++) LG[i][c] = BG[i][c];
    while (now() < tend) {
      if (!f_setup(q, Ld, LG)) break;
      Fm = Lm; memcpy(SH, LS, sizeof(SH[0]) * Lm); f_build();
      if (Fnu != 0) break;
      int more = can_split(q);
      double frac = !more ? 1.0 : Ld == 3 ? 0.15 : Ld == 2 ? 0.4 : 0.5;
      f_descend(q, k, Ld, LG, now() + (tend - now()) * frac);
      if (!more) break;
      for (int i = 0; i < Lm; i++) for (int a = 0; a < q; a++) for (int c = 0; c < k; c++)
        EXS[i * q + a][c] = (unsigned char)((LS[i][c] + a * LG[Ld - 1][c]) % q);
      Lm *= q; Ld--;
      memcpy(LS, EXS, sizeof(EXS[0]) * Lm);
    }
  }
  out->d = Fbd; out->m = Fbm;
  for (int i = 0; i < Fbd; i++) for (int c = 0; c < k; c++) out->G[i][c] = BG[i][c];
  for (int i = 0; i < Fbm; i++) for (int c = 0; c < k; c++) out->S[i][c] = BS[i][c];
  fprintf(stderr, "  chose d=%d m=%d est=%ld at %.1f s\n", Fbd, Fbm, Fbest, now() - t0);
}

__attribute__((unused)) static void old_construct(int q, int k, int table_N, double seconds, Construction *out) {
  double t0 = now(), tend = t0 + seconds - 0.3;
  int GLSmax = q * q;
  mulc = calloc(MAXNT, sizeof *mulc);
  gcnt = malloc(sizeof(int) * NT * GLSmax); gul = malloc(sizeof(int) * NT * GLSmax); gup = malloc(sizeof(int) * NT * GLSmax);
  gW = malloc(sizeof(unsigned) * NT * GLSmax);
  SH = calloc(MAXM, sizeof *SH); glab = malloc(sizeof *glab * MAXM); gtabu = malloc(sizeof *gtabu * MAXM);
  long best_est = 1L << 40; int bestd = 2, bestm = 0;
  for (int pi = 0; pi < 2; pi++) {
    int d = pi == 0 ? 3 : 2;
    double pend = pi == 0 ? t0 + seconds * 0.35 : tend;
    if (now() >= pend) continue;
    GQ = q; GD = d; GE = 4 - d; GLS = ipow(q, GE);
    memset(mulc, 0, MAXNT * sizeof *mulc);
    int P = ipow(q, d), mmax = 40000 / P; if (mmax > MAXM) mmax = MAXM;
    int m, last = -1, step = 1;
    int m0 = structured(q, k, d, out, pend, best_est);
    if (m0 > 0) {
      if ((long)m0 * P < best_est) {
        best_est = (long)m0 * P; bestd = d; bestm = m0; memcpy(BS, SH, sizeof(SH[0]) * m0);
        for (int i = 0; i < d; i++) for (int c = 0; c < k; c++) BG[i][c] = out->G[i][c];
      }
      fprintf(stderr, "  d=%d structured m=%d (N=%d) at %.1f s\n", d, m0, m0 * P, now() - t0);
      if (m0 - 1 < GLS) continue;
      last = m0; step = 1; memcpy(LS, SH, sizeof(SH[0]) * m0);
      drop_rows(m0 - 1); m = m0 - 1;
    } else {
      if (!setup_code(q, k, d, out, now() + (pend - now()) * 0.3)) { fprintf(stderr, "  no code d=%d\n", d); continue; }
      m = (int)(1.3 * table_N / P) + 1; if (m < GLS) m = GLS; if (m > mmax) m = mmax;
      step = m / 20 > 1 ? m / 20 : 1;
      GM_ = m;
      for (int i = 0; i < m; i++) for (int c = 0; c < k; c++) SH[i][c] = (unsigned char)(rnd() % q);
    }
    while (now() < pend) {
      GM_ = m; gbuild();
      double te = now() + (pend - now()) * 0.4;
      if (step == 1 && last >= 0) te = pend;
      if (te - now() < 0.5) te = pend;
      gsearch(te);
      if (gnu == 0) {
        last = m; memcpy(LS, SH, sizeof(SH[0]) * m);
        if ((long)m * P < best_est) {
          best_est = (long)m * P; bestd = d; bestm = m; memcpy(BS, SH, sizeof(SH[0]) * m);
          for (int i = 0; i < d; i++) for (int c = 0; c < k; c++) BG[i][c] = out->G[i][c];
        }
        fprintf(stderr, "  d=%d solved m=%d (N=%d) at %.1f s\n", d, m, m * P, now() - t0);
        if (m - step < GLS) { if (step > 1) { step = 1; } if (m - 1 < GLS) break; }
        drop_rows(m - step); m -= step;
      } else {
        long est = (long)m * P + (long)gnu * P;
        if (est < best_est) {
          best_est = est; bestd = d; bestm = m; memcpy(BS, SH, sizeof(SH[0]) * m);
          for (int i = 0; i < d; i++) for (int c = 0; c < k; c++) BG[i][c] = out->G[i][c];
        }
        if (last < 0) {
          if (m >= mmax) break;
          m = m * 3 / 2 + 1; if (m > mmax) m = mmax;
          for (int i = 0; i < m; i++) for (int c = 0; c < k; c++) SH[i][c] = (unsigned char)(rnd() % q);
          step = m / 20 > 1 ? m / 20 : 1;
        } else if (step > 1) {
          step /= 2; m = last; memcpy(SH, LS, sizeof(SH[0]) * m);
          GM_ = m; gbuild(); drop_rows(m - step); m -= step;
        } else break;
      }
    }
  }
  out->d = bestd; out->m = bestm;
  for (int i = 0; i < bestd; i++) for (int c = 0; c < k; c++) out->G[i][c] = BG[i][c];
  for (int i = 0; i < bestm; i++) for (int c = 0; c < k; c++) out->S[i][c] = BS[i][c];
  fprintf(stderr, "  chose d=%d m=%d est=%ld\n", bestd, bestm, best_est);
}
// EVOLVE-BLOCK-END

/* ---- fixed evaluation: expand, complete, prune, print ---- */
static unsigned char R[MAXROWS][MAXK];
static int nr;

static inline int tidx(const unsigned char *r, int t) {
  return ((r[cols[t][0]] * Q + r[cols[t][1]]) * Q + r[cols[t][2]]) * Q + r[cols[t][3]];
}

static int expand(const Construction *c) {
  if (c->d < 0 || c->d > MAXD || c->m < 0 || c->m > MAXM) return -1;
  int P = ipow(Q, c->d);
  if ((long)c->m * P > MAXROWS) return -1;
  nr = 0;
  for (int i = 0; i < c->m; i++)
    for (int u = 0; u < P; u++) {
      int uu[MAXD], x = u;
      for (int j = 0; j < c->d; j++) { uu[j] = x % Q; x /= Q; }
      for (int col = 0; col < K; col++) {
        int v = c->S[i][col] % Q;
        for (int j = 0; j < c->d; j++) v += uu[j] * (c->G[j][col] % Q);
        R[nr][col] = (unsigned char)(((v % Q) + Q) % Q);
      }
      nr++;
    }
  return nr;
}

static int complete_and_prune(void) {
  int Q4 = Q * Q * Q * Q, TT = NT * Q4;
  unsigned short *cnt = calloc(TT, sizeof *cnt);
  for (int r = 0; r < nr; r++) for (int t = 0; t < NT; t++) { int k = t * Q4 + tidx(R[r], t); if (cnt[k] < 65535) cnt[k]++; }
  int *miss = malloc(sizeof(int) * TT), nm = 0;
  for (int k = 0; k < TT; k++) if (!cnt[k]) miss[nm++] = k;
  int added = 0, start = 0;
  while (start < nm) {
    if (nr >= MAXROWS) { free(cnt); free(miss); return -1; }
    int fixed[MAXK]; for (int c = 0; c < K; c++) fixed[c] = -1;
    for (int i = start; i < nm; i++) {           /* greedy merge of compatible missing tuples */
      int k = miss[i]; if (cnt[k]) continue;
      int t = k / Q4, x = k % Q4, v[4] = { x / (Q * Q * Q), (x / (Q * Q)) % Q, (x / Q) % Q, x % Q }, ok = 1;
      for (int j = 0; j < 4 && ok; j++) if (fixed[cols[t][j]] >= 0 && fixed[cols[t][j]] != v[j]) ok = 0;
      if (!ok) continue;
      for (int j = 0; j < 4; j++) fixed[cols[t][j]] = v[j];
    }
    for (int c = 0; c < K; c++) R[nr][c] = (unsigned char)(fixed[c] >= 0 ? fixed[c] : (int)(rnd() % Q));
    for (int t = 0; t < NT; t++) { int k = t * Q4 + tidx(R[nr], t); if (cnt[k] < 65535) cnt[k]++; }
    nr++; added++;
    while (start < nm && cnt[miss[start]]) start++;
  }
  int pruned = 0;                                 /* drop rows whose every tuple is covered twice */
  for (int r = nr - 1; r >= 0; r--) {
    int red = 1;
    for (int t = 0; t < NT && red; t++) if (cnt[t * Q4 + tidx(R[r], t)] < 2) red = 0;
    if (!red) continue;
    for (int t = 0; t < NT; t++) cnt[t * Q4 + tidx(R[r], t)]--;
    nr--; if (r != nr) memcpy(R[r], R[nr], K); pruned++;
  }
  fprintf(stderr, "  completion added %d rows, pruning removed %d, N=%d\n", added, pruned, nr);
  free(cnt); free(miss);
  return nr;
}

static void setup(int q, int k) {
  Q = q; K = k; NT = 0;
  for (int a = 0; a < k; a++) for (int b = a + 1; b < k; b++) for (int c = b + 1; c < k; c++) for (int d = c + 1; d < k; d++) {
    cols[NT][0] = a; cols[NT][1] = b; cols[NT][2] = c; cols[NT][3] = d;
    smask[NT] = (1 << a) | (1 << b) | (1 << c) | (1 << d); NT++;
  }
}

static Construction CON;

int main(void) {
  const char *grid = getenv("GRID"), *ss = getenv("SIZE_SECONDS");
  double secs = ss ? atof(ss) : 30;
  if (!grid) { fprintf(stderr, "GRID not set\n"); return 1; }
  char *g = strdup(grid);
  for (char *tok = strtok(g, ";"); tok; tok = strtok(NULL, ";")) {
    int q, k, tn;
    if (sscanf(tok, "%d,%d,%d", &q, &k, &tn) != 3 || k > MAXK || k < 4) continue;
    printf("# size %d %d\n", q, k); fflush(stdout);
    fprintf(stderr, "size q=%d k=%d table=%d\n", q, k, tn);
    pid_t pid = fork();
    if (pid == 0) {
      alarm((unsigned)(secs * 1.5 + 30));
      rs ^= (unsigned long long)time(0) * 2654435761ULL ^ ((unsigned long long)q << 32) ^ (unsigned long long)k;
      for (int i = 0; i < 10; i++) rnd();
      setup(q, k);
      memset(&CON, 0, sizeof CON);
      double t0 = now();
      construct(q, k, tn, secs, &CON);
      fprintf(stderr, "  construct: d=%d m=%d (%d rows before fix-up) in %.1f s\n", CON.d, CON.m,
              CON.m * ipow(q, CON.d < 0 ? 0 : CON.d), now() - t0);
      if (expand(&CON) < 0 || complete_and_prune() < 0) { fprintf(stderr, "  invalid construction\n"); _exit(2); }
      size_t len = (size_t)nr * (2 * K + 1);
      char *buf = malloc(len + 1), *p = buf;
      for (int r = 0; r < nr; r++) { for (int c = 0; c < K; c++) { *p++ = (char)('0' + R[r][c]); *p++ = c + 1 < K ? ' ' : '\n'; } }
      fwrite(buf, 1, (size_t)(p - buf), stdout); fflush(stdout);
      _exit(0);
    }
    int st = 0; waitpid(pid, &st, 0);
    if (!(WIFEXITED(st) && WEXITSTATUS(st) == 0)) fprintf(stderr, "  size q=%d k=%d failed (status %d)\n", q, k, st);
    printf("# end %d %d\n", q, k); fflush(stdout);
  }
  return 0;
}
