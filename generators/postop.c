/* Classic post-optimisation of a covering array CA(N; t, k, v), any t <= 6, no LLM.
 *
 *   postop FILE t v SECONDS [OUT]
 *
 * Prints one JSON line: the array's weakness profile (how many t-tuples each row is the only one
 * to cover) and the result of a greedy delete-and-patch loop run for SECONDS:
 * pick a row with few uniquely covered tuples, delete it, and re-home each of its holes by
 * re-assigning entries of another row, accepted only if that uncovers nothing. Each success
 * removes one row; OUT receives the smallest verified array found. SECONDS = 0 profiles only.
 * This is the baseline the LLM-written searches are measured against (same core, same clock).
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

static int N, K, V, T;
static unsigned char **A;
static unsigned short *cnt;
static long long B[2048][8];
static long long VT, NSETS;

static long long rank_set(const int *c) { long long r = 0; for (int i = 0; i < T; i++) r += B[c[i]][i + 1]; return r; }
static long long idx(const unsigned char *row, const int *c) {
  long long x = 0; for (int i = 0; i < T; i++) x = x * V + row[c[i]];
  return rank_set(c) * VT + x;
}

/* enumerate every t-subset of columns, calling f */
static void each_set(void (*f)(const int *, void *), void *ctx) {
  int c[8]; for (int i = 0; i < T; i++) c[i] = i;
  while (1) {
    f(c, ctx);
    int i = T - 1; while (i >= 0 && c[i] == K - T + i) i--;
    if (i < 0) break;
    c[i]++; for (int j = i + 1; j < T; j++) c[j] = c[j - 1] + 1;
  }
}

typedef struct { unsigned char *row; int s; long u; } RowCtx;
static void f_add(const int *c, void *p) { RowCtx *r = p; cnt[idx(r->row, c)] += r->s; }
static void f_uniq(const int *c, void *p) { RowCtx *r = p; if (cnt[idx(r->row, c)] == 1) r->u++; }
static void addrow(unsigned char *row, int s) { RowCtx r = {row, s, 0}; each_set(f_add, &r); }
static long uniq(unsigned char *row) { RowCtx r = {row, 0, 0}; each_set(f_uniq, &r); return r.u; }

/* holes of a deleted row = its tuples whose count dropped to 0 */
typedef struct { int c[8]; unsigned char v[8]; } Hole;
static Hole *holes; static int nh, hcap;
typedef struct { unsigned char *row; } HCtx;
static void f_hole(const int *c, void *p) {
  HCtx *h = p;
  if (cnt[idx(h->row, c)] == 0) {
    if (nh == hcap) { hcap = hcap ? 2 * hcap : 256; holes = realloc(holes, hcap * sizeof(Hole)); }
    for (int i = 0; i < T; i++) { holes[nh].c[i] = c[i]; holes[nh].v[i] = h->row[c[i]]; }
    nh++;
  }
}

/* would changing columns ch[0..m) of row x to y uncover anything? enumerate the t-sets that
   contain at least one changed column, each once (its smallest changed column leads). */
static unsigned char *Y; static int ischg[2048];
static unsigned char **logrow; static int *logidx, nlog; static long long lcap;
static int safe_change(unsigned char *x, int *ch, int m) {
  int rest[8];
  for (int q = 0; q < m; q++) {
    int c0 = ch[q];
    /* choose T-1 other columns from K\{c0}, skipping sets whose smaller changed column exists */
    int o[8]; int r = T - 1;
    for (int i = 0; i < r; i++) o[i] = i;
    int others = K - 1;
    while (1) {
      int s[8], n = 0, ok = 1;
      for (int i = 0; i < r; i++) { int col = o[i] >= c0 ? o[i] + 1 : o[i]; if (ischg[col] && col < c0) { ok = 0; break; } rest[i] = col; }
      if (ok) {
        /* sorted set */
        int ins = 0; for (int i = 0; i < r; i++) { if (!ins && c0 < rest[i]) { s[n++] = c0; ins = 1; } s[n++] = rest[i]; }
        if (!ins) s[n++] = c0;
        long long io = idx(x, s), in = idx(Y, s);
        if (io != in && cnt[io] <= 1) return 0;
      }
      if (r == 0) break;
      int i = r - 1; while (i >= 0 && o[i] == others - r + i) i--;
      if (i < 0) break;
      o[i]++; for (int j = i + 1; j < r; j++) o[j] = o[j - 1] + 1;
    }
  }
  return 1;
}

static int patch(int ri, Hole *h) {
  unsigned char *x = A[ri]; int ch[8], m = 0;
  for (int i = 0; i < T; i++) if (x[h->c[i]] != h->v[i]) ch[m++] = h->c[i];
  if (m == 0) return 1;
  memcpy(Y, x, K); for (int i = 0; i < T; i++) Y[h->c[i]] = h->v[i];
  for (int i = 0; i < m; i++) ischg[ch[i]] = 1;
  int ok = safe_change(x, ch, m);
  for (int i = 0; i < m; i++) ischg[ch[i]] = 0;
  if (ok) { addrow(x, -1); memcpy(x, Y, K); addrow(x, 1); }
  return ok;
}

static void write_array(const char *path) {
  FILE *f = fopen(path, "w");
  for (int i = 0; i < N; i++) { for (int j = 0; j < K; j++) fprintf(f, j ? " %d" : "%d", A[i][j]); fputc('\n', f); }
  fclose(f);
}

int main(int argc, char **argv) {
  if (argc < 5) { fprintf(stderr, "postop FILE t v SECONDS [OUT]\n"); return 2; }
  T = atoi(argv[2]); V = atoi(argv[3]); double lim = atof(argv[4]); const char *out = argc > 5 ? argv[5] : NULL;
  FILE *f = fopen(argv[1], "r"); if (!f) { perror(argv[1]); return 2; }
  /* rows: whitespace/comma separated ints, one row per line; symbols may start at 1 */
  int cap = 1024; A = malloc(cap * sizeof(*A)); char *line = NULL; size_t ln = 0; int minsym = 1 << 30;
  int *tmp = malloc(4096 * sizeof(int));
  while (getline(&line, &ln, f) > 0) {
    int k = 0; char *p = line, *e;
    while (1) { while (*p == ' ' || *p == ',' || *p == '\t' || *p == '[' || *p == ']') p++; long val = strtol(p, &e, 10); if (e == p) break; tmp[k++] = val; p = e; }
    if (!k) continue;
    if (N == 0) K = k; else if (k != K) { printf("{\"error\":\"ragged rows\"}\n"); return 1; }
    if (N == cap) { cap *= 2; A = realloc(A, cap * sizeof(*A)); }
    A[N] = malloc(K); for (int i = 0; i < K; i++) { if (tmp[i] < minsym) minsym = tmp[i]; A[N][i] = tmp[i]; }
    N++;
  }
  if (N == 0) { printf("{\"error\":\"no rows parsed\"}\n"); return 1; }
  if (minsym > 0) for (int i = 0; i < N; i++) for (int j = 0; j < K; j++) A[i][j] -= minsym;
  if (K >= 2048 || T > 7) { printf("{\"error\":\"too large\"}\n"); return 1; }
  for (int n = 0; n < 2048; n++) { B[n][0] = 1; for (int j = 1; j < 8; j++) B[n][j] = n == 0 ? 0 : B[n - 1][j - 1] + B[n - 1][j]; }
  VT = 1; for (int i = 0; i < T; i++) VT *= V; NSETS = B[K][T];
  double cells = (double)NSETS * VT;
  if (cells > 4e9) { printf("{\"N\":%d,\"k\":%d,\"error\":\"too large (%.2g counters)\"}\n", N, K, cells); return 1; }
  cnt = calloc((size_t)NSETS * VT, sizeof(unsigned short)); Y = malloc(K);
  if (!cnt) { printf("{\"error\":\"out of memory\"}\n"); return 1; }
  for (int i = 0; i < N; i++) addrow(A[i], 1);
  long long missing = 0; for (long long i = 0; i < NSETS * VT; i++) if (!cnt[i]) missing++;
  long *u = malloc(N * sizeof(long)); long umin = 1L << 60; int nmin = 0, n_le4 = 0;
  for (int i = 0; i < N; i++) { u[i] = uniq(A[i]); if (u[i] < umin) { umin = u[i]; nmin = 0; } if (u[i] == umin) nmin++; if (u[i] <= 4) n_le4++; }
  long long hmax = 1; for (int i = 0; i < N; i++) if (u[i] > hmax) hmax = u[i];
  lcap = hmax;
  logrow = malloc(hmax * sizeof(*logrow)); logidx = malloc(hmax * sizeof(int));
  for (long long i = 0; i < hmax; i++) logrow[i] = malloc(K);
  int N0 = N; clock_t t0 = clock(); long tries = 0;
  srand(12345);
  while (!missing && lim > 0 && (double)(clock() - t0) / CLOCKS_PER_SEC < lim) {
    /* candidate rows by increasing unique count, random tie-break; retry the list until time runs out */
    int *ord = malloc(N * sizeof(int)); for (int i = 0; i < N; i++) { ord[i] = i; u[i] = uniq(A[i]); }
    for (int i = N - 1; i > 0; i--) { int j = rand() % (i + 1); int t = ord[i]; ord[i] = ord[j]; ord[j] = t; }
    for (int i = 1; i < N; i++) { int x = ord[i], j = i; while (j > 0 && u[ord[j - 1]] > u[x]) { ord[j] = ord[j - 1]; j--; } ord[j] = x; }
    int removed = 0;
    for (int oi = 0; oi < N && (double)(clock() - t0) / CLOCKS_PER_SEC < lim; oi++) {
      int r = ord[oi]; tries++;
      addrow(A[r], -1); nh = 0; HCtx hc = {A[r]}; each_set(f_hole, &hc);
      int ok = 1; nlog = 0;
      for (int h = 0; h < nh && ok; h++) {
        int done = 0, start = rand() % N;
        for (int q = 0; q < N && !done; q++) {
          int y = (start + q) % N; if (y == r) continue;
          if (nlog == lcap) { lcap *= 2; logrow = realloc(logrow, lcap * sizeof(*logrow)); logidx = realloc(logidx, lcap * sizeof(int)); for (long long i = lcap / 2; i < lcap; i++) logrow[i] = malloc(K); }
          memcpy(logrow[nlog], A[y], K); logidx[nlog] = y;
          done = patch(y, &holes[h]);
          if (done && memcmp(logrow[nlog], A[y], K)) nlog++;
        }
        ok = done;
      }
      if (ok) {
        free(A[r]); A[r] = A[N - 1]; N--;
        removed = 1;
        if (out) write_array(out);
        break;
      }
      for (int i = nlog - 1; i >= 0; i--) { int y = logidx[i]; addrow(A[y], -1); memcpy(A[y], logrow[i], K); addrow(A[y], 1); }
      addrow(A[r], 1);
    }
    free(ord);
    if (!removed) {
      /* plateau move (the NCK idea): randomise entries whose change uncovers nothing, then retry */
      for (int it = 0; it < 20 * N; it++) {
        int y = rand() % N, c = rand() % K, val = rand() % V; if (A[y][c] == val) continue;
        memcpy(Y, A[y], K); Y[c] = val; ischg[c] = 1; int ch[1] = {c};
        int ok = safe_change(A[y], ch, 1); ischg[c] = 0;
        if (ok) { addrow(A[y], -1); memcpy(A[y], Y, K); addrow(A[y], 1); }
      }
    }
  }
  printf("{\"N\":%d,\"t\":%d,\"k\":%d,\"v\":%d,\"missing\":%lld,\"min_unique\":%ld,\"rows_at_min\":%d,\"rows_unique_le4\":%d,"
         "\"postop_N\":%d,\"postop_removed\":%d,\"tries\":%ld,\"seconds\":%.0f}\n",
         N0, T, K, V, missing, umin, nmin, n_le4, N, N0 - N, tries, (double)(clock() - t0) / CLOCKS_PER_SEC);
  return 0;
}
