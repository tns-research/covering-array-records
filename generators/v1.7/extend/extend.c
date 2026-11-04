/* Add one column to a strength-4 covering array: extend in.txt (N x k, symbols 0..v-1) to k+1 columns.
   1. local search on the new column's values: every 3-set S of old columns and every pattern p on S must see all v values;
   2. greedy extra rows for what is still missing (old 4-sets stay covered: old rows are unchanged).
   usage: extend in.txt k v seconds seed > out.txt */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
int N, K, V, P3; unsigned char *A; int *S3; int NS; unsigned short *cnt; long missing;
static inline int pat(int r, int s) { int *c = S3 + 3 * s; unsigned char *x = A + (long)r * (K + 1); return (x[c[0]] * V + x[c[1]]) * V + x[c[2]]; }
static inline unsigned short *C(int s, int p) { return cnt + ((long)s * P3 + p) * V; }
int main(int ac, char **av) {
  FILE *f = fopen(av[1], "r"); K = atoi(av[2]); V = atoi(av[3]); double secs = atof(av[4]); srand(atoi(av[5]));
  int cap = 1 << 16; A = malloc((long)cap * (K + 1)); N = 0; int x;
  while (fscanf(f, "%d", &x) == 1) { A[(long)N * (K + 1)] = x; for (int j = 1; j < K; j++) { fscanf(f, "%d", &x); A[(long)N * (K + 1) + j] = x; } N++; if (N + 8 >= cap) { cap *= 2; A = realloc(A, (long)cap * (K + 1)); } }
  for (int j = K - 1; j >= 0; j--) ; 
  P3 = V * V * V; NS = K * (K - 1) * (K - 2) / 6; S3 = malloc(sizeof(int) * 3 * NS); int s = 0;
  for (int a = 0; a < K; a++) for (int b = a + 1; b < K; b++) for (int c = b + 1; c < K; c++) { S3[3*s] = a; S3[3*s+1] = b; S3[3*s+2] = c; s++; }
  cnt = calloc((long)NS * P3 * V, sizeof(unsigned short));
  for (int r = 0; r < N; r++) A[(long)r * (K + 1) + K] = rand() % V;
  for (int r = 0; r < N; r++) for (s = 0; s < NS; s++) C(s, pat(r, s))[A[(long)r * (K + 1) + K]]++;
  missing = 0; for (long i = 0; i < (long)NS * P3 * V; i++) missing += cnt[i] == 0;
  fprintf(stderr, "start missing %ld\n", missing);
  time_t t0 = time(0); long it = 0;
  while (missing && difftime(time(0), t0) < secs) {
    int r = rand() % N, old = A[(long)r * (K + 1) + K], best = -1; long bd = 1L << 60;
    for (int nv = 0; nv < V; nv++) if (nv != old) {
      long d = 0; for (s = 0; s < NS; s++) { unsigned short *c = C(s, pat(r, s)); d += (c[old] == 1) - (c[nv] == 0); }
      if (d < bd || (d == bd && rand() & 1)) { bd = d; best = nv; }
    }
    if (bd <= 0 || rand() % 50 == 0) {
      for (s = 0; s < NS; s++) { unsigned short *c = C(s, pat(r, s)); if (--c[old] == 0) missing++; if (c[best]++ == 0) missing--; }
      A[(long)r * (K + 1) + K] = best;
    }
    if (++it % 20000 == 0) fprintf(stderr, "it %ld missing %ld\n", it, missing);
  }
  fprintf(stderr, "after search missing %ld\n", missing);
  int K1 = K + 1, added = 0; unsigned char *row = malloc(K1); signed char *fix = malloc(K1);
  while (missing) {
    memset(fix, -1, K1);
    for (s = 0; s < NS; s++) for (int p = 0; p < P3; p++) for (int v = 0; v < V; v++) if (!C(s, p)[v]) {
      int *c = S3 + 3 * s, pv[3] = {p / (V * V), (p / V) % V, p % V}, ok = fix[K] < 0 || fix[K] == v;
      for (int i = 0; i < 3 && ok; i++) ok = fix[c[i]] < 0 || fix[c[i]] == pv[i];
      if (ok) { for (int i = 0; i < 3; i++) fix[c[i]] = pv[i]; fix[K] = v; }
    }
    for (int j = 0; j < K1; j++) row[j] = fix[j] < 0 ? rand() % V : fix[j];
    if (N + 1 >= cap) { cap *= 2; A = realloc(A, (long)cap * K1); }
    memcpy(A + (long)N * K1, row, K1);
    for (s = 0; s < NS; s++) if (C(s, pat(N, s))[row[K]]++ == 0) missing--;
    N++; added++;
  }
  fprintf(stderr, "added %d rows, N = %d\n", added, N);
  for (int r = 0; r < N; r++) { for (int j = 0; j < K1; j++) printf(j ? " %d" : "%d", A[(long)r * K1 + j]); printf("\n"); }
  return 0;
}
