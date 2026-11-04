"""Constructions of the thesis of K. Shokri (2026, doi:10.20381/ruor-31663) for strength 3, and an exact triple-coverage check.
Symbols are field elements coded 0..q-1 (prime q: residues; q = 2^n: bit vectors of polynomials in w)."""
import hashlib, itertools, json, os, time
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "out")
MODPOLY = {4: 0b111, 8: 0b1011, 16: 0b10011}        # w^2+w+1, w^3+w+1, w^4+w+1 (report, section 0)
PRIM = {2: 1, 3: 2, 5: 2, 7: 3}                       # primitive e for prime q; w (= 2) for q = 2^n


class GF:
    def __init__(self, q):
        self.q = q
        r = np.arange(q)
        if q in PRIM:
            self.add = (r[:, None] + r[None, :]) % q
            self.mul = (r[:, None] * r[None, :]) % q
            self.e = PRIM[q]
        else:
            m, n = MODPOLY[q], q.bit_length() - 1
            def pm(a, b):
                p = 0
                while b:
                    if b & 1:
                        p ^= a
                    b >>= 1
                    a <<= 1
                    if a >> n:
                        a ^= m
                return p
            self.add = r[:, None] ^ r[None, :]
            self.mul = np.array([[pm(a, b) for b in range(q)] for a in range(q)])
            self.e = 2
        self.neg = np.array([int(np.where(self.add[a] == 0)[0][0]) for a in range(q)])
        self.inv = np.array([0] + [int(np.where(self.mul[a] == 1)[0][0]) for a in range(1, q)])
        pw, x = [], 1
        for _ in range(q - 1):
            pw.append(x)
            x = int(self.mul[x, self.e])
        if sorted(pw) != list(range(1, q)):
            raise ValueError(f"e = {self.e} is not primitive in GF({q})")
        self.epow = pw                                  # e^i, i = 0..q-2
        self.elt = [0] + pw                             # elt(0) = 0, elt(h) = e^(h-1)

    def sub(self, a, b):
        return int(self.add[a, self.neg[b]])

    def dot(self, U, G):
        """U: n x 3, G: 3 x m  ->  n x m matrix of u . G[:, col]."""
        U, G = np.asarray(U), np.asarray(G)
        acc = self.mul[U[:, 0:1], G[0][None, :]]
        for i in range(1, U.shape[1]):
            acc = self.add[acc, self.mul[U[:, i:i + 1], G[i][None, :]]]
        return acc


def allvec(q, n=3):
    return np.array(list(itertools.product(range(q), repeat=n)), dtype=np.int64)


def cubic_powers(F, f):
    """f = (a0, a1, a2) for x^3 + a2 x^2 + a1 x + a0. Returns list of L(x^j), j = 0..q^3-2, or None if x is not primitive."""
    q, (a0, a1, a2) = F.q, f
    cur, out, seen = (1, 0, 0), [], set()
    for _ in range(q ** 3 - 1):
        if cur in seen:
            return None
        seen.add(cur)
        out.append(cur)
        c0, c1, c2 = cur
        cur = (F.sub(0, int(F.mul[c2, a0])), F.sub(c0, int(F.mul[c2, a1])), F.sub(c1, int(F.mul[c2, a2])))
    return out if cur == (1, 0, 0) else None


def primitive_cubic(F):
    if F.q == 3:
        f = (1, 2, 1)                                   # x^3 + x^2 + 2x + 1, thesis Fig. 5.13
        if cubic_powers(F, f) is None:
            raise ValueError("Fig. 5.13 cubic x^3+x^2+2x+1 is not primitive over GF(3)")
        return f
    for a2, a1, a0 in itertools.product(range(F.q), range(F.q), range(1, F.q)):
        if cubic_powers(F, (a0, a1, a2)) is not None:
            return (a0, a1, a2)
    raise ValueError("no primitive cubic")


def normalize(F, v):
    i = next(i for i in range(3) if v[i])
    s = int(F.inv[v[i]])
    return tuple(int(F.mul[s, x]) for x in v)


class Ingredients:
    def __init__(self, q, cubic=None):
        F = self.F = GF(q)
        self.q, self.k = q, q * q + q + 1
        self.f = cubic or primitive_cubic(F)
        pw = self.pow = cubic_powers(F, self.f)
        M, k = q ** 3 - 1, self.k
        self.G1 = np.array([pw[j] for j in range(k)]).T                 # column j = L(alpha^j)
        self.Gm1 = np.array([pw[(-j) % M] for j in range(k)]).T          # column j = L(alpha^-j)
        self.U = allvec(q)
        self.A1, self.Am1 = F.dot(self.U, self.G1), F.dot(self.U, self.Gm1)
        jof = {normalize(F, pw[j]): j for j in range(k)}
        if len(jof) != k:
            raise ValueError("alpha^0..alpha^(k-1) are not k distinct projective points")
        self.jof = jof
        self.P = self._P()

    def rms(self):
        """Thm 3.2.2: A(G^1) over A(G^-1) minus one zero row = CA(2q^3-1; 3, q^2+q+1, q)."""
        B = np.vstack([self.A1, self.Am1])
        z = [i for i in range(len(B)) if not B[i].any()]
        return np.delete(B, z[1], axis=0)

    def _P(self):
        """Cor 5.1.9: CA(2q^2-q; 2, q^2+2q, q); first q^2+q+1 columns: all (w,z), all II_w, III_0."""
        F, q = self.F, self.q
        rho = [(a1, a0) for a1 in range(1, q) for a0 in range(q)]
        Fm = [[int(F.add[F.mul[a1, F.elt[z]], a0]) for z in range(q)] for a1, a0 in rho]
        cv = [a1 for a1, _ in rho]
        cols = []
        for w in range(q):
            for z in range(q):
                cols.append([r[z] for r in Fm] + [r[w] for r in Fm] + list(range(q)))
        for w in range(q):
            cols.append(cv + [r[w] for r in Fm] + [0] * q)
        z = 0
        cols.append([r[z] for r in Fm] + cv + [0] * q)
        return np.array(cols[:self.k]).T


A2_Q3 = [[1, 1, 1, 1, 1, 1, 1, 1, 1], [2, 1, 1, 1, 2, 0, 2, 0, 0], [1, 0, 1, 2, 0, 0, 2, 1, 2]]   # Fig. 3.2 A2, report


def orthogoval_bc(F):
    """b, c with x^3 + b x + c irreducible (no root), report 3.3.7 step 1."""
    q = F.q
    for b, c in itertools.product(range(q), range(1, q)):
        if all(F.add[F.add[F.mul[F.mul[x, x], x], F.mul[b, x]], c] for x in range(q)):
            return b, c
    raise ValueError("no irreducible x^3+bx+c")


def g2star(F, x, y, b, c):
    """Cor 5.2.12, q = 2^n: (1, x^2+y, y^2+by+cx)."""
    m = F.mul
    return (1, int(F.add[m[x, x], y]), int(F.add[F.add[m[y, y], m[b, y]], m[c, x]]))


def z_cubic(q):
    """Thm 5.2.13 needs |Z| = q - 1. Over GF(q^3) the zero set of the first entry of L(alpha^(-j)) on the affine points
    has q - 1 or q + 1 elements depending on the primitive cubic; the first primitive cubic gives q + 1 at q = 2, 4."""
    F = GF(q)
    for a2, a1, a0 in itertools.product(range(q), range(q), range(1, q)):
        f = (a0, a1, a2)
        if q == 3 and f != (1, 2, 1):
            continue
        if cubic_powers(F, f) is None:
            continue
        try:
            I = Ingredients(q, f)
        except ValueError:
            continue
        if sum(1 for x in range(q) for y in range(q) if I.pow[(-I.jof[(1, F.elt[y], F.elt[x])]) % (q ** 3 - 1)][0] == 0) == q - 1:
            return I
    raise ValueError("no primitive cubic with |Z| = q - 1")


def build5(thm, q, swap=False, noscale=False, ing=None):
    """Thm 5.2.11 / Cor 5.2.12 / Thm 5.2.13. Returns (array, info). Columns are (c, j) at index (block position)*k + j."""
    I = ing or (z_cubic(q) if thm == "5.2.13" else Ingredients(q))
    F, k, U = I.F, I.k, I.U
    elt, M = F.elt, q ** 3 - 1
    gh = [(c // q, c % q) for c in range(q * q)]
    idx = [(h * q + g) if swap else c for c, (g, h) in enumerate(gh)]       # swap: last block uses point (elt(g), elt(h))
    lp = [(elt[gh[i][1]], elt[gh[i][0]]) for i in idx]                     # last-block point (x, y) of block c
    pc = [(1, x, y) for x, y in lp]
    if thm in ("5.2.11", "5.2.13"):
        g2 = [I.pow[(-I.jof[p]) % M] for p in pc]                          # L(alpha^(-j_c))
    elif thm == "5.2.12":
        if q == 3:
            g2 = [tuple(A2_Q3[r][i] for r in range(3)) for i in idx]
        elif q & (q - 1) == 0:
            b, cc = orthogoval_bc(F)
            g2 = [g2star(F, x, y, b, cc) for x, y in lp]
        else:
            raise ValueError("Cor 5.2.12 needs q = 2^n or 3")
    info = dict(theorem=thm, q=q, swap=swap, noscale=noscale, cubic=I.f)
    blocks = list(range(q * q))
    scaled = thm in ("5.2.12",) or (thm == "5.2.13" and not noscale)
    if thm == "5.2.13":
        Z = [c for c in blocks if g2[c][0] == 0]
        info["Z"] = Z
        if len(Z) != q - 1:
            info["Z_size_wrong"] = True
        blocks = [c for c in blocks if c not in Z]
        if not noscale:
            g2 = [tuple(int(F.mul[F.inv[v[0]], t]) for t in v) if v[0] else v for v in g2]
    nb = len(blocks)
    xcol = np.repeat([elt[gh[c][1]] for c in blocks], k)
    ycol = np.repeat([elt[gh[c][0]] for c in blocks], k)
    tile = lambda A: np.tile(A, (1, nb))
    parts, names = [tile(I.A1), F.add[tile(I.Am1), xcol[None, :]]], ["R0", "R1"]
    Pt = tile(I.P)
    for t in range(2, q):
        parts.append(F.add[Pt, F.mul[F.epow[t - 1], xcol][None, :]]); names.append(f"R{t}")
    for s in range(1, q):
        parts.append(F.add[Pt, F.mul[F.epow[s - 1], ycol][None, :]]); names.append(f"R'{s}")
    U1 = U[(U[:, 1] != 0) & (U[:, 2] != 0)]
    parts.append(np.repeat(F.dot(U1, np.array([pc[c] for c in blocks]).T), k, axis=1)); names.append("last-G1")
    U2 = U[(U[:, 1] != 0) | (U[:, 2] != 0)] if scaled else U[U.any(axis=1)]
    parts.append(np.repeat(F.dot(U2, np.array([g2[c] for c in blocks]).T), k, axis=1)); names.append("last-G2")
    info["row_blocks"] = [[n, len(p)] for n, p in zip(names, parts)]
    info["col_blocks"] = blocks
    return np.vstack(parts).astype(np.int64), info


def build337(q):
    """Thm 3.3.7 via [17] Cor 3.9 / Thm 4.5, q = 2^n: CA(2q^3-q; 3, q^2+2, q)."""
    F = GF(q)
    b, c = orthogoval_bc(F)
    m, a = F.mul, F.add
    h1 = [(x, y, 1) for x in range(q) for y in range(q)] + [(1, 0, 0), (0, 1, 0)]
    h2 = [(int(a[m[x, x], y]), int(a[a[m[y, y], m[b, y]], m[c, x]]), 1) for x in range(q) for y in range(q)] + [(1, 0, 0), (0, 1, 0)]
    U = allvec(q)
    U2 = U[(U[:, 0] != 0) | (U[:, 1] != 0)]                                  # delete s = 2, u = (0, 0, c)
    A = np.vstack([F.dot(U, np.array(h1).T), F.dot(U2, np.array(h2).T)]).astype(np.int64)
    return A, dict(theorem="3.3.7", q=q, b=b, c=c)


def check2(A, q):
    N, k = A.shape
    for a in range(k):
        for b in range(a + 1, k):
            cnt = np.bincount(A[:, a] * q + A[:, b], minlength=q * q)
            if (cnt == 0).any():
                v = int(np.argmin(cnt))
                return dict(cols=[a, b], missing=[v // q, v % q])
    return None


def sha(A):
    return hashlib.sha1(np.ascontiguousarray(A, dtype=np.int64).tobytes()).hexdigest()


def check3(A, q, ckpt=None, save_every=60.0):
    """Exact strength-3 check. Returns None if every triple is covered, else dict with the first failing triple.
    ckpt: path of a resume file {"sha", "a"}: columns a' < a have been fully checked as first column."""
    N, k = A.shape
    dt = np.uint8 if q <= 8 else np.uint16
    S = np.left_shift(np.ones(A.shape, dtype=dt), A.astype(dt))
    full = (1 << q) - 1
    h, a0 = sha(A), 0
    if ckpt and os.path.exists(ckpt):
        d = json.load(open(ckpt))
        if d.get("sha") == h:
            a0 = d["a"]
    last = time.time()
    for a in range(a0, k - 2):
        for b in range(a + 1, k - 1):
            v = A[:, a] * q + A[:, b]
            cnt = np.bincount(v, minlength=q * q)
            if (cnt == 0).any():
                i = int(np.argmin(cnt))
                return dict(cols=[a, b, b + 1], missing=[[i // q, i % q, "any"]])
            order = np.argsort(v, kind="stable")
            starts = np.concatenate([[0], np.cumsum(cnt)[:-1]])
            R = np.bitwise_or.reduceat(S[order, b + 1:], starts, axis=0)
            bad = R != full
            if bad.any():
                i, jj = (int(x) for x in np.argwhere(bad)[0])
                miss = [[i // q, i % q, s] for s in range(q) if not (int(R[i, jj]) >> s) & 1]
                return dict(cols=[a, b, b + 1 + jj], missing=miss)
        if ckpt and time.time() - last > save_every:
            json.dump(dict(sha=h, a=a + 1, k=k), open(ckpt + ".tmp", "w"))
            os.replace(ckpt + ".tmp", ckpt)
            last = time.time()
    if ckpt:
        json.dump(dict(sha=h, a=k, k=k, done=True), open(ckpt + ".tmp", "w"))
        os.replace(ckpt + ".tmp", ckpt)
    return None


def colname(info, col):
    """(c, j) of a column index for build5 arrays."""
    if "col_blocks" not in info:
        return col
    k = info["q"] ** 2 + info["q"] + 1
    return dict(col=col, c=info["col_blocks"][col // k], j=col % k)


def write_array(A, q):
    os.makedirs(OUT, exist_ok=True)
    N, k = A.shape
    path = os.path.join(OUT, f"CA_{N}_3_{k}_{q}.txt")
    with open(path + ".tmp", "w") as f:
        for r in A:
            f.write(" ".join(map(str, r.tolist())) + "\n")
    os.replace(path + ".tmp", path)
    return path


def formula(thm, q):
    """(N, k) claimed by the report; for 5.2.13 noscale the N of 5.2.11."""
    K = q * q + q + 1
    return {"5.2.11": (8 * q ** 3 - 10 * q ** 2 + 4 * q - 1, q * q * K),
            "5.2.12": (8 * q ** 3 - 10 * q ** 2 + 3 * q, q * q * K),
            "5.2.13": (8 * q ** 3 - 10 * q ** 2 + 3 * q, (q * q - q + 1) * K),
            "5.2.13-noscale": (8 * q ** 3 - 10 * q ** 2 + 4 * q - 1, (q * q - q + 1) * K),
            "3.3.7": (2 * q ** 3 - q, q * q + 2)}[thm]
