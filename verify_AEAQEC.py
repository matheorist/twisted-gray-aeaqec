# -*- coding: utf-8 -*-
"""
Machine verification of all matrix identities, NSC properties and code
parameters appearing in AEAQEC.tex.  Pure Python, no dependencies.
Run:  python verify_AEAQEC.py
"""
import itertools, sys

# ---------------------------------------------------------------- finite fields
def find_irreducible(p, e):
    """Brute-force a monic irreducible polynomial of degree e over F_p.
    Coefficient lists are low-to-high, length e+1, leading coeff 1."""
    def polymulmod_none(a, b):  # plain poly multiply over F_p
        r = [0]*(len(a)+len(b)-1)
        for i, x in enumerate(a):
            for j, y in enumerate(b):
                r[i+j] = (r[i+j] + x*y) % p
        return r
    def divides(d, f):  # does monic d divide f over F_p?
        f = f[:]
        while len(f) >= len(d) and any(f):
            if f[-1] == 0:
                f.pop(); continue
            c = f[-1]
            shift = len(f)-len(d)
            for i, x in enumerate(d):
                f[shift+i] = (f[shift+i] - c*x) % p
            while f and f[-1] == 0:
                f.pop()
        return not any(f)
    for tail in itertools.product(range(p), repeat=e):
        f = list(tail)+[1]
        if f[0] == 0:
            continue
        ok = True
        for de in range(1, e//2+1):
            for dt in itertools.product(range(p), repeat=de):
                d = list(dt)+[1]
                if divides(d, f):
                    ok = False; break
            if not ok:
                break
        if ok:
            return f
    raise RuntimeError

class GF:
    """GF(p^e); elements are integers 0..p^e-1 encoding base-p coeff vectors."""
    def __init__(self, p, e=1):
        self.p, self.e, self.q = p, e, p**e
        if e == 1:
            self.add = lambda a, b: (a+b) % p
            self.mul = lambda a, b: (a*b) % p
            self.neg = lambda a: (-a) % p
        else:
            mod = find_irreducible(p, e)
            def enc(v): return sum(c*(p**i) for i, c in enumerate(v))
            def dec(a):
                v = []
                for _ in range(e):
                    v.append(a % p); a //= p
                return v
            def add(a, b):
                va, vb = dec(a), dec(b)
                return enc([(x+y) % p for x, y in zip(va, vb)])
            def neg(a):
                return enc([(-x) % p for x in dec(a)])
            def mul(a, b):
                va, vb = dec(a), dec(b)
                r = [0]*(2*e-1)
                for i, x in enumerate(va):
                    for j, y in enumerate(vb):
                        r[i+j] = (r[i+j] + x*y) % p
                for k in range(2*e-2, e-1, -1):  # reduce by mod (monic)
                    c = r[k]
                    if c:
                        for i in range(e):
                            r[k-e+i] = (r[k-e+i] - c*mod[i]) % p
                        r[k] = 0
                return enc(r[:e])
            self.add, self.mul, self.neg = add, mul, neg
        self.zero, self.one = 0, 1 if e == 1 else 1
        self.elements = list(range(self.q))
    def sub(self, a, b): return self.add(a, self.neg(b))
    def inv(self, a):
        for b in range(1, self.q):
            if self.mul(a, b) == self.one:
                return b
        raise ZeroDivisionError
    def pow(self, a, n):
        r = self.one
        for _ in range(n):
            r = self.mul(r, a)
        return r
    def order(self, a):
        r, k = a, 1
        while r != self.one:
            r = self.mul(r, a); k += 1
            if k > self.q: raise RuntimeError
        return k
    def element_of_order(self, s):
        for a in range(2, self.q):
            try:
                if self.order(a) == s:
                    return a
            except RuntimeError:
                pass
        raise RuntimeError(f"no element of order {s}")
    def sqrt(self, a):
        for b in range(self.q):
            if self.mul(b, b) == a:
                return b
        return None

# ---------------------------------------------------------------- matrices
def mat_mul(F, A, B):
    n, m, r = len(A), len(B), len(B[0])
    return [[
        __import__('functools').reduce(F.add, (F.mul(A[i][k], B[k][j]) for k in range(m)))
        for j in range(r)] for i in range(n)]

def transpose(A): return [list(r) for r in zip(*A)]

def det(F, M):
    M = [row[:] for row in M]; n = len(M); d = F.one
    for c in range(n):
        piv = next((r for r in range(c, n) if M[r][c] != 0), None)
        if piv is None: return 0
        if piv != c:
            M[c], M[piv] = M[piv], M[c]; d = F.neg(d)
        d = F.mul(d, M[c][c]); inv = F.inv(M[c][c])
        for r in range(c+1, n):
            f = F.mul(M[r][c], inv)
            for j in range(c, n):
                M[r][j] = F.sub(M[r][j], F.mul(f, M[c][j]))
    return d

def is_NSC(F, A):
    s = len(A)
    for k in range(1, s+1):
        for cols in itertools.combinations(range(s), k):
            sub = [[A[i][j] for j in cols] for i in range(k)]
            if det(F, sub) == 0:
                return False
    return True

def gram(F, A): return mat_mul(F, A, transpose(A))

def check_monomial(F, G, rho, Dexp=None):
    """G == diag(a)*I_rho with rho a 1-based involution list; entry (rho(i),i)."""
    s = len(G)
    for i in range(s):
        for j in range(s):
            expect_nonzero = (rho[j] - 1 == i)
            if expect_nonzero and G[i][j] == 0: return False
            if not expect_nonzero and G[i][j] != 0: return False
    if Dexp is not None:
        for j in range(s):
            if G[rho[j]-1][j] != Dexp[j]: return False
    return True

def fourier(F, s, w):
    return [[F.pow(w, i*j) for j in range(s)] for i in range(s)]

def rho_fourier(s):  # rho_s: 1->1, i->s+2-i
    return [1] + [s+2-i for i in range(2, s+1)]

OK = True
def report(name, cond):
    global OK
    print(("PASS  " if cond else "FAIL  ") + name)
    if not cond: OK = False

# ================================================================ matrix checks
# --- Fourier matrices over F_5, F_7, F_9, F_49, F_81
for (p, e, s) in [(5,1,4), (7,1,3), (3,2,4), (7,2,4), (3,4,4)]:
    F = GF(p, e); q = F.q
    assert (q-1) % s == 0
    w = F.element_of_order(s)
    A = fourier(F, s, w)
    G = gram(F, A)
    sF = F.one
    for _ in range(s-1): sF = F.add(sF, F.one)   # s in the field
    ok = check_monomial(F, G, rho_fourier(s), [sF]*s) and is_NSC(F, A)
    report(f"Fourier s={s} over F_{q}: FF^T = {s}*I_rho_s, NSC", ok)

# --- explicit F_5 Fourier matrix printed in the paper
F5 = GF(5)
Apaper = [[1,1,1,1],[1,2,4,3],[1,4,1,4],[1,3,4,2]]
ok = (fourier(F5, 4, 2) == Apaper and
      check_monomial(F5, gram(F5, Apaper), [1,4,3,2], [4,4,4,4]) and is_NSC(F5, Apaper))
report("paper's F_5 matrix equals F_4(omega=2), Gram=4*I_(2,4), NSC", ok)

# --- A(gamma) over F_4 and F_8, rho=(1,3)
for (p, e) in [(2,2),(2,3)]:
    F = GF(p, e)
    g = next(a for a in F.elements if a not in (0, F.one))
    one = F.one; gp1 = F.add(one, g)
    A = [[one, g, gp1],[one, one, one],[one, gp1, g]]
    ok = check_monomial(F, gram(F, A), [3,2,1], [one]*3) and is_NSC(F, A)
    report(f"A(gamma) over F_{F.q}: AA^T = I_(1,3), NSC", ok)

# --- explicit id-OD 4x4 over F_7 from the elimination (Example)
F7 = GF(7)
A7 = [[1,1,1,1],[2,3,4,5],[1,6,6,1],[6,3,4,1]]
ok = (check_monomial(F7, gram(F7, A7), [1,2,3,4], [4,5,4,6]) and is_NSC(F7, A7))
report("F_7 example: AA^T = diag(4,5,4,6), NSC", ok)

# --- elimination over F_7 and F_11 from Vandermonde nodes 1,2,3,4 (auto-run)
def sym_eliminate(F, N):
    """Run the Case-1-only elimination; return (A=LN, D) or None if a pivot dies."""
    s = len(N); A = [row[:] for row in N]
    G = gram(F, A); D = []
    for c in range(s):
        if G[c][c] == 0: return None
        D.append(G[c][c]); inv = F.inv(G[c][c])
        for r in range(c+1, s):
            f = F.neg(F.mul(G[r][c], inv))
            for j in range(s):
                A[r][j] = F.add(A[r][j], F.mul(f, A[c][j]))
        G = gram(F, A)
    return A, D

for (q, Dexp) in [(7, [4,5,4,6]), (11, [4,5,4,4])]:
    F = GF(q)
    N = [[F.pow(x, i) for x in (1,2,3,4)] for i in range(4)]
    res = sym_eliminate(F, N)
    ok = res is not None and res[1] == Dexp and is_NSC(F, res[0]) \
         and check_monomial(F, gram(F, res[0]), [1,2,3,4], Dexp)
    report(f"elimination over F_{q} (nodes 1..4): rho=id, D=diag{tuple(Dexp)}, NSC", ok)
    if q == 7 and res is not None:
        report("F_7 elimination reproduces the printed matrix", res[0] == A7)

# --- dihedral twist over F_17: A = F_4*diag(1,8,13,2), Gram = 4*I_(2,1)(4,3)
F17 = GF(17)
w = 4; assert F17.order(w) == 4
sigma = 8; assert F17.mul(sigma, sigma) == F17.inv(w)   # sigma^2 = w^{-1}
Lam = [1, sigma, F17.mul(sigma, sigma), F17.mul(F17.mul(sigma, sigma), sigma)]
report("F_17 dihedral: Lambda = diag(1,8,13,2)", Lam == [1,8,13,2])
A = [[F17.mul(fourier(F17,4,w)[i][j], Lam[j]) for j in range(4)] for i in range(4)]
ok = check_monomial(F17, gram(F17, A), [2,1,4,3], [4,4,4,4]) and is_NSC(F17, A)
report("F_17 dihedral c=1: AA^T = 4*I_(1,2)(3,4), NSC (fixed-point-free)", ok)

# --- s=2 classification witnesses
for q in (3,5,7,11):
    F = GF(q)
    b = 1
    A = [[1,b],[b,F.neg(1)]]
    ok = check_monomial(F, gram(F, A), [1,2]) and is_NSC(F, A)
    report(f"s=2 id-OD over F_{q}: [[1,1],[1,-1]]", ok)
for q in (5,13,17):
    F = GF(q)
    a = F.sqrt(F.neg(1)); assert a is not None
    A = [[1,a],[a,1]]
    ok = check_monomial(F, gram(F, A), [2,1]) and is_NSC(F, A)
    report(f"s=2 (1,2)-OD over F_{q} (q=1 mod 4): [[1,a],[a,1]]", ok)
F4 = GF(2,2); g = 2
A = [[1,g],[g,1]]
report("s=2 id-OD over F_4: [[1,gamma],[gamma,1]]",
       check_monomial(F4, gram(F4, A), [1,2]) and is_NSC(F4, A))
# q=2 exhaustive non-existence
F2 = GF(2)
found = False
for A in itertools.product([0,1], repeat=4):
    M = [[A[0],A[1]],[A[2],A[3]]]
    if det(F2, M) == 0 or not is_NSC(F2, M):
        continue
    G = gram(F2, M)
    if check_monomial(F2, G, [1,2]) or check_monomial(F2, G, [2,1]):
        found = True
report("q=2: exhaustive search confirms NO 2x2 Euclidean rho-OD matrix", not found)

# --- full reversal at s = q  (new theorem): A = (I - (1/2)E_{q,1}) * V(all of F_q)
for q in (3,5,7,11):
    F = GF(q)
    nodes = list(range(q))
    V = [[F.pow(x, i) for x in nodes] for i in range(q)]
    c = F.neg(F.inv(2 % q))          # c = -1/2
    A = [row[:] for row in V]
    A[q-1] = [F.add(A[q-1][j], F.mul(c, A[0][j])) for j in range(q)]
    rho = [q+1-i for i in range(1, q+1)]          # full reversal
    mJ = F.neg(F.one)
    ok = check_monomial(F, gram(F, A), rho, [mJ]*q) and is_NSC(F, A)
    report(f"full reversal at s=q over F_{q}: AA^T = -J, NSC", ok)

# ================================================================ table checks
def dzdx(s, rho, k, l):
    dz = min((s-i+1)*(k[rho[i-1]-1]+1) for i in range(1, s+1))
    dx = min((s-i+1)*(l[rho[s-i]-1]+1) for i in range(1, s+1))
    return dz, dx

def check_row(name, q, s, n, k, l, t, rho, code, net):
    N, K, dZ, dX, c = code
    ok = True
    for i in range(s):
        lo, hi = max(0, l[i]-k[i]), min(l[i], n-k[i])
        ok &= lo <= t[i] <= hi
        if n == q+1:
            ok &= (n-k[i], l[i], t[i]) not in [(2,1,1),(1,2,1)]
    cc = sum(l) - sum(t)                    # c = sum(l_i - t_i)
    KK = s*n - (sum(k)+sum(t))              # K = sn - sum(k_i + t_i)
    dz, dx = dzdx(s, rho, k, l)
    ok &= (N == s*n and c == cc and K == KK and dz >= dZ and dx >= dX)
    ok &= abs((K-c)/N - net) < 5e-4
    report(f"table row {name}: [[{N},{K},>={dZ}/>={dX};{c}]]_{q}, net {net}", ok)

rid4 = [1,2,3,4]; rf4 = rho_fourier(4); rf3 = rho_fourier(3); rg3 = [3,2,1]
check_row("q5 s4 F",  5,4,4,  [1,2,1,1],[2,1,1,1],[1,0,0,0], rf4, (16,10,3,3,4), 0.375)
check_row("q7 s3 F",  7,3,8,  [1,2,1],[2,0,1],[1,0,0],       rf3, (24,19,3,3,2), 0.708)
check_row("q7 s3 F'", 7,3,8,  [1,3,1],[2,0,1],[1,0,0],       rf3, (24,18,4,3,2), 0.667)
check_row("q7 s4 id", 7,4,8,  [1,1,1,2],[2,1,1,1],[1,0,0,0], rid4,(32,26,3,3,4), 0.688)
check_row("q9 s4 F",  9,4,8,  [0,3,1,1],[3,0,1,1],[3,0,0,0], rf4, (32,24,4,4,2), 0.688)
check_row("q11 s4 id",11,4,12,[1,1,2,3],[3,2,1,1],[2,1,0,0], rid4,(48,38,4,4,4), 0.708)
check_row("q49 s4 F", 49,4,12,[0,3,1,1],[3,0,1,1],[3,0,0,0], rf4, (48,40,4,4,2), 0.792)
check_row("q81 s4 F", 81,4,20,[1,4,2,1],[4,1,1,2],[3,0,0,1], rf4, (80,68,5,5,4), 0.800)
check_row("q4 s3 Ag", 4,3,5,  [2,1,0],[0,1,2],[0,0,2],       rg3, (15,10,3,3,1), 0.600)
check_row("q8 s3 Ag", 8,3,9,  [3,1,1],[1,1,3],[0,0,2],       rg3, (27,20,4,4,3), 0.630)

# ================================================================ exhaustive
# classification of realizable involutions for s = 3 over small fields
def classify_s3(F):
    """Return the set of involutions rho (as tuples) realizable by 3x3
    Euclidean rho-OD matrices over F, by exhaustive search."""
    q = F.q
    rows_all = list(itertools.product(F.elements, repeat=3))
    rows_nz  = [r for r in rows_all if all(x != 0 for x in r)]
    def dot(u, v):
        return F.add(F.add(F.mul(u[0],v[0]), F.mul(u[1],v[1])), F.mul(u[2],v[2]))
    found = set()
    targets = {(1,2,3),(2,1,3),(3,2,1),(1,3,2)}
    for r1 in rows_nz:
        for r2 in rows_all:
            # NSC 2x2 minors from rows 1,2 (all three column pairs)
            m01 = F.sub(F.mul(r1[0],r2[1]), F.mul(r1[1],r2[0]))
            m02 = F.sub(F.mul(r1[0],r2[2]), F.mul(r1[2],r2[0]))
            m12 = F.sub(F.mul(r1[1],r2[2]), F.mul(r1[2],r2[1]))
            if m01 == 0 or m02 == 0 or m12 == 0:
                continue
            d11, d12 = dot(r1,r1), dot(r1,r2)
            for r3 in rows_all:
                A = [list(r1), list(r2), list(r3)]
                if det(F, A) == 0:
                    continue
                G = [[d11,d12,dot(r1,r3)],[d12,dot(r2,r2),dot(r2,r3)],
                     [dot(r1,r3),dot(r2,r3),dot(r3,r3)]]
                # monomial test -> involution
                rho = []
                good = True
                for i in range(3):
                    nz = [j for j in range(3) if G[i][j] != 0]
                    if len(nz) != 1:
                        good = False; break
                    rho.append(nz[0]+1)
                if not good:
                    continue
                t = tuple(rho)
                if t in targets and t not in found:
                    found.add(t)
                    if found == targets:
                        return found
    return found

names = {(1,2,3):"id",(2,1,3):"(1,2)",(3,2,1):"(1,3)",(1,3,2):"(2,3)"}
# F_3: rho(1)=1 excluded (all-nonzero rows are isotropic since 3|s); (1,2) fails deeper.
# F_5: rho(1)!=1 excluded (no all-nonzero isotropic vector in F_5^3).
# F_4: everything with f>=1; F_7: expected unobstructed.
expected = {3: {"(1,3)"},
            4: {"id","(1,2)","(1,3)","(2,3)"},
            5: {"id"},
            7: {"id","(1,2)","(1,3)","(2,3)"}}
for q, e in [(3,1),(4,2),(5,1),(7,1)]:
    F = GF(q if e == 1 else 2, e)
    got = {names[t] for t in classify_s3(F)}
    print(f"      s=3 over F_{F.q}: realizable involutions = {sorted(got)}")
    report(f"s=3 over F_{F.q}: matches expectation {sorted(expected[F.q])}",
           got == expected[F.q])

print()
print("ALL CHECKS PASSED" if OK else "*** SOME CHECKS FAILED ***")
sys.exit(0 if OK else 1)
