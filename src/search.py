"""Tabu local search for mixed ternary/binary covering codes K(t, b, R).

Words are integers in mixed radix (t ternary digits first, then b binary digits).
The objective is the number of uncovered words; a move relocates one codeword.

Move selection (per iteration):
  * pick a random uncovered word u;
  * for every non-codeword x in the ball B(u) (so the move always covers u), and every
    codeword c, compute the EXACT change in the number of uncovered words for c -> x:
        delta = loss[c] - corr(c, x) - gain(x)
    where gain(x) = #{w in B(x): cnt[w] == 0}, loss[c] = #{w in B(c): cnt[w] == 1}, and
    corr(c, x) = #{w in B(x) ∩ B(c): cnt[w] == 1} (found via sumcov[w], which equals the
    unique covering codeword when cnt[w] == 1);
  * take the best non-tabu move (ties broken at random); aspiration if it beats the best.
Tabu: a removed word may not be re-added, and an added word may not be removed, for
`tenure` iterations.

Descending mode: whenever coverage is complete at size M, the code is saved, the
codeword whose removal uncovers the fewest words is deleted, and the search continues
at M - 1, until the target size is reached or time runs out.

Results only count after src/checker.py (independent) validates the saved file.
"""
import argparse
import csv
import datetime as dt
import os
import sys
import time

import numpy as np
from numba import njit

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)


# ---------------------------------------------------------------- space
def radices(t, b):
    return np.array([3] * t + [2] * b, dtype=np.int64)


def all_digits(t, b):
    r = radices(t, b)
    n = len(r)
    N = int(np.prod(r)) if n else 1
    idx = np.arange(N, dtype=np.int64)
    digs = np.zeros((N, n), dtype=np.int64)
    stride = 1
    for i in range(n - 1, -1, -1):  # last coordinate least significant
        digs[:, i] = (idx // stride) % r[i]
        stride *= r[i]
    return digs


def to_index(digs, t, b):
    r = radices(t, b)
    out = np.zeros(digs.shape[0], dtype=np.int64)
    for i in range(len(r)):
        out = out * r[i] + digs[:, i]
    return out


def build_balls(t, b, R):
    """balls[w] = indices of all words within distance R of w (group translate of D)."""
    digs = all_digits(t, b)
    r = radices(t, b)
    wt = (digs != 0).sum(axis=1)
    D = digs[wt <= R]  # difference patterns of weight <= R
    N = digs.shape[0]
    balls = np.empty((N, D.shape[0]), dtype=np.int64)
    for j in range(D.shape[0]):
        balls[:, j] = to_index((digs + D[j]) % r, t, b)
    return balls


def word_str(i, t, b):
    r = radices(t, b)
    d = []
    for k in range(len(r) - 1, -1, -1):
        d.append(int(i % r[k]))
        i //= r[k]
    return "".join(str(x) for x in reversed(d))


def parse_word(s, t, b):
    s = s.strip().replace(" ", "").replace("|", "").replace(",", "")
    r = radices(t, b)
    assert len(s) == len(r), s
    v = 0
    for k, ch in enumerate(s):
        dgt = int(ch)
        assert 0 <= dgt < r[k]
        v = v * int(r[k]) + dgt
    return v


# ---------------------------------------------------------------- core (numba)
@njit(cache=True)
def _add(x, s, balls, code, slot, cnt, sumcov, loss, unc_list, unc_pos, nunc):
    code[s] = x
    slot[x] = s
    loss[s] = 0
    for j in range(balls.shape[1]):
        w = balls[x, j]
        if cnt[w] == 0:
            loss[s] += 1
            # remove w from uncovered list
            p = unc_pos[w]
            last = unc_list[nunc - 1]
            unc_list[p] = last
            unc_pos[last] = p
            unc_pos[w] = -1
            nunc -= 1
        elif cnt[w] == 1:
            loss[slot[sumcov[w]]] -= 1
        cnt[w] += 1
        sumcov[w] += x
    return nunc


@njit(cache=True)
def _remove(s, balls, code, slot, cnt, sumcov, loss, unc_list, unc_pos, nunc):
    c = code[s]
    slot[c] = -1
    for j in range(balls.shape[1]):
        w = balls[c, j]
        cnt[w] -= 1
        sumcov[w] -= c
        if cnt[w] == 0:
            unc_list[nunc] = w
            unc_pos[w] = nunc
            nunc += 1
        elif cnt[w] == 1:
            loss[slot[sumcov[w]]] += 1
    loss[s] = 0
    code[s] = -1
    return nunc


@njit(cache=True)
def _seed(sd):
    np.random.seed(sd)


@njit(cache=True)
def run_chunk(balls, code, slot, cnt, sumcov, loss, unc_list, unc_pos, nunc_arr,
              M, tabu_add, tabu_rem, it0, iters, tenure, best_arr, corr):
    """Run up to `iters` tabu iterations. Returns the iteration counter.
    Stops early when nunc == 0. best_arr[0] tracks the best uncovered count."""
    nunc = nunc_arr[0]
    B = balls.shape[1]
    it = it0
    while it < it0 + iters and nunc > 0:
        it += 1
        u = unc_list[np.random.randint(nunc)]
        bestd = 1 << 30
        bx = -1
        bs = -1
        nties = 0
        for jj in range(B):
            x = balls[u, jj]
            if slot[x] >= 0:
                continue
            gain = 0
            for j in range(B):
                w = balls[x, j]
                if cnt[w] == 0:
                    gain += 1
                elif cnt[w] == 1:
                    corr[slot[sumcov[w]]] += 1
            xtabu = tabu_add[x] > it
            for s in range(M):
                d = loss[s] - corr[s] - gain
                if xtabu or tabu_rem[code[s]] > it:
                    if nunc + d >= best_arr[0]:  # aspiration only if new best
                        continue
                if d < bestd:
                    bestd = d
                    bx = x
                    bs = s
                    nties = 1
                elif d == bestd:
                    nties += 1
                    if np.random.randint(nties) == 0:
                        bx = x
                        bs = s
            for j in range(B):
                w = balls[x, j]
                if cnt[w] == 1:
                    corr[slot[sumcov[w]]] = 0
        if bx < 0:
            continue
        c = code[bs]
        nunc = _remove(bs, balls, code, slot, cnt, sumcov, loss, unc_list, unc_pos, nunc)
        nunc = _add(bx, bs, balls, code, slot, cnt, sumcov, loss, unc_list, unc_pos, nunc)
        tabu_add[c] = it + tenure
        tabu_rem[bx] = it + tenure
        if nunc < best_arr[0]:
            best_arr[0] = nunc
    nunc_arr[0] = nunc
    return it


# ---------------------------------------------------------------- driver
class State:
    def __init__(self, balls, words):
        N = balls.shape[0]
        self.balls = balls
        self.N = N
        self.M = len(words)
        self.code = np.full(self.M, -1, dtype=np.int64)
        self.slot = np.full(N, -1, dtype=np.int64)
        self.cnt = np.zeros(N, dtype=np.int64)
        self.sumcov = np.zeros(N, dtype=np.int64)
        self.loss = np.zeros(self.M, dtype=np.int64)
        self.unc_list = np.arange(N, dtype=np.int64)
        self.unc_pos = np.arange(N, dtype=np.int64)
        self.nunc = np.array([N], dtype=np.int64)
        self.corr = np.zeros(self.M, dtype=np.int64)
        for s, x in enumerate(words):
            self.nunc[0] = _add(x, s, balls, self.code, self.slot, self.cnt, self.sumcov,
                                self.loss, self.unc_list, self.unc_pos, self.nunc[0])

    def words(self):
        return [int(x) for x in self.code if x >= 0]

    def drop_cheapest(self, rng):
        # delete the codeword that uncovers the fewest words (ties random)
        m = self.loss.min()
        cands = np.flatnonzero(self.loss == m)
        s = int(rng.choice(cands))
        self.nunc[0] = _remove(s, self.balls, self.code, self.slot, self.cnt, self.sumcov,
                               self.loss, self.unc_list, self.unc_pos, self.nunc[0])
        # compact: move last slot into s
        last = self.M - 1
        if s != last:
            x = self.code[last]
            self.nunc[0] = _remove(last, self.balls, self.code, self.slot, self.cnt, self.sumcov,
                                   self.loss, self.unc_list, self.unc_pos, self.nunc[0])
            self.nunc[0] = _add(x, s, self.balls, self.code, self.slot, self.cnt, self.sumcov,
                                self.loss, self.unc_list, self.unc_pos, self.nunc[0])
        self.M -= 1
        self.code = self.code[:self.M].copy()
        self.loss = self.loss[:self.M].copy()
        self.corr = np.zeros(self.M, dtype=np.int64)


def read_code(path, t, b):
    out = []
    with open(path) as f:
        for line in f:
            s = line.strip()
            if s and not s.startswith("#"):
                out.append(parse_word(s, t, b))
    return out


def write_code(path, words, t, b, R, note=""):
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, "w") as f:
        f.write(f"# K(t={t}, b={b}, R={R}) covering, M={len(words)} {note}\n")
        for x in sorted(words):
            f.write(word_str(x, t, b) + "\n")


def log_run(row):
    path = os.path.join(ROOT, "results", "runs.csv")
    new = not os.path.exists(path)
    with open(path, "a", newline="") as f:
        w = csv.writer(f)
        if new:
            w.writerow(["timestamp", "target", "start", "M_start", "M_target", "seed", "tenure",
                        "iterations", "wall_s", "smallest_M_covered", "best_uncovered_at_final_M",
                        "final_M", "saved_file"])
        w.writerow(row)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("t", type=int, help="number of ternary coordinates")
    ap.add_argument("b", type=int, help="number of binary coordinates")
    ap.add_argument("R", type=int)
    ap.add_argument("--start", type=int, help="random start size (default = target)")
    ap.add_argument("--init", help="file with a starting code (e.g. a record code)")
    ap.add_argument("--delete", type=int, default=0, help="codewords to delete from --init before searching")
    ap.add_argument("--target", type=int, required=True, help="stop when a covering of this size is found")
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--tenure", type=int, default=0, help="tabu tenure (0 = auto)")
    ap.add_argument("--time", type=float, default=60.0, help="wall-time limit in seconds (hard cap 600)")
    ap.add_argument("--chunk", type=int, default=20000)
    ap.add_argument("--label", default="")
    a = ap.parse_args()
    tlim = min(a.time, 600.0)

    t, b, R = a.t, a.b, a.R
    target_name = f"K({t},{b},{R})"
    balls = build_balls(t, b, R)
    N, B = balls.shape
    rng = np.random.default_rng(a.seed)
    _seed(a.seed)

    if a.init:
        words = read_code(a.init, t, b)
        words = list(dict.fromkeys(words))
        rng.shuffle(words)
        words = words[:len(words) - a.delete] if a.delete else words
        start_desc = f"file:{os.path.basename(a.init)}-{a.delete}"
    else:
        M0 = a.start or a.target
        words = [int(x) for x in rng.choice(N, size=M0, replace=False)]
        start_desc = "random"
    st = State(balls, words)
    M_start = st.M
    tenure = a.tenure or max(2, min(15, st.M // 10 + 1))
    print(f"{target_name}: N={N} |B|={B} start={start_desc} M={st.M} target={a.target} "
          f"seed={a.seed} tenure={tenure} uncovered={st.nunc[0]}", flush=True)

    t0 = time.time()
    it = 0
    total_it = 0
    smallest_ok = None
    saved = ""
    best = np.array([st.nunc[0]], dtype=np.int64)
    tabu_add = np.zeros(N, dtype=np.int64)
    tabu_rem = np.zeros(N, dtype=np.int64)
    last_print = t0
    while True:
        if st.nunc[0] == 0:
            smallest_ok = st.M
            fn = os.path.join(ROOT, "results", "codes", f"K_{t}_{b}_{R}_M{st.M}_seed{a.seed}.txt")
            write_code(fn, st.words(), t, b, R, note=f"seed={a.seed} found by src/search.py")
            if saved and os.path.join(ROOT, saved) != fn:
                os.remove(os.path.join(ROOT, saved))  # keep only the smallest covering of this run
            saved = os.path.relpath(fn, ROOT)
            print(f"  covered with M={st.M} after {total_it} it, {time.time() - t0:.1f}s -> {saved}", flush=True)
            if st.M <= a.target:
                break
            st.drop_cheapest(rng)
            best[0] = st.nunc[0]
            tabu_add[:] = 0
            tabu_rem[:] = 0
            it = 0
        if time.time() - t0 > tlim:
            break
        new_it = run_chunk(balls, st.code, st.slot, st.cnt, st.sumcov, st.loss, st.unc_list,
                           st.unc_pos, st.nunc, st.M, tabu_add, tabu_rem, it, a.chunk, tenure,
                           best, st.corr)
        total_it += new_it - it
        it = new_it
        if time.time() - last_print > 15:
            print(f"  M={st.M} it={total_it} uncovered={st.nunc[0]} best={best[0]} "
                  f"t={time.time() - t0:.0f}s", flush=True)
            last_print = time.time()
    wall = time.time() - t0
    print(f"done: smallest covered M={smallest_ok}, final M={st.M}, best uncovered at final M={best[0]}, "
          f"{total_it} it, {wall:.1f}s", flush=True)
    log_run([dt.datetime.now().isoformat(timespec="seconds"), target_name + (f" {a.label}" if a.label else ""),
             start_desc, M_start, a.target, a.seed, tenure, total_it, f"{wall:.1f}", smallest_ok,
             int(best[0]), st.M, saved])
    return 0


if __name__ == "__main__":
    sys.exit(main())
