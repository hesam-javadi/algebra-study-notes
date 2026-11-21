(function () {
"use strict";

/* ===================== Exact rational arithmetic ===================== */
const gcd = (a, b) => { a = a < 0n ? -a : a; b = b < 0n ? -b : b; while (b) { [a, b] = [b, a % b]; } return a; };

class Frac {
  constructor(n, d = 1n) {
    n = BigInt(n); d = BigInt(d);
    if (d === 0n) throw new Error("division by zero");
    if (d < 0n) { n = -n; d = -d; }
    const g = gcd(n, d) || 1n;
    this.n = n / g; this.d = d / g;
  }
  add(o) { return new Frac(this.n * o.d + o.n * this.d, this.d * o.d); }
  sub(o) { return new Frac(this.n * o.d - o.n * this.d, this.d * o.d); }
  mul(o) { return new Frac(this.n * o.n, this.d * o.d); }
  div(o) { return new Frac(this.n * o.d, this.d * o.n); }
  neg() { return new Frac(-this.n, this.d); }
  abs() { return this.n < 0n ? this.neg() : this; }
  isZero() { return this.n === 0n; }
  isInt() { return this.d === 1n; }
  cmp(o) { const x = this.n * o.d - o.n * this.d; return x < 0n ? -1 : x > 0n ? 1 : 0; }
  eq(o) { return this.n === o.n && this.d === o.d; }
  toNumber() { return Number(this.n) / Number(this.d); }
  toTeX(mode = "fraction") {
    const neg = this.n < 0n, an = neg ? -this.n : this.n, s = neg ? "-" : "";
    if (this.d === 1n) return s + an;
    if (mode === "decimal") {
      const k = 10000n;
      if (k % this.d === 0n) {
        const sc = an * (k / this.d);
        const fp = (sc % k).toString().padStart(4, "0").replace(/0+$/, "");
        return s + (sc / k) + (fp ? "." + fp : "");
      }
    }
    return s + "\\frac{" + an + "}{" + this.d + "}";
  }
  toText() { return this.d === 1n ? String(this.n) : this.n + "/" + this.d; }
}

const F = (n, d = 1) => new Frac(n, d);
const ZERO = F(0), ONE = F(1);

function normDigits(s) {
  return String(s)
    .replace(/[۰-۹]/g, c => String(c.charCodeAt(0) - 0x06F0))
    .replace(/[٠-٩]/g, c => String(c.charCodeAt(0) - 0x0660))
    .replace(/[٫,]/g, ".")
    .replace(/[−–—]/g, "-")
    .replace(/\s+/g, "");
}
function parseDec(s) {
  const m = /^([+-]?)(\d*)(?:\.(\d*))?$/.exec(s);
  if (!m) return null;
  const ip = m[2] || "", fp = m[3] || "";
  if (!ip && !fp) return null;
  let n = BigInt((ip + fp) || "0");
  if (m[1] === "-") n = -n;
  return new Frac(n, 10n ** BigInt(fp.length));
}
function parseFrac(str) {
  if (str == null) return null;
  const s = normDigits(str);
  if (!s) return null;
  if (s.includes("/")) {
    const parts = s.split("/");
    if (parts.length !== 2) return null;
    const a = parseDec(parts[0]), b = parseDec(parts[1]);
    if (!a || !b || b.isZero()) return null;
    return a.div(b);
  }
  return parseDec(s);
}

/* ===================== Matrix helpers ===================== */
const clone = M => M.map(r => r.slice());
const eye = n => Array.from({ length: n }, (_, i) => Array.from({ length: n }, (_, j) => (i === j ? ONE : ZERO)));
const mmul = (A, B) => A.map(r => B[0].map((_, j) => r.reduce((s, a, k) => s.add(a.mul(B[k][j])), ZERO)));
const meq = (A, B) => A.length === B.length && A.every((r, i) => r.length === B[i].length && r.every((v, j) => v.eq(B[i][j])));
const mat = rows => rows.map(r => r.map(v => (v instanceof Frac ? v : parseFrac(String(v)))));
const augment = (A, B) => A.map((r, i) => r.concat(B[i]));
const augI = A => augment(A, eye(A.length));
const sub = (M, r0, r1, c0, c1) => M.slice(r0, r1).map(r => r.slice(c0, c1));

function lcmBig(a, b) { return (a / gcd(a, b)) * b; }
function texMat(M, mode) {
  return "\\begin{bmatrix}" + M.map(r => r.map(v => v.toTeX(mode)).join("&")).join("\\\\") + "\\end{bmatrix}";
}
function texMatScaled(M, mode) {
  let D = 1n;
  M.forEach(r => r.forEach(v => { D = lcmBig(D, v.d); }));
  if (D === 1n) return texMat(M, mode);
  const Z = M.map(r => r.map(v => v.mul(new Frac(D))));
  return "\\frac{1}{" + D + "}" + texMat(Z, "fraction");
}

function det(A) {
  const n = A.length, M = clone(A);
  let d = ONE;
  for (let c = 0; c < n; c++) {
    let p = -1;
    for (let r = c; r < n; r++) if (!M[r][c].isZero()) { p = r; break; }
    if (p < 0) return ZERO;
    if (p !== c) { [M[p], M[c]] = [M[c], M[p]]; d = d.neg(); }
    d = d.mul(M[c][c]);
    for (let r = c + 1; r < n; r++) {
      const f = M[r][c].div(M[c][c]);
      if (f.isZero()) continue;
      for (let k = c; k < n; k++) M[r][k] = M[r][k].sub(f.mul(M[c][k]));
    }
  }
  return d;
}

function inverse(A) {
  const n = A.length, M = augI(A);
  for (let c = 0; c < n; c++) {
    let p = -1;
    for (let r = c; r < n; r++) if (!M[r][c].isZero()) { p = r; break; }
    if (p < 0) return null;
    [M[p], M[c]] = [M[c], M[p]];
    const pv = M[c][c];
    M[c] = M[c].map(v => v.div(pv));
    for (let r = 0; r < n; r++) {
      if (r === c || M[r][c].isZero()) continue;
      const f = M[r][c];
      M[r] = M[r].map((v, k) => v.sub(f.mul(M[c][k])));
    }
  }
  return sub(M, 0, n, n, 2 * n);
}

function normInf(A) {
  let best = ZERO;
  A.forEach(r => { const s = r.reduce((a, v) => a.add(v.abs()), ZERO); if (s.cmp(best) > 0) best = s; });
  return best;
}
function leadingMinors(A) {
  return A.map((_, k) => det(sub(A, 0, k + 1, 0, k + 1)));
}
function diagDominance(A) {
  const rows = A.map((r, i) => {
    const d = r[i].abs();
    const off = r.reduce((s, v, j) => (j === i ? s : s.add(v.abs())), ZERO);
    return { diag: d, off, cmp: d.cmp(off) };
  });
  const verdict = rows.every(r => r.cmp > 0) ? "strict" : rows.every(r => r.cmp >= 0) ? "weak" : "none";
  return { rows, verdict };
}
function kappaInf(A) {
  const inv = inverse(A);
  if (!inv) return null;
  const na = normInf(A), ni = normInf(inv);
  return { normA: na, normInv: ni, kappa: na.mul(ni), inv };
}

/* ===================== Step builders ===================== */
const ORD = ["اول", "دوم", "سوم", "چهارم", "پنجم", "ششم"];
const Rr = i => "R_{" + (i + 1) + "}";
const Aij = (i, j) => "a_{" + (i + 1) + (j + 1) + "}";

function makeCtx(S, opts) {
  opts = opts || {};
  return {
    M: clone(S), n: S.length, m: S[0].length, steps: [], mode: opts.mode || "fraction",
    names: opts.names || null, stopMsg: opts.stopMsg || "دستگاه جواب یکتا ندارد",
    L: null, perm: null, lu: false, plu: false, lenient: !!opts.lenient, x: null,
  };
}
function T(c, f) { return f.toTeX(c.mode); }
function Pn(c, f) { return f.n < 0n ? "\\left(" + T(c, f) + "\\right)" : T(c, f); }
function push(c, o) {
  c.steps.push(Object.assign({
    M: clone(c.M), L: c.L ? clone(c.L) : null, perm: c.perm ? c.perm.slice() : null,
    x: c.x ? c.x.slice() : null,
  }, o));
}
function xName(c, i) { return c.names ? c.names[i] : "x_{" + (i + 1) + "}"; }

function elim(c, i, j, phase, label, extra) {
  const M = c.M, piv = M[j][j], before = M[i].slice(), l = before[j].div(piv);
  if (!l.isZero()) for (let k = 0; k < c.m; k++) M[i][k] = M[i][k].sub(l.mul(M[j][k]));
  const changed = [], zeros = [];
  for (let k = 0; k < c.m; k++) {
    if (!before[k].eq(M[i][k])) changed.push([i, k]);
    if (k === j && !before[k].isZero() && M[i][k].isZero()) zeros.push([i, k]);
  }
  const op = l.isZero()
    ? "$" + Rr(i) + " \\leftarrow " + Rr(i) + "$ (تغییری نمی‌کند)"
    : "$" + Rr(i) + " \\leftarrow " + Rr(i) + " - " + Pn(c, l) + "\\,{" + Rr(j) + "}$";
  const text = label + "$l_{" + (i + 1) + (j + 1) + "}=\\frac{" + Aij(i, j) + "}{" + Aij(j, j) + "}=\\frac{" +
    T(c, before[j]) + "}{" + T(c, piv) + "}=" + T(c, l) + "$ ⇒ " + op;
  push(c, Object.assign({ kind: "elim", phase, j: j + 1, i: i + 1, l, pivot: [j, j], pivotRow: j, target: i, changed, zeros, text }, extra || {}));
  return l;
}

function forward(c) {
  const n = c.n;
  for (let j = 0; j < n - 1; j++) {
    let r = -1;
    for (let k = j; k < n; k++) if (!c.M[k][j].isZero()) { r = k; break; }
    if (r < 0) {
      if (c.lenient) {
        push(c, { kind: "pivot", phase: "forward", j: j + 1, pivot: null, pivotRow: null,
          text: "تکرار $j=" + (j + 1) + "$ — ستون " + (j + 1) + " از سطر " + (j + 1) + " به پایین صفر است؛ از این ستون می‌گذریم." });
        continue;
      }
      push(c, { kind: "stop", phase: "forward", j: j + 1,
        text: "ستون " + (j + 1) + " از سطر " + (j + 1) + " به پایین تمام صفر است ⇒ " + c.stopMsg });
      return false;
    }
    const head = "تکرار $j=" + (j + 1) + "$ — محورگیری: ";
    if (r === j) {
      push(c, { kind: "pivot", phase: "forward", j: j + 1, pivot: [j, j], pivotRow: j,
        text: head + "$" + Aij(j, j) + " = " + T(c, c.M[j][j]) + " \\neq 0$ پس $r = " + (j + 1) + "$ و جابه‌جایی لازم نیست." });
    } else {
      if (c.lu && !c.plu) {
        push(c, { kind: "stop", phase: "forward", j: j + 1, pivot: [r, j], pivotRow: r,
          text: head + "$" + Aij(j, j) + " = 0$ و باید سطرها جابه‌جا شوند ⇒ تجزیهٔ $LU$ بدون جابه‌جایی وجود ندارد؛ از $PA=LU$ استفاده کنید." });
        return false;
      }
      push(c, { kind: "pivot", phase: "forward", j: j + 1, pivot: [r, j], pivotRow: r,
        text: head + "$" + Aij(j, j) + " = 0$ ولی $" + Aij(r, j) + " = " + T(c, c.M[r][j]) + " \\neq 0$ ⇒ $r = " + (r + 1) +
          "$ و $" + Rr(j) + " \\leftrightarrow " + Rr(r) + "$" });
      [c.M[j], c.M[r]] = [c.M[r], c.M[j]];
      if (c.plu) {
        for (let k = 0; k < j; k++) { const t = c.L[j][k]; c.L[j][k] = c.L[r][k]; c.L[r][k] = t; }
        [c.perm[j], c.perm[r]] = [c.perm[r], c.perm[j]];
      }
      push(c, { kind: "swap", phase: "forward", j: j + 1, swap: [j, r], pivot: [j, j], pivotRow: j,
        text: "جابه‌جایی: $" + Rr(j) + " \\leftrightarrow " + Rr(r) + "$" +
          (c.plu && j > 0 ? " (ضرایب محاسبه‌شدهٔ این دو سطر در $L$ هم جابه‌جا می‌شوند)" : "") });
    }
    const ls = [];
    for (let i = j + 1; i < n; i++) {
      const l = elim(c, i, j, "forward", "حذف: ", c.lu ? { lcell: [i, j] } : null);
      if (c.lu) {
        c.L[i][j] = l;
        c.steps[c.steps.length - 1].L = clone(c.L);
      }
      ls.push(l);
    }
    if (c.lu) {
      const Pj = eye(n);
      if (r !== j) { [Pj[j], Pj[r]] = [Pj[r], Pj[j]]; }
      const Mj = eye(n);
      for (let i = j + 1; i < n; i++) Mj[i][j] = ls[i - j - 1].neg();
      push(c, { kind: "iter", phase: "forward", j: j + 1, Pj, Mj,
        text: "پایان تکرار $j=" + (j + 1) + "$: ماتریس جایگشت $P^{(" + (j + 1) + ")}$ و ماتریس تبدیل گاوس $M^{(" + (j + 1) +
          ")}=I-\\ell^{(" + (j + 1) + ")}e_{" + (j + 1) + "}^{T}$؛ داریم $A^{(" + (j + 2) + ")}=M^{(" + (j + 1) + ")}P^{(" + (j + 1) + ")}A^{(" + (j + 1) + ")}$" });
    }
  }
  return true;
}

function backSub(c) {
  const n = c.n, M = c.M;
  if (M[n - 1][n - 1].isZero()) {
    push(c, { kind: "stop", phase: "back", pivot: [n - 1, n - 1],
      text: "بررسی: $" + Aij(n - 1, n - 1) + " = 0$ ⇒ دستگاه جواب یکتا ندارد" });
    return false;
  }
  c.x = Array(n).fill(null);
  push(c, { kind: "check", phase: "back", pivot: [n - 1, n - 1],
    text: "بررسی: $" + Aij(n - 1, n - 1) + " = " + T(c, M[n - 1][n - 1]) + " \\neq 0$ ⇒ دستگاه جواب یکتا دارد؛ جایگذاری پسرو را از $" + xName(c, n - 1) + "$ شروع می‌کنیم." });
  for (let i = n - 1; i >= 0; i--) {
    let s = M[i][n]; const terms = [];
    for (let j = i + 1; j < n; j++) {
      s = s.sub(M[i][j].mul(c.x[j]));
      terms.push(" - " + Pn(c, M[i][j]) + "\\cdot " + Pn(c, c.x[j]));
    }
    const xi = s.div(M[i][i]);
    c.x[i] = xi;
    let expr;
    if (!terms.length) expr = "\\frac{" + T(c, M[i][n]) + "}{" + T(c, M[i][i]) + "}";
    else if (M[i][i].isInt() && M[i][i].n > 0n) expr = "\\frac{1}{" + T(c, M[i][i]) + "}\\left(" + T(c, M[i][n]) + terms.join("") + "\\right)";
    else expr = "\\frac{" + T(c, M[i][n]) + terms.join("") + "}{" + T(c, M[i][i]) + "}";
    push(c, { kind: "back", phase: "back", i: i + 1, pivot: [i, i], target: i,
      text: "جایگذاری پسرو: $" + xName(c, i) + " = " + expr + " = " + T(c, xi) + "$" });
  }
  return true;
}

function initText(kind) {
  if (kind === "b") return "ماتریس افزوده اولیه $S = [A \\mid b]$";
  if (kind === "I") return "ماتریس افزوده اولیه $S = [A \\mid I]$";
  return "ماتریس اولیه $A$";
}

function gaussSteps(S, opts) {
  opts = opts || {};
  const c = makeCtx(S, opts);
  const n = c.n;
  const hasB = c.m === n + 1;
  push(c, { kind: "init", phase: "init", text: initText(opts.init || (hasB ? "b" : "I")) });
  const ok = forward(c);
  const result = { ok: false, M: null, U: null, c: null, x: null };
  if (!ok) return { steps: c.steps, result };
  if (opts.forwardOnly) {
    const rows = [];
    for (let i = 0; i < n; i++) if (c.M[i].slice(0, n).every(v => v.isZero())) rows.push(i);
    const rank = n - rows.length;
    const text = rows.length
      ? "سطر " + ORD[rows[0]] + " بخش $A$ صفر شد ⇒ $\\operatorname{rank}(A)=" + rank + "$ ⇒ $A$ وارون‌پذیر نیست"
      : "هیچ سطر صفری در بخش $A$ نیست ⇒ $\\operatorname{rank}(A)=" + n + "$ ⇒ $A$ وارون‌پذیر است";
    push(c, { kind: "final", phase: "final", text });
    result.ok = true; result.M = clone(c.M); result.zeroRows = rows; result.rank = rank;
    return { steps: c.steps, result };
  }
  if (hasB) {
    const ok2 = backSub(c);
    if (!ok2) return { steps: c.steps, result };
    result.ok = true; result.x = c.x.slice();
  } else result.ok = true;
  result.M = clone(c.M);
  result.U = sub(c.M, 0, n, 0, n);
  if (hasB) result.c = c.M.map(r => r[n]);
  return { steps: c.steps, result };
}

function gaussJordanSteps(S, opts) {
  opts = opts || {};
  const c = makeCtx(S, opts);
  const n = c.n;
  push(c, { kind: "init", phase: "init", text: initText("b") });
  const result = { ok: false, M: null, x: null };
  if (!forward(c)) return { steps: c.steps, result };
  if (c.M[n - 1][n - 1].isZero()) {
    push(c, { kind: "stop", phase: "up", pivot: [n - 1, n - 1],
      text: "بررسی: $" + Aij(n - 1, n - 1) + " = 0$ ⇒ دستگاه جواب یکتا ندارد" });
    return { steps: c.steps, result };
  }
  push(c, { kind: "check", phase: "up", pivot: [n - 1, n - 1],
    text: "فاز رو‌به‌جلو تمام شد و $" + Aij(n - 1, n - 1) + " = " + T(c, c.M[n - 1][n - 1]) + " \\neq 0$؛ حالا از ستون آخر به سمت چپ، بالای قطر را صفر می‌کنیم." });
  for (let j = n - 1; j >= 1; j--)
    for (let i = 0; i < j; i++) elim(c, i, j, "up", "حذف رو‌به‌بالا: ");
  c.x = c.M.map((r, i) => r[n].div(r[i]));
  const parts = c.x.map((xi, i) => "$" + xName(c, i) + "=\\frac{" + T(c, c.M[i][n]) + "}{" + T(c, c.M[i][i]) + "}=" + T(c, xi) + "$");
  push(c, { kind: "final", phase: "final", text: "جواب: " + parts.join("، ") });
  result.ok = true; result.M = clone(c.M); result.x = c.x.slice();
  return { steps: c.steps, result };
}

function gjInverseSteps(A, opts) {
  opts = opts || {};
  const c = makeCtx(augI(A), Object.assign({}, opts, { stopMsg: "ماتریس وارون‌پذیر نیست" }));
  const n = c.n;
  push(c, { kind: "init", phase: "init", text: initText("I") });
  const result = { ok: false, M: null, inverse: null };
  if (!forward(c)) return { steps: c.steps, result };
  const afterForward = clone(c.M);
  if (c.M[n - 1][n - 1].isZero()) {
    push(c, { kind: "stop", phase: "forward", pivot: [n - 1, n - 1],
      text: "بررسی: $" + Aij(n - 1, n - 1) + " = 0$ ⇒ ماتریس وارون‌پذیر نیست" });
    return { steps: c.steps, result };
  }
  push(c, { kind: "check", phase: "forward", pivot: [n - 1, n - 1],
    text: "فاز رو‌به‌جلو تمام شد و $" + Aij(n - 1, n - 1) + " \\neq 0$؛ حالا از آخر به اول قطر را ۱ و بالای قطر را صفر می‌کنیم." });
  for (let j = n - 1; j >= 0; j--) {
    if (!c.M[j][j].eq(ONE)) {
      const inv = ONE.div(c.M[j][j]), before = c.M[j].slice();
      c.M[j] = c.M[j].map(v => v.mul(inv));
      const changed = [];
      c.M[j].forEach((v, k) => { if (!v.eq(before[k])) changed.push([j, k]); });
      push(c, { kind: "scale", phase: "scale", j: j + 1, pivot: [j, j], target: j, changed,
        text: "قطر یک نیست: $" + Aij(j, j) + " = " + T(c, before[j]) + " \\neq 1$ ⇒ $" + Rr(j) + " \\leftarrow " + Pn(c, inv) + "\\,{" + Rr(j) + "}$" });
    }
    for (let i = 0; i < j; i++) elim(c, i, j, "up", "حذف رو‌به‌بالا: ");
  }
  const inv = sub(c.M, 0, n, n, 2 * n);
  push(c, { kind: "final", phase: "final", inverse: inv,
    text: "پایان: $[A \\mid I] \\to [I \\mid A^{-1}]$ ⇒ $A^{-1} = " + texMatScaled(inv, c.mode) + "$" });
  result.ok = true; result.M = clone(c.M); result.inverse = inv; result.afterForward = afterForward;
  return { steps: c.steps, result };
}

function luCore(A, plu, opts) {
  opts = opts || {};
  const c = makeCtx(A, Object.assign({}, opts, { stopMsg: "$A$ وارون‌پذیر نیست؛ الگوریتم متوقف شد" }));
  const n = c.n;
  c.lu = true; c.plu = !!plu;
  c.perm = Array.from({ length: n }, (_, i) => i);
  c.L = Array.from({ length: n }, (_, i) => Array.from({ length: n }, (_, j) => (i === j ? ONE : i < j ? ZERO : null)));
  push(c, { kind: "init", phase: "init",
    text: "ماتریس اولیه $A^{(1)}=A$؛ هدف: " + (plu ? "$PA=LU$ (با جابه‌جایی سطر)" : "$A=LU$ (بدون جابه‌جایی)") });
  const result = { ok: false };
  if (!forward(c)) return { steps: c.steps, result };
  const L = c.L.map(r => r.map(v => v || ZERO)), U = clone(c.M);
  const P = Array.from({ length: n }, (_, i) => Array.from({ length: n }, (_, j) => (c.perm[i] === j ? ONE : ZERO)));
  const verified = meq(mmul(P, A), mmul(L, U));
  push(c, { kind: "final", phase: "final", final: { L, U, P: plu ? P : null }, verified,
    text: plu ? "پایان: $U=A^{(n)}$ و $PA=LU$ برقرار است (درستی با ضرب ماتریس‌ها بررسی شد)."
      : "پایان: $U=A^{(n)}$ و $A=LU$ برقرار است (درستی با ضرب ماتریس‌ها بررسی شد)." });
  Object.assign(result, { ok: true, L, U, P, verified, perm: c.perm.slice() });
  return { steps: c.steps, result };
}
const luSteps = (A, opts) => luCore(A, false, opts);
const pluSteps = (A, opts) => luCore(A, true, opts);

/* ===================== Presets ===================== */
const EX43 = [[1, 1, -1, 1, 1], [1, 1, 1, -1, 1], [0, 1, 1, -1, -1], [1, -1, -1, -1, -1]];
const PRESETS = {
  ex43: { title: "مثال ۴.۳ — حذف گاوسی و جایگذاری پسرو", sub: "حذف رو‌به‌جلو با محورگیری «اولین درایهٔ ناصفر»، سپس جایگذاری پسرو", alg: "gauss", S: EX43 },
  ex53: { title: "مثال ۵.۳ — سهمی گذرنده از سه نقطه", sub: "مجهولات α، β، γ ضرایب سهمی‌اند", alg: "gauss", S: [[1, 1, 1, 1], [1, 2, 4, 2], [1, 3, 9, 0]], names: ["\\alpha", "\\beta", "\\gamma"] },
  ex63: { title: "مثال ۶.۳ — زنبورهای عسل", sub: "همان حذف گاوسی با نمایش اعشاری", alg: "gauss", display: "decimal",
    S: [["0.4", 0, 0, "0.2", 12], [0, "0.4", "0.3", "0.2", 25], [0, "0.3", "0.4", "0.2", 26], ["0.6", "0.3", "0.3", "0.4", 37]] },
  ex73: { title: "مثال ۷.۳ — روش حذفی گاوس-جردن", sub: "الگوریتم ۲.۳: بعد از فاز رو‌به‌جلو، بالای قطر هم صفر می‌شود", alg: "gj", S: EX43 },
  ex83: { title: "مثال ۸.۳ — آیا وارون‌پذیر است؟", sub: "حذف رو‌به‌جلو روی [A | I]؛ سطر صفر یعنی وارون‌ناپذیر", alg: "forward", init: "I", S: [[1, 2, 3, 1, 0, 0], [2, 3, 4, 0, 1, 0], [3, 4, 5, 0, 0, 1]] },
  ex93: { title: "مثال ۹.۳ — محاسبه وارون با گاوس-جردن", sub: "[A | I] به [I | A⁻¹] تبدیل می‌شود", alg: "inv", S: [[1, 2, -1], [2, 1, 0], [-1, 1, 2]] },
  ex113: { title: "مثال ۱۱.۳ — تجزیه LU", sub: "ضرایب حذف در ماتریس L جمع می‌شوند", alg: "lu", S: [[1, 1, -1, 1], [1, -1, -1, -1], [1, 1, 1, -1], [0, 1, 1, -1]] },
  ex133: { title: "مثال ۱۳.۳ — تجزیه PA = LU", sub: "با جابه‌جایی سطرها؛ ضرایب قبلی هم جابه‌جا می‌شوند", alg: "plu", S: [[0, 1, -1, 1], [1, 1, -1, 2], [-1, -1, 1, 0], [1, 2, 0, 2]] },
};

function buildSteps(cfg) {
  const S = mat(cfg.S);
  const o = { mode: cfg.display || "fraction", names: cfg.names || null, init: cfg.init };
  switch (cfg.alg) {
    case "gauss": return gaussSteps(S, o);
    case "forward": return gaussSteps(S, Object.assign(o, { forwardOnly: true }));
    case "gj": return gaussJordanSteps(S, o);
    case "inv": return gjInverseSteps(S, o);
    case "lu": return luSteps(S, o);
    case "plu": return pluSteps(S, o);
  }
  return gaussSteps(S, o);
}

/* ===================== DOM helpers (browser only) ===================== */
const toFa = s => String(s).replace(/[0-9]/g, d => "۰۱۲۳۴۵۶۷۸۹"[d]);
const esc = s => String(s).replace(/[&<>"]/g, c => ({ "&": "&amp;", "<": "&lt;", ">": "&gt;", '"': "&quot;" }[c]));

let uidCounter = 0;
const uid = () => ++uidCounter;
const kcache = new Map();
function ktx(t) {
  let v = kcache.get(t);
  if (v === undefined) {
    try { v = typeof katex !== "undefined" ? katex.renderToString(t, { throwOnError: false, displayMode: false }) : esc(t); }
    catch (e) { v = esc(t); }
    kcache.set(t, v);
  }
  return v;
}
function rich(s) {
  return String(s).split("$").map((p, i) => (i % 2 ? ktx(p) : esc(toFa(p)))).join("");
}
const reduceMotion = () => typeof matchMedia === "function" && matchMedia("(prefers-reduced-motion: reduce)").matches;

function h(tag, attrs, ...kids) {
  const e = document.createElement(tag);
  if (attrs) for (const k in attrs) {
    const v = attrs[k];
    if (v == null || v === false) continue;
    if (k === "class") e.className = v;
    else if (k === "html") e.innerHTML = v;
    else if (k.startsWith("on")) e.addEventListener(k.slice(2), v);
    else e.setAttribute(k, v === true ? "" : v);
  }
  kids.flat(Infinity).forEach(c => { if (c == null || c === false) return; e.append(c.nodeType ? c : document.createTextNode(c)); });
  return e;
}
const SVGNS = "http://www.w3.org/2000/svg";
function sv(tag, attrs, ...kids) {
  const e = document.createElementNS(SVGNS, tag);
  if (attrs) for (const k in attrs) { if (attrs[k] != null) e.setAttribute(k, attrs[k]); }
  kids.flat(Infinity).forEach(c => { if (c == null || c === false) return; e.append(c.nodeType ? c : document.createTextNode(c)); });
  return e;
}
function card(title, sub, extraClass) {
  const body = h("div", { class: "w-body" });
  const root = h("div", { class: "w-card " + (extraClass || "") },
    h("div", { class: "w-head" }, h("div", { class: "w-title" }, title), sub ? h("div", { class: "w-sub" }, sub) : null), body);
  return { root, body };
}
function cellHtml(v, mode) {
  if (v == null) return '<span class="w-dim">·</span>';
  return ktx(v.toTeX(mode));
}
/* Matrix grid with CSS brackets. o: {labels, sep, mode, cls(i,j), cell(i,j)} */
function matEl(M, o) {
  o = o || {};
  const n = M.length, m = M[0].length;
  const g = h("div", { class: "w-mx" + (o.labels ? " w-lab" : "") });
  g.style.gridTemplateColumns = (o.labels ? "auto " : "") + "repeat(" + m + ", minmax(2.6em, auto))";
  for (let i = 0; i < n; i++) {
    if (o.labels) g.append(h("span", { class: "w-rl", "data-r": i }, "R" + (i + 1)));
    for (let j = 0; j < m; j++) {
      let cl = "w-c";
      if (j === 0) cl += " w-l";
      if (j === m - 1) cl += " w-r";
      if (i === 0) cl += " w-t";
      if (i === n - 1) cl += " w-b";
      if (o.sep != null && j === o.sep) cl += " w-sep";
      if (o.cls) cl += " " + o.cls(i, j);
      g.append(h("span", { class: cl, "data-r": i, html: o.cell ? o.cell(i, j) : cellHtml(M[i][j], o.mode) }));
    }
  }
  return h("div", { class: "w-scroll", dir: "ltr" }, g);
}
const fig = (capTeX, el, cls) => h("figure", { class: "w-fig " + (cls || ""), dir: "ltr" }, h("figcaption", { class: "w-cap", html: ktx(capTeX) }), el);

/* ===================== Widget: stepper ===================== */
function mountStepper(host, cfg, mo) {
  mo = mo || {};
  const built = buildSteps(cfg), steps = built.steps, mode = cfg.display || "fraction";
  const A0 = mat(cfg.S);
  const total = steps.length - 1;
  let idx = 0, timer = null;

  const textEl = h("div", { class: "w-text", "aria-live": "polite" });
  const matBox = h("div", { class: "w-matbox" });
  const side = h("div", { class: "w-side", dir: "ltr" });
  const info = h("span", { class: "w-count" });
  const barFill = h("i");
  const bar = h("div", { class: "w-bar", role: "progressbar", "aria-valuemin": "0", "aria-valuemax": String(total) }, barFill);
  const bPrev = h("button", { class: "btn", type: "button", onclick: () => { stop(); go(idx - 1); } }, "◀ قبلی");
  const bPlay = h("button", { class: "btn", type: "button", "aria-pressed": "false", onclick: () => (timer ? stop() : play()) }, "پخش خودکار");
  const bReset = h("button", { class: "btn", type: "button", onclick: () => { stop(); go(0); } }, "از اول");
  const bNext = h("button", { class: "btn primary", type: "button", onclick: () => { stop(); go(idx + 1); } }, "بعدی ▶");
  const ctrl = h("div", { class: "w-ctrl" }, bPrev, bPlay, bReset, bNext);
  const status = h("div", { class: "w-status" }, info, bar);

  const root = mo.bare ? h("div", { class: "w-stepper", tabindex: "0", "aria-label": cfg.title || "گام‌به‌گام" }) :
    h("div", { class: "w-card w-stepper", tabindex: "0", "aria-label": cfg.title });
  if (!mo.bare) root.append(h("div", { class: "w-head" }, h("div", { class: "w-title" }, cfg.title), cfg.sub ? h("div", { class: "w-sub" }, cfg.sub) : null));
  root.append(h("div", { class: "w-body" }, textEl, matBox, side, ctrl, status));
  host.textContent = "";
  host.append(root);

  function stop() { if (timer) { clearInterval(timer); timer = null; } bPlay.setAttribute("aria-pressed", "false"); bPlay.textContent = "پخش خودکار"; }
  function play() {
    if (idx >= total) go(0);
    bPlay.setAttribute("aria-pressed", "true"); bPlay.textContent = "توقف";
    timer = setInterval(() => { if (idx >= total) { stop(); return; } go(idx + 1); if (idx >= total) stop(); }, 1600);
  }
  function go(i) {
    i = Math.max(0, Math.min(total, i));
    const prev = idx; idx = i; render(prev);
  }
  function xPanel(st) {
    const list = h("div", { class: "w-xp", dir: "ltr" });
    st.x.forEach((v, i) => {
      const name = cfg.names ? cfg.names[i] : "x_{" + (i + 1) + "}";
      list.append(h("div", { class: "w-xr" + (st.kind === "back" && st.target === i ? " w-new" : ""),
        html: ktx(name + (v ? "=" + v.toTeX(mode) : "=\\,?")) }));
    });
    return h("figure", { class: "w-fig", dir: "ltr" }, h("figcaption", { class: "w-cap", html: rich("بردار جواب") }), list);
  }
  function render(prev) {
    const st = steps[idx];
    textEl.innerHTML = rich(st.text);
    const n = st.M.length, m = st.M[0].length;
    const ch = new Set((st.changed || []).map(p => p[0] + "," + p[1]));
    const zs = new Set((st.zeros || []).map(p => p[0] + "," + p[1]));
    const cls = (i, j) => {
      let c = "";
      if (st.pivotRow === i) c += " w-pr";
      if (st.target === i || (st.swap && st.swap.includes(i))) c += " w-tr";
      if (st.pivot && st.pivot[0] === i && st.pivot[1] === j) c += " w-pv";
      if (ch.has(i + "," + j)) c += " w-ch";
      if (zs.has(i + "," + j)) c += " w-z";
      return c;
    };
    const mx = matEl(st.M, { labels: true, sep: m > n ? n : null, mode, cls });
    matBox.textContent = "";
    matBox.append(mx);
    mx.setAttribute("aria-label", "ماتریس در گام " + toFa(idx));

    side.textContent = "";
    if (st.L) {
      const lc = st.lcell;
      side.append(fig("L", matEl(st.L, { mode, cls: (i, j) => (lc && lc[0] === i && lc[1] === j ? "w-new" : "") })));
    }
    if (st.perm && cfg.alg === "plu") {
      const P = st.perm.map(pi => st.perm.map((_, j) => (j === pi ? ONE : ZERO)));
      side.append(fig("P", matEl(P, { mode })));
    }
    if (st.Pj) {
      side.append(fig("P^{(" + st.j + ")}", matEl(st.Pj, { mode })));
      side.append(fig("M^{(" + st.j + ")}", matEl(st.Mj, { mode })));
    }
    if (st.x) side.append(xPanel(st));
    if (st.final) {
      const f = st.final, eq = (t) => h("span", { class: "w-op", html: ktx(t) });
      const row = h("div", { class: "w-eq", dir: "ltr" });
      if (f.P) row.append(fig("P", matEl(f.P, { mode })), eq("\\times"), fig("A", matEl(A0, { mode })));
      else row.append(fig("A", matEl(A0, { mode })));
      row.append(eq("="), fig("L", matEl(f.L, { mode })), eq("\\times"), fig("U", matEl(f.U, { mode })));
      side.append(row);
    }
    if (st.inverse) side.append(fig("A^{-1}", matEl(st.inverse, { mode })));

    info.textContent = "گام " + toFa(idx) + " از " + toFa(total);
    barFill.style.width = (total ? (idx / total) * 100 : 100) + "%";
    bar.setAttribute("aria-valuenow", String(idx));
    bPrev.disabled = idx === 0; bNext.disabled = idx === total;

    if (st.swap && idx > prev && !reduceMotion()) {
      try {
        const tops = st.swap.map(r => mx.querySelector('.w-c[data-r="' + r + '"]').getBoundingClientRect().top);
        st.swap.forEach((r, k) => {
          const d = tops[1 - k] - tops[k];
          mx.querySelectorAll('[data-r="' + r + '"]').forEach(el => el.animate([{ transform: "translateY(" + d + "px)" }, { transform: "translateY(0)" }], { duration: 480, easing: "ease" }));
        });
      } catch (e) { /* animation is optional */ }
    }
  }
  root.addEventListener("keydown", e => {
    if (e.target !== root && /^(INPUT|SELECT|TEXTAREA)$/.test(e.target.tagName)) return;
    if (e.key === "ArrowLeft") { e.preventDefault(); stop(); go(idx + 1); }
    else if (e.key === "ArrowRight") { e.preventDefault(); stop(); go(idx - 1); }
  });
  render(0);
}

/* ===================== Widget: custom (editor + stepper) ===================== */
function mountCustom(host) {
  const id = uid();
  const { root, body } = card("حل دلخواه — گام‌به‌گام", "ماتریس، اندازه و الگوریتم را خودتان انتخاب کنید");
  const vals = { A: [["1", "1", "1"], ["1", "2", "4"], ["1", "3", "9"]], b: ["1", "2", "0"] };
  let n = 3;
  const selN = h("select", { id: "w-cu" + id + "-n" }, [2, 3, 4, 5].map(k => h("option", { value: k, selected: k === n }, toFa(k))));
  const ALGS = [["gauss", "حذف گاوسی + جایگذاری پسرو"], ["gj", "گاوس-جردن"], ["lu", "تجزیه LU"], ["plu", "تجزیه PA=LU"], ["inv", "وارون با گاوس-جردن"]];
  const selA = h("select", { id: "w-cu" + id + "-alg" }, ALGS.map(a => h("option", { value: a[0] }, a[1])));
  const selD = h("select", { id: "w-cu" + id + "-disp" }, h("option", { value: "fraction" }, "کسری"), h("option", { value: "decimal" }, "اعشاری"));
  const grid = h("div", { class: "w-gridin", dir: "ltr" });
  const msg = h("div", { class: "w-msg", role: "alert" });
  const out = h("div", { class: "w-out" });
  const inputs = { A: [], b: [] };

  function readVals() {
    inputs.A.forEach((r, i) => r.forEach((el, j) => { vals.A[i] = vals.A[i] || []; vals.A[i][j] = el.value; }));
    inputs.b.forEach((el, i) => { vals.b[i] = el.value; });
  }
  function build() {
    grid.textContent = ""; inputs.A = []; inputs.b = [];
    grid.style.gridTemplateColumns = "repeat(" + n + ", auto) auto auto";
    for (let i = 0; i < n; i++) {
      inputs.A[i] = [];
      for (let j = 0; j < n; j++) {
        const el = h("input", { id: "w-cu" + id + "-a" + i + j, class: "w-in", type: "text", inputmode: "text", dir: "ltr", autocomplete: "off", "aria-label": "a" + (i + 1) + (j + 1), value: (vals.A[i] && vals.A[i][j]) != null ? vals.A[i][j] : "0" });
        inputs.A[i][j] = el; grid.append(el);
      }
      grid.append(h("span", { class: "w-vsep", "aria-hidden": "true" }));
      const eb = h("input", { id: "w-cu" + id + "-b" + i, class: "w-in", type: "text", inputmode: "text", dir: "ltr", autocomplete: "off", "aria-label": "b" + (i + 1), value: vals.b[i] != null ? vals.b[i] : "0" });
      inputs.b[i] = eb; grid.append(eb);
    }
    syncB();
  }
  function syncB() { const on = selA.value === "gauss" || selA.value === "gj"; inputs.b.forEach(el => { el.disabled = !on; }); }
  function run() {
    readVals();
    let bad = 0;
    const parse = el => { const v = parseFrac(el.value); el.classList.toggle("w-invalid", !v); if (!v) bad++; return v; };
    const A = inputs.A.map(r => r.map(parse));
    const useB = selA.value === "gauss" || selA.value === "gj";
    inputs.b.forEach(el => el.classList.remove("w-invalid"));
    const b = useB ? inputs.b.map(parse) : null;
    if (bad) { msg.textContent = "ورودی نامعتبر است: فقط عدد صحیح، کسر مثل 1/3 یا اعشار مثل 0.25 مجاز است (خانه‌های قرمز را اصلاح کنید)."; return; }
    msg.textContent = "";
    const S = useB ? A.map((r, i) => r.concat([b[i]])) : A;
    const alg = selA.value;
    mountStepper(out, { title: "حل دلخواه", alg, S, display: selD.value, init: alg === "inv" ? "I" : undefined }, { bare: true });
  }
  selN.addEventListener("change", () => { readVals(); n = +selN.value; build(); });
  selA.addEventListener("change", syncB);
  body.append(
    h("div", { class: "w-row" },
      h("label", { for: selN.id }, "اندازه n"), selN,
      h("label", { for: selA.id }, "الگوریتم"), selA,
      h("label", { for: selD.id }, "نمایش"), selD),
    h("div", { class: "w-scroll", dir: "ltr" }, grid),
    h("div", { class: "w-row" }, h("button", { class: "btn primary", type: "button", onclick: run }, "اجرا کن")),
    msg, out);
  host.textContent = "";
  host.append(root);
  build();
  run();
}

/* ===================== Number formatting ===================== */
function fmt(x) {
  if (!isFinite(x)) return "∞";
  const neg = x < 0 ? "−" : "";
  x = Math.abs(x);
  let s = x >= 1e6 ? Math.round(x).toString() : String(+x.toPrecision(6));
  if (/e/.test(s)) s = x.toFixed(10).replace(/0+$/, "").replace(/\.$/, "");
  const parts = s.split(".");
  const ii = parts[0].replace(/\B(?=(\d{3})+(?!\d))/g, "٬");
  return neg + toFa(ii + (parts[1] ? "٫" + parts[1] : ""));
}
const num3 = x => (x >= 10 ? Math.round(x) : +x.toPrecision(3));
const SUP = "⁰¹²³⁴⁵⁶⁷⁸⁹";
const sup = k => String(k).replace(/[0-9]/g, d => SUP[d]);

function mkSlider(id, min, max, step, value, label, onInput) {
  const out = h("output", { class: "w-num", for: id });
  const input = h("input", { id, type: "range", min, max, step, value, class: "w-range" });
  input.addEventListener("input", () => onInput(+input.value));
  return { input, out, row: h("div", { class: "w-srow" }, h("label", { for: id }, label), input, out) };
}

/* ===================== Widget: illcond ===================== */
function mountIllcond(host) {
  const id = uid();
  const { root, body } = card("آزمایشگاه بدوضعی — مثال ۱.۳", "دو خط تقریباً موازی؛ باقی‌مانده کوچک ولی خطا بزرگ");
  const W = 520, H = 372, L = 46, R = 504, Tp = 12, Bt = 338;
  const X = v => L + ((v + 1) / 5) * (R - L), Y = v => Bt - ((v + 1) / 3.5) * (Bt - Tp);
  let sE = -4, t = 1;

  const svg = sv("svg", { viewBox: "0 0 " + W + " " + H, class: "w-svg", role: "img", "aria-label": "نمودار دو خط و جواب دقیق و تقریبی" });
  const clipId = "w-ic" + id + "-clip";
  svg.append(sv("defs", null, sv("clipPath", { id: clipId }, sv("rect", { x: L, y: Tp, width: R - L, height: Bt - Tp }))));
  for (let v = -1; v <= 4; v++) {
    svg.append(sv("line", { x1: X(v), y1: Tp, x2: X(v), y2: Bt, class: "w-gridl" }));
    svg.append(sv("text", { x: X(v), y: Bt + 16, "text-anchor": "middle", class: "w-tick" }, v < 0 ? "−" + Math.abs(v) : String(v)));
  }
  [-1, 0, 1, 2].forEach(v => {
    svg.append(sv("line", { x1: L, y1: Y(v), x2: R, y2: Y(v), class: "w-gridl" }));
    svg.append(sv("text", { x: L - 7, y: Y(v) + 4, "text-anchor": "end", class: "w-tick" }, v < 0 ? "−" + Math.abs(v) : String(v)));
  });
  svg.append(sv("line", { x1: X(0), y1: Tp, x2: X(0), y2: Bt, class: "w-axisl" }));
  svg.append(sv("line", { x1: L, y1: Y(0), x2: R, y2: Y(0), class: "w-axisl" }));
  svg.append(sv("text", { x: R, y: Bt + 32, "text-anchor": "end", class: "w-tick" }, "x₁"));
  svg.append(sv("text", { x: L - 7, y: Tp + 8, "text-anchor": "end", class: "w-tick" }, "x₂"));
  const g = sv("g", { "clip-path": "url(#" + clipId + ")" });
  const l1 = sv("line", { class: "w-ln w-s1" }), l2 = sv("line", { class: "w-ln w-s2" });
  const seg = sv("line", { class: "w-seg" });
  g.append(l1, l2, seg);
  svg.append(g);
  const lab1 = sv("text", { class: "w-dl w-t1", "text-anchor": "end" }, "خط ۱");
  const lab2 = sv("text", { class: "w-dl w-t2", "text-anchor": "end" }, "خط ۲");
  const labX = sv("text", { class: "w-pl" }, "x (جواب دقیق)");
  const labT = sv("text", { class: "w-pl" }, "x̃ (تقریب)");
  const labE = sv("text", { class: "w-pl w-errl", "text-anchor": "middle" }, "خطا");
  const dotX = sv("circle", { r: 6, class: "w-dot-fill" }), dotT = sv("circle", { r: 6, class: "w-dot-hollow" });
  const hitX = sv("circle", { r: 16, class: "w-hit", tabindex: "0", "aria-label": "جواب دقیق" });
  const hitT = sv("circle", { r: 16, class: "w-hit", tabindex: "0", "aria-label": "جواب تقریبی" });
  svg.append(lab1, lab2, labE, dotX, dotT, labX, labT, hitX, hitT);

  const tip = h("div", { class: "w-tip", hidden: true, dir: "ltr" });
  const plot = h("div", { class: "w-plot", dir: "ltr" }, svg, tip);
  function showTip(pt, name) {
    const p = pt();
    tip.hidden = false;
    tip.textContent = name + " = (" + (+p[0].toFixed(3)) + ", " + (+p[1].toFixed(3)) + ")";
    tip.style.left = (X(p[0]) / W) * 100 + "%"; tip.style.top = (Y(p[1]) / H) * 100 + "%";
  }
  const pX = () => [1, 1], pT = () => [1 + 2 * t, 1 - t];
  [[hitX, pX, "x"], [hitT, pT, "x̃"]].forEach(([el, pt, nm]) => {
    el.addEventListener("pointerenter", () => showTip(pt, nm));
    el.addEventListener("pointerdown", () => showTip(pt, nm));
    el.addEventListener("focus", () => showTip(pt, nm));
    el.addEventListener("pointerleave", e => { if (e.pointerType === "mouse") tip.hidden = true; });
    el.addEventListener("blur", () => { tip.hidden = true; });
  });
  svg.addEventListener("pointerdown", e => { if (!e.target.classList.contains("w-hit")) tip.hidden = true; });

  const cR = h("div", { class: "w-num" }), cE = h("div", { class: "w-num" }), cK = h("div", { class: "w-num" });
  const mk = (lbl, tex, val) => h("div", { class: "w-stat" }, h("div", { class: "w-statl" }, lbl, " ", h("span", { html: ktx(tex), dir: "ltr" })), val);
  const cards = h("div", { class: "w-cards" },
    mk("باقی‌مانده", "\\|r\\|_\\infty", cR), mk("خطا", "\\|x-\\tilde{x}\\|_\\infty", cE), mk("عدد وضعیت", "\\kappa_\\infty(A)", cK));
  const sentence = h("p", { class: "w-sentence", "aria-live": "polite" });
  const slopes = h("p", { class: "w-note" });
  const sE_ = mkSlider("w-ic" + id + "-e", -4, 0, 0.01, -4, "مقدار ε (مقیاس لگاریتمی)", v => { sE = v; update(); });
  const sT_ = mkSlider("w-ic" + id + "-t", -1, 1, 0.05, 1, "پارامتر t در x̃", v => { t = v; update(); });
  const legend = h("div", { class: "w-legend", dir: "ltr" },
    h("span", { class: "w-lg" }, h("i", { class: "w-sw w-b1" }), h("span", { html: ktx("x_1+2x_2=3") })),
    h("span", { class: "w-lg" }, h("i", { class: "w-sw w-b2" }), h("span", { html: ktx("(1+\\varepsilon)x_1+2x_2=3+\\varepsilon") })));

  function update() {
    const eps = +Math.pow(10, sE).toPrecision(3);
    const yEnd = x => (3 + eps - (1 + eps) * x) / 2;
    l1.setAttribute("x1", X(-1)); l1.setAttribute("y1", Y(2)); l1.setAttribute("x2", X(4)); l1.setAttribute("y2", Y(-0.5));
    l2.setAttribute("x1", X(-1)); l2.setAttribute("y1", Y(yEnd(-1))); l2.setAttribute("x2", X(4)); l2.setAttribute("y2", Y(yEnd(4)));
    const y1l = Y(-0.5) - 10;
    const y2l = Math.min(Bt - 4, Math.max(Y(Math.max(yEnd(4), -0.98)) + 18, y1l + 22));
    lab1.setAttribute("x", R - 4); lab1.setAttribute("y", y1l);
    lab2.setAttribute("x", R - 4); lab2.setAttribute("y", y2l);
    const [tx, ty] = pT();
    dotX.setAttribute("cx", X(1)); dotX.setAttribute("cy", Y(1));
    hitX.setAttribute("cx", X(1)); hitX.setAttribute("cy", Y(1));
    dotT.setAttribute("cx", X(tx)); dotT.setAttribute("cy", Y(ty));
    hitT.setAttribute("cx", X(tx)); hitT.setAttribute("cy", Y(ty));
    seg.setAttribute("x1", X(1)); seg.setAttribute("y1", Y(1)); seg.setAttribute("x2", X(tx)); seg.setAttribute("y2", Y(ty));
    const flip = X(tx) > R - 110;
    labX.setAttribute("x", X(1) + 11); labX.setAttribute("y", Y(1) - 11);
    labT.setAttribute("x", X(tx) + (flip ? -11 : 11)); labT.setAttribute("y", Y(ty) + 20);
    labT.setAttribute("text-anchor", flip ? "end" : "start");
    labE.setAttribute("x", (X(1) + X(tx)) / 2); labE.setAttribute("y", (Y(1) + Y(ty)) / 2 - 10);
    labE.style.display = Math.abs(t) < 0.05 ? "none" : "";
    const r = 2 * eps * Math.abs(t), err = 2 * Math.abs(t), kap = (2 * (3 + eps)) / eps;
    cR.textContent = fmt(r); cE.textContent = fmt(err); cK.textContent = fmt(kap);
    sE_.out.textContent = "ε = " + fmt(eps); sT_.out.textContent = "t = " + fmt(t);
    sentence.textContent = Math.abs(t) < 1e-9
      ? "x̃ همان جواب دقیق است؛ باقی‌مانده و خطا هر دو صفرند."
      : "باقی‌مانده " + fmt(r) + " است ولی خطا " + fmt(err) + " است؛ یعنی خطا " + fmt(err / r) + " برابر باقی‌مانده است.";
    slopes.innerHTML = "شیب دو خط: " + ktx("-0.5") + " و " + ktx("-\\frac{1+\\varepsilon}{2}=-" + (+((1 + eps) / 2).toFixed(6)));
    if (!tip.hidden) tip.hidden = true;
  }
  body.append(legend, plot, sE_.row, sT_.row, cards, sentence, slopes);
  host.textContent = ""; host.append(root);
  update();
}

/* ===================== Widget: cost ===================== */
function factBig(k) { let r = 1n; for (let i = 2n; i <= BigInt(k); i++) r *= i; return r; }
function log10Big(b) {
  const s = b.toString();
  if (s.length <= 15) return Math.log10(Number(s));
  return Math.log10(Number(s.slice(0, 15))) + s.length - 15;
}
function sciHtml(b) {
  const s = b.toString();
  if (s.length <= 9) return '<span dir="ltr">' + toFa(s.replace(/\B(?=(\d{3})+(?!\d))/g, "٬")) + "</span>";
  let exp = s.length - 1, mant = Math.round(Number(s.slice(0, 4)) / 10) / 100;
  if (mant >= 10) { mant = 1; exp++; }
  return '<span dir="ltr">' + toFa(mant.toFixed(2).replace(".", "٫")) + " × " + toFa(10) + "<sup>" + toFa(exp) + "</sup></span>";
}
function timeStr(sec) {
  if (sec < 1e-9) return "کمتر از یک میلیاردم ثانیه";
  if (sec < 60) return fmt(+sec.toPrecision(3)) + " ثانیه";
  if (sec < 3600) return fmt(num3(sec / 60)) + " دقیقه";
  if (sec < 86400) return fmt(num3(sec / 3600)) + " ساعت";
  if (sec < 365.25 * 86400) return fmt(num3(sec / 86400)) + " روز";
  return fmt(num3(sec / (365.25 * 86400))) + " سال";
}
const costData = n => {
  const N = BigInt(n);
  return { n, cramer: factBig(n + 1), gauss: (N ** 3n + 3n * N * N - N) / 3n, gj: (N * N * (N + 1n)) / 2n };
};

function mountCost(host) {
  const id = uid();
  const { root, body } = card("هزینهٔ کرامر در برابر حذف گاوسی", "تعداد ضرب و تقسیم بر حسب n (مقیاس لگاریتمی)");
  const W = 560, H = 350, L = 52, R = 452, Tp = 12, Bt = 306, N0 = 2, N1 = 30, YM = 35;
  const X = n => L + ((n - N0) / (N1 - N0)) * (R - L), Y = v => Bt - (v / YM) * (Bt - Tp);
  const data = []; for (let n = N0; n <= N1; n++) data.push(costData(n));
  const SER = [
    { key: "cramer", name: "کرامر", cls: "s2", f: d => d.cramer },
    { key: "gauss", name: "حذف گاوسی", cls: "s1", f: d => d.gauss },
    { key: "gj", name: "گاوس-جردن", cls: "s3", f: d => d.gj },
  ];
  let sel = 28;

  const svg = sv("svg", { viewBox: "0 0 " + W + " " + H, class: "w-svg", role: "img", "aria-label": "نمودار تعداد عمل بر حسب n با مقیاس لگاریتمی" });
  for (let k = 0; k <= YM; k += 5) {
    svg.append(sv("line", { x1: L, y1: Y(k), x2: R, y2: Y(k), class: "w-gridl" }));
    svg.append(sv("text", { x: L - 6, y: Y(k) + 4, "text-anchor": "end", class: "w-tick" }, "10" + sup(k)));
  }
  [2, 5, 10, 15, 20, 25, 30].forEach(n => {
    svg.append(sv("line", { x1: X(n), y1: Bt, x2: X(n), y2: Bt + 4, class: "w-axisl" }));
    svg.append(sv("text", { x: X(n), y: Bt + 18, "text-anchor": "middle", class: "w-tick" }, String(n)));
  });
  svg.append(sv("line", { x1: L, y1: Bt, x2: R, y2: Bt, class: "w-axisl" }), sv("line", { x1: L, y1: Tp, x2: L, y2: Bt, class: "w-axisl" }));
  svg.append(sv("text", { x: R, y: Bt + 34, "text-anchor": "end", class: "w-tick" }, "n"));
  const markers = {};
  SER.forEach(s => {
    const pts = data.map(d => X(d.n) + "," + Y(log10Big(s.f(d)))).join(" ");
    svg.append(sv("polyline", { points: pts, class: "w-ln w-" + s.cls, fill: "none" }));
    markers[s.key] = sv("circle", { r: 5, class: "w-mk w-f" + s.cls });
  });
  const ends = SER.map(s => ({ s, y: Y(log10Big(s.f(data[data.length - 1]))) })).sort((a, b) => a.y - b.y);
  for (let i = 1; i < ends.length; i++) ends[i].ly = Math.max(ends[i].y, (ends[i - 1].ly ?? ends[i - 1].y) + 15);
  ends[0].ly = ends[0].y;
  ends.forEach(e => svg.append(sv("text", { x: R + 8, y: e.ly + 4, class: "w-dl w-t" + e.s.cls.slice(1) }, e.s.name)));
  SER.forEach(s => svg.append(markers[s.key]));
  const cross = sv("line", { y1: Tp, y2: Bt, class: "w-cross", style: "display:none" });
  const hov = {}; SER.forEach(s => { hov[s.key] = sv("circle", { r: 4, class: "w-hv w-f" + s.cls, style: "display:none" }); });
  svg.append(cross, ...SER.map(s => hov[s.key]));
  const catcher = sv("rect", { x: L, y: Tp, width: R - L, height: Bt - Tp, class: "w-catch" });
  svg.append(catcher);
  const tip = h("div", { class: "w-tip w-tip-wide", hidden: true });
  const plot = h("div", { class: "w-plot", dir: "ltr" }, svg, tip);

  function hover(e) {
    const rc = svg.getBoundingClientRect();
    const vx = ((e.clientX - rc.left) / rc.width) * W;
    const n = Math.max(N0, Math.min(N1, Math.round(N0 + ((vx - L) / (R - L)) * (N1 - N0))));
    const d = data[n - N0];
    cross.style.display = ""; cross.setAttribute("x1", X(n)); cross.setAttribute("x2", X(n));
    SER.forEach(s => { const c = hov[s.key]; c.style.display = ""; c.setAttribute("cx", X(n)); c.setAttribute("cy", Y(log10Big(s.f(d)))); });
    tip.hidden = false;
    tip.innerHTML = "<b>n = " + toFa(n) + "</b>" + SER.map(s => '<div class="w-tl"><i class="w-sw w-b' + s.cls.slice(1) + '"></i>' + s.name + ": " + sciHtml(s.f(d)) + "</div>").join("");
    const frac = X(n) / W;
    tip.style.left = frac * 100 + "%"; tip.style.top = "6%";
    tip.style.transform = frac > 0.55 ? "translateX(calc(-100% - 10px))" : "translateX(10px)";
  }
  function unhover() { cross.style.display = "none"; SER.forEach(s => { hov[s.key].style.display = "none"; }); tip.hidden = true; }
  catcher.addEventListener("pointermove", hover);
  catcher.addEventListener("pointerdown", hover);
  catcher.addEventListener("pointerleave", unhover);

  const legend = h("div", { class: "w-legend" }, SER.map(s => h("span", { class: "w-lg" }, h("i", { class: "w-sw w-b" + s.cls.slice(1) }), s.name === "کرامر" ? "کرامر ≈ (n+1)! ضرب" : s.name === "حذف گاوسی" ? "حذف گاوسی + جایگذاری پسرو" : "گاوس-جردن")));
  const tbody = h("tbody");
  const table = h("table", { class: "w-table" }, h("thead", null, h("tr", null, h("th", null, "روش"), h("th", null, "تعداد ضرب/تقسیم"), h("th", null, "زمان با ۱۰²⁰ عمل در ثانیه"))), tbody);
  const callout = h("p", { class: "w-callout", "aria-live": "polite" });
  const sl = mkSlider("w-co" + id + "-n", N0, N1, 1, sel, "مقدار n", v => { sel = v; update(); });
  const fullBody = h("tbody");
  const fullWrap = h("div", { class: "w-scroll w-fulltable", hidden: true },
    h("table", { class: "w-table" }, h("thead", null, h("tr", null, h("th", null, "n"), SER.map(s => h("th", null, s.name)))), fullBody));
  data.forEach(d => fullBody.append(h("tr", null, h("td", null, toFa(d.n)), SER.map(s => h("td", { html: sciHtml(s.f(d)) })))));
  const tog = h("button", { class: "btn", type: "button", "aria-expanded": "false", "aria-controls": "w-co" + id + "-full", onclick: () => {
    const open = fullWrap.hidden; fullWrap.hidden = !open; tog.setAttribute("aria-expanded", String(open));
    tog.textContent = open ? "پنهان کردن جدول" : "نمایش جدول";
  } }, "نمایش جدول");
  fullWrap.id = "w-co" + id + "-full";

  function update() {
    const d = data[sel - N0];
    SER.forEach(s => { const c = markers[s.key]; c.setAttribute("cx", X(sel)); c.setAttribute("cy", Y(log10Big(s.f(d)))); });
    sl.out.textContent = "n = " + toFa(sel);
    tbody.textContent = "";
    SER.forEach(s => tbody.append(h("tr", null,
      h("td", null, h("i", { class: "w-sw w-b" + s.cls.slice(1) }), " " + s.name),
      h("td", { html: sciHtml(s.f(d)) }),
      h("td", null, timeStr(Number(s.f(d)) / 1e20)))));
    const nf = factBig(sel);
    callout.innerHTML = "طبق کتاب، محاسبهٔ یک دترمینان " + toFa(sel) + "×" + toFa(sel) + " با بسط کهاد حدود <span dir=\"ltr\">" + toFa(sel) + "!</span> ضرب لازم دارد: <span dir=\"ltr\">" + toFa(sel) + "! ≈</span> " + sciHtml(nf) + " ⇒ حدود " + timeStr(Number(nf) / 1e20) + " با ۱۰²⁰ عمل در ثانیه.";
  }
  body.append(legend, plot, sl.row, h("div", { class: "w-scroll" }, table), callout, h("div", { class: "w-row" }, tog), fullWrap);
  host.textContent = ""; host.append(root);
  update();
}

/* ===================== Widget: perm ===================== */
function mountPerm(host) {
  const id = uid();
  const { root, body } = card("ماتریس جایگشت مقدماتی P(r,s)", "ضرب از چپ سطرها را جابه‌جا می‌کند و ضرب از راست ستون‌ها را");
  let n = 3, r = 1, s = 3;
  const selN = h("select", { id: "w-pm" + id + "-n" }, [3, 4].map(k => h("option", { value: k }, toFa(k))));
  const selR = h("select", { id: "w-pm" + id + "-r" }), selS = h("select", { id: "w-pm" + id + "-s" });
  const figs = h("div", { class: "w-figs", dir: "ltr" });
  const facts = h("div", { class: "w-facts" });
  function fillSel(el, val) {
    el.textContent = "";
    for (let k = 1; k <= n; k++) el.append(h("option", { value: k, selected: k === val }, toFa(k)));
  }
  function update() {
    const sg = k => (k === r - 1 ? s - 1 : k === s - 1 ? r - 1 : k);
    const P = eye(n).map((row, i) => eye(n)[sg(i)]);
    const sw = k => k === r - 1 || k === s - 1;
    const a = (i, j) => ktx("a_{" + (i + 1) + (j + 1) + "}");
    const dummy = eye(n);
    const hl = (i, j, byRow) => (r !== s && sw(byRow ? i : j) ? "w-tr w-ch" : "");
    figs.textContent = "";
    figs.append(
      fig("P_{(" + r + "," + s + ")}", matEl(P, { cls: (i, j) => hl(i, j, true) })),
      fig("A", matEl(dummy, { cell: a })),
      fig("P\\,A", matEl(dummy, { cell: (i, j) => a(sg(i), j), cls: (i, j) => hl(i, j, true) })),
      fig("A\\,P", matEl(dummy, { cell: (i, j) => a(i, sg(j)), cls: (i, j) => hl(i, j, false) })));
    const sq = meq(mmul(P, P), eye(n));
    facts.innerHTML = "";
    facts.append(
      h("p", { html: ktx("P^2=I") + " " + (sq ? "برقرار است ✓" : "برقرار نیست") + " — دو بار جابه‌جا کردن همان دو سطر، همه‌چیز را به حالت اول برمی‌گرداند." }),
      h("p", { html: ktx(r === s ? "\\det P=1" : "\\det P=-1") + (r === s ? " — r = s یعنی همان I (هیچ جابه‌جایی‌ای انجام نمی‌شود)." : " — یک جابه‌جایی سطر، علامت دترمینان را عوض می‌کند.") }));
  }
  selN.addEventListener("change", () => { n = +selN.value; r = Math.min(r, n); s = Math.min(s, n); fillSel(selR, r); fillSel(selS, s); update(); });
  selR.addEventListener("change", () => { r = +selR.value; update(); });
  selS.addEventListener("change", () => { s = +selS.value; update(); });
  body.append(h("div", { class: "w-row" },
    h("label", { for: selN.id }, "اندازه n"), selN, h("label", { for: selR.id }, "r"), selR, h("label", { for: selS.id }, "s"), selS), figs, facts);
  host.textContent = ""; host.append(root);
  fillSel(selR, r); fillSel(selS, s);
  update();
}

/* ===================== Widget: lab ===================== */
const LAB_PRESETS = [
  ["مثال ۱۲.۳ (A)", [[2, 1, 0, 0], [1, -5, 1, -2], [-1, 0, 6, 3], [1, -3, -4, 10]]],
  ["مثال ۱۲.۳ (B)", [[2, 1, -1, 0], [3, -5, -1, 0], [0, 0, 2, -1], [1, -2, 0, 4]]],
  ["مثال ۱۳.۳", [[0, 1, -1, 1], [1, 1, -1, 2], [-1, -1, 1, 0], [1, 2, 0, 2]]],
  ["مثال ۱۱.۳", [[1, 1, -1, 1], [1, -1, -1, -1], [1, 1, 1, -1], [0, 1, 1, -1]]],
  ["مثال ۱.۳", [["1", "2"], ["1.0001", "2"]]],
];
function decimalOf(f) {
  const x = f.toNumber();
  return String(+x.toFixed(6));
}
function mountLab(host) {
  const id = uid();
  const { root, body } = card("آزمایشگاه ماتریس", "ماتریس را وارد کنید؛ ویژگی‌ها به‌صورت دقیق و زنده محاسبه می‌شوند");
  let n = 4, vals = LAB_PRESETS[0][1].map(r => r.map(String));
  const selN = h("select", { id: "w-lb" + id + "-n" }, [2, 3, 4].map(k => h("option", { value: k, selected: k === n }, toFa(k))));
  const grid = h("div", { class: "w-gridin", dir: "ltr" });
  const msg = h("div", { class: "w-msg", role: "alert" });
  const results = h("div", { class: "w-lcards" });
  let inputs = [];

  function readVals() { inputs.forEach((r, i) => r.forEach((el, j) => { vals[i] = vals[i] || []; vals[i][j] = el.value; })); }
  function build() {
    grid.textContent = ""; inputs = [];
    grid.style.gridTemplateColumns = "repeat(" + n + ", auto)";
    for (let i = 0; i < n; i++) {
      inputs[i] = [];
      for (let j = 0; j < n; j++) {
        const el = h("input", { id: "w-lb" + id + "-a" + i + j, class: "w-in", type: "text", inputmode: "text", dir: "ltr", autocomplete: "off",
          "aria-label": "a" + (i + 1) + (j + 1), value: vals[i] && vals[i][j] != null ? vals[i][j] : "0" });
        el.addEventListener("input", compute);
        inputs[i][j] = el; grid.append(el);
      }
    }
  }
  const chip = (ok, txt) => h("span", { class: "w-chip " + (ok ? "w-ok" : "w-no") }, txt);
  const lcard = (title, ...kids) => h("section", { class: "w-lcard" }, h("h4", { class: "w-lt" }, title), kids);
  const tex = (t, cls) => h("div", { class: "w-tex " + (cls || ""), dir: "ltr", html: ktx(t) });

  function icon(k) {
    const cs = 22 / n;
    const s = sv("svg", { viewBox: "0 0 22 22", width: 22, height: 22, class: "w-ico", "aria-hidden": "true" });
    for (let i = 0; i < n; i++) for (let j = 0; j < n; j++)
      s.append(sv("rect", { x: j * cs + 0.5, y: i * cs + 0.5, width: cs - 1, height: cs - 1, rx: 1, class: i < k && j < k ? "w-ico-on" : "w-ico-off" }));
    return s;
  }

  function compute() {
    readVals();
    let bad = 0;
    const A = inputs.map(r => r.map(el => { const v = parseFrac(el.value); el.classList.toggle("w-invalid", !v); if (!v) bad++; return v; }));
    results.textContent = "";
    if (bad) { msg.textContent = "ورودی نامعتبر است: عدد صحیح، کسر (مثل 1/3) یا اعشار (مثل 0.25) بنویسید."; return; }
    msg.textContent = "";

    const d = det(A), invertible = !d.isZero();
    results.append(lcard("دترمینان",
      tex("\\det(A)=" + d.toTeX("fraction") + (d.isInt() ? "" : "\\approx " + decimalOf(d))),
      h("div", null, chip(invertible, invertible ? "وارون‌پذیر" : "وارون‌ناپذیر"))));

    const minors = leadingMinors(A);
    const firstZero = minors.findIndex(v => v.isZero());
    const allNZ = firstZero < 0;
    const ml = h("ul", { class: "w-minors" }, minors.map((v, k) => h("li", { class: v.isZero() ? "w-mno" : "w-mok" }, icon(k + 1),
      h("span", { dir: "ltr", html: ktx("\\det(A_{" + (k + 1) + "})=" + v.toTeX("fraction")) }), chip(!v.isZero(), v.isZero() ? "صفر" : "ناصفر"))));
    results.append(lcard("زیرماتریس‌های اصلی پیشرو", ml, h("p", { class: "w-verdict " + (allNZ ? "w-vok" : "w-vbad") },
      allNZ ? "همهٔ زیرماتریس‌های اصلی پیشرو وارون‌پذیرند ⇒ طبق قضیهٔ ۵.۳ تجزیهٔ LU (بدون جابه‌جایی) وجود دارد"
        : "$A_" + (firstZero + 1) + "$ وارون‌ناپذیر است ⇒ تجزیهٔ LU وجود ندارد؛ با جابه‌جایی سطر $PA=LU$ (قضیهٔ ۷.۳)")));
    const vp = results.lastChild.querySelector(".w-verdict");
    if (!allNZ) vp.innerHTML = rich(vp.textContent);

    const dd = diagDominance(A);
    const mx = Math.max(1e-9, ...dd.rows.map(r => Math.max(r.diag.toNumber(), r.off.toNumber())));
    const dl = h("ul", { class: "w-dom" }, dd.rows.map((r, i) => h("li", null,
      h("span", { class: "w-dn", dir: "ltr", html: ktx("|a_{" + (i + 1) + (i + 1) + "}|=" + r.diag.toTeX("fraction") + ",\\ \\sum_{j\\neq " + (i + 1) + "}|a_{" + (i + 1) + "j}|=" + r.off.toTeX("fraction")) }),
      h("span", { class: "w-bars", dir: "ltr" },
        h("span", { class: "w-hb" }, h("i", { class: "w-hb-d", style: "width:" + (r.diag.toNumber() / mx) * 100 + "%" })),
        h("span", { class: "w-hb" }, h("i", { class: "w-hb-o", style: "width:" + (r.off.toNumber() / mx) * 100 + "%" }))),
      h("span", { class: "w-chip " + (r.cmp > 0 ? "w-ok" : r.cmp === 0 ? "w-eq0" : "w-no"), dir: "ltr" }, r.cmp > 0 ? ">" : r.cmp === 0 ? "=" : "<"))));
    const vtxt = dd.verdict === "strict" ? "غالب قطری اکید" : dd.verdict === "weak" ? "غالب قطری (نه اکید)" : "غالب قطری نیست";
    results.append(lcard("غلبهٔ قطری", dl,
      h("p", { class: "w-verdict " + (dd.verdict === "none" ? "w-vbad" : "w-vok") }, vtxt),
      dd.verdict === "strict" ? h("p", { class: "w-note" }, "طبق قضیهٔ ۶.۳: وارون‌پذیر است و حذف گاوسی بدون جابه‌جایی سطر انجام می‌شود") : null));

    if (invertible) {
      const k = kappaInf(A);
      results.append(lcard("وارون و عدد وضعیت",
        fig("A^{-1}", matEl(k.inv, { mode: "fraction" })),
        tex("\\|A\\|_\\infty=" + k.normA.toTeX("fraction") + (k.normA.isInt() ? "" : "=" + decimalOf(k.normA))),
        tex("\\|A^{-1}\\|_\\infty=" + k.normInv.toTeX("fraction") + (k.normInv.isInt() ? "" : "\\approx " + decimalOf(k.normInv))),
        tex("\\kappa_\\infty(A)=" + k.kappa.toTeX("fraction") + (k.kappa.isInt() ? "" : "\\approx " + decimalOf(k.kappa)))));
    } else {
      results.append(lcard("وارون و عدد وضعیت", h("p", { class: "w-verdict w-vbad" }, "A وارون‌ناپذیر است؛ وارون و عدد وضعیت تعریف نمی‌شوند.")));
    }

    const res = allNZ ? luSteps(A).result : pluSteps(A, { lenient: true }).result;
    if (res && res.ok) {
      const row = h("div", { class: "w-figs", dir: "ltr" });
      if (!allNZ) row.append(fig("P", matEl(res.P, { mode: "fraction" })));
      row.append(fig("L", matEl(res.L, { mode: "fraction" })), fig("U", matEl(res.U, { mode: "fraction" })));
      results.append(lcard(allNZ ? "تجزیهٔ LU" : "تجزیهٔ PA = LU", row,
        h("p", { class: "w-note", html: rich(allNZ ? "$A=LU$ (درستی بررسی شد)" : "$PA=LU$ " + (res.verified ? "(درستی با ضرب دقیق بررسی شد)" : "")) })));
    }
  }
  selN.addEventListener("change", () => { readVals(); n = +selN.value; build(); compute(); });
  const presetBtns = h("div", { class: "w-row" }, LAB_PRESETS.map(p => h("button", { class: "btn", type: "button", onclick: () => {
    vals = p[1].map(r => r.map(String)); n = vals.length; selN.value = String(n); build(); compute();
  } }, p[0])));
  body.append(h("div", { class: "w-row" }, h("label", { for: selN.id }, "اندازه n"), selN), presetBtns,
    h("div", { class: "w-scroll", dir: "ltr" }, grid), msg, results);
  host.textContent = ""; host.append(root);
  build(); compute();
}

/* ===================== Widget: flashcards ===================== */
const LEITNER_KEY = ((typeof document !== "undefined" && document.documentElement.getAttribute("data-page")) || "ch3") + "-leitner";
const BOX_DAYS = [0, 1, 3, 7, 16];
function mountFlashcards(host) {
  let cards = [];
  try { cards = JSON.parse(document.getElementById("cards-data").textContent) || []; } catch (e) { cards = []; }
  let prog = {};
  try { prog = JSON.parse(localStorage.getItem(LEITNER_KEY)) || {}; } catch (e) { prog = {}; }
  const save = () => { try { localStorage.setItem(LEITNER_KEY, JSON.stringify(prog)); } catch (e) { /* storage unavailable */ } };
  const today = () => { const d = new Date(); return Math.floor((d.getTime() - d.getTimezoneOffset() * 60000) / 86400000); };
  const stOf = c => { const p = prog[c.id]; return p && p.box ? p : { box: 1, due: 0 }; };
  const tags = ["همه"].concat(cards.map(c => c.tag).filter((t, i, a) => t && a.indexOf(t) === i));
  let tag = "همه", mode = "today", flipped = false, last = null, confirming = false;

  const { root, body } = card("کارت‌های مرور", "روش لایتنر: کارت‌هایی که بلد نیستید زودتر برمی‌گردند");
  root.setAttribute("tabindex", "0");
  host.textContent = ""; host.append(root);

  function queue() {
    const t = today();
    let list = cards.map((c, i) => ({ c, i, s: stOf(c) })).filter(x => tag === "همه" || x.c.tag === tag);
    if (mode === "today") list = list.filter(x => x.s.due <= t);
    list.sort((a, b) => a.s.due - b.s.due || a.s.box - b.s.box || a.i - b.i);
    if (list.length > 1 && last && list[0].c.id === last) list.push(list.shift());
    return list;
  }
  function math(el) {
    if (typeof renderMathInElement === "function") {
      try { renderMathInElement(el, { delimiters: [{ left: "$$", right: "$$", display: true }, { left: "$", right: "$", display: false }], throwOnError: false }); } catch (e) { /* ignore */ }
    }
  }
  function answer(known) {
    const q = queue()[0]; if (!q) return;
    const box = known ? Math.min(5, q.s.box + 1) : 1;
    prog[q.c.id] = { box, due: today() + BOX_DAYS[box - 1] };
    save(); last = q.c.id; flipped = false; render();
  }
  function flip() { flipped = !flipped; const el = body.querySelector(".w-fc"); if (!el) return; setFlip(el); }
  function setFlip(el) {
    el.classList.toggle("w-flipped", flipped);
    el.setAttribute("aria-pressed", String(flipped));
    const f = el.querySelector(".w-front"), b = el.querySelector(".w-back");
    f.setAttribute("aria-hidden", String(flipped)); b.setAttribute("aria-hidden", String(!flipped));
    const act = body.querySelector(".w-answers"), hint = body.querySelector(".w-hint");
    if (act) act.hidden = !flipped; if (hint) hint.hidden = flipped;
  }

  function render() {
    body.textContent = "";
    const t = today();
    const dueAll = cards.filter(c => stOf(c).due <= t).length;
    const counts = [1, 2, 3, 4, 5].map(b => cards.filter(c => stOf(c).box === b).length);
    const boxes = h("div", { class: "w-boxes", role: "img", "aria-label": "تعداد کارت‌ها در هر جعبه" }, counts.map((k, i) =>
      h("div", { class: "w-boxseg" }, h("i", { style: "width:" + (cards.length ? (k / cards.length) * 100 : 0) + "%" }), h("span", null, "جعبه " + toFa(i + 1) + ": " + toFa(k)))));
    const chips = h("div", { class: "w-chips", role: "group", "aria-label": "فیلتر برچسب" }, tags.map(tg => h("button", { class: "btn w-tagbtn", type: "button", "aria-pressed": String(tg === tag), onclick: () => { tag = tg; flipped = false; render(); } }, tg)));
    const modeBtn = h("button", { class: "btn", type: "button", "aria-pressed": String(mode === "today"), onclick: () => { mode = mode === "today" ? "all" : "today"; flipped = false; render(); } },
      mode === "today" ? "فقط کارت‌های امروز" : "همهٔ کارت‌ها");
    let resetUi;
    if (confirming) resetUi = h("span", { class: "w-confirm", role: "alert" }, "همهٔ پیشرفت پاک می‌شود. مطمئنید؟ ",
      h("button", { class: "btn", type: "button", onclick: () => { prog = {}; save(); confirming = false; last = null; flipped = false; render(); } }, "بله، پاک کن"),
      h("button", { class: "btn", type: "button", onclick: () => { confirming = false; render(); } }, "انصراف"));
    else resetUi = h("button", { class: "btn", type: "button", onclick: () => { confirming = true; render(); } }, "شروع دوباره");
    body.append(h("div", { class: "w-fc-top" }, h("div", { class: "w-fc-today" }, "امروز: " + toFa(dueAll) + " کارت"), boxes),
      chips, h("div", { class: "w-row" }, modeBtn, resetUi));

    const q = queue();
    if (!q.length) {
      const pool = cards.filter(c => tag === "همه" || c.tag === tag);
      let txt;
      if (!pool.length) txt = "در این دسته کارتی وجود ندارد.";
      else {
        const next = Math.min(...pool.map(c => stOf(c).due)) - t;
        txt = "همهٔ کارت‌های امروز را مرور کردید. کارت بعدی " + (next <= 1 ? "فردا" : toFa(next) + " روز دیگر") + " آمادهٔ مرور است.";
      }
      body.append(h("div", { class: "w-done" }, h("p", null, txt),
        mode === "today" && pool.length ? h("button", { class: "btn", type: "button", onclick: () => { mode = "all"; render(); } }, "نمایش همهٔ کارت‌ها") : null));
      return;
    }
    const c = q[0].c;
    const face = (cls, html) => { const f = h("div", { class: "w-face " + cls }, h("div", { class: "w-face-tag" }, c.tag || ""), h("div", { class: "w-face-c", html })); math(f); return f; };
    const fc = h("div", { class: "w-fc", role: "button", tabindex: "0", "aria-label": "کارت — برای دیدن پشت یا روی کارت بزنید", onclick: flip },
      h("div", { class: "w-fc-in" }, face("w-front", c.front), face("w-back", c.back)));
    const answers = h("div", { class: "w-answers w-row", hidden: true },
      h("button", { class: "btn primary", type: "button", onclick: () => answer(true) }, "بلد بودم"),
      h("button", { class: "btn", type: "button", onclick: () => answer(false) }, "بلد نبودم"));
    body.append(fc, h("p", { class: "w-hint" }, "کارت را بزنید یا Space را فشار دهید تا پاسخ را ببینید."), answers,
      h("p", { class: "w-note" }, "جعبهٔ این کارت: " + toFa(q[0].s.box)));
    setFlip(fc);
  }
  root.addEventListener("keydown", e => {
    if (e.key === " " && !/^(BUTTON|INPUT|SELECT|TEXTAREA)$/.test(e.target.tagName)) { e.preventDefault(); flip(); }
    else if (e.key === "Enter" && e.target.classList && e.target.classList.contains("w-fc")) { e.preventDefault(); flip(); }
  });
  render();
}

/* ===================== Mounting & exports ===================== */
function mountStepperPreset(el) {
  const p = el.getAttribute("data-preset");
  if (p === "custom") mountCustom(el);
  else if (PRESETS[p]) mountStepper(el, PRESETS[p]);
}
const REGISTRY = { stepper: mountStepperPreset, illcond: mountIllcond, cost: mountCost, perm: mountPerm, lab: mountLab, flashcards: mountFlashcards };
function mountAll() {
  document.querySelectorAll("[data-widget]").forEach(el => {
    const f = REGISTRY[el.getAttribute("data-widget")];
    if (!f) return;
    try { f(el); } catch (e) { if (typeof console !== "undefined") console.error(e); }
  });
}

globalThis.LAW = { Frac, parseFrac, gaussSteps, gaussJordanSteps, gjInverseSteps, luSteps, pluSteps, det, inverse,
  mat, leadingMinors, diagDominance, normInf, kappaInf, buildSteps, PRESETS, meq, mmul, eye, toFa };
if (typeof module !== "undefined" && module.exports) module.exports = globalThis.LAW;

if (typeof document !== "undefined") {
  if (document.readyState === "loading") document.addEventListener("DOMContentLoaded", mountAll);
  else mountAll();
}
})();
