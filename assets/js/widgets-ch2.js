(function () {
"use strict";

/* ===================== Pure math ===================== */
const EPS = 1e-12;
const le = (a, b) => a <= b + EPS * Math.max(1, Math.abs(b));
const SQ2 = Math.SQRT2;

function pnorm(x, p) {
  const a = x.map(Math.abs), m = Math.max(...a);
  if (p === Infinity || m === 0) return m;
  if (p === 1) return a.reduce((s, v) => s + v, 0);
  if (p === 2) return Math.hypot(...a);
  return m * Math.pow(a.reduce((s, v) => s + Math.pow(v / m, p), 0), 1 / p);
}
function pCurve(p, n) {
  const out = [], e = 2 / p;
  for (let i = 0; i < n; i++) {
    const t = (2 * Math.PI * i) / n;
    let c = Math.cos(t), s = Math.sin(t);
    if (Math.abs(c) < 1e-15) c = 0;
    if (Math.abs(s) < 1e-15) s = 0;
    if (p === Infinity) { const m = Math.max(Math.abs(c), Math.abs(s)); out.push([c / m, s / m]); }
    else out.push([Math.sign(c) * Math.pow(Math.abs(c), e), Math.sign(s) * Math.pow(Math.abs(s), e)]);
  }
  return out;
}
function exercise5(x) {
  const a = pnorm(x, 1), b = pnorm(x, 2), c = pnorm(x, Infinity);
  return [
    { ok: le(c, b) && le(b, a), vals: [c, b, a] },
    { ok: le(a, SQ2 * b), vals: [a, SQ2 * b] },
    { ok: le(b, SQ2 * c), vals: [b, SQ2 * c] },
    { ok: le(a, 2 * c), vals: [a, 2 * c] },
  ];
}

function normDigits(s) {
  return String(s)
    .replace(/[۰-۹]/g, c => String(c.charCodeAt(0) - 0x06F0))
    .replace(/[٠-٩]/g, c => String(c.charCodeAt(0) - 0x0660))
    .replace(/[٫,]/g, ".")
    .replace(/[−–—]/g, "-")
    .replace(/[÷⁄]/g, "/")
    .replace(/\s+/g, "");
}
function parseDec(s) {
  if (!/^[+-]?(\d+\.?\d*|\.\d+)$/.test(s)) return null;
  const v = Number(s);
  return Number.isFinite(v) ? v : null;
}
function parseNum(str) {
  if (str == null) return null;
  const s = normDigits(str);
  if (!s) return null;
  if (s.includes("/")) {
    const parts = s.split("/");
    if (parts.length !== 2) return null;
    const a = parseDec(parts[0]), b = parseDec(parts[1]);
    if (a == null || b == null || b === 0) return null;
    return a / b;
  }
  return parseDec(s);
}

const mmul = (A, B) => [
  [A[0][0] * B[0][0] + A[0][1] * B[1][0], A[0][0] * B[0][1] + A[0][1] * B[1][1]],
  [A[1][0] * B[0][0] + A[1][1] * B[1][0], A[1][0] * B[0][1] + A[1][1] * B[1][1]],
];
const I2 = () => [[1, 0], [0, 1]];
function mpow(A, k) { let R = I2(); for (let i = 0; i < k; i++) R = mmul(R, A); return R; }
function powers(A, K) {
  const out = [I2()];
  for (let k = 1; k <= K; k++) out.push(mmul(out[k - 1], A));
  return out;
}
const norm1 = A => Math.max(Math.abs(A[0][0]) + Math.abs(A[1][0]), Math.abs(A[0][1]) + Math.abs(A[1][1]));
const normInf = A => Math.max(Math.abs(A[0][0]) + Math.abs(A[0][1]), Math.abs(A[1][0]) + Math.abs(A[1][1]));
const normFro = A => Math.hypot(A[0][0], A[0][1], A[1][0], A[1][1]);
function gram(A) {
  return {
    p: A[0][0] * A[0][0] + A[1][0] * A[1][0],
    q: A[0][0] * A[0][1] + A[1][0] * A[1][1],
    r: A[0][1] * A[0][1] + A[1][1] * A[1][1],
  };
}
function norm2(A) {
  const { p, q, r } = gram(A);
  return Math.sqrt(Math.max(0, (p + r) / 2 + Math.sqrt(((p - r) / 2) ** 2 + q * q)));
}
function topSingular(A) {
  const { p, q, r } = gram(A);
  const th = 0.5 * Math.atan2(2 * q, p - r);
  const v = [Math.cos(th), Math.sin(th)];
  return { sigma: norm2(A), v, Av: [A[0][0] * v[0] + A[0][1] * v[1], A[1][0] * v[0] + A[1][1] * v[1]] };
}
function rho(A) {
  const tr = A[0][0] + A[1][1], det = A[0][0] * A[1][1] - A[0][1] * A[1][0], disc = tr * tr - 4 * det;
  if (disc < 0) return Math.sqrt(det);
  const s = Math.sqrt(disc);
  return Math.max(Math.abs((tr + s) / 2), Math.abs((tr - s) / 2));
}
function matNorms(A) {
  const n1 = norm1(A), ni = normInf(A);
  return { n1, ni, nF: normFro(A), n2: norm2(A), rho: rho(A), geo: Math.sqrt(n1 * ni) };
}
function exercise19(A) {
  const N = matNorms(A);
  const two = (lo, mid, hi) => ({ ok: le(lo, mid) && le(mid, hi), vals: [lo, mid, hi] });
  return [
    two(N.n1 / SQ2, N.n2, SQ2 * N.n1),
    two(N.ni / SQ2, N.n2, SQ2 * N.ni),
    two(N.n2, N.nF, SQ2 * N.n2),
    { ok: le(N.n2, N.geo), vals: [N.n2, N.geo] },
    { ok: le(N.rho, N.n1), vals: [N.rho, N.n1] },
    { ok: le(N.rho, N.ni), vals: [N.rho, N.ni] },
    { ok: le(N.rho, N.n2), vals: [N.rho, N.n2] },
  ];
}
function csData(x, y) {
  const nx = Math.hypot(x[0], x[1]), ny = Math.hypot(y[0], y[1]);
  let ip = x[0] * y[0] + x[1] * y[1];
  if (Math.abs(ip) < 1e-12) ip = 0;
  const prod = nx * ny;
  const ratio = prod > 0 ? Math.min(1, Math.abs(ip) / prod) : 0;
  return { ip, nx, ny, prod, ratio, alpha: ny > 0 ? ip / (ny * ny) : null, equal: prod > 0 && ratio >= 0.9995 };
}
function snapParallel(x, y) {
  const ny = Math.hypot(y[0], y[1]);
  if (ny === 0) return x.slice();
  let len = Math.hypot(x[0], x[1]);
  if (len === 0) len = 1.5;
  const sg = x[0] * y[0] + x[1] * y[1] < 0 ? -1 : 1;
  return fitBox([sg * len * y[0] / ny, sg * len * y[1] / ny], 3);
}
function snapPerp(x, y) {
  const ny = Math.hypot(y[0], y[1]);
  if (ny === 0) return x.slice();
  let len = Math.hypot(x[0], x[1]);
  if (len === 0) len = 1.5;
  const u = [-y[1] / ny, y[0] / ny];
  const sg = x[0] * u[0] + x[1] * u[1] < 0 ? -1 : 1;
  return fitBox([sg * len * u[0], sg * len * u[1]], 3);
}
function fitBox(v, lim) {
  const m = Math.max(Math.abs(v[0]), Math.abs(v[1]));
  return m > lim ? [v[0] * lim / m, v[1] * lim / m] : v;
}
function ellipsePts(A, n) {
  const out = [];
  for (let i = 0; i < n; i++) {
    const t = (2 * Math.PI * i) / n, c = Math.cos(t), s = Math.sin(t);
    out.push([A[0][0] * c + A[0][1] * s, A[1][0] * c + A[1][1] * s]);
  }
  return out;
}
function squareImage(A) {
  return [[1, 1], [-1, 1], [-1, -1], [1, -1]].map(v => [A[0][0] * v[0] + A[0][1] * v[1], A[1][0] * v[0] + A[1][1] * v[1]]);
}
function niceStep(E) {
  const base = Math.pow(10, Math.floor(Math.log10(E / 4)));
  for (const c of [1, 2, 5, 10]) if (c * base >= E / 4) return c * base;
  return 10 * base;
}

const LAW2 = { pnorm, pCurve, exercise5, parseNum, mmul, mpow, powers, norm1, normInf, normFro, norm2, topSingular, rho,
  matNorms, exercise19, csData, snapParallel, snapPerp, ellipsePts, squareImage, niceStep, le };

/* ===================== DOM helpers (browser only) ===================== */
const toFa = s => String(s).replace(/[0-9]/g, d => "۰۱۲۳۴۵۶۷۸۹"[d]);
const esc = s => String(s).replace(/[&<>"]/g, c => ({ "&": "&amp;", "<": "&lt;", ">": "&gt;", '"': "&quot;" }[c]));
function fmtN(x) {
  if (x === 0) return "0";
  const a = Math.abs(x);
  if (a >= 1e5 || a < 1e-3) return x.toExponential(3).replace(/\.?0+e/, "e").replace("e+", "e");
  return String(parseFloat(x.toPrecision(4)));
}
const fa = (x, d) => toFa(d == null ? fmtN(x) : x.toFixed(d)).replace(/\./g, "٫").replace(/-/g, "−");
function latin(x, d) {
  const s = x.toFixed(d == null ? 3 : d);
  return /^-0\.?0*$/.test(s) ? s.slice(1) : s;
}
function sig4tex(x) {
  if (x === 0) return "0";
  const a = Math.abs(x);
  if (a >= 1e6 || a < 1e-4) {
    const parts = x.toExponential(3).split("e");
    return parts[0].replace(/\.?0+$/, "") + "\\times10^{" + (+parts[1]) + "}";
  }
  return String(parseFloat(x.toPrecision(4)));
}
const SUP = { "0": "⁰", "1": "¹", "2": "²", "3": "³", "4": "⁴", "5": "⁵", "6": "⁶", "7": "⁷", "8": "⁸", "9": "⁹", "-": "⁻" };
const sup = n => String(n).replace(/[0-9-]/g, c => SUP[c]);

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
function card(title, sub) {
  const body = h("div", { class: "w-body" });
  const root = h("div", { class: "w-card w2" },
    h("div", { class: "w-head" }, h("div", { class: "w-title" }, title), sub ? h("div", { class: "w-sub" }, sub) : null), body);
  return { root, body };
}
function mkSlider(id, min, max, step, value, label, onInput) {
  const out = h("output", { class: "w-num", for: id, dir: "ltr" });
  const input = h("input", { id, type: "range", min, max, step, value, class: "w-range" });
  input.addEventListener("input", () => onInput(+input.value));
  return { input, out, row: h("div", { class: "w-srow" }, h("label", { for: id }, label), input, out) };
}
const chip = ok => h("span", { class: "w-chip " + (ok ? "w-ok" : "w-no"), role: "img", "aria-label": ok ? "برقرار" : "نقض شد" }, ok ? "✓" : "✗");
const stat = (tex, id) => {
  const val = h("div", { class: "w-num", dir: "ltr", id });
  return { el: h("div", { class: "w-stat" }, h("div", { class: "w-statl", dir: "ltr", html: ktx(tex) }), val), val };
};

function dragHandle(svg, hit, VB, move, step) {
  const toVB = e => {
    const rc = svg.getBoundingClientRect();
    return [((e.clientX - rc.left) / rc.width) * VB, ((e.clientY - rc.top) / rc.height) * VB];
  };
  let active = false;
  hit.addEventListener("pointerdown", e => {
    active = true;
    if (hit.setPointerCapture) { try { hit.setPointerCapture(e.pointerId); } catch (err) { /* ignore */ } }
    if (hit.focus) hit.focus({ preventScroll: true });
    e.preventDefault();
    const v = toVB(e); move(v[0], v[1]);
  });
  hit.addEventListener("pointermove", e => { if (active) { const v = toVB(e); move(v[0], v[1]); } });
  const end = e => {
    active = false;
    if (hit.releasePointerCapture) { try { hit.releasePointerCapture(e.pointerId); } catch (err) { /* ignore */ } }
  };
  hit.addEventListener("pointerup", end);
  hit.addEventListener("pointercancel", end);
  hit.addEventListener("keydown", e => {
    const d = { ArrowLeft: [-1, 0], ArrowRight: [1, 0], ArrowUp: [0, 1], ArrowDown: [0, -1] }[e.key];
    if (!d) return;
    e.preventDefault();
    step(d[0], d[1], e.shiftKey);
  });
}

function arrow(cls) {
  const line = sv("line", { class: "w-ln w-" + cls, "stroke-width": 3, "stroke-linecap": "round" });
  const head = sv("polygon", { class: "w-f" + cls });
  const g = sv("g", null, line, head);
  function set(x1, y1, x2, y2) {
    const dx = x2 - x1, dy = y2 - y1, L = Math.hypot(dx, dy);
    if (L < 1) { g.style.display = "none"; return; }
    g.style.display = "";
    const ux = dx / L, uy = dy / L, hl = Math.min(14, L * 0.6), hw = 5.5;
    const bx = x2 - ux * hl, by = y2 - uy * hl;
    line.setAttribute("x1", x1); line.setAttribute("y1", y1);
    line.setAttribute("x2", bx + ux * 1); line.setAttribute("y2", by + uy * 1);
    head.setAttribute("points", x2 + "," + y2 + " " + (bx - uy * hw) + "," + (by + ux * hw) + " " + (bx + uy * hw) + "," + (by - ux * hw));
  }
  return { g, set };
}

const poly = pts => pts.map(p => p[0].toFixed(2) + "," + p[1].toFixed(2)).join(" ");
function gridPlane(svg, VB, C, S, lim, step, labelAt) {
  const n = Math.floor(lim / step + 1e-9);
  for (let i = -n; i <= n; i++) {
    const v = i * step;
    if (i === 0) continue;
    svg.append(sv("line", { x1: C + S * v, y1: C - S * lim, x2: C + S * v, y2: C + S * lim, class: "w-gridl" }),
      sv("line", { x1: C - S * lim, y1: C + S * v, x2: C + S * lim, y2: C + S * v, class: "w-gridl" }));
  }
  svg.append(sv("line", { x1: C - S * lim, y1: C, x2: C + S * lim, y2: C, class: "w-axisl" }),
    sv("line", { x1: C, y1: C - S * lim, x2: C, y2: C + S * lim, class: "w-axisl" }));
  labelAt.forEach(v => {
    const t = String(parseFloat(v.toPrecision(6)));
    svg.append(sv("text", { x: C + S * v, y: C + 16, "text-anchor": "middle", class: "w-tick" }, t),
      sv("text", { x: C - 6, y: C - S * v + 4, "text-anchor": "end", class: "w-tick" }, t));
  });
  return VB;
}

/* ===================== Widget: unitball ===================== */
function mountUnitBall(host) {
  const id = uid();
  const { root, body } = card("گوی واحد در نرم‌های مختلف", "نقطه را بکش؛ نرم‌ها و نابرابری‌های تمرین ۵ زنده حساب می‌شوند");
  const VB = 480, LIM = 1.6, S = (VB - 72) / (2 * LIM), C = VB / 2;
  const X = v => C + S * v, Y = v => C - S * v;
  const st = { p: 4, inf: false, P: [0.6, 0.8] };

  const svg = sv("svg", { viewBox: "0 0 " + VB + " " + VB, class: "w-svg", role: "img", "aria-label": "گوی‌های واحد نرم‌های ۱، ۲ و بی‌نهایت در صفحه و نقطهٔ قابل کشیدن P" });
  gridPlane(svg, VB, C, S, 1.5, 0.5, [-1, 1]);
  svg.append(sv("text", { x: X(LIM) - 2, y: C - 8, "text-anchor": "end", class: "w-tick" }, "x₁"),
    sv("text", { x: C + 8, y: Y(LIM) + 10, class: "w-tick" }, "x₂"));
  svg.append(
    sv("polygon", { points: poly([[X(1), Y(0)], [X(0), Y(1)], [X(-1), Y(0)], [X(0), Y(-1)]]), class: "w-ln w-s2" }),
    sv("circle", { cx: C, cy: C, r: S, class: "w-ln w-s1" }),
    sv("rect", { x: X(-1), y: Y(1), width: 2 * S, height: 2 * S, class: "w-ln w-s3" }));
  const pcurve = sv("polygon", { class: "w2-pcurve", "stroke-dasharray": "7 5" });
  svg.append(pcurve);
  const lab = (t, x, y) => svg.append(sv("text", { x: X(x), y: Y(y) + 4, "text-anchor": "middle", class: "w-pl" }, t));
  lab("‖·‖₁", -0.5, -0.5); lab("‖·‖₂", -0.7071, -0.7071); lab("‖·‖∞", -1, -1);
  const plab = sv("text", { "text-anchor": "middle", class: "w-pl" }, "‖·‖p");
  svg.append(plab);

  const dot = sv("circle", { r: 7, class: "w-dot-fill" });
  const pl = sv("text", { class: "w-pl", "font-weight": 800 }, "P");
  const hit = sv("circle", { r: 22, class: "w2-hit", tabindex: "0", "aria-roledescription": "نقطهٔ قابل کشیدن" });
  svg.append(dot, pl, hit);
  const plot = h("div", { class: "w-plot", dir: "ltr" }, svg);

  const legendP = h("span", { class: "w-lg" });
  const legend = h("div", { class: "w-legend", dir: "ltr" },
    h("span", { class: "w-lg" }, h("i", { class: "w-sw w-b2" }), h("span", { html: ktx("\\|x\\|_1=1") })),
    h("span", { class: "w-lg" }, h("i", { class: "w-sw w-b1" }), h("span", { html: ktx("\\|x\\|_2=1") })),
    h("span", { class: "w-lg" }, h("i", { class: "w-sw w-b3" }), h("span", { html: ktx("\\|x\\|_\\infty=1") })),
    legendP);
  legendP.append(h("i", { class: "w-sw w2-swd w2-swd-p" }), h("span"));

  const sl = mkSlider("w2-ub" + id + "-p", 1, 30, 0.1, st.p, "پارامتر p", v => { st.p = v; update(); });
  const tog = h("button", { class: "btn", type: "button", "aria-pressed": "false", id: "w2-ub" + id + "-inf", onclick: () => {
    st.inf = !st.inf; update();
  } }, "⁦p = ∞⁩");
  const ctl = h("div", { class: "w-row" }, tog);

  const sN1 = stat("\\|P\\|_1"), sN2 = stat("\\|P\\|_2"), sNi = stat("\\|P\\|_\\infty"), sNp = stat("\\|P\\|_p");
  const cards = h("div", { class: "w-cards" }, sN1.el, sN2.el, sNi.el, sNp.el);
  const coord = h("p", { class: "w-note", dir: "ltr" });
  const checks = h("ul", { class: "w2-checks", dir: "ltr" });
  const sentence = h("p", { class: "w-sentence" }, "هرچه p بزرگ‌تر شود، منحنی به مربع نرم بی‌نهایت نزدیک‌تر می‌شود (تمرین ۴).");
  body.append(legend, plot, sentence, sl.row, ctl, cards, coord, checks);

  function update() {
    const P = st.P, pv = st.inf ? Infinity : st.p;
    const pts = pCurve(pv, 480).map(q => [X(q[0]), Y(q[1])]);
    pcurve.setAttribute("points", poly(pts));
    const qq = pCurve(pv, 8)[7];
    plab.setAttribute("x", X(qq[0]) + 6); plab.setAttribute("y", Y(qq[1]) + 4);
    const pTxt = st.inf ? "\\infty" : String(+st.p.toFixed(1));
    legendP.lastChild.innerHTML = ktx("\\|x\\|_{p}=1,\\ p=" + pTxt);
    sNp.el.firstChild.innerHTML = ktx("\\|P\\|_{p},\\ p=" + pTxt);
    sl.out.textContent = "p = " + (st.inf ? "∞" : fa(st.p, 1));
    sl.input.disabled = st.inf;
    tog.setAttribute("aria-pressed", String(st.inf));
    dot.setAttribute("cx", X(P[0])); dot.setAttribute("cy", Y(P[1]));
    hit.setAttribute("cx", X(P[0])); hit.setAttribute("cy", Y(P[1]));
    pl.setAttribute("x", X(P[0]) + 11); pl.setAttribute("y", Y(P[1]) - 9);
    hit.setAttribute("aria-label", "نقطهٔ P با مختصات " + fa(P[0], 2) + " و " + fa(P[1], 2) + "؛ با کلیدهای جهت‌دار جابه‌جا کنید");
    sN1.val.textContent = fa(pnorm(P, 1), 3);
    sN2.val.textContent = fa(pnorm(P, 2), 3);
    sNi.val.textContent = fa(pnorm(P, Infinity), 3);
    sNp.val.textContent = fa(pnorm(P, pv), 3);
    coord.innerHTML = ktx("P=(" + latin(P[0], 2) + ",\\ " + latin(P[1], 2) + ")");
    const f = n => latin(n, 3);
    const ex = exercise5(P), v = ex.map(e => e.vals);
    const tex = [
      "\\|x\\|_\\infty=" + f(v[0][0]) + "\\ \\le\\ \\|x\\|_2=" + f(v[0][1]) + "\\ \\le\\ \\|x\\|_1=" + f(v[0][2]),
      "\\|x\\|_1=" + f(v[1][0]) + "\\ \\le\\ \\sqrt2\\,\\|x\\|_2=" + f(v[1][1]),
      "\\|x\\|_2=" + f(v[2][0]) + "\\ \\le\\ \\sqrt2\\,\\|x\\|_\\infty=" + f(v[2][1]),
      "\\|x\\|_1=" + f(v[3][0]) + "\\ \\le\\ 2\\,\\|x\\|_\\infty=" + f(v[3][1]),
    ];
    checks.textContent = "";
    tex.forEach((t, i) => checks.append(h("li", { class: "w2-chk" }, h("span", { class: "w-tex", html: ktx(t) }), chip(ex[i].ok))));
  }

  dragHandle(svg, hit, VB, (vx, vy) => {
    const m = LIM - 0.05;
    st.P = [Math.max(-m, Math.min(m, (vx - C) / S)), Math.max(-m, Math.min(m, (C - vy) / S))];
    update();
  }, (dx, dy, shift) => {
    const m = LIM - 0.05, s = shift ? 0.01 : 0.05;
    st.P = [Math.max(-m, Math.min(m, st.P[0] + dx * s)), Math.max(-m, Math.min(m, st.P[1] + dy * s))];
    update();
  });

  host.textContent = ""; host.append(root);
  update();
}

/* ===================== Widget: cs ===================== */
function mountCS(host) {
  const { root, body } = card("نابرابری کوشی-شوارتز", "⁦|(x,y)|⁩ هیچ‌وقت از ⁦‖x‖₂‖y‖₂⁩ بیشتر نمی‌شود");
  const VB = 480, LIM = 3, S = (VB - 60) / (2 * LIM), C = VB / 2;
  const X = v => C + S * v, Y = v => C - S * v;
  const st = { x: [2, 1], y: [0.5, 2] };

  const svg = sv("svg", { viewBox: "0 0 " + VB + " " + VB, class: "w-svg", role: "img", "aria-label": "صفحه با دو بردار قابل کشیدن x و y" });
  gridPlane(svg, VB, C, S, LIM, 1, [-3, -2, -1, 1, 2, 3]);
  const ax = arrow("s1"), ay = arrow("s2");
  const lx = sv("text", { class: "w-pl w2-vl" }, "x"), ly = sv("text", { class: "w-pl w2-vl" }, "y");
  svg.append(ax.g, ay.g, lx, ly);
  const mkHandle = (cls, name) => {
    const dot = sv("circle", { r: 7, class: "w-mk w-f" + cls });
    const hit = sv("circle", { r: 22, class: "w2-hit", tabindex: "0", "aria-roledescription": "نوک بردار قابل کشیدن" });
    svg.append(dot, hit);
    return { dot, hit, name };
  };
  const hx = mkHandle("s1", "x"), hy = mkHandle("s2", "y");
  const plot = h("div", { class: "w-plot", dir: "ltr" }, svg);

  const legend = h("div", { class: "w-legend", dir: "ltr" },
    h("span", { class: "w-lg" }, h("i", { class: "w-sw w-b1" }), h("span", { html: ktx("x") })),
    h("span", { class: "w-lg" }, h("i", { class: "w-sw w-b2" }), h("span", { html: ktx("y") })));

  const sIp = stat("(x,y)=x_1y_1+x_2y_2"), sPr = stat("\\|x\\|_2\\,\\|y\\|_2");
  const meterFill = h("i");
  const meter = h("div", { class: "w2-meter", role: "meter", "aria-valuemin": "0", "aria-valuemax": "1", "aria-label": "نسبت |(x,y)| به ‖x‖₂‖y‖₂" }, meterFill);
  const ratioVal = h("span", { class: "w-num", dir: "ltr" });
  const ratioBox = h("div", { class: "w-stat w2-ratio" },
    h("div", { class: "w-statl", dir: "ltr", html: ktx("\\dfrac{|(x,y)|}{\\|x\\|_2\\|y\\|_2}=|\\cos\\theta|") }),
    h("div", { class: "w2-ratio-row", dir: "ltr" }, meter, ratioVal));
  const cards = h("div", { class: "w-cards" }, sIp.el, sPr.el);
  const eqMsg = h("div", { class: "w2-eq", "aria-live": "polite" });
  const btnPar = h("button", { class: "btn primary", type: "button", onclick: () => { st.x = snapParallel(st.x, st.y); update(); } }, "x را موازی y کن");
  const btnPer = h("button", { class: "btn", type: "button", onclick: () => { st.x = snapPerp(st.x, st.y); update(); } }, "عمود کن");
  const note = h("p", { class: "w-note" }, "عمودکردن یعنی ", h("span", { html: ktx("(x,y)=0") }), " (تعریف ۳.۲)؛ در این حالت نسبت صفر است.");
  body.append(legend, plot, h("div", { class: "w-row" }, btnPar, btnPer), cards, ratioBox, eqMsg, note);

  function place(v, a, lbl, hd) {
    const cx = X(v[0]), cy = Y(v[1]);
    a.set(C, C, cx, cy);
    hd.dot.setAttribute("cx", cx); hd.dot.setAttribute("cy", cy);
    hd.hit.setAttribute("cx", cx); hd.hit.setAttribute("cy", cy);
    hd.hit.setAttribute("aria-label", "نوک بردار " + hd.name + " در " + fa(v[0], 2) + " و " + fa(v[1], 2) + "؛ با کلیدهای جهت‌دار جابه‌جا کنید");
    const L = Math.hypot(v[0], v[1]) || 1;
    lbl.setAttribute("x", Math.max(10, Math.min(VB - 10, cx + (v[0] / L) * 18)));
    lbl.setAttribute("y", Math.max(14, Math.min(VB - 6, cy - (v[1] / L) * 18 + 5)));
    lbl.setAttribute("text-anchor", "middle");
  }
  function update() {
    place(st.x, ax, lx, hx); place(st.y, ay, ly, hy);
    const d = csData(st.x, st.y);
    sIp.val.textContent = fa(d.ip, 3);
    sPr.val.textContent = fa(d.prod, 3);
    ratioVal.textContent = fa(d.ratio, 4);
    meterFill.style.width = (d.ratio * 100).toFixed(2) + "%";
    meter.setAttribute("aria-valuenow", d.ratio.toFixed(4));
    eqMsg.textContent = "";
    if (d.equal) eqMsg.append(h("span", { class: "w-chip w-ok" }, "برابری برقرار است: x مضربی از y است (x = αy)"),
      h("span", { class: "w-num w2-alpha", dir: "ltr" }, " α = " + fa(d.alpha, 3)));
  }
  const clamp = v => Math.max(-LIM, Math.min(LIM, v));
  const wire = (hd, key) => dragHandle(svg, hd.hit, VB, (vx, vy) => {
    st[key] = [clamp((vx - C) / S), clamp((C - vy) / S)]; update();
  }, (dx, dy, shift) => {
    const s = shift ? 0.05 : 0.25;
    st[key] = [clamp(st[key][0] + dx * s), clamp(st[key][1] + dy * s)]; update();
  });
  wire(hx, "x"); wire(hy, "y");

  host.textContent = ""; host.append(root);
  update();
}

/* ===================== Matrix input ===================== */
function fmtEntry(v) { return String(+v.toPrecision(12)); }
function mkMatInput(pfx, onChange) {
  const inputs = [[], []];
  const msg = h("div", { class: "w-msg", role: "alert" });
  const grid = h("div", { class: "w-gridin w2-mbr", dir: "ltr" });
  for (let i = 0; i < 2; i++) for (let j = 0; j < 2; j++) {
    const inp = h("input", { id: pfx + "-a" + (i + 1) + (j + 1), class: "w-in", type: "text", inputmode: "text", autocomplete: "off", spellcheck: "false",
      "aria-label": "درایهٔ سطر " + toFa(i + 1) + " ستون " + toFa(j + 1), dir: "ltr" });
    inp.addEventListener("input", read);
    inputs[i][j] = inp; grid.append(inp);
  }
  function read() {
    let ok = true;
    const A = [[0, 0], [0, 0]];
    for (let i = 0; i < 2; i++) for (let j = 0; j < 2; j++) {
      const v = parseNum(inputs[i][j].value);
      inputs[i][j].classList.toggle("w-invalid", v == null);
      inputs[i][j].setAttribute("aria-invalid", String(v == null));
      if (v == null) ok = false; else A[i][j] = v;
    }
    msg.textContent = ok ? "" : "ورودی نامعتبر است؛ عدد، اعشار یا کسر (مثل 1/2) وارد کنید.";
    if (ok) onChange(A);
  }
  function set(A) {
    for (let i = 0; i < 2; i++) for (let j = 0; j < 2; j++) {
      inputs[i][j].value = fmtEntry(A[i][j]);
      inputs[i][j].classList.remove("w-invalid"); inputs[i][j].setAttribute("aria-invalid", "false");
    }
    msg.textContent = "";
  }
  return { grid, msg, set };
}
const sameMat = (A, B) => A.every((r, i) => r.every((v, j) => Math.abs(v - B[i][j]) < 1e-12));
function presetRow(presets, apply) {
  const btns = presets.map(p => h("button", { class: "btn", type: "button", "aria-pressed": "false", onclick: () => apply(p.A) }, p.label));
  return {
    row: h("div", { class: "w-row" }, btns),
    mark(A) { presets.forEach((p, i) => btns[i].setAttribute("aria-pressed", String(sameMat(p.A, A)))); },
  };
}

/* ===================== Widget: matnorm ===================== */
const MN_PRESETS = [
  { label: "نمونه", A: [[1, -2], [3, 4]] },
  { label: "متقارن", A: [[2, 1], [1, 2]] },
  { label: "مثال ۸.۲", A: [[0.5, 0], [0.25, 0.5]] },
  { label: "مثال ۷.۲", A: [[0, 2], [0, 0]] },
  { label: "مثال ۹.۲ (B)", A: [[0.5, 0], [100, 0.5]] },
];

function mountMatNorm(host) {
  const id = uid();
  const { root, body } = card("نرم‌های یک ماتریس ۲×۲", "دایرهٔ واحد زیر A به بیضی تبدیل می‌شود؛ بلندترین قطرش ⁦‖A‖₂⁩ است");
  const VB = 480, C = VB / 2, HALF = 212;
  const st = { A: MN_PRESETS[0].A.map(r => r.slice()), sq: false };

  const mi = mkMatInput("w2-mn" + id, A => { st.A = A; pr.mark(A); update(); });
  const pr = presetRow(MN_PRESETS, A => { st.A = A.map(r => r.slice()); mi.set(st.A); pr.mark(st.A); update(); });
  const inRow = h("div", { class: "w-row" }, h("span", { class: "w-note", html: ktx("A=") }), mi.grid);
  body.append(inRow, mi.msg, pr.row);

  const svg = sv("svg", { viewBox: "0 0 " + VB + " " + VB, class: "w-svg", role: "img", "aria-label": "دایرهٔ واحد و تصویر آن زیر ماتریس A" });
  const gBack = sv("g"), gCurves = sv("g");
  svg.append(gBack, gCurves);
  const plot = h("div", { class: "w-plot", dir: "ltr" }, svg);

  const swSq = h("span", { class: "w-lg", hidden: true }, h("i", { class: "w-sw w-b3" }), h("span", { html: ktx("\\{Ax:\\|x\\|_\\infty=1\\}") }));
  const legend = h("div", { class: "w-legend", dir: "ltr" },
    h("span", { class: "w-lg" }, h("i", { class: "w-sw w2-swd w2-swd-ink" }), h("span", { html: ktx("\\|x\\|_2=1") })),
    h("span", { class: "w-lg" }, h("i", { class: "w-sw w-b1" }), h("span", { html: ktx("\\{Ax:\\|x\\|_2=1\\}") })),
    h("span", { class: "w-lg" }, h("i", { class: "w-sw w-b2" }), h("span", { html: ktx("Av,\\ \\|v\\|_2=1") })),
    swSq);
  const tog = h("button", { class: "btn", type: "button", "aria-pressed": "false", onclick: () => { st.sq = !st.sq; update(); } }, "نمایش تصویر مربع ⁦‖x‖∞ = 1⁩");
  const bars = h("div", { class: "w2-bars", dir: "ltr" });
  const checks = h("ul", { class: "w2-checks", dir: "ltr" });
  body.append(legend, plot, h("div", { class: "w-row" }, tog),
    h("div", { class: "w-cap" }, "مقایسهٔ نرم‌ها روی یک مقیاس"), bars,
    h("div", { class: "w-cap" }, "تمرین ۱۹ و قضیهٔ ۱۷.۲ برای ماتریس ۲×۲"), checks);

  const BAR_TEX = {
    rho: "\\rho(A)", n2: "\\|A\\|_2", nF: "\\|A\\|_F", n1: "\\|A\\|_1", ni: "\\|A\\|_\\infty", geo: "\\sqrt{\\|A\\|_1\\|A\\|_\\infty}",
  };

  function drawPlot() {
    const A = st.A;
    const ell = ellipsePts(A, 240), sq = squareImage(A), top = topSingular(A);
    let R = 1;
    ell.forEach(q => { R = Math.max(R, Math.abs(q[0]), Math.abs(q[1])); });
    if (st.sq) sq.forEach(q => { R = Math.max(R, Math.abs(q[0]), Math.abs(q[1])); });
    const E = R * 1.08, S = HALF / E, X = v => C + S * v, Y = v => C - S * v;
    gBack.textContent = ""; gCurves.textContent = "";
    const step = niceStep(E), n = Math.floor(E / step + 1e-9);
    const labels = [];
    for (let i = -n; i <= n; i++) if (i !== 0 && (n <= 4 || i % 2 === 0)) labels.push(i * step);
    const tmp = sv("g");
    gridPlane(tmp, VB, C, S, n * step, step, labels);
    gBack.append(tmp);
    gBack.append(sv("text", { x: C + HALF - 2, y: C - 8, "text-anchor": "end", class: "w-tick" }, "x₁"),
      sv("text", { x: C + 8, y: C - HALF + 10, class: "w-tick" }, "x₂"));
    gCurves.append(sv("circle", { cx: C, cy: C, r: S, class: "w2-unit" }));
    if (st.sq) {
      gCurves.append(sv("polygon", { points: poly(sq.map(q => [X(q[0]), Y(q[1])])), class: "w-ln w-s3" }));
      let best = sq[0];
      sq.forEach(q => { if (Math.max(Math.abs(q[0]), Math.abs(q[1])) > Math.max(Math.abs(best[0]), Math.abs(best[1])) + 1e-12) best = q; });
      const bn = Math.max(Math.abs(best[0]), Math.abs(best[1]));
      gCurves.append(sv("circle", { cx: X(best[0]), cy: Y(best[1]), r: 6, class: "w-mk w-fs3" }),
        sv("text", { x: Math.max(8, Math.min(VB - 8, X(best[0]))), y: Y(best[1]) + (best[1] >= 0 ? -12 : 22), "text-anchor": "middle", class: "w-pl" }, "‖A‖∞ = " + fmtN(bn)));
    }
    gCurves.append(sv("polygon", { points: poly(ell.map(q => [X(q[0]), Y(q[1])])), class: "w-ln w-s1" }));
    const ar = arrow("s2");
    gCurves.append(ar.g);
    ar.set(C, C, X(top.Av[0]), Y(top.Av[1]));
    const tx = X(top.Av[0]), ty = Y(top.Av[1]);
    gCurves.append(sv("text", { x: tx + (tx > C ? -8 : 8), y: ty + (ty > C ? 20 : -10), "text-anchor": tx > C ? "end" : "start", class: "w-pl" }, "‖A‖₂ = " + fmtN(top.sigma)));
  }

  function update() {
    const A = st.A, N = matNorms(A);
    pr.mark(A);
    tog.setAttribute("aria-pressed", String(st.sq));
    swSq.hidden = !st.sq;
    drawPlot();
    const rows = Object.keys(BAR_TEX).map(k => ({ k, v: N[k] })).sort((a, b) => b.v - a.v);
    const mx = rows[0].v > 0 ? rows[0].v : 1;
    bars.textContent = "";
    rows.forEach(r => bars.append(h("div", { class: "w2-brow" },
      h("span", { class: "w2-blab", html: ktx(BAR_TEX[r.k]) }),
      h("span", { class: "w2-btrack" }, h("i", { class: r.k === "rho" ? "w2-bo" : "w2-bd", style: "width:" + (r.v / mx * 100).toFixed(2) + "%" })),
      h("span", { class: "w-num w2-bval", dir: "ltr" }, fa(r.v)))));
    const ex = exercise19(A), f = x => latin(x, 3);
    const tex = [
      ["\\tfrac{1}{\\sqrt2}\\|A\\|_1=", "\\ \\le\\ \\|A\\|_2=", "\\ \\le\\ \\sqrt2\\,\\|A\\|_1="],
      ["\\tfrac{1}{\\sqrt2}\\|A\\|_\\infty=", "\\ \\le\\ \\|A\\|_2=", "\\ \\le\\ \\sqrt2\\,\\|A\\|_\\infty="],
      ["\\|A\\|_2=", "\\ \\le\\ \\|A\\|_F=", "\\ \\le\\ \\sqrt2\\,\\|A\\|_2="],
      ["\\|A\\|_2=", "\\ \\le\\ \\sqrt{\\|A\\|_1\\|A\\|_\\infty}="],
      ["\\rho(A)=", "\\ \\le\\ \\|A\\|_1="],
      ["\\rho(A)=", "\\ \\le\\ \\|A\\|_\\infty="],
      ["\\rho(A)=", "\\ \\le\\ \\|A\\|_2="],
    ];
    checks.textContent = "";
    ex.forEach((e, i) => {
      const t = tex[i].map((s, j) => s + f(e.vals[j])).join("");
      checks.append(h("li", { class: "w2-chk" }, h("span", { class: "w-tex", html: ktx(t) }), chip(e.ok)));
    });
  }

  host.textContent = ""; host.append(root);
  mi.set(st.A);
  update();
}

/* ===================== Widget: powers ===================== */
const PW_PRESETS = [
  { label: "مثال ۸.۲", A: [[0.5, 0], [0.25, 0.5]] },
  { label: "مثال ۹.۲ (A)", A: [[1, 0], [0.25, 0.5]] },
  { label: "مثال ۹.۲ (B)", A: [[0.5, 0], [100, 0.5]] },
  { label: "مثال ۷.۲ (C)", A: [[0, 2], [0, 0]] },
  { label: "چرخشی همگرا", A: [[0.6, 0.5], [-0.5, 0.6]] },
  { label: "چرخشی واگرا", A: [[0.9, -0.5], [0.5, 0.9]] },
];
const PW_K = 40, PW_FLOOR = 1e-16;
function powerSeries(A) {
  const seq = powers(A, PW_K), r = rho(A);
  const nk = seq.map(normInf), rk = seq.map((_, k) => Math.pow(r, k));
  return { seq, rho: r, nk, rk, n1: norm1(A), ni: normInf(A) };
}

function mountPowers(host) {
  const id = uid();
  const { root, body } = card("توان‌های ماتریس و همگرایی", "قضیهٔ ۱۸.۲: ⁦A^k → 0⁩ اگر و فقط اگر ⁦ρ(A) < 1⁩");
  const W = 520, H = 340, L = 64, R = 506, Tp = 14, Bt = 296;
  const X = k => L + (k / PW_K) * (R - L);
  const st = { A: PW_PRESETS[0].A.map(r => r.slice()), sel: 5, hover: null, D: null, lo: 0, hi: 1 };
  const Y = v => Bt - ((Math.log10(Math.max(v, PW_FLOOR)) - st.lo) / (st.hi - st.lo)) * (Bt - Tp);

  const mi = mkMatInput("w2-pw" + id, A => { st.A = A; pr.mark(A); recompute(); });
  const pr = presetRow(PW_PRESETS, A => { st.A = A.map(r => r.slice()); mi.set(st.A); pr.mark(st.A); recompute(); });
  body.append(h("div", { class: "w-row" }, h("span", { class: "w-note", html: ktx("A=") }), mi.grid), mi.msg, pr.row);

  const sRho = stat("\\rho(A)"), sN1 = stat("\\|A\\|_1"), sNi = stat("\\|A\\|_\\infty");
  const cards = h("div", { class: "w-cards" }, sRho.el, sN1.el, sNi.el);
  const verdict = h("p", { class: "w-verdict", "aria-live": "polite" });
  const note = h("p", { class: "w-note", "aria-live": "polite" });
  body.append(cards, verdict, note);

  const svg = sv("svg", { viewBox: "0 0 " + W + " " + H, class: "w-svg", role: "img", "aria-label": "نمودار لگاریتمی نرم A^k و ρ(A)^k برای k از ۰ تا ۴۰" });
  const gStatic = sv("g");
  const cross = sv("line", { y1: Tp, y2: Bt, class: "w-cross", style: "display:none" });
  const m1 = sv("circle", { r: 7, class: "w-mk w-fs1" }), m2 = sv("circle", { r: 7, class: "w-mk w-fs2" });
  const catcher = sv("rect", { x: L, y: Tp, width: R - L, height: Bt - Tp, class: "w-catch" });
  svg.append(gStatic, cross, m2, m1, catcher);
  const tip = h("div", { class: "w-tip w-tip-wide", hidden: true });
  const plot = h("div", { class: "w-plot", dir: "ltr" }, svg, tip);
  const legend = h("div", { class: "w-legend", dir: "ltr" },
    h("span", { class: "w-lg" }, h("i", { class: "w-sw w-b1" }), h("span", { html: ktx("\\|A^k\\|_\\infty") })),
    h("span", { class: "w-lg" }, h("i", { class: "w-sw w2-swd w2-swd-s2" }), h("span", { html: ktx("\\rho(A)^k") })));
  const floorNote = h("p", { class: "w-note", hidden: true }, "مقدارهای صفر (یا کمتر از ۱۰⁻¹۶) روی مقیاس لگاریتمی قابل نمایش نیستند؛ در نمودار با کف ۱۰⁻¹۶ جایگزین شده‌اند.");

  const sl = mkSlider("w2-pw" + id + "-k", 0, PW_K, 1, st.sel, "توان k", v => { st.sel = v; updateK(); });
  const matOut = h("div", { class: "w-tex", dir: "ltr", "aria-live": "polite" });
  body.append(legend, plot, floorNote, sl.row, matOut);

  function drawChart() {
    const D = st.D;
    const logs = D.nk.concat(D.rk).map(v => Math.log10(Math.max(v, PW_FLOOR)));
    let lo = Math.floor(Math.min(...logs) - 0.05), hi = Math.ceil(Math.max(...logs) + 0.05);
    if (hi - lo < 2) hi = lo + 2;
    st.lo = lo; st.hi = hi;
    gStatic.textContent = "";
    const stepD = Math.ceil((hi - lo) / 8);
    for (let e = lo; e <= hi; e += stepD) {
      const y = Bt - ((e - lo) / (hi - lo)) * (Bt - Tp);
      gStatic.append(sv("line", { x1: L, y1: y, x2: R, y2: y, class: "w-gridl" }),
        sv("text", { x: L - 6, y: y + 4, "text-anchor": "end", class: "w-tick" }, "10" + sup(e)));
    }
    for (let k = 0; k <= PW_K; k += 5) {
      gStatic.append(sv("line", { x1: X(k), y1: Bt, x2: X(k), y2: Bt + 4, class: "w-axisl" }),
        sv("text", { x: X(k), y: Bt + 20, "text-anchor": "middle", class: "w-tick" }, String(k)));
    }
    gStatic.append(sv("line", { x1: L, y1: Bt, x2: R, y2: Bt, class: "w-axisl" }), sv("line", { x1: L, y1: Tp, x2: L, y2: Bt, class: "w-axisl" }),
      sv("text", { x: R, y: Bt + 38, "text-anchor": "end", class: "w-tick" }, "k"),
      sv("text", { x: L + 6, y: Tp + 12, class: "w-tick" }, "log₁₀"));
    gStatic.append(
      sv("polyline", { points: poly(D.rk.map((v, k) => [X(k), Y(v)])), class: "w-ln w-s2", "stroke-dasharray": "7 5" }),
      sv("polyline", { points: poly(D.nk.map((v, k) => [X(k), Y(v)])), class: "w-ln w-s1" }));
  }
  function marker() {
    const k = st.hover != null ? st.hover : st.sel, D = st.D;
    m1.setAttribute("cx", X(k)); m1.setAttribute("cy", Y(D.nk[k]));
    m2.setAttribute("cx", X(k)); m2.setAttribute("cy", Y(D.rk[k]));
  }
  const fv = v => (v === 0 ? "0 (کف نمودار)" : fmtN(v));
  function hover(e) {
    const rc = svg.getBoundingClientRect();
    const vx = ((e.clientX - rc.left) / rc.width) * W;
    const k = Math.max(0, Math.min(PW_K, Math.round(((vx - L) / (R - L)) * PW_K)));
    st.hover = k; marker();
    cross.style.display = ""; cross.setAttribute("x1", X(k)); cross.setAttribute("x2", X(k));
    tip.hidden = false;
    tip.innerHTML = "<b>k = " + k + "</b><div class=\"w-tl\"><i class=\"w-sw w-b1\"></i>‖A^k‖∞ = " + esc(fv(st.D.nk[k])) +
      "</div><div class=\"w-tl\"><i class=\"w-sw w2-swd w2-swd-s2\"></i>ρ^k = " + esc(fv(st.D.rk[k])) + "</div>";
    const frac = X(k) / W;
    tip.style.left = frac * 100 + "%"; tip.style.top = "6%";
    tip.style.transform = frac > 0.55 ? "translateX(calc(-100% - 10px))" : "translateX(10px)";
    if (e.type === "pointerdown") { st.sel = k; sl.input.value = k; updateK(); }
  }
  function unhover() { st.hover = null; cross.style.display = "none"; tip.hidden = true; marker(); }
  catcher.addEventListener("pointermove", hover);
  catcher.addEventListener("pointerdown", hover);
  catcher.addEventListener("pointerleave", unhover);

  function updateK() {
    const k = st.sel, M = st.D.seq[k];
    sl.out.textContent = "k = " + k;
    const t = "A^{" + k + "}=\\begin{bmatrix}" + sig4tex(M[0][0]) + "&" + sig4tex(M[0][1]) + "\\\\" + sig4tex(M[1][0]) + "&" + sig4tex(M[1][1]) +
      "\\end{bmatrix},\\qquad \\|A^{" + k + "}\\|_\\infty=" + sig4tex(st.D.nk[k]);
    matOut.innerHTML = ktx(t);
    marker();
  }
  function recompute() {
    const A = st.A, D = powerSeries(A);
    st.D = D; st.hover = null; cross.style.display = "none"; tip.hidden = true;
    pr.mark(A);
    sRho.val.textContent = fa(D.rho); sN1.val.textContent = fa(D.n1); sNi.val.textContent = fa(D.ni);
    const conv = D.rho < 1;
    verdict.className = "w-verdict " + (conv ? "w-vok" : "w-vbad");
    verdict.innerHTML = conv ? rich("$\\rho(A)<1$ ⇒ A همگراست (قضیهٔ ۱۸.۲)") : rich("$\\rho(A)\\ge1$ ⇒ A همگرا نیست");
    floorNote.hidden = !(D.nk.some(v => v < PW_FLOOR) || D.rk.some(v => v < PW_FLOOR));
    let kmax = 0;
    D.nk.forEach((v, k) => { if (v > D.nk[kmax]) kmax = k; });
    const rlt = "$\\rho(A)=" + latin(D.rho, 3) + "\\lt1$", rge = "$\\rho(A)=" + latin(D.rho, 3) + "\\ge1$";
    let msg;
    const zeroAt = D.seq.findIndex(M => M[0][0] === 0 && M[0][1] === 0 && M[1][0] === 0 && M[1][1] === 0);
    if (conv && zeroAt > 0) msg = "$A^k$ برای $k\\ge" + zeroAt + "$ دقیقاً صفر است (ماتریس پوچ‌توان)؛ چون " + rlt + " همگرا است.";
    else if (conv && D.nk[kmax] > 1 + 1e-9) msg = "نرم $A^k$ اول تا حدود $k=" + kmax + "$ بزرگ می‌شود (تا حدود " + fa(D.nk[kmax], 3) + ") ولی چون " + rlt + " در نهایت به صفر می‌رود.";
    else if (conv) msg = "نرم $A^k$ از همان ابتدا کاهش می‌یابد و چون " + rlt + " به صفر می‌رود.";
    else msg = "چون " + rge + " نرم $A^k$ به صفر نمی‌رود؛ در $k=40$ برابر " + fa(D.nk[PW_K], 3) + " است.";
    note.innerHTML = rich(msg);
    drawChart();
    updateK();
  }

  host.textContent = ""; host.append(root);
  mi.set(st.A);
  recompute();
}

/* ===================== Mounting & exports ===================== */
const REGISTRY = { unitball: mountUnitBall, cs: mountCS, matnorm: mountMatNorm, powers: mountPowers };
function mountAll() {
  document.querySelectorAll("[data-widget]").forEach(el => {
    const f = REGISTRY[el.getAttribute("data-widget")];
    if (!f) return;
    try { f(el); } catch (e) { if (typeof console !== "undefined") console.error(e); }
  });
}

globalThis.LAW2 = LAW2;
if (typeof module !== "undefined" && module.exports) module.exports = LAW2;

if (typeof document !== "undefined") {
  if (document.readyState === "loading") document.addEventListener("DOMContentLoaded", mountAll);
  else mountAll();
}
})();
