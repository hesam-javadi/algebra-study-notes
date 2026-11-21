/* جزوه جبر خطی عددی — page infrastructure
   theme · math rendering · table of contents · reading progress ·
   step-by-step proofs · quizzes · "learned" marks */
(function () {
  "use strict";

  var FA = "۰۱۲۳۴۵۶۷۸۹";
  function toFa(x) { return String(x).replace(/[0-9]/g, function (d) { return FA[+d]; }); }
  window.toFa = toFa;

  var store = {
    get: function (k, def) { try { var v = localStorage.getItem(k); return v === null ? def : JSON.parse(v); } catch (e) { return def; } },
    set: function (k, v) { try { localStorage.setItem(k, JSON.stringify(v)); } catch (e) { /* storage unavailable */ } }
  };
  var PAGE = document.documentElement.getAttribute("data-page") || "page";

  /* ---------- theme ---------- */
  var root = document.documentElement;
  var savedTheme = store.get("theme", null);
  if (savedTheme === "dark" || savedTheme === "light") root.setAttribute("data-theme", savedTheme);
  function isDark() {
    var t = root.getAttribute("data-theme");
    if (t) return t === "dark";
    return window.matchMedia && window.matchMedia("(prefers-color-scheme: dark)").matches;
  }

  document.addEventListener("DOMContentLoaded", function () {
    var themeBtn = document.getElementById("themeBtn");
    if (themeBtn) {
      var paint = function () {
        themeBtn.setAttribute("aria-label", isDark() ? "حالت روشن" : "حالت تیره");
        themeBtn.innerHTML = isDark()
          ? '<svg viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="2" stroke-linecap="round"><circle cx="12" cy="12" r="4.5"/><path d="M12 2v2M12 20v2M4.9 4.9l1.4 1.4M17.7 17.7l1.4 1.4M2 12h2M20 12h2M4.9 19.1l1.4-1.4M17.7 6.3l1.4-1.4"/></svg>'
          : '<svg viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="2" stroke-linecap="round" stroke-linejoin="round"><path d="M20 14.5A8 8 0 1 1 9.5 4a6.5 6.5 0 0 0 10.5 10.5z"/></svg>';
      };
      paint();
      themeBtn.addEventListener("click", function () {
        var next = isDark() ? "light" : "dark";
        root.setAttribute("data-theme", next);
        store.set("theme", next);
        paint();
      });
    }

    /* ---------- math ---------- */
    if (window.renderMathInElement) {
      renderMathInElement(document.body, {
        delimiters: [
          { left: "$$", right: "$$", display: true },
          { left: "\\[", right: "\\]", display: true },
          { left: "$", right: "$", display: false }
        ],
        ignoredClasses: ["widget", "no-math"],
        throwOnError: false
      });
    }

    fitInlineMath();
    var fitT;
    window.addEventListener("resize", function () { clearTimeout(fitT); fitT = setTimeout(fitInlineMath, 150); });
    if (document.fonts && document.fonts.ready) document.fonts.ready.then(fitInlineMath);

    buildToc();
    readingProgress();
    proofs();
    quizzes();
    doneMarks();
  });

  /* ---------- inline math: one line unless it is wider than its line box ---------- */
  function fitInlineMath() {
    document.querySelectorAll(".content .katex").forEach(function (k) {
      if (k.closest(".katex-display") || k.closest(".widget")) return;
      k.classList.remove("k-wrap");
      var host = k.parentElement && k.parentElement.closest("p, li, td, th, summary, button, h2, h3, h4, h5, div");
      if (!host || !host.clientWidth) return;
      var cs = getComputedStyle(host);
      var avail = host.clientWidth - parseFloat(cs.paddingLeft) - parseFloat(cs.paddingRight);
      if (k.getBoundingClientRect().width > avail - 2) k.classList.add("k-wrap");
    });
  }

  /* ---------- table of contents ---------- */
  function buildToc() {
    var list = document.getElementById("tocList");
    var toc = document.getElementById("toc");
    if (!list || !toc) return;
    var heads = document.querySelectorAll(".content h2[id], .content h3[id]");
    var links = [];
    heads.forEach(function (h) {
      var li = document.createElement("li");
      li.className = h.tagName === "H2" ? "lvl2" : "lvl3";
      var a = document.createElement("a");
      a.href = "#" + h.id;
      var num = h.querySelector(".num");
      var title = h.getAttribute("data-toc") || Array.prototype.filter.call(h.childNodes, function (n) {
        return !(n.nodeType === 1 && n.classList.contains("num"));
      }).map(function (n) { return n.textContent; }).join("").trim();
      a.innerHTML = (num ? '<span class="n">' + num.textContent + "</span>" : "") +
        "<span>" + title + "</span>" + '<span class="chk" aria-hidden="true">✓</span>';
      a.setAttribute("data-target", h.id);
      li.appendChild(a);
      list.appendChild(li);
      links.push(a);
      a.addEventListener("click", function () { closeToc(); });
    });

    var btn = document.getElementById("tocBtn");
    var close = toc.querySelector(".toc-close");
    var backdrop = document.getElementById("backdrop");
    function openToc() { toc.classList.add("open"); if (backdrop) backdrop.classList.add("show"); if (btn) btn.setAttribute("aria-expanded", "true"); }
    function closeToc() { toc.classList.remove("open"); if (backdrop) backdrop.classList.remove("show"); if (btn) btn.setAttribute("aria-expanded", "false"); }
    if (btn) btn.addEventListener("click", function () { toc.classList.contains("open") ? closeToc() : openToc(); });
    if (close) close.addEventListener("click", closeToc);
    if (backdrop) backdrop.addEventListener("click", closeToc);
    document.addEventListener("keydown", function (e) { if (e.key === "Escape") closeToc(); });

    if ("IntersectionObserver" in window) {
      var current = null;
      var io = new IntersectionObserver(function (entries) {
        entries.forEach(function (en) {
          if (en.isIntersecting) {
            var id = en.target.id;
            if (current) current.classList.remove("active");
            current = list.querySelector('a[data-target="' + id + '"]');
            if (current) {
              current.classList.add("active");
              if (window.innerWidth >= 1100 && current.scrollIntoView) {
                var r = current.getBoundingClientRect(), tr = toc.getBoundingClientRect();
                if (r.top < tr.top || r.bottom > tr.bottom) current.scrollIntoView({ block: "nearest" });
              }
            }
          }
        });
      }, { rootMargin: "-70px 0px -70% 0px" });
      heads.forEach(function (h) { io.observe(h); });
    }
  }

  /* ---------- reading progress ---------- */
  function readingProgress() {
    var bar = document.querySelector(".read-progress span");
    if (!bar) return;
    var ticking = false;
    function upd() {
      var h = document.documentElement;
      var max = h.scrollHeight - h.clientHeight;
      bar.style.width = (max > 0 ? (h.scrollTop / max) * 100 : 0) + "%";
      ticking = false;
    }
    window.addEventListener("scroll", function () { if (!ticking) { ticking = true; requestAnimationFrame(upd); } }, { passive: true });
    upd();
  }

  /* ---------- proofs: full view ⇄ step-by-step practice ---------- */
  function proofs() {
    var all = document.querySelectorAll(".proof");
    all.forEach(function (p) {
      var steps = p.querySelectorAll(".step");
      if (!steps.length) return;
      var head = p.querySelector(".proof-h");
      var toggle = document.createElement("button");
      toggle.type = "button";
      toggle.className = "btn";
      toggle.setAttribute("aria-pressed", "false");
      toggle.textContent = "تمرین گام‌به‌گام";
      head.appendChild(toggle);

      var bar = document.createElement("div");
      bar.className = "practice-bar";
      bar.innerHTML = '<button type="button" class="btn primary nx">گام بعدی را نشان بده</button>' +
        '<button type="button" class="btn all">نمایش کامل</button>' +
        '<span class="practice-hint">قبل از دیدن هر گام، روی کاغذ حدس بزن گام بعد چیست.</span>';
      var after = p.querySelector(".steps");
      after.parentNode.insertBefore(bar, after.nextSibling);
      var nx = bar.querySelector(".nx");
      var extras = p.querySelectorAll(".skeleton, .qed");
      var shown = 1;

      function render() {
        steps.forEach(function (s, i) {
          var hide = p.classList.contains("practice") && i >= shown;
          s.classList.toggle("concealed", hide);
        });
        var finished = shown >= steps.length;
        extras.forEach(function (x) { x.classList.toggle("concealed", p.classList.contains("practice") && !finished); });
        nx.disabled = finished;
        nx.textContent = finished ? "اثبات کامل شد ✓" : "گام بعدی را نشان بده (" + toFa(shown) + " از " + toFa(steps.length) + ")";
      }
      function setPractice(on) {
        p.classList.toggle("practice", on);
        toggle.setAttribute("aria-pressed", on ? "true" : "false");
        shown = 1;
        render();
      }
      p._setPractice = setPractice;
      toggle.addEventListener("click", function () { setPractice(!p.classList.contains("practice")); });
      nx.addEventListener("click", function () {
        if (shown < steps.length) {
          steps[shown].classList.add("revealed");
          shown++;
          render();
          steps[shown - 1].scrollIntoView && steps[shown - 1].scrollIntoView({ block: "nearest", behavior: "smooth" });
        }
      });
      bar.querySelector(".all").addEventListener("click", function () { shown = steps.length; render(); });
      render();
    });

    var global = document.getElementById("allPractice");
    if (global) {
      global.addEventListener("click", function () {
        var on = global.getAttribute("aria-pressed") !== "true";
        global.setAttribute("aria-pressed", on ? "true" : "false");
        global.textContent = on ? "همهٔ اثبات‌ها: حالت تمرین (روشن)" : "همهٔ اثبات‌ها را در حالت تمرین ببین";
        all.forEach(function (p) { if (p._setPractice) p._setPractice(on); });
      });
    }
  }

  /* ---------- quizzes ---------- */
  function quizzes() {
    var KEYS = ["آ", "ب", "پ", "ت", "ث"];
    document.querySelectorAll(".quiz").forEach(function (quiz) {
      var qs = quiz.querySelectorAll(".q");
      var scoreEl = quiz.querySelector(".score");
      qs.forEach(function (q, qi) {
        var t = q.querySelector(".q-t");
        if (t && !t.querySelector(".qn")) t.insertAdjacentHTML("afterbegin", '<span class="qn">' + toFa(qi + 1) + ".</span>");
        q.querySelectorAll(".opt").forEach(function (o, i) {
          o.type = "button";
          // keep the option text in one wrapping box (the button is a flex row)
          var label = document.createElement("span");
          label.className = "opt-t";
          while (o.firstChild) label.appendChild(o.firstChild);
          o.appendChild(label);
          o.insertAdjacentHTML("afterbegin", '<span class="k">' + (KEYS[i] || i + 1) + "</span>");
          o.addEventListener("click", function () {
            if (q.classList.contains("answered")) return;
            q.classList.add("answered");
            var ok = o.hasAttribute("data-c");
            o.classList.add(ok ? "right" : "wrong");
            q.setAttribute("data-ok", ok ? "1" : "0");
            q.querySelectorAll(".opt").forEach(function (x) {
              if (x.hasAttribute("data-c")) x.classList.add("right");
              x.disabled = true;
            });
            score();
          });
        });
      });
      function score() {
        var done = quiz.querySelectorAll(".q.answered").length;
        var ok = quiz.querySelectorAll('.q[data-ok="1"]').length;
        if (scoreEl) scoreEl.textContent = done ? toFa(ok) + " درست از " + toFa(done) + " (کل: " + toFa(qs.length) + ")" : toFa(qs.length) + " پرسش";
      }
      var reset = quiz.querySelector(".quiz-reset");
      if (reset) reset.addEventListener("click", function () {
        qs.forEach(function (q) {
          q.classList.remove("answered"); q.removeAttribute("data-ok");
          q.querySelectorAll(".opt").forEach(function (x) { x.classList.remove("right", "wrong"); x.disabled = false; });
        });
        score();
      });
      score();
    });
  }

  /* ---------- "I learned this section" marks ---------- */
  function doneMarks() {
    var key = PAGE + "-done";
    var done = store.get(key, {});
    var btns = document.querySelectorAll(".done-btn[data-sec]");
    var prog = document.getElementById("tocProg");
    function paint() {
      var n = 0;
      btns.forEach(function (b) {
        var id = b.getAttribute("data-sec");
        var on = !!done[id];
        if (on) n++;
        b.setAttribute("aria-pressed", on ? "true" : "false");
        b.textContent = on ? "✓ یاد گرفتم" : "این بخش را یاد گرفتم";
        var link = document.querySelector('#tocList a[data-target="' + id + '"]');
        if (link) link.classList.toggle("done", on);
      });
      if (prog) {
        prog.querySelector(".t").textContent = toFa(n) + " از " + toFa(btns.length) + " بخش";
        prog.querySelector(".bar span").style.width = (btns.length ? (n / btns.length) * 100 : 0) + "%";
      }
    }
    btns.forEach(function (b) {
      b.type = "button";
      b.addEventListener("click", function () {
        var id = b.getAttribute("data-sec");
        done[id] = !done[id];
        store.set(key, done);
        paint();
      });
    });
    paint();
  }
})();
