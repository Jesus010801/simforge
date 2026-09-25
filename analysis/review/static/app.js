// SimForge review dashboard — presentation only.
// Every scientific statement shown here comes from the ReviewDataset served by
// the local runtime; nothing is computed, interpolated or reinterpreted.
"use strict";

const TOKEN = new URLSearchParams(location.search).get("token") || "";
const PALETTE = ["#1f6feb", "#d1242f", "#1a7f37", "#9a6700", "#8250df", "#0a7ea4", "#bf3989"];
const state = { session: null, frame: 0, time: null, times: [], charts: [], lastSeq: 0 };

const $ = (id) => document.getElementById(id);
const el = (tag, cls, text) => {
  const e = document.createElement(tag);
  if (cls) e.className = cls;
  if (text !== undefined && text !== null) e.textContent = String(text);
  return e;
};
const badge = (s) => el("span", "state s-" + String(s).toLowerCase(), s);
const fmt = (v) => (v === null || v === undefined) ? "—" :
  (Math.abs(v) >= 1e4 || (Math.abs(v) < 1e-3 && v !== 0)) ? v.toExponential(4) : (+v).toPrecision(6);

async function api(path, body) {
  const opt = { headers: { "X-SimForge-Token": TOKEN } };
  if (body !== undefined) {
    opt.method = "POST";
    opt.headers["Content-Type"] = "application/json";
    opt.body = JSON.stringify(body);
  }
  const r = await fetch(path, opt);
  const data = await r.json();
  if (!r.ok) throw new Error(data.error || r.statusText);
  return data;
}

// ── sync ────────────────────────────────────────────────────────────────────
let clientSeq = 0, pendingFrame = null, throttle = null;
function requestFrame(frame, final) {
  pendingFrame = frame;
  if (final) { flushFrame(); return; }
  if (!throttle) throttle = setTimeout(flushFrame, 100);   // UI request rate only
}
function flushFrame() {
  clearTimeout(throttle); throttle = null;
  if (pendingFrame === null) return;
  const f = pendingFrame; pendingFrame = null;
  api("/api/sync/frame", { frame: f, client_event_id: "b" + (++clientSeq) })
    .then((r) => r.event && showMapping(r.event)).catch((e) => showError(e));
}
function requestTime(t) {
  api("/api/sync/time", { time_ps: t, client_event_id: "b" + (++clientSeq) })
    .then((r) => r.event && showMapping(r.event)).catch((e) => showError(e));
}
function showError(e) { $("mapping").textContent = "⚠ " + e.message; }

function showMapping(ev) {
  let t = `${ev.origin} → frame ${ev.frame} at ${ev.time_ps} ps (${ev.mapping_status})`;
  if (ev.requested_time_ps !== null && ev.requested_time_ps !== undefined)
    t += `; requested ${ev.requested_time_ps} ps, Δ=${fmt(ev.delta_ps)} ps${ev.tie ? ", tie → earlier frame" : ""}`;
  if (ev.duplicate_frames && ev.duplicate_frames.length > 1)
    t += `; frames sharing this time: ${ev.duplicate_frames.join(", ")}`;
  $("mapping").textContent = t;
}

function applyState(frame, time, ev) {
  if (ev && ev.sequence <= state.lastSeq) return;           // stale / duplicate delivery
  if (ev) state.lastSeq = ev.sequence;
  state.frame = frame; state.time = time;
  $("slider").value = frame; $("frame").value = frame; $("time").value = time;
  if (ev) showMapping(ev);
  for (const c of state.charts) c.cursor(frame, time);
}

function connect() {
  const es = new EventSource("/api/sync/events?token=" + encodeURIComponent(TOKEN));
  es.addEventListener("open", () => { $("connection").textContent = "live"; });
  es.addEventListener("error", () => { $("connection").textContent = "reconnecting…"; });
  es.addEventListener("state", (m) => {
    const s = JSON.parse(m.data);
    state.lastSeq = Math.max(0, s.sequence - 1);
    applyState(s.frame, s.time_ps, s.last_event);
  });
  es.addEventListener("sync", (m) => {
    const ev = JSON.parse(m.data);
    applyState(ev.frame, ev.time_ps, ev);
  });
}

// ── charts (canvas; one base draw, cursor on an overlay) ─────────────────────
function sampleForFrame(cursor, frame) {
  if (!cursor) return null;
  if (cursor.mode === "identity") return frame;
  if (cursor.mode === "frames") { const s = cursor.frame_samples[frame]; return s >= 0 ? s : null; }
  return null;
}

function makeCanvas(parent) {
  const stack = el("div", "stack"), base = el("canvas"), over = el("canvas", "overlay");
  stack.append(base, over); parent.append(stack);
  const size = () => {
    const w = stack.clientWidth || 600, h = 190, r = window.devicePixelRatio || 1;
    for (const c of [base, over]) { c.width = w * r; c.height = h * r; c.style.width = w + "px"; c.style.height = h + "px"; }
    return { w, h, r };
  };
  return { base, over, size };
}

function extent(arrays) {
  let lo = Infinity, hi = -Infinity;
  for (const a of arrays) for (const v of a) if (v !== null) { if (v < lo) lo = v; if (v > hi) hi = v; }
  if (lo === Infinity) return [0, 1];
  if (lo === hi) { lo -= 0.5; hi += 0.5; }
  return [lo, hi];
}

function drawAxes(ctx, box, xr, yr, xl, yl) {
  const css = getComputedStyle(document.body);
  ctx.strokeStyle = css.getPropertyValue("--line"); ctx.fillStyle = css.getPropertyValue("--muted");
  ctx.font = "11px system-ui"; ctx.lineWidth = 1;
  ctx.strokeRect(box.x, box.y, box.w, box.h);
  ctx.fillText(fmt(yr[1]), 4, box.y + 10); ctx.fillText(fmt(yr[0]), 4, box.y + box.h);
  ctx.fillText(fmt(xr[0]), box.x, box.y + box.h + 14);
  const hi = fmt(xr[1]); ctx.fillText(hi, box.x + box.w - ctx.measureText(hi).width, box.y + box.h + 14);
  ctx.fillText(xl, box.x + box.w / 2 - 20, box.y + box.h + 14);
  ctx.save(); ctx.translate(12, box.y + box.h / 2); ctx.rotate(-Math.PI / 2);
  ctx.textAlign = "center"; ctx.fillText(yl.length > 28 ? yl.slice(0, 27) + "…" : yl, 0, 0); ctx.restore();
}

function lineChart(parent, xs, series, xl, yl, withCursor, payload) {
  const cv = makeCanvas(parent);
  const legend = el("div", "legend"); parent.append(legend);
  series.forEach((s, i) => { const sp = el("span", null, "■ " + s.label); sp.style.color = PALETTE[i % PALETTE.length]; legend.append(sp); });
  const xr = extent([xs]), yr = extent(series.map((s) => s.values));
  let geo = null;
  const draw = () => {
    const { w, h, r } = cv.size();
    const ctx = cv.base.getContext("2d"); ctx.setTransform(r, 0, 0, r, 0, 0);
    const box = { x: 70, y: 8, w: w - 78, h: h - 30 };
    geo = { box, r, X: (x) => box.x + (x - xr[0]) / (xr[1] - xr[0]) * box.w,
            Y: (y) => box.y + box.h - (y - yr[0]) / (yr[1] - yr[0]) * box.h };
    drawAxes(ctx, box, xr, yr, xl, yl);
    series.forEach((s, i) => {
      ctx.strokeStyle = PALETTE[i % PALETTE.length]; ctx.beginPath(); let pen = false;
      for (let k = 0; k < xs.length; k++) {           // every stored sample, as stored
        const v = s.values[k];
        if (v === null || xs[k] === null) { pen = false; continue; }
        const px = geo.X(xs[k]), py = geo.Y(v);
        if (pen) ctx.lineTo(px, py); else { ctx.moveTo(px, py); pen = true; }
      }
      ctx.stroke();
    });
  };
  draw();
  window.addEventListener("resize", () => { draw(); api_cursor(); });
  const readout = el("div", "value muted"); parent.append(readout);
  let last = null;
  function api_cursor() { if (last) cursor(last[0], last[1]); }
  function cursor(frame, time) {
    last = [frame, time];
    const ctx = cv.over.getContext("2d"); ctx.setTransform(geo.r, 0, 0, geo.r, 0, 0);
    ctx.clearRect(0, 0, cv.over.width, cv.over.height);
    if (!withCursor || time === null) return;
    const css = getComputedStyle(document.body);
    ctx.strokeStyle = css.getPropertyValue("--cursor"); ctx.lineWidth = 1.5;
    const px = geo.X(time); ctx.beginPath(); ctx.moveTo(px, geo.box.y); ctx.lineTo(px, geo.box.y + geo.box.h); ctx.stroke();
    const s = sampleForFrame(payload.cursor, frame);
    if (s === null || s >= xs.length) { readout.textContent = `t = ${time} ps: value unavailable at this frame (no exact sample)`; return; }
    readout.textContent = `t = ${xs[s]} ps (exact sample ${s}): ` +
      series.map((sr) => `${sr.label} = ${fmt(sr.values[s])} ${payload.unit || ""}`).join(", ");
    series.forEach((sr, i) => { const v = sr.values[s]; if (v === null) return;
      ctx.fillStyle = PALETTE[i % PALETTE.length]; ctx.beginPath(); ctx.arc(geo.X(xs[s]), geo.Y(v), 3, 0, 7); ctx.fill(); });
  }
  return { cursor };
}

function panelFor(desc) {
  const p = el("section", "card panel");
  const head = el("div", "head");
  head.append(el("h3", null, `${desc.display_name || desc.observable} — ${desc.instance_id}`));
  const tags = el("span", "muted", `${desc.quantity} [${desc.unit || "–"}] · axes (${desc.axis_signature.join(", ")}) · ${desc.availability}` +
    (desc.syncable ? " · time-synchronized" : " · not time-syncable"));
  head.append(tags); p.append(head);
  const v = desc.view;
  p.append(el("div", "muted", `computed on: ${v.analysis_kind || "?"}` +
    (v.analysis_operations.length ? ` (${v.analysis_operations.join(" + ")})` : "") +
    ` · displayed coordinates: ${v.display_kind}`));
  if (!v.same_coordinates) p.append(el("div", "note", v.note));
  if (desc.externally_supplied) p.append(el("div", "note", "EXPLICIT IMPORT — externally supplied, not independently verified"));
  const body = el("div"); p.append(body);
  $("panels").append(p);
  if (desc.renderer === "unsupported") {
    body.append(el("div", null, `Scientific result available. Renderer not implemented for axis signature: (${desc.axis_signature.join(", ")})`));
    body.append(el("div", "muted", `axes: ${desc.axes.map((a) => `${a.name}:${a.kind}${a.unit ? " [" + a.unit + "]" : ""}`).join(", ")} · storage: ${desc.storage_format}`));
    return;
  }
  body.append(el("div", "muted", "loading data…"));
  api("/api/results/" + desc.id).then((d) => {
    body.textContent = "";
    if (d.renderer === "profile") {
      body.append(el("div", "muted", `static result over ${d.data.x_label} (${d.data.x_kind}) — no trajectory cursor`));
      const numeric = d.data.x.every((x) => typeof x === "number");
      if (numeric) lineChart(body, d.data.x, [{ label: d.name, values: d.data.values }], d.data.x_label, d.unit || "", false, d);
      const t = el("table"), wrap = el("div", "scroll"); wrap.append(t); body.append(wrap);
      d.data.x.forEach((x, i) => { const tr = el("tr"); tr.append(el("td", null, x), el("td", "value", fmt(d.data.values[i]))); t.append(tr); });
      return;
    }
    const chart = lineChart(body, d.data.time_ps, d.data.series, "time (ps)", `${d.quantity} [${d.unit || "–"}]`, d.syncable, d);
    if (!d.syncable) body.append(el("div", "note", "not time-syncable: " + (d.cursor.reason || d.sync.reason)));
    if (d.data.non_finite) body.append(el("div", "note", `${d.data.non_finite} non-finite stored values shown as gaps`));
    state.charts.push(chart); chart.cursor(state.frame, state.time);
  }).catch((e) => { body.textContent = "result unavailable: " + e.message; });
}

// ── side panels (stored text only) ──────────────────────────────────────────
function renderSide(s) {
  const b = $("blocked");
  if (!s.blocked.length) b.append(el("div", "muted", "none — every requested observable is available"));
  for (const o of s.blocked) {
    const d = el("div"); d.append(el("b", null, o.instance_id + " "), badge(o.state));
    d.append(el("div", "muted", o.reason)); b.append(d);
  }
  const g = s.diagnostics, dg = $("diagnostics");
  dg.append(el("div", null, `${g.execution} / ${g.status}` + (g.reason ? ` — ${g.reason}` : "")));
  if (g.counts) dg.append(el("div", "muted", `warnings ${g.counts.warnings} · review required ${g.counts.review_required} · errors ${g.counts.errors} · info ${g.counts.info}`));
  for (const f of g.findings || []) {
    const d = el("div"); d.append(badge(f.severity), el("span", null, ` ${f.code}: ${f.message}`)); dg.append(d);
  }
  const labels = { ran_with_findings: "findings", ran_found_nothing: "ran / no finding",
                   not_applicable: "not applicable", deferred: "deferred", failed: "failed" };
  for (const [k, lab] of Object.entries(labels)) {
    const list = (g.detectors || {})[k] || [];
    if (list.length) dg.append(el("div", "muted", `${lab}: ${list.map((x) => x.detector).join(", ")}`));
  }
  const t = el("table"); t.append(el("tr"));
  t.firstChild.append(el("th", null, "id"), el("th", null, "kind"), el("th", null, "state"), el("th", null, "origin"));
  for (const a of s.annotations) {
    const tr = el("tr"), st = el("td"); st.append(badge(a.state));
    tr.append(el("td", null, a.annotation_id), el("td", null, a.kind || ""), st, el("td", null, a.origin || ""));
    tr.title = (a.reasons || []).join("; "); t.append(tr);
  }
  const wrap = el("div", "scroll"); wrap.append(t); $("annotations").append(wrap);
}

async function main() {
  const s = await api("/api/session");
  state.session = s; state.times = s.timeline.times_ps;
  $("identity").textContent = `session ${s.session_id} · system ${s.system.system_id} · ${s.sources.trajectory || ""}`;
  const tl = s.timeline;
  $("timeline").textContent = `${tl.n_frames} frames · ${tl.start_time_ps} → ${tl.end_time_ps} ps · timeline ${tl.state}` +
    (tl.n_duplicate_groups ? ` · ${tl.n_duplicate_groups} duplicate-time groups (kept)` : "");
  $("viewer").textContent = `viewer: ${s.viewer.kind}${s.viewer.connected ? "" : " (not connected)"}`;
  const dv = s.display_view;
  $("display").textContent = `Displayed coordinates: ${dv.kind} (${dv.status}; purpose ${dv.purpose}; frames ${tl.display_relation})` +
    (dv.reason ? ` — ${dv.reason}` : "");
  const sl = $("slider"); sl.max = tl.n_frames - 1;
  $("frame").max = tl.n_frames - 1;
  sl.addEventListener("input", () => requestFrame(+sl.value, false));
  sl.addEventListener("change", () => requestFrame(+sl.value, true));
  $("frame").addEventListener("change", () => requestFrame(parseInt($("frame").value, 10), true));
  $("time").addEventListener("change", () => requestTime(parseFloat($("time").value)));
  renderSide(s);
  if (!s.results.length) $("panels").append(el("section", "card muted", "no available results in this session"));
  for (const d of s.results) panelFor(d);
  applyState(s.sync.frame, s.sync.time_ps, null);
  connect();
}
main().catch((e) => { $("identity").textContent = "cannot load session: " + e.message; });
