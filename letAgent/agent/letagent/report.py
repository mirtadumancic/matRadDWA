"""Writes the local report: report.html, runs.csv, pareto.csv, agent_log.md."""
import csv
import html
import json
import os

from .pareto import LET_METRICS, OBJECTIVE_LABELS
from .privacy import banner


def _as_list(x):
    if x is None:
        return []
    return x if isinstance(x, list) else [x]


def _row(r):
    p = r["params"]
    row = {"runId": r["runId"], "round": r["round"], "status": r["status"], "beamSet": p["beamSet"],
           "letMode": p["mode"], "letPenalty": round(p["penalty"], 3) if p["mode"] != "none" else "",
           "oarScale": round(p["oarScale"], 3), "robust": p["robust"], "coverageMode": p["coverageMode"],
           "seconds": round(r.get("seconds") or 0, 1)}
    if r["status"] != "ok":
        row["issues"] = r.get("error", "")
        return row
    for k, v in r["obj"].items():
        row[k] = round(v, 3)
    m = r["metrics"]
    for t in _as_list(m["targets"]):
        n = t["name"]
        row[n + " D95 nom"] = round(t["D95"]["nom"], 2)
        row[n + " D95 wc"] = round(t["D95"]["wc"], 2)
        row[n + " D2 nom"] = round(t["D2"]["nom"], 2)
        row[n + " V95 %"] = round(t["V95pct"]["nom"], 1)
    for o in _as_list(m["oars"]):
        n = o["name"]
        row[n + " D2"] = round(o["D2"]["nom"], 2)
        row[n + " Dmean"] = round(o["Dmean"]["nom"], 2)
        if o.get("let"):
            row[n + " dirty D2"] = round(o["dirtyDoseD2"]["nom"], 2)
            row[n + " LETd p98"] = round(o["LETdP98"]["nom"], 2)
            row[n + " VDL %"] = round(o["VDLpct"]["nom"], 2)
    row["issues"] = "; ".join(r.get("issues") or [])
    return row


def write(case_dir, case, state, front, settings):
    rows = [_row(r) for r in state["history"]]
    cols = []
    for r in rows:
        for k in r:
            if k not in cols:
                cols.append(k)
    for name, sub in (("runs.csv", rows), ("pareto.csv", [r for r in rows if r["runId"] in front])):
        with open(os.path.join(case_dir, name), "w", newline="") as f:
            w = csv.DictWriter(f, fieldnames=cols)
            w.writeheader()
            w.writerows(sub)

    with open(os.path.join(case_dir, "agent_log.md"), "w") as f:
        f.write("# Agent log: {}\n\n> {}\n\n".format(case["alias"], banner(case)))
        for e in state["log"]:
            f.write("- **r{} {}** ({}): {}\n".format(e["round"], e["rule"], e["time"], e["text"]))

    labels = dict(OBJECTIVE_LABELS)
    labels["letBurden"] = LET_METRICS[settings["letMetric"]][1]
    data = {"alias": case["alias"], "rx": case["prescription"]["totalDoseGy"], "front": front,
            "labels": labels, "rows": rows, "log": state["log"], "cols": cols,
            "letOars": [o["name"] for o in case["structures"].get("oars", []) if o.get("let")],
            "target": case["structures"]["targets"][0]["name"]}
    page = _TEMPLATE.replace("__DATA__", json.dumps(data).replace("</", "<\\/")) \
                    .replace("__TITLE__", html.escape("LET plan search " + case["alias"])) \
                    .replace("__BANNER__", html.escape(banner(case)))
    with open(os.path.join(case_dir, "report.html"), "w", encoding="utf-8") as f:
        f.write(page)


_TEMPLATE = r"""<!doctype html>
<html lang="en"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1">
<title>__TITLE__</title>
<style>
:root{--bg:#fbfaf8;--fg:#1d2433;--muted:#5d6577;--line:#d9dde5;--card:#fff;--accent:#1f5fbf;
--m0:#8a8f99;--m1:#c2410c;--m2:#0f766e;--ok:#e8f1ff}
@media (prefers-color-scheme:dark){:root{--bg:#12151c;--fg:#e6e9ef;--muted:#9aa3b5;--line:#2c3240;--card:#1a1f29;
--accent:#7aa7ff;--m0:#9aa0aa;--m1:#fb923c;--m2:#2dd4bf;--ok:#1e2c45}}
*{box-sizing:border-box}body{margin:0;background:var(--bg);color:var(--fg);font:14px/1.45 system-ui,-apple-system,Segoe UI,sans-serif}
main{max-width:1180px;margin:0 auto;padding:20px 16px 60px}
h1{font-size:22px;margin:0 0 4px}h2{font-size:16px;margin:28px 0 10px}
.banner{border:1px solid var(--line);border-left:4px solid var(--m1);background:var(--card);padding:8px 12px;margin:12px 0;color:var(--muted)}
.grid{display:grid;grid-template-columns:repeat(auto-fit,minmax(340px,1fr));gap:16px}
.card{background:var(--card);border:1px solid var(--line);border-radius:8px;padding:12px}
.stats{display:flex;gap:24px;flex-wrap:wrap;color:var(--muted)}.stats b{display:block;font-size:20px;color:var(--fg)}
svg{max-width:640px;display:block;margin:0 auto}svg text{fill:var(--muted);font-size:11px}.axis{stroke:var(--line)}
.legend span{display:inline-flex;align-items:center;gap:6px;margin-right:14px;color:var(--muted)}
.dot{width:10px;height:10px;border-radius:50%;display:inline-block}
label{display:block;color:var(--muted);margin:6px 0 2px}input[type=range]{width:100%}
.tbl{overflow-x:auto}table{border-collapse:collapse;width:100%;font-size:12.5px}
th,td{padding:5px 8px;border-bottom:1px solid var(--line);text-align:right;white-space:nowrap}
th:first-child,td:first-child{text-align:left}th{color:var(--muted);font-weight:600;position:sticky;top:0;background:var(--card)}
tr.meets td{background:var(--ok)}tr:hover td{outline:1px solid var(--line)}tr{cursor:pointer}
ol.log{padding-left:18px}ol.log li{margin:4px 0}.rule{font-weight:600;color:var(--accent)}
#detail pre{white-space:pre-wrap;font-size:12px;margin:0}
</style></head><body><main>
<h1 id="title"></h1>
<div class="banner">__BANNER__</div>
<div class="stats" id="stats"></div>
<h2>Trade-offs</h2>
<div class="legend" id="legend"></div>
<div class="grid" id="plots"></div>
<h2>Set your limits</h2>
<div class="grid"><div class="card" id="sliders"></div><div class="card" id="pick"></div></div>
<h2>Pareto-optimal plans</h2>
<div class="card tbl"><table id="front"></table></div>
<div class="card" id="detail" style="margin-top:12px;display:none"></div>
<h2>What the agent did</h2>
<div class="card"><ol class="log" id="log"></ol></div>
<h2>All plans</h2>
<div class="card tbl"><table id="all"></table></div>
</main>
<script>
const D=__DATA__;
const MODE_COL={none:'var(--m0)',dirtyDose:'var(--m1)',LETxDose:'var(--m2)'};
const MODE_NAME={none:'dose only',dirtyDose:'dirty dose objective',LETxDose:'LET×dose objective'};
const ok=D.rows.filter(r=>r.status==='ok');const front=new Set(D.front);
const keys=['coverageLoss','letBurden','robustLoss'].filter(k=>ok.length&&k in ok[0]);
document.getElementById('title').textContent='LET-aware plan search · '+D.alias;
const st=document.getElementById('stats');
[['Plans evaluated',D.rows.length],['On the Pareto front',D.front.length],['Failed',D.rows.length-ok.length],['Prescription',D.rx+' Gy']]
.forEach(([a,b])=>{st.insertAdjacentHTML('beforeend',`<div>${a}<b>${b}</b></div>`)});
document.getElementById('legend').innerHTML=Object.keys(MODE_COL).map(m=>`<span><i class="dot" style="background:${MODE_COL[m]}"></i>${MODE_NAME[m]}</span>`).join('')+'<span>large = Pareto-optimal</span>';
const lim={};keys.forEach(k=>{lim[k]=Math.max(...ok.map(r=>r[k]))});
function plot(xk,yk){
 const W=520,H=320,L=56,B=44,T=12,R=12;const xs=ok.map(r=>r[xk]),ys=ok.map(r=>r[yk]);
 const x0=Math.min(...xs),x1=Math.max(...xs),y0=Math.min(...ys),y1=Math.max(...ys);
 const sx=v=>L+(v-x0)/((x1-x0)||1)*(W-L-R),sy=v=>H-B-(v-y0)/((y1-y0)||1)*(H-B-T);
 let s=`<svg viewBox="0 0 ${W} ${H}" width="100%" role="img" aria-label="${D.labels[yk]} versus ${D.labels[xk]}">`;
 s+=`<line class="axis" x1="${L}" y1="${H-B}" x2="${W-R}" y2="${H-B}"/><line class="axis" x1="${L}" y1="${T}" x2="${L}" y2="${H-B}"/>`;
 for(let i=0;i<=4;i++){const vx=x0+(x1-x0)*i/4,vy=y0+(y1-y0)*i/4;
  s+=`<text x="${sx(vx)}" y="${H-B+14}" text-anchor="middle">${vx.toFixed(1)}</text><text x="${L-6}" y="${sy(vy)+4}" text-anchor="end">${vy.toFixed(1)}</text>`;}
 s+=`<text x="${(L+W-R)/2}" y="${H-6}" text-anchor="middle">${D.labels[xk]}</text>`;
 s+=`<text transform="translate(12 ${(T+H-B)/2}) rotate(-90)" text-anchor="middle">${D.labels[yk]}</text>`;
 const fp=ok.filter(r=>front.has(r.runId)).sort((a,b)=>a[xk]-b[xk]);
 ok.forEach(r=>{const f=front.has(r.runId);s+=`<circle data-id="${r.runId}" cx="${sx(r[xk])}" cy="${sy(r[yk])}" r="${f?6:3.5}" fill="${MODE_COL[r.letMode]}" fill-opacity="${f?0.95:0.45}" stroke="${f?'var(--fg)':'none'}" stroke-width="1"><title>${r.runId} · ${MODE_NAME[r.letMode]} · ${r.beamSet}</title></circle>`;});
 return s+'</svg>';}
const plots=document.getElementById('plots');
plots.insertAdjacentHTML('beforeend',`<div class="card">${plot('letBurden','coverageLoss')}</div>`);
if(keys.includes('robustLoss'))plots.insertAdjacentHTML('beforeend',`<div class="card">${plot('letBurden','robustLoss')}</div>`);
plots.addEventListener('click',e=>{const id=e.target.getAttribute&&e.target.getAttribute('data-id');if(id)detail(id)});
const sl=document.getElementById('sliders');
keys.forEach(k=>{const mx=Math.max(...ok.map(r=>r[k])),mn=Math.min(...ok.map(r=>r[k]));
 sl.insertAdjacentHTML('beforeend',`<label for="s_${k}">${D.labels[k]}: at most <b id="v_${k}">${mx.toFixed(2)}</b></label><input id="s_${k}" type="range" min="${mn}" max="${mx}" step="${(mx-mn)/200||0.01}" value="${mx}">`);
 sl.querySelector('#s_'+k).addEventListener('input',e=>{lim[k]=+e.target.value;document.getElementById('v_'+k).textContent=lim[k].toFixed(2);render()});});
const fcols=['runId','beamSet','letMode','letPenalty','oarScale','robust','coverageMode',...keys,
 D.target+' D95 nom',D.target+' D95 wc',D.target+' D2 nom',...D.letOars.flatMap(o=>[o+' D2',o+' dirty D2',o+' LETd p98']),'issues'].filter(c=>D.cols.includes(c));
function meets(r){return keys.every(k=>r[k]<=lim[k]+1e-9)}
function table(id,rows,cols){const t=document.getElementById(id);
 t.innerHTML='<tr>'+cols.map(c=>`<th>${c}</th>`).join('')+'</tr>'+rows.map(r=>`<tr data-id="${r.runId}" class="${r.status==='ok'&&meets(r)?'meets':''}">`+cols.map(c=>`<td>${r[c]===undefined?'':(c==='letMode'?MODE_NAME[r[c]]:r[c])}</td>`).join('')+'</tr>').join('');}
function render(){const fr=ok.filter(r=>front.has(r.runId)).sort((a,b)=>a.letBurden-b.letBurden);table('front',fr,fcols);
 const m=fr.filter(meets);const pick=document.getElementById('pick');
 if(!m.length){pick.innerHTML='<b>No Pareto plan meets all limits.</b><br>Relax one of the sliders.';return;}
 const best=m.reduce((a,b)=>a.letBurden<b.letBurden?a:b);
 pick.innerHTML=`<b>${m.length} Pareto plan(s) meet your limits.</b><br>Lowest OAR LET burden among them: <b>${best.runId}</b> (${MODE_NAME[best.letMode]}, beam set ${best.beamSet}`+
 (best.letPenalty!==''?`, penalty ${best.letPenalty}`:'')+`): coverage loss ${best.coverageLoss.toFixed(1)}%, LET burden ${best.letBurden.toFixed(2)}`+
 (keys.includes('robustLoss')?`, worst-case coverage loss ${best.robustLoss.toFixed(1)}%`:'')+`.<br><span style="color:var(--muted)">Spot weights: plans/${best.runId}_w.mat in the case folder.</span>`;}
function detail(id){const r=D.rows.find(x=>x.runId===id);const d=document.getElementById('detail');d.style.display='block';
 d.innerHTML=`<b>${id}</b><pre>${Object.entries(r).map(([k,v])=>k+': '+v).join('\n')}</pre>`;d.scrollIntoView({behavior:'smooth',block:'nearest'});}
['front','all'].forEach(id=>document.getElementById(id).addEventListener('click',e=>{const tr=e.target.closest('tr');if(tr&&tr.dataset.id)detail(tr.dataset.id)}));
table('all',D.rows,D.cols);render();
document.getElementById('log').innerHTML=D.log.map(e=>`<li><span class="rule">r${e.round} ${e.rule}</span> ${e.text.replace(/</g,'&lt;')}</li>`).join('');
</script></body></html>
"""
