#!/usr/bin/env python3
"""
Parse Fortran .F90 files in pmd/ and extract subroutine/function call
relationships.  Outputs:
  pmd_structure.json -- nodes (routines) and links (call relations)
  diagram.html       -- interactive network diagram using vis.js
"""

import re
import json
import sys
from pathlib import Path
from collections import defaultdict

SCRIPT_DIR = Path(__file__).parent
PMD_DIR    = SCRIPT_DIR.parent
JSON_FILE  = SCRIPT_DIR / 'pmd_structure.json'
HTML_FILE  = SCRIPT_DIR / 'diagram.html'

# ---------------------------------------------------------------------------
# Regex patterns (matched against stripped, lowercased lines)
# ---------------------------------------------------------------------------
RE_MODULE_START  = re.compile(r'^\s*module\s+(\w+)\s*(?:!.*)?$', re.I)
RE_MODULE_END    = re.compile(r'^\s*end\s+module\b', re.I)
RE_PROGRAM_START = re.compile(r'^\s*program\s+(\w+)', re.I)
RE_PROGRAM_END   = re.compile(r'^\s*end\s+program\b', re.I)
RE_SUB_START     = re.compile(
    r'^\s*(?:(?:pure|elemental|recursive|impure)\s+)*subroutine\s+(\w+)', re.I)
RE_SUB_END       = re.compile(r'^\s*end\s+subroutine\b', re.I)
RE_FUNC_START    = re.compile(
    r'^\s*(?:(?:pure|elemental|recursive|impure|\w[\w\s(),:=*]*?)\s+)?'
    r'function\s+(\w+)\s*\(', re.I)
RE_FUNC_END      = re.compile(r'^\s*end\s+function\b', re.I)
RE_INTERFACE_S   = re.compile(r'^\s*(?:abstract\s+)?interface\b', re.I)
RE_INTERFACE_E   = re.compile(r'^\s*end\s+interface\b', re.I)
RE_CALL          = re.compile(r'^\s*call\s+(\w+)\s*[\s(]?', re.I)
RE_USE           = re.compile(r'^\s*use\s+(\w+)', re.I)
RE_COMMENT       = re.compile(r'^\s*!')


def strip_comment(line: str) -> str:
    """Remove everything from the first unquoted '!' to end of line."""
    in_str, sc = False, None
    for i, c in enumerate(line):
        if in_str:
            if c == sc:
                in_str = False
        else:
            if c in ('"', "'"):
                in_str, sc = True, c
            elif c == '!':
                return line[:i]
    return line


def join_continuations(raw_lines: list[str]) -> list[str]:
    """Join Fortran continuation lines (trailing '&')."""
    result, i = [], 0
    while i < len(raw_lines):
        line = raw_lines[i].rstrip('\n')
        body = strip_comment(line).rstrip()
        while body.endswith('&'):
            body = body[:-1].rstrip()
            i += 1
            if i < len(raw_lines):
                nxt = strip_comment(raw_lines[i].rstrip('\n')).strip()
                # leading & on continuation line
                if nxt.startswith('&'):
                    nxt = nxt[1:].lstrip()
                body = body + ' ' + nxt
            else:
                break
        result.append(body)
        i += 1
    return result


# ---------------------------------------------------------------------------
# Parser
# ---------------------------------------------------------------------------

def parse_files(pmd_dir: Path):
    """Return (nodes, calls, use_edges, name_to_ids)."""
    nodes      = []  # list of dicts
    calls      = []  # (caller_id, callee_name_lower)
    use_edges  = []  # (file_or_unit_id, module_name_lower)
    name_to_ids = defaultdict(list)  # name_lower -> [node_id, ...]

    for fpath in sorted(pmd_dir.glob('*.F90')):
        _parse_one(fpath.name, fpath, nodes, calls, use_edges, name_to_ids)

    return nodes, calls, use_edges, name_to_ids


def _make_id(fname, name, existing_ids):
    base = f"{fname}::{name}"
    if base not in existing_ids:
        return base
    n = 2
    while f"{base}_{n}" in existing_ids:
        n += 1
    return f"{base}_{n}"


def _parse_one(fname, fpath, nodes, calls, use_edges, name_to_ids):
    raw = fpath.read_text(encoding='utf-8', errors='replace').splitlines()
    lines = join_continuations(raw)

    existing_ids = {n['id'] for n in nodes}

    # Stack entries: {'type': str, 'name': str, 'id': str|None}
    stack          = []
    interface_depth = 0   # track interface blocks (skip definitions inside)

    def current_module():
        for ctx in reversed(stack):
            if ctx['type'] == 'module':
                return ctx['name']
        return None

    def current_callable_id():
        for ctx in reversed(stack):
            if ctx['type'] in ('program', 'subroutine', 'function'):
                return ctx['id']
        return None

    def push_node(name, ntype):
        mod = current_module()
        nid = _make_id(fname, name, existing_ids)
        existing_ids.add(nid)
        node = {'id': nid, 'name': name, 'type': ntype,
                'file': fname, 'module': mod}
        nodes.append(node)
        name_to_ids[name.lower()].append(nid)
        stack.append({'type': ntype, 'name': name, 'id': nid})
        return nid

    for line in lines:
        clean = strip_comment(line).strip()
        if not clean:
            continue

        # ---- interface blocks: skip definitions inside ----
        if RE_INTERFACE_S.match(clean):
            interface_depth += 1
            continue
        if RE_INTERFACE_E.match(clean):
            interface_depth = max(0, interface_depth - 1)
            continue
        if interface_depth > 0:
            continue

        # ---- end statements ----
        if RE_SUB_END.match(clean):
            if stack and stack[-1]['type'] == 'subroutine':
                stack.pop()
            continue
        if RE_FUNC_END.match(clean):
            if stack and stack[-1]['type'] == 'function':
                stack.pop()
            continue
        if RE_PROGRAM_END.match(clean):
            if stack and stack[-1]['type'] == 'program':
                stack.pop()
            continue
        if RE_MODULE_END.match(clean):
            if stack and stack[-1]['type'] == 'module':
                stack.pop()
            continue

        # ---- definitions ----
        m = RE_MODULE_START.match(clean)
        if m:
            mod_name = m.group(1)
            stack.append({'type': 'module', 'name': mod_name, 'id': None})
            continue

        m = RE_PROGRAM_START.match(clean)
        if m:
            push_node(m.group(1), 'program')
            continue

        m = RE_SUB_START.match(clean)
        if m:
            push_node(m.group(1), 'subroutine')
            continue

        # function: avoid matching 'end function' / 'use ..., only: func'
        lo = clean.lower()
        if 'function' in lo and not lo.startswith('end ') and not lo.startswith('use '):
            m = RE_FUNC_START.match(clean)
            if m:
                push_node(m.group(1), 'function')
                continue

        # ---- call statements ----
        m = RE_CALL.match(clean)
        if m:
            callee = m.group(1).lower()
            caller = current_callable_id()
            if caller:
                calls.append((caller, callee))
            continue

        # ---- use statements (module-level dependencies) ----
        m = RE_USE.match(clean)
        if m:
            mod_name = m.group(1).lower()
            caller = current_callable_id() or fname
            use_edges.append((caller, mod_name))

    return


# ---------------------------------------------------------------------------
# Graph builder
# ---------------------------------------------------------------------------

def build_graph(nodes, calls, name_to_ids, include_external=True,
                include_mpi=False):
    """
    Returns dict with 'nodes' and 'links'.
    include_external: add external (non-pmd) callee nodes
    include_mpi:      include mpi_* calls (very noisy; default False)
    """
    # program ノードは pmd のみ残す
    nodes = [n for n in nodes
             if not (n['type'] == 'program' and n['name'].lower() != 'pmd')]
    node_set = {n['id'] for n in nodes}
    ext_nodes = {}
    links = []
    seen = set()

    for caller_id, callee_name in calls:
        if caller_id not in node_set:
            continue
        targets = name_to_ids.get(callee_name, [])
        if targets:
            for tid in targets:
                key = (caller_id, tid)
                if key not in seen:
                    links.append({'source': caller_id, 'target': tid,
                                  'kind': 'call'})
                    seen.add(key)
        elif include_external:
            if callee_name.startswith('mpi_') and not include_mpi:
                continue
            ext_id = f"external::{callee_name}"
            if ext_id not in ext_nodes:
                etype = 'external_mpi' if callee_name.startswith('mpi_') \
                    else 'external'
                ext_nodes[ext_id] = {
                    'id': ext_id, 'name': callee_name,
                    'type': etype, 'file': None, 'module': None
                }
            key = (caller_id, ext_id)
            if key not in seen:
                links.append({'source': caller_id, 'target': ext_id,
                              'kind': 'call'})
                seen.add(key)

    all_nodes = nodes + list(ext_nodes.values())
    return {'nodes': all_nodes, 'links': links}


# ---------------------------------------------------------------------------
# HTML generation  (dark + warm-gold, CSS Grid, custom SVG)
# ---------------------------------------------------------------------------

HTML_TEMPLATE = '''\
<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="UTF-8">
<title>pmd Call Graph</title>
<style>
:root {
  --bg:#0c1117; --bg-2:#11171f; --panel:#161d27; --line:#2a3340;
  --text:#d6dde6; --muted:#8a94a3;
  --edge:#2c3441; --edge-active:#ffd166; --edge-active-glow:rgba(255,209,102,0.55);
  --anno-h:260px;
}
*{box-sizing:border-box;margin:0;padding:0}
body{background:var(--bg);color:var(--text);font-family:system-ui,sans-serif;overflow:hidden}
.layout{
  display:grid;
  grid-template-columns:320px 1fr;
  grid-template-rows:48px 1fr 6px var(--anno-h);
  grid-template-areas:"header header" "sidebar canvas" "sidebar resizer" "sidebar annotations";
  height:100vh;
}
header{
  grid-area:header;background:var(--bg-2);border-bottom:1px solid var(--line);
  display:flex;align-items:center;padding:0 20px;gap:14px;
}
header h1{font-size:14px;font-weight:700;letter-spacing:.03em;white-space:nowrap}
header .path{font-size:11px;color:var(--muted);font-family:monospace}
header .hint{margin-left:auto;font-size:10px;color:var(--muted);white-space:nowrap}
.sidebar{
  grid-area:sidebar;background:var(--panel);border-right:1px solid var(--line);
  padding:14px 12px;overflow-y:auto;display:flex;flex-direction:column;gap:18px;
}
.s-title{font-size:9px;font-weight:800;letter-spacing:.12em;text-transform:uppercase;
          color:var(--muted);margin-bottom:6px}
.legend{display:flex;flex-direction:column;gap:5px}
.li{display:flex;align-items:center;gap:8px;font-size:11px;padding:3px 6px;
    border-radius:4px;cursor:pointer;transition:background .15s}
.li:hover{background:rgba(255,255,255,.05)}
.swatch{width:22px;height:12px;border-radius:3px;border:1.5px solid;flex-shrink:0}
#stats{font-size:11px;color:var(--muted);line-height:1.9}
.btn{width:100%;background:transparent;border:1px solid var(--line);color:var(--text);
     padding:7px 12px;border-radius:5px;cursor:pointer;font-size:11px;
     transition:border-color .2s,color .2s;margin-top:auto}
.btn:hover{border-color:var(--edge-active);color:var(--edge-active)}
.canvas{grid-area:canvas;overflow:hidden;position:relative}
#diagram{width:100%;height:100%;display:block;cursor:grab}
#diagram.panning{cursor:grabbing}
.resizer{grid-area:resizer;background:var(--line);cursor:ns-resize;
          transition:background .2s}
.resizer:hover{background:var(--edge-active)}
.annotations{
  grid-area:annotations;background:var(--bg-2);border-top:1px solid var(--line);
  padding:12px 16px;overflow-y:auto;font-size:11px;
}
.a-title{font-size:9px;font-weight:800;letter-spacing:.1em;text-transform:uppercase;
          color:var(--muted);margin-bottom:8px}
.a-empty{color:var(--muted)}
.a-head{margin-bottom:10px}
.a-name{font-size:13px;font-weight:700}
.a-sub{color:var(--muted);margin-top:3px;font-size:10px}
.a-section{margin-top:10px}
.a-list{list-style:none;display:flex;flex-direction:column;gap:3px;margin-top:5px}
.a-item{
  display:flex;align-items:baseline;gap:8px;padding:5px 10px;border-radius:4px;
  cursor:pointer;border-left:2px solid transparent;transition:background .15s,border-color .15s;
}
.a-item:hover{background:rgba(255,209,102,.07);border-left-color:var(--edge-active)}
.a-item.focused{background:rgba(255,209,102,.13);border-left-color:var(--edge-active);
                 animation:anno-flash 1.2s ease-out}
.a-ft{font-family:monospace;font-size:10px;flex:1;color:var(--text)}
.a-cnt{font-size:10px;color:var(--muted)}
#tooltip{
  position:fixed;display:none;background:rgba(14,20,28,.97);
  border:1px solid var(--edge-active);border-radius:6px;padding:10px 14px;
  font-size:11px;max-width:250px;pointer-events:none;z-index:9999;
  backdrop-filter:blur(4px);line-height:1.75;
}
#tooltip .tt-name{font-size:12px;font-weight:700;margin-bottom:5px}
#tooltip .tt-row{color:var(--muted)}
#tooltip .tt-row span{color:var(--text)}
@keyframes anno-flash{0%{background:rgba(255,209,102,.38)}100%{background:rgba(255,209,102,.13)}}
@keyframes node-glow{
  0%,100%{filter:drop-shadow(0 0 5px rgba(255,209,102,.7))}
  50%{filter:drop-shadow(0 0 14px rgba(255,255,255,.85))}
}
@keyframes edge-pulse{
  0%,100%{stroke-width:2;stroke:#ffd166}
  50%{stroke-width:3.5;stroke:#fff}
}
.node-glow rect.body{animation:node-glow 1.1s ease-in-out infinite}
.edge-active-anim{animation:edge-pulse 1.1s ease-in-out infinite}
</style>
</head>
<body>
<div class="layout" id="layout">
  <header>
    <h1>pmd&#8202;&#x2215;&#8202;Call Graph</h1>
    <span class="path">pmd/*.F90 &rarr; file-level dependencies</span>
    <span class="hint">F=fit &nbsp; Esc=deselect &nbsp; scroll=zoom &nbsp; drag=pan</span>
  </header>
  <nav class="sidebar">
    <div><div class="s-title">Legend</div><div class="legend" id="legend"></div></div>
    <div><div class="s-title">Stats</div><div id="stats"></div></div>
    <button class="btn" id="btn-reset">&#x27F2; Reset view</button>
  </nav>
  <main class="canvas">
    <svg id="diagram" xmlns="http://www.w3.org/2000/svg">
      <defs id="svg-defs"></defs>
      <g id="col-layer"></g>
      <g id="edge-layer"></g>
      <g id="node-layer"></g>
      <g id="badge-layer"></g>
    </svg>
  </main>
  <div class="resizer" id="resizer"></div>
  <div class="annotations">
    <div class="a-title">Connections</div>
    <div id="anno"></div>
  </div>
</div>
<div id="tooltip"></div>

<script type="application/json" id="workflow-data">
__JSON_DATA__
</script>
<script>
/* =====================================================================
   Data loading
   ===================================================================== */
async function loadData(){
  try{const r=await fetch('./pmd_structure.json');if(!r.ok)throw 0;return r.json()}
  catch(_){return JSON.parse(document.getElementById('workflow-data').textContent)}
}

/* =====================================================================
   File classification & column definitions
   ===================================================================== */
const COLS=[
  {id:'program', label:'PROGRAMS',  x:0,    stroke:'#e96060', fill:'#1e0c0c'},
  {id:'core',    label:'CORE',      x:260,  stroke:'#4cc9f0', fill:'#0b1f28'},
  {id:'module',  label:'MODULES',   x:520,  stroke:'#a29bfe', fill:'#16132a'},
  {id:'force',   label:'FORCES',    x:780,  stroke:'#55efc4', fill:'#0a1e18'},
  {id:'utility', label:'UTILITIES', x:1040, stroke:'#fdcb6e', fill:'#1e1a0a'},
  {id:'external',label:'EXTERNAL',  x:1300, stroke:'#636e72', fill:'#141819'},
];
const COL_W=180, COL_H=54, ROW_GAP=18;
const CORE_FILES=new Set([
  'pmd_core.F90','read_input.F90','descriptor.F90','common_parallel.F90'
]);
// JSONにprogramノードを持つファイルだけがprogramカラムになる（固定リストは持たない）
function classifyFile(f){
  if(!f||f==='(external)')return'external';
  if(CORE_FILES.has(f))return'core';
  const lo=f.toLowerCase();
  if(lo.startsWith('force_'))return'force';
  if(lo.startsWith('mod_'))return'module';
  return'utility';
}

/* =====================================================================
   Build file-level graph from routine-level JSON
   ===================================================================== */
function buildFileGraph(raw){
  const nById=Object.fromEntries(raw.nodes.map(n=>[n.id,n]));
  const fmap={};
  for(const n of raw.nodes){
    const k=n.file||'(external)';
    if(!fmap[k]){
      fmap[k]={id:k,name:k.replace(/\\.F90$/i,''),file:k,col:null,
               sub:0,fn:0,prog:0,routines:[]};
    }
    fmap[k].routines.push(n.id);
    if(n.type==='subroutine')fmap[k].sub++;
    else if(n.type==='function')fmap[k].fn++;
    else if(n.type==='program')fmap[k].prog++;
  }
  // prog>0 のファイルだけ 'program' カラム（JSONに実際にprogramノードがある場合のみ）
  for(const fm of Object.values(fmap)){
    fm.col = fm.prog>0 ? 'program' : classifyFile(fm.file);
  }
  const eMap={};
  for(const lk of raw.links){
    const sn=nById[lk.source], tn=nById[lk.target];
    if(!sn||!tn)continue;
    const sf=sn.file||'(external)', tf=tn.file||'(external)';
    if(sf===tf)continue;
    const k=sf+'\\x00'+tf;
    eMap[k]=(eMap[k]||0)+1;
  }
  const nodes=Object.values(fmap);
  const links=Object.entries(eMap).map(([k,cnt])=>{
    const[src,tgt]=k.split('\\x00');
    return{source:src,target:tgt,count:cnt};
  });
  return{nodes,links};
}

/* =====================================================================
   Column layout  (vertically centred per column)
   ===================================================================== */
function computeLayout(fg){
  const groups={};
  for(const n of fg.nodes){
    if(!groups[n.col])groups[n.col]=[];
    groups[n.col].push(n);
  }
  for(const g of Object.values(groups))g.sort((a,b)=>a.name.localeCompare(b.name));
  let maxH=0;
  for(const g of Object.values(groups))
    maxH=Math.max(maxH,g.length*(COL_H+ROW_GAP)-ROW_GAP);
  const pos={};
  for(const c of COLS){
    const g=groups[c.id]||[];
    const h=g.length*(COL_H+ROW_GAP)-ROW_GAP;
    const startY=(maxH-h)/2;
    g.forEach((n,i)=>{ pos[n.id]={x:c.x, y:startY+i*(COL_H+ROW_GAP), col:c.id}; });
  }
  return{pos, maxH};
}

/* =====================================================================
   SVG helpers
   ===================================================================== */
const NS='http://www.w3.org/2000/svg';
function mk(tag,attrs,text){
  const e=document.createElementNS(NS,tag);
  for(const[k,v]of Object.entries(attrs||{}))e.setAttribute(k,v);
  if(text!==undefined)e.textContent=text;
  return e;
}

/* =====================================================================
   Render defs
   ===================================================================== */
function renderDefs(defs){
  function arrow(id,color,w=8,h=8){
    const m=mk('marker',{id,'markerWidth':w,'markerHeight':h,
      'refX':w-1,'refY':h/2,'orient':'auto'});
    m.appendChild(mk('path',{d:`M0,0 L0,${h} L${w},${h/2} z`,fill:color}));
    defs.appendChild(m);
  }
  arrow('arr','#2c3441');
  arrow('arr-act','#ffd166',10,10);

  const sh=mk('filter',{id:'sh',x:'-20%',y:'-20%',width:'140%',height:'140%'});
  sh.appendChild(mk('feDropShadow',{dx:'0',dy:'2',stdDeviation:'3','flood-color':'rgba(0,0,0,.65)'}));
  defs.appendChild(sh);

  const glow=mk('filter',{id:'glow',x:'-40%',y:'-40%',width:'180%',height:'180%'});
  const blur=mk('feGaussianBlur',{stdDeviation:'5',result:'b'});
  const flood=mk('feFlood',{'flood-color':'rgba(255,209,102,.65)',result:'c'});
  const comp=mk('feComposite',{in:'c',in2:'b',operator:'in',result:'g'});
  const merge=mk('feMerge');
  merge.appendChild(mk('feMergeNode',{in:'g'}));
  merge.appendChild(mk('feMergeNode',{in:'SourceGraphic'}));
  [blur,flood,comp,merge].forEach(e=>glow.appendChild(e));
  defs.appendChild(glow);
}

/* =====================================================================
   Render col-layer  (headers + dividers)
   ===================================================================== */
function renderCols(layer,fg,pos,offX,offY,svgH){
  const colBounds={};
  for(const[id,p]of Object.entries(pos)){
    const col=p.col;
    if(!colBounds[col])colBounds[col]={minY:Infinity,maxY:-Infinity};
    const ay=p.y+offY;
    colBounds[col].minY=Math.min(colBounds[col].minY,ay);
    colBounds[col].maxY=Math.max(colBounds[col].maxY,ay+COL_H);
  }
  for(const c of COLS){
    const b=colBounds[c.id]; if(!b)continue;
    const cx=c.x+offX+COL_W/2;
    layer.appendChild(mk('text',{
      x:cx, y:b.minY-26, 'text-anchor':'middle',
      'font-size':'9','letter-spacing':'.15em',
      fill:'#2e3d50','font-weight':'800'
    },c.label));
    if(c.x>0){
      const lx=c.x+offX-30;
      layer.appendChild(mk('line',{
        x1:lx,y1:20,x2:lx,y2:svgH-20,
        stroke:'#1a2535','stroke-width':'1','stroke-dasharray':'3,5'
      }));
    }
  }
}

/* =====================================================================
   Bezier midpoint at t=0.5
   ===================================================================== */
function bezMid(x1,y1,cx1,cy1,cx2,cy2,x2,y2){
  return{
    x:.125*x1+.375*cx1+.375*cx2+.125*x2,
    y:.125*y1+.375*cy1+.375*cy2+.125*y2
  };
}

/* =====================================================================
   Render edge-layer + badge-layer
   ===================================================================== */
function renderEdges(eLayer,bLayer,fg,pos,offX,offY){
  const pairSet=new Set(fg.links.map(l=>l.source+'\\x00'+l.target));
  const colOf=id=>{const c=COLS.find(c=>c.id===(pos[id]?.col));return c?COLS.indexOf(c):-1};

  for(const lk of fg.links){
    const sp=pos[lk.source], tp=pos[lk.target]; if(!sp||!tp)continue;
    const sx=sp.x+offX, sy=sp.y+offY, tx=tp.x+offX, ty=tp.y+offY;
    const ci=colOf(lk.source), cj=colOf(lk.target);
    const isBidi=pairSet.has(lk.target+'\\x00'+lk.source);
    const lane=isBidi?(lk.source<lk.target?-16:16):0;

    let x1,y1,cx1,cy1,cx2,cy2,x2,y2;
    if(ci===cj){
      // same column: vertical
      x1=sx+COL_W/2+lane; y1=sy+COL_H;
      x2=tx+COL_W/2+lane; y2=ty;
      const my=(y1+y2)/2;
      cx1=x1;cy1=my;cx2=x2;cy2=my;
    } else if(sx<tx){
      // forward (left→right)
      x1=sx+COL_W; y1=sy+COL_H/2+lane;
      x2=tx;       y2=ty+COL_H/2+lane;
      const mx=(x1+x2)/2;
      cx1=mx;cy1=y1;cx2=mx;cy2=y2;
    } else {
      // backward (right→left) — arc over top
      x1=sx+COL_W/2; y1=sy;
      x2=tx+COL_W/2; y2=ty;
      const topY=Math.min(y1,y2)-60+lane*1.5;
      cx1=x1;cy1=topY;cx2=x2;cy2=topY;
    }

    const d=`M${x1},${y1} C${cx1},${cy1} ${cx2},${cy2} ${x2},${y2}`;
    const path=mk('path',{
      d, stroke:'var(--edge)', 'stroke-width':'1.5', fill:'none',
      'marker-end':'url(#arr)',opacity:'0.8',
      'data-src':lk.source,'data-tgt':lk.target,class:'epath'
    });
    eLayer.appendChild(path);

    // Badge (pill) with call count
    const mid=bezMid(x1,y1,cx1,cy1,cx2,cy2,x2,y2);
    const label=String(lk.count);
    const pw=Math.max(18,label.length*6+10);
    const bg=mk('g',{transform:`translate(${mid.x},${mid.y})`,
                      class:'ebadge','data-src':lk.source,'data-tgt':lk.target});
    bg.appendChild(mk('rect',{x:-pw/2,y:-8,width:pw,height:16,rx:8,
      fill:'#111a24',stroke:'#2a3a4a','stroke-width':'1'}));
    bg.appendChild(mk('text',{
      'text-anchor':'middle','dominant-baseline':'central',
      'font-size':'9',fill:'#5a6a7a','font-family':'monospace'
    },label));
    bLayer.appendChild(bg);
  }
}

/* =====================================================================
   Render node-layer
   ===================================================================== */
function renderNodes(layer,fg,pos,offX,offY){
  const colMap=Object.fromEntries(COLS.map(c=>[c.id,c]));
  for(const n of fg.nodes){
    const p=pos[n.id]; if(!p)continue;
    const x=p.x+offX, y=p.y+offY;
    const c=colMap[n.col];

    const g=mk('g',{class:'fnode','data-id':n.id,
                    transform:`translate(${x},${y})`,cursor:'pointer'});

    const rect=mk('rect',{class:'body',width:COL_W,height:COL_H,rx:8,
      fill:c.fill,stroke:c.stroke,'stroke-width':'1.5',filter:'url(#sh)'});
    g.appendChild(rect);

    // Title
    const nm=n.name.length>22?n.name.slice(0,21)+'…':n.name;
    g.appendChild(mk('text',{x:11,y:21,'font-size':'11','font-weight':'700',
      fill:'#d6dde6'},nm));

    // Subtitle
    const parts=[];
    if(n.prog)parts.push(n.prog+'prog');
    if(n.sub)parts.push(n.sub+'sub');
    if(n.fn)parts.push(n.fn+'fn');
    g.appendChild(mk('text',{x:11,y:38,'font-size':'9.5',fill:'#6a7a8a'},
      parts.join('  ')||n.col));

    // Column colour dot (right edge)
    g.appendChild(mk('circle',{cx:COL_W-13,cy:COL_H/2,r:4,fill:c.stroke,opacity:'.65'}));
    layer.appendChild(g);
  }
}

/* =====================================================================
   Pan / zoom via SVG viewBox
   ===================================================================== */
let VB={x:0,y:0,w:0,h:0};
const svg=document.getElementById('diagram');

function setVB(){svg.setAttribute('viewBox',`${VB.x} ${VB.y} ${VB.w} ${VB.h}`)}

function fitView(svgW,svgH){
  const r=svg.getBoundingClientRect();
  const a=r.width/r.height, da=svgW/svgH;
  const P=60;
  if(a>da){VB.h=svgH+P*2;VB.w=VB.h*a}
  else{VB.w=svgW+P*2;VB.h=VB.w/a}
  VB.x=(svgW-VB.w)/2; VB.y=(svgH-VB.h)/2;
  setVB();
}

function setupPanZoom(){
  let drag=false,sx,sy,svb;
  svg.addEventListener('mousedown',e=>{
    if(e.target.closest('.fnode'))return;
    drag=true;sx=e.clientX;sy=e.clientY;svb={...VB};
    svg.classList.add('panning');
  });
  window.addEventListener('mousemove',e=>{
    if(!drag)return;
    const r=svg.getBoundingClientRect();
    VB.x=svb.x-(e.clientX-sx)*VB.w/r.width;
    VB.y=svb.y-(e.clientY-sy)*VB.h/r.height;
    setVB();
  });
  window.addEventListener('mouseup',()=>{drag=false;svg.classList.remove('panning')});
  svg.addEventListener('wheel',e=>{
    e.preventDefault();
    const r=svg.getBoundingClientRect();
    const mx=(e.clientX-r.left)/r.width, my=(e.clientY-r.top)/r.height;
    const f=e.deltaY>0?1.1:.909;
    const nw=VB.w*f, nh=VB.h*f;
    VB.x-=(nw-VB.w)*mx; VB.y-=(nh-VB.h)*my;
    VB.w=nw; VB.h=nh; setVB();
  },{passive:false});
}

/* =====================================================================
   Focus system
   ===================================================================== */
let focusId=null;

function setFocus(id,fg){
  focusId=id;
  const adj=new Set([id]);
  fg.links.forEach(l=>{if(l.source===id)adj.add(l.target);
                        if(l.target===id)adj.add(l.source);});

  document.querySelectorAll('.fnode').forEach(g=>{
    const nid=g.getAttribute('data-id');
    const dim=!adj.has(nid);
    g.style.opacity=dim?'0.15':'1';
    const rect=g.querySelector('rect.body');
    if(nid===id){
      rect.setAttribute('stroke-width','2.5');
      rect.setAttribute('filter','url(#glow)');
      g.classList.add('node-glow');
    } else {
      rect.setAttribute('stroke-width','1.5');
      rect.setAttribute('filter','url(#sh)');
      g.classList.remove('node-glow');
    }
  });

  document.querySelectorAll('.epath').forEach(p=>{
    const src=p.getAttribute('data-src'), tgt=p.getAttribute('data-tgt');
    const active=src===id||tgt===id;
    p.setAttribute('stroke',active?'var(--edge-active)':'var(--edge)');
    p.setAttribute('stroke-width',active?'2.5':'1.5');
    p.setAttribute('marker-end',active?'url(#arr-act)':'url(#arr)');
    p.setAttribute('opacity',active?'1':(adj.has(src)&&adj.has(tgt)?'0.7':'0.1'));
    if(active)p.classList.add('edge-active-anim');
    else p.classList.remove('edge-active-anim');
  });

  document.querySelectorAll('.ebadge').forEach(b=>{
    const src=b.getAttribute('data-src'), tgt=b.getAttribute('data-tgt');
    b.style.opacity=(src===id||tgt===id)?'1':'0.15';
  });

  updateAnno(id,fg);
}

function clearFocus(fg){
  focusId=null;
  document.querySelectorAll('.fnode').forEach(g=>{
    g.style.opacity='1'; g.classList.remove('node-glow');
    const r=g.querySelector('rect.body');
    r.setAttribute('stroke-width','1.5'); r.setAttribute('filter','url(#sh)');
  });
  document.querySelectorAll('.epath').forEach(p=>{
    p.setAttribute('stroke','var(--edge)');p.setAttribute('stroke-width','1.5');
    p.setAttribute('marker-end','url(#arr)');p.setAttribute('opacity','0.8');
    p.classList.remove('edge-active-anim');
  });
  document.querySelectorAll('.ebadge').forEach(b=>b.style.opacity='1');
  document.getElementById('anno').innerHTML=
    '<p class="a-empty">Click a node to see its connections.</p>';
}

/* =====================================================================
   Annotation panel
   ===================================================================== */
function updateAnno(id,fg){
  const n=fg.nodes.find(n=>n.id===id); if(!n)return;
  const nmap=Object.fromEntries(fg.nodes.map(n=>[n.id,n]));
  const callers=fg.links.filter(l=>l.target===id);
  const callees=fg.links.filter(l=>l.source===id);

  const parts=[];
  if(n.prog)parts.push(n.prog+' program');
  if(n.sub)parts.push(n.sub+' subroutine');
  if(n.fn)parts.push(n.fn+' function');

  let h=`<div class="a-head">
    <div class="a-name">${n.name}</div>
    <div class="a-sub">${parts.join(' &nbsp;·&nbsp; ')||n.col}</div>
  </div>`;

  function section(title,list,dir){
    if(!list.length)return'';
    let s=`<div class="a-section">
      <div class="a-title">${title} (${list.length})</div><ol class="a-list">`;
    for(const lk of list){
      const peer=nmap[dir==='from'?lk.source:lk.target];
      s+=`<li class="a-item" data-peer="${peer?.id||''}">
        <span class="a-ft">${peer?.name||'?'}</span>
        <span class="a-cnt">${lk.count}&nbsp;call${lk.count>1?'s':''}</span>
      </li>`;
    }
    return s+'</ol></div>';
  }
  h+=section('Called by',callers,'from');
  h+=section('Calls into',callees,'to');

  document.getElementById('anno').innerHTML=h;

  document.querySelectorAll('.a-item').forEach(li=>{
    li.addEventListener('click',()=>{
      li.classList.remove('focused');
      void li.offsetWidth; // reflow for re-trigger
      li.classList.add('focused');
      const pid=li.getAttribute('data-peer');
      if(pid)setFocus(pid,fg);
    });
  });
}

/* =====================================================================
   Tooltip
   ===================================================================== */
function showTip(n,fg,e){
  const tip=document.getElementById('tooltip');
  const inc=fg.links.filter(l=>l.target===n.id).length;
  const out=fg.links.filter(l=>l.source===n.id).length;
  tip.innerHTML=`<div class="tt-name">${n.name}</div>
    <div class="tt-row">type: <span>${n.col}</span></div>
    <div class="tt-row">routines: <span>${(n.sub+n.fn+n.prog)||0}</span></div>
    <div class="tt-row">callers: <span>${inc} file${inc!==1?'s':''}</span></div>
    <div class="tt-row">callees: <span>${out} file${out!==1?'s':''}</span></div>`;
  const pw=240,ph=100;
  let lx=e.clientX+16,ly=e.clientY+12;
  if(lx+pw>innerWidth)lx=e.clientX-pw-8;
  if(ly+ph>innerHeight)ly=e.clientY-ph-8;
  tip.style.left=lx+'px'; tip.style.top=ly+'px'; tip.style.display='block';
}

/* =====================================================================
   Resizer
   ===================================================================== */
function setupResizer(){
  const rs=document.getElementById('resizer');
  const lay=document.getElementById('layout');
  let drag=false,sy,sh;
  rs.addEventListener('mousedown',e=>{
    drag=true;sy=e.clientY;
    sh=parseFloat(getComputedStyle(lay).getPropertyValue('--anno-h'))||260;
  });
  window.addEventListener('mousemove',e=>{
    if(!drag)return;
    const nh=Math.max(80,Math.min(600,sh+(sy-e.clientY)));
    lay.style.setProperty('--anno-h',nh+'px');
    localStorage.setItem('pmd-anno-h',nh);
  });
  window.addEventListener('mouseup',()=>drag=false);
  const saved=localStorage.getItem('pmd-anno-h');
  if(saved)lay.style.setProperty('--anno-h',saved+'px');
}

/* =====================================================================
   Legend
   ===================================================================== */
function renderLegend(fg){
  const leg=document.getElementById('legend');
  const cnt={};
  for(const n of fg.nodes)cnt[n.col]=(cnt[n.col]||0)+1;
  for(const c of COLS){
    if(!cnt[c.id])continue;
    const d=document.createElement('div');
    d.className='li';
    d.innerHTML=`<div class="swatch" style="background:${c.fill};border-color:${c.stroke}"></div>
      <span>${c.label.toLowerCase()} (${cnt[c.id]})</span>`;
    leg.appendChild(d);
  }
}

/* =====================================================================
   Main
   ===================================================================== */
async function main(){
  const raw=await loadData();
  const fg=buildFileGraph(raw);
  const{pos,maxH}=computeLayout(fg);

  // Compute SVG dimensions
  let minX=Infinity,minY=Infinity,maxX=-Infinity,maxY=-Infinity;
  for(const p of Object.values(pos)){
    minX=Math.min(minX,p.x); minY=Math.min(minY,p.y);
    maxX=Math.max(maxX,p.x+COL_W); maxY=Math.max(maxY,p.y+COL_H);
  }
  const PAD=90;
  const svgW=maxX-minX+PAD*2, svgH=maxY-minY+PAD*2;
  const offX=-minX+PAD, offY=-minY+PAD;

  const defs=document.getElementById('svg-defs');
  const colL=document.getElementById('col-layer');
  const edgeL=document.getElementById('edge-layer');
  const nodeL=document.getElementById('node-layer');
  const badgeL=document.getElementById('badge-layer');

  renderDefs(defs);
  renderCols(colL,fg,pos,offX,offY,svgH);
  renderEdges(edgeL,badgeL,fg,pos,offX,offY);
  renderNodes(nodeL,fg,pos,offX,offY);

  // Initial viewBox
  VB={x:0,y:0,w:svgW,h:svgH};
  setVB();
  setTimeout(()=>fitView(svgW,svgH),80);

  setupPanZoom();
  setupResizer();
  renderLegend(fg);

  // Stats
  const rCnt=raw.nodes.filter(n=>!n.type.startsWith('external')).length;
  document.getElementById('stats').innerHTML=
    `<div>${fg.nodes.length} files</div>
     <div>${rCnt} routines</div>
     <div>${fg.links.length} inter-file edges</div>`;

  // Node events
  const tip=document.getElementById('tooltip');
  svg.addEventListener('click',e=>{
    const nd=e.target.closest('.fnode');
    if(!nd){clearFocus(fg);return}
    setFocus(nd.getAttribute('data-id'),fg);
  });
  svg.addEventListener('mousemove',e=>{
    const nd=e.target.closest('.fnode');
    if(nd&&focusId===null){
      const n=fg.nodes.find(n=>n.id===nd.getAttribute('data-id'));
      if(n)showTip(n,fg,e);
    } else tip.style.display='none';
  });
  svg.addEventListener('mouseleave',()=>tip.style.display='none');

  // Keyboard
  document.addEventListener('keydown',e=>{
    if(e.key==='Escape'){clearFocus(fg);}
    else if(e.key==='f'||e.key==='F'){fitView(svgW,svgH);}
    else if(e.key==='+'||e.key==='='){
      VB.x+=VB.w*.1;VB.y+=VB.h*.1;VB.w*=.8;VB.h*=.8;setVB();
    } else if(e.key==='-'){
      VB.x-=VB.w*.125;VB.y-=VB.h*.125;VB.w*=1.25;VB.h*=1.25;setVB();
    }
  });

  // Reset button
  document.getElementById('btn-reset').addEventListener('click',()=>{
    clearFocus(fg); fitView(svgW,svgH);
  });
}

main().catch(console.error);
</script>
</body>
</html>
'''


def generate_html(graph, html_file: Path):
    json_str = json.dumps(graph, indent=2, ensure_ascii=False)
    html = HTML_TEMPLATE.replace('__JSON_DATA__', json_str)
    html_file.write_text(html, encoding='utf-8')
    print(f"Written: {html_file}")


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    print(f"Parsing .F90 files in {PMD_DIR} ...")
    nodes, calls, use_edges, name_to_ids = parse_files(PMD_DIR)
    print(f"  {len(nodes)} routines  |  {len(calls)} call statements")

    graph = build_graph(nodes, calls, name_to_ids,
                        include_external=True, include_mpi=False)
    print(f"  Graph: {len(graph['nodes'])} nodes, {len(graph['links'])} links")

    JSON_FILE.write_text(
        json.dumps(graph, indent=2, ensure_ascii=False), encoding='utf-8')
    print(f"Written: {JSON_FILE}")

    generate_html(graph, HTML_FILE)
    print("Done.")


if __name__ == '__main__':
    main()
