import numpy as np, re
def parse(fn):
    """Return {energy: (depth_mm, dose)} from TOPAS/RayStation commissioning csv."""
    out, E, cur, meta = {}, None, None, {}
    for l in open(fn, encoding='utf-8', errors='ignore'):
        l = l.strip()
        if not l: continue
        if l.startswith('Nominal beam energy'):
            if E is not None and cur: out[E] = np.array(cur).T
            E = float(l.split(';')[-1]); cur = []; continue
        m = re.match(r'^(-?[\d.eE+-]+);\s*(-?[\d.eE+-]+)$', l)
        if m and E is not None: cur.append([float(m.group(1)), float(m.group(2))])
        elif ':' in l and E is None: meta[l.split(':')[0]] = l.split(';')[-1].strip()
    if E is not None and cur: out[E] = np.array(cur).T
    return out, meta
