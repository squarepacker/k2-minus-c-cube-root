# -*- coding: utf-8 -*-
"""summary of the cross-check files out/xcheck_*.json"""
import json, glob, sys
tot = dict(files=0, flows=0, points=0, agree=0, disagree=0, merges=0, merge_ok=0, merge_bad=0, quant_bad=0,
           degenerate=0, asmx_maxdiff=0.0, asmx_missing=0, anom=0, Mmeasure=0.0, Tmeasure=0.0, Hmeasure=0.0)
rows = []
for fn in sorted(glob.glob(sys.argv[1] if len(sys.argv) > 1 else 'out/xcheck_*.json')):
    L = json.load(open(fn))
    tot['files'] += 1
    for r in L:
        for side in ('floor', 'ceil'):
            s = r[side]
            pt = s['point']
            tot['flows'] += 1
            tot['points'] += pt['points']; tot['agree'] += pt['agree']; tot['disagree'] += len(pt['disagree'])
            tot['merges'] += pt['merges']; tot['merge_ok'] += pt['merge_ok']; tot['merge_bad'] += len(pt['merge_bad'])
            tot['quant_bad'] += len(pt['quant_bad']); tot['degenerate'] += pt['degenerate']
            tot['anom'] += len(s['anom'])
            if s['asmx_maxdiff'] is None:
                tot['asmx_missing'] += 1
            else:
                tot['asmx_maxdiff'] = max(tot['asmx_maxdiff'], s['asmx_maxdiff'])
            m = s['measures']
            tot['Mmeasure'] += m.get('M', 0.0); tot['Tmeasure'] += m.get('T', 0.0); tot['Hmeasure'] += m.get('H', 0.0)
            rows.append((fn.split('xcheck_')[-1], r['idx'], side, s['nrec'], pt['points'], pt['agree'],
                         len(pt['disagree']), pt['merges'], pt['merge_ok'], len(pt['quant_bad']), s['asmx_maxdiff'],
                         {kk: round(v, 4) for kk, v in m.items()}))
for row in rows:
    print(row)
print(tot)
