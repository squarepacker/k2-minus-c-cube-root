# -*- coding: utf-8 -*-
"""tbx_summary.py -- aggregate the campaign records out/camp*.jsonl into the tables of the report.

usage: python tbx_summary.py "out/campA_*.jsonl" [more patterns ...]  > out/summary.txt
"""
import sys, json, glob, math
from collections import defaultdict


def fmax(a, b):
    if b is None:
        return a
    if a is None:
        return b
    return max(a, b)


def fmin(a, b):
    if b is None:
        return a
    if a is None:
        return b
    return min(a, b)


recs = []
hdr_only = []
for pat in sys.argv[1:]:
    for fn in sorted(glob.glob(pat)):
        for line in open(fn, encoding='utf-8'):
            line = line.strip()
            if not line:
                continue
            r = json.loads(line)
            if 'paths' in r:
                recs.append(r)
            else:
                hdr_only.append(r)

print('=' * 100)
print('records (packing x theta):', len(recs), '  header-only lines (errors / budget / unverified):', len(hdr_only))
for h in hdr_only[:20]:
    print('  HDR', {k: (v if k != 'error' else v[-400:]) for k, v in h.items()})
packs = set((r['fam'], r['k'], r['seed'], r['idx']) for r in recs)
print('distinct packings:', len(packs), ' all verified exactly:', all(r['verified'] for r in recs))
byfam = defaultdict(set)
for r in recs:
    byfam[r['fam']].add((r['k'], r['seed'], r['idx']))
print('packings per family:', {f: len(v) for f, v in sorted(byfam.items())})
print('N range', min(r['N'] for r in recs), max(r['N'] for r in recs), ' W range', min(r['W'] for r in recs), max(r['W'] for r in recs))
print('max time per record %.2f s, max rss %.1f MB' % (max(r['t_total'] for r in recs), max(r['rss'][1] for r in recs)))
anoms = [(r['fam'], r['idx'], r['anom']) for r in recs if r['anom'][0] or r['anom'][1]]
print('tracer anomalies:', len(anoms), anoms[:5])
viol = [(r['fam'], r['k'], r['seed'], r['idx'], r['tmode'], r['viol']) for r in recs if r['viol']]
print('records with violations:', len(viol))
for v in viol[:30]:
    print('  VIOL', v)

# ------------------------------------------------------------------ item 1 and Lemma 2
print('=' * 100)
print('ITEM 1 (Lemma 1) and Lemma 2, per theta mode')
modes = sorted(set(r['tmode'] for r in recs))
for mode in modes + ['ALL']:
    R = [r for r in recs if mode == 'ALL' or r['tmode'] == mode]
    agg = defaultdict(lambda: None)
    s = dict(flows=0, recs=0, passes=0, entries=0, halfA=0, noY2=0, ext=0, halfrhs=0, noD=0, nogap=0, halfDbound=0)
    for r in R:
        for p in r['paths']:
            s['flows'] += 1; s['recs'] += p['nrec']; s['passes'] += p['npass']; s['entries'] += p['nentry']
            s['halfA'] += p['L1_neg_halfA']; s['ext'] += p['L1_neg_ext_overlap']; s['halfrhs'] += p['L1_neg_capped']
            s['noD'] += p['L3_neg_noD']; s['nogap'] += p['L3_neg_nogap']
            s['noY2'] += p.get('L1_neg_noY2', 0); s['halfDbound'] += p.get('L3_neg_halfbound', 0)
            for key in ('L1_ratio_path_max', 'L1_ratio_glob_max', 'L1_ratio_draft_max', 'L2_ratio_max', 'L2W_ratio_max',
                        'L3_Dmax', 'L3_ratio_D_bound_max', 'L3_gapfrac_max', 'L3_distq_max', 'maxpass', 'L1_lhs_max'):
                agg[key] = fmax(agg[key], p[key])
            if p['L3_D_chain_bound'] and p['L3_D_chain_bound'] > 0:
                agg['D_over_chainbound'] = fmax(agg['D_over_chainbound'], p['L3_Dmax'] / p['L3_D_chain_bound'])
    print(' mode %-6s' % mode, s)
    print('   worst ratios:', {k: (round(v, 6) if isinstance(v, float) else v) for k, v in agg.items()})

# ------------------------------------------------------------------ item 2 by family (deficits)
print('=' * 100)
print('ITEM 2 (Lemma 3): max deficit D per family and theta mode; chain mode must have D <= tau1 = 0.15')
for fam in sorted(byfam):
    line = []
    for mode in modes:
        R = [r for r in recs if r['fam'] == fam and r['tmode'] == mode]
        if R:
            line.append('%s: maxD %.4f (D/chainbound %.3f)' % (mode, max(r['maxD'] for r in R),
                        max(max(p['L3_Dmax'] / p['L3_D_chain_bound'] if p['L3_D_chain_bound'] else 0 for p in r['paths']) for r in R)))
    print(' %-9s' % fam, ' | '.join(line))
chainR = [r for r in recs if r['tmode'] == 'chain']
if chainR:
    th = [r['sinT'] for r in chainR]
    print(' chain theta: sin theta in [%.4f, %.4f]; max D in chain mode %.5f' % (min(th), max(th), max(r['maxD'] for r in chainR)))

# ------------------------------------------------------------------ item 3
print('=' * 100)
print('ITEM 3 (Lemma 4) per tau set: type-Z squares and entries on them')
keys = sorted(set(k for r in recs for k in r['L4']))
for key in keys:
    st = dict(flows_cond_smallD=0, flows_cond_largeD=0, flows_nocond=0, nZ_sum=0, entries_smallD=0,
              entries_largeD=0, measure_largeD=0.0, minD_largeD=None, viol_impl=0, geo_min=None, packings_withZ=0)
    for r in recs:
        l4 = r['L4'][key]
        cond = l4['cond']; small = l4['maxD_le_tau1']
        n_ent = l4['floor']['n_entry_pieces'] + l4['ceil']['n_entry_pieces']
        st['nZ_sum'] += l4['floor']['nZ']
        if l4['floor']['nZ'] > 0:
            st['packings_withZ'] += 1
        st['geo_min'] = fmin(st['geo_min'], fmin(l4['floor']['geo_min_margin'], l4['ceil']['geo_min_margin']))
        if not cond:
            st['flows_nocond'] += 1
            continue
        if small:
            st['flows_cond_smallD'] += 1
            st['entries_smallD'] += n_ent
        else:
            st['flows_cond_largeD'] += 1
            st['entries_largeD'] += n_ent
            st['measure_largeD'] += l4['floor']['entry_measure'] + l4['ceil']['entry_measure']
            st['minD_largeD'] = fmin(st['minD_largeD'], fmin(l4['floor']['minD_on_Z'], l4['ceil']['minD_on_Z']))
        st['viol_impl'] += l4['floor']['n_entries_D_le_tau1'] + l4['ceil']['n_entries_D_le_tau1']
    print(' %-22s' % key, st)

# ------------------------------------------------------------------ item 4
print('=' * 100)
print('ITEM 4 (Lemma 6) per Z set')
zkeys = sorted(set(k for r in recs for k in r['L6'][0] if k.startswith('Z_')))
tot = dict(samples=0, exc_hits=0, tracer_bad=0)
for r in recs:
    for l6 in r['L6']:
        tot['samples'] += l6['nsamples']; tot['exc_hits'] += l6['exc_hits']; tot['tracer_bad'] += l6['tracer_bad']
print(' line samples:', tot)
for zk in zkeys:
    st = dict(flows=0, flows_nonempty=0, nZ_used_sum=0, nZ_passed_sum=0, ratio_global_max=None, ratio_int_max=None,
              ptw_ratio_max=None, ptw_margin_min=None, neg_with_passed=0, neg_LTonly=0)
    for r in recs:
        for l6 in r['L6']:
            z = l6[zk]
            st['flows'] += 1
            if z['nZ_used'] > 0:
                st['flows_nonempty'] += 1
            st['nZ_used_sum'] += z['nZ_used']; st['nZ_passed_sum'] += z['nZ_passed']
            st['ratio_global_max'] = fmax(st['ratio_global_max'], z['ratio_global'])
            st['ratio_int_max'] = fmax(st['ratio_int_max'], z['ratio_int'])
            st['ptw_ratio_max'] = fmax(st['ptw_ratio_max'], z['ptw_max_ratio'])
            st['ptw_margin_min'] = fmin(st['ptw_margin_min'], z['ptw_min_margin'])
            st['neg_with_passed'] += z['neg_with_passed_fail']; st['neg_LTonly'] += z['neg_LTonly_fail']
    print(' %-26s' % zk, st)
Lam = defaultdict(float)
for r in recs:
    for l6 in r['L6']:
        for t in ('LT', 'D', 'M', 'E', 'W', 'H', 'supOvgap'):
            Lam[t] += l6[t]
print(' summed loss measures over all flows:', {k: round(v, 3) for k, v in Lam.items()})
print(' measured A\' = theta (L_T + L~_T)/W : max %.4f' % max(r['Aprime_meas'] for r in recs))

# ------------------------------------------------------------------ item 5
print('=' * 100)
print('ITEM 5 (Theorem A with measured constants) per (tau0, tau1, beta, nu0)')
akeys = sorted(set(k for r in recs for k in r['TA']))
for ak in akeys:
    st = dict(n=0, cond=0, calH_min=None, I_max_ratio=None, II_max_ratio=None, K_pos=0, K_max=0.0, den_pos=0,
              ratio_K_max=None, Emax_on_K_max=None, L5_viol=0, TA_ratio_max=None, negII=0, K2_pos=0, K2_meas_max=0.0,
              K2_chain_ratio_max=None, K2_L6_max=None, ZK_passed_cond=0, ZK_passed_any=0, int_K_max=0.0)
    for r in recs:
        if ak not in r['TA']:
            continue          # campaign A records predate the last two tau sets
        a = r['TA'][ak]
        st['n'] += 1
        if a['cond'] and a['maxD_le_tau1']:
            st['cond'] += 1
            st['ZK_passed_cond'] += a['ZKb_passed'] + a['ZKt_passed'] + a['ZK2b_passed'] + a['ZK2t_passed']
        st['ZK_passed_any'] += a['ZKb_passed'] + a['ZKt_passed']
        st['calH_min'] = fmin(st['calH_min'], a['calH'] - a['calH_lower'])
        st['I_max_ratio'] = fmax(st['I_max_ratio'], a['ratio_I'])
        st['II_max_ratio'] = fmax(st['II_max_ratio'], a['ratio_II'])
        if a['K'] > 0:
            st['K_pos'] += 1; st['K_max'] = max(st['K_max'], a['K'])
            st['Emax_on_K_max'] = fmax(st['Emax_on_K_max'], a['Emax'])
            st['int_K_max'] = max(st['int_K_max'], a['intK'])
            if a['den'] > 0:
                st['den_pos'] += 1
        st['ratio_K_max'] = fmax(st['ratio_K_max'], a['ratio_K'])
        st['TA_ratio_max'] = fmax(st['TA_ratio_max'], a['TA_ratio'])
        st['negII'] += 1 if a['neg_II_halfA'] else 0
        if a['K2'] > 0:
            st['K2_pos'] += 1; st['K2_meas_max'] = max(st['K2_meas_max'], a['K2'])
        st['K2_chain_ratio_max'] = fmax(st['K2_chain_ratio_max'], a['K2_chain_ratio'])
        st['K2_L6_max'] = fmax(st['K2_L6_max'], fmax(a['K2_L6_ratio_b'], a['K2_L6_ratio_t']))
        if a['L5_int_K'] > a['L5_sum_int_c'] + 1e-12 or a['L5_sum_int_c'] > a['L5_sum_sinacosa'] + 1e-12:
            st['L5_viol'] += 1
    print(' %-28s' % ak, {k: (round(v, 6) if isinstance(v, float) else v) for k, v in st.items()})
