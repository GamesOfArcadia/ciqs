"""
From the published paper:
"Page_Curve_Tested_on_Hardware_Entropy_rises_Correlations_do_not_v4"
Jessy Pensédent - July 12, 2026
---
Get CIQS and CIQA on GitHub: https://github.com/GamesOfArcadia/ciqs
---

Base setup:
IBM Heron r2 (156 qubits) | N_BH=6 (Page time t=3) | SCRAMBLER_DEPTH = 4

54 / 156 physical qubits:
  q[0..29]  BH data     (6 blocks x 5, CIQA n=5)
  q[30..41] BH syndrome (6 blocks x 2)
  q[42..47] Reference R (6 qubits, unprotected, never scrambled)
  q[48..53] Radiation D (up to 6 qubits)

24 classical bits per shot:
  cr[0..11]  syndrome
  cr[12..17] R
  cr[18..23] D

300 circuits x 1024 shots = 307,200 total shots.
"""

import math, sys, time, json
from collections import defaultdict, deque

import numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt

try:
    from CIQA import CIQA
except ImportError:
    print("ERROR: CIQA.py not found (downloadable for free though)"); sys.exit(1)
try:
    from CIQS_pipeline import CIQS_pipeline  # noqa: F401
    from CIQM import CouplingMap
except ImportError:
    print("ERROR: CIQS_pipeline.py / CIQM.py not found (also downloadable for free)"); sys.exit(1)
try:
    from qiskit import QuantumCircuit, QuantumRegister, ClassicalRegister
    from qiskit.compiler import transpile as qtranspile
    from qiskit_ibm_runtime import QiskitRuntimeService, SamplerV2 as IBMSampler
    from qiskit.transpiler.preset_passmanagers import generate_preset_pass_manager
except ImportError:
    print("ERROR: pip install qiskit qiskit-ibm-runtime"); sys.exit(1)


# config

# USE_SIMULATOR = False (submit tp QPU)
USE_SIMULATOR  = False

# Set an IBM API token to submit to IBM QPU. Setup another method to submit to another platform.
IBM_API_TOKEN  = ""

# Set the backend name to display in any log and output.
BACKEND_NAME   = ""

# Circuit size (default 6 logical qubits - Page time at 3 | higher count isn't consistantly reliable, it will depend on the hardware you're using)
N_BH           = 6

N_SHADOWS      = 50
SHOTS_PER_CIRC = 1024

# Set a job ID to skip execution and replay post-processing only without QPU usage. (USE_SIMULATOR must be false)
# Set job ID back to "" to run a fresh job.
REPLAY_JOB_ID  = ""

# Base depth = 4
SCRAMBLER_DEPTH = 4

SCRAMBLER_SEED  = 42
SHADOW_SEED     = 1337

BLK           = 5
SYN_PER_BLOCK = 2
CLEAN_RATE    = 0.942

# Simulation variables
SIM_N_BH       = 4
SIM_N_SHADOWS  = 10
SIM_SHOTS      = 512



# index functions
# N_BH=6: anchors at 0,5,10,15,20,25 | syn at 30-41 | R at 42-47 | D at 48-53
# cr: syn 0-11 | R 12-17 | D 18-23

def bh_anchor(i):       return BLK * i
def bh_syn1(i, n):      return BLK*n + SYN_PER_BLOCK*i
def bh_syn2(i, n):      return BLK*n + SYN_PER_BLOCK*i + 1
def ref_q(i, n):        return (BLK+SYN_PER_BLOCK)*n + i
def rad_q(r, n):        return (BLK+SYN_PER_BLOCK)*n + n + r
def cr_syn1(i, n):           return 2*i
def cr_syn2(i, n):           return 2*i+1
def cr_ref(i, n, syn=True):  return (2*n if syn else 0) + i
def cr_rad(r, n, syn=True):  return (3*n if syn else n) + r


# random parameters

def make_scrambler(n, depth, seed):
    rng = np.random.default_rng(seed)
    return {l: {q: float(rng.uniform(0, 2*math.pi)) for q in range(n)}
            for l in range(depth)}

def make_shadows(n_u, n_bh, seed):
    """Shape (n_u, 2*n_bh, 3). Cols 0..n_bh-1: D. Cols n_bh..2*n_bh-1: R."""
    rng = np.random.default_rng(seed)
    n = 2*n_bh
    a = rng.uniform(0, 2*math.pi, (n_u, n))
    u = rng.uniform(0, 1, (n_u, n))
    b = 2*np.arccos(np.sqrt(np.clip(1-u, 0, 1)))
    g = rng.uniform(0, 2*math.pi, (n_u, n))
    return np.stack([a, b, g], axis=-1)


# CIQA encoding

def encode_bh(n_bh):
    ciqa = CIQA(n=5, d=2)
    enc  = ciqa.encode([{"type":"Ry","qudit":i,"theta":0.0,"d":2}
                         for i in range(n_bh)])
    n_p  = enc["n_physical"]
    blk  = n_p // n_bh
    assert blk == BLK
    ps   = CLEAN_RATE**n_bh
    print(f"  CIQA: {n_bh} logical -> {n_p} data + {SYN_PER_BLOCK*n_bh} syn"
          f" + {n_bh} R + {n_bh} D = {n_p+SYN_PER_BLOCK*n_bh+2*n_bh} qubits"
          f" ({n_p+SYN_PER_BLOCK*n_bh+2*n_bh}/156)")
    print(f"  Post-sel yield: {CLEAN_RATE:.3f}^{n_bh} = {ps:.3f}"
          f"  (~{int(ps*SHOTS_PER_CIRC)} clean shots/{SHOTS_PER_CIRC})")
    return enc["gates"], n_p, blk, ciqa


# circuit

def add_scrambler(qc, n_bh, params, depth):
    anc = [bh_anchor(i) for i in range(n_bh)]
    for l in range(depth):
        for qi, ph in enumerate(anc): qc.ry(params[l][qi], ph)
        for i in range(l%2, n_bh-1, 2): qc.cx(anc[i], anc[i+1])

def add_syndrome(qc, ciqa_obj, n_bh, n_phys_bh, cr):
    """
    CIQA syndrome round. Ancilla layout:
    block i -> bh[n_phys_bh + 2*i] and bh[n_phys_bh + 2*i + 1]
    maps to cr[0..2*n_bh-1].
    """
    dummy  = {"gates": [], "n_logical": n_bh,
              "n_physical": n_phys_bh, "d": 2, "report": None}
    result = ciqa_obj.syndrome(dummy)
    for g in result["syndrome_gates"]:
        gt = g["type"]
        if   gt == "Ry": qc.ry(g["theta"], g["qudit"])
        elif gt == "CX": qc.cx(g["control"], g["target"])
        elif gt == "M":
            qc.measure(g["qudit"], cr[g["qudit"] - n_phys_bh])

def build_circuit(t, sidx, bh_gates, n_phys_bh, blk, sp, su, ciqa_obj,
                   n_bh=N_BH, depth=SCRAMBLER_DEPTH, syn=True):
    nd  = blk*n_bh
    ns  = SYN_PER_BLOCK*n_bh if syn else 0
    ncr = (2*n_bh if syn else 0) + n_bh + t

    bh  = QuantumRegister(nd+ns,  "bh")
    ref = QuantumRegister(n_bh,   "ref")
    rad = QuantumRegister(t,      "r")
    cr  = ClassicalRegister(ncr,  "c")
    qc  = QuantumCircuit(bh, ref, rad, cr, name=f"hp_t{t}_s{sidx}")

    # CIQA BH encoding prefix
    for g in bh_gates:
        if g["type"]=="Ry":  qc.ry(g["theta"], bh[g["qudit"]])
        elif g["type"]=="CX": qc.cx(bh[g["control"]], bh[g["target"]])
    qc.barrier()

    # Bell init: H(R[i]), CX(R[i] -> anchor(i))
    for i in range(n_bh):
        qc.h(ref[i])
        qc.cx(ref[i], bh[bh_anchor(i)])
    qc.barrier()

    # Evaporation
    for s in range(t):
        add_scrambler(qc, n_bh, sp, depth); qc.barrier()
        qc.cx(bh[bh_anchor(s%n_bh)], rad[s]); qc.barrier()
        if syn: add_syndrome(qc, ciqa_obj, n_bh, blk*n_bh, cr); qc.barrier()

    # Shadow unitaries on D and R
    u = sidx % su.shape[0]
    for r in range(t):
        a,b,g = su[u, r];       qc.rz(float(a),rad[r]); qc.ry(float(b),rad[r]); qc.rz(float(g),rad[r])
    for i in range(n_bh):
        a,b,g = su[u, n_bh+i]; qc.rz(float(a),ref[i]); qc.ry(float(b),ref[i]); qc.rz(float(g),ref[i])
    qc.barrier()

    for i in range(n_bh):   qc.measure(ref[i],    cr[cr_ref(i,n_bh,syn)])
    for r in range(t):      qc.measure(rad[r],     cr[cr_rad(r,n_bh,syn)])
    return qc


# layout fix: place R[i] adjacent to anchor(i)

def fix_layout_for_R(pinned, n_bh, cm_edges):
    """
    Reassign R[i] qubits to physical qubits adjacent to anchor(i).
    Reduces Bell-init SWAPs.
    """
    adj = defaultdict(set)
    for (a,b) in cm_edges:
        adj[a].add(b); adj[b].add(a)

    n_bh_total = (BLK+SYN_PER_BLOCK)*n_bh  # 42
    bh_phys    = {pinned[l] for l in range(n_bh_total)}

    layout  = {l: pinned[l] for l in range(n_bh_total)}
    in_use  = set(bh_phys)

    for i in range(n_bh):
        anc_phys = layout[bh_anchor(i)]
        free     = sorted(adj[anc_phys] - in_use)
        if free:
            r_phys = free[0]
        else:
            # BFS to nearest available
            q, seen, r_phys = deque([anc_phys]), {anc_phys}, None
            while q and r_phys is None:
                cur = q.popleft()
                for nb in sorted(adj[cur]):
                    if nb not in seen:
                        seen.add(nb)
                        if nb not in in_use: r_phys = nb; break
                        q.append(nb)
            if r_phys is None:
                r_phys = next(p for p in range(156) if p not in in_use)
        layout[n_bh_total + i] = r_phys
        in_use.add(r_phys)

    rad_start = n_bh_total + n_bh  # 48
    avail = sorted(p for p in range(156) if p not in in_use)
    for r in range(n_bh):
        layout[rad_start + r] = avail[r]
        in_use.add(avail[r])

    return [layout[i] for i in range(len(layout))]


# backend helpers

def get_ibm_cm(backend):
    n     = backend.configuration().n_qubits
    edges = {tuple(sorted(e)) for e in backend.coupling_map.get_edges()}
    return CouplingMap(n_qubits=n, edges=edges, topology=BACKEND_NAME), list(edges)

def get_sim_cm(n): return CouplingMap.line(n), [(i,i+1) for i in range(n-1)]

def read_cal(backend):
    props = backend.properties()
    n     = backend.configuration().n_qubits
    cx_err, ro_err, t1_us = {}, {}, {}
    for q in range(n):
        try:    ro_err[q] = props.readout_error(q)
        except: ro_err[q] = 0.05
        try:    t1_us[q]  = props.t1(q)*1e6
        except: t1_us[q]  = 0.0
    for a,b in backend.coupling_map.get_edges():
        err = None
        for gn in ("cx","ecr","rzx"):
            try: err = props.gate_error(gn,[a,b]);
            except: pass
            if err is not None: break
        if err is None: err = 0.01
        k = (min(a,b),max(a,b))
        if k not in cx_err or err<cx_err[k]: cx_err[k]=err
    return cx_err, ro_err, t1_us

def ciqs_layout(circuit, cm, verbose=True):
    if verbose:
        print(f"    {circuit.num_qubits}q  depth={circuit.depth()}  gates={circuit.size()}")
        print("    routing ...")
    r = CIQS_pipeline(circuit, cm)
    if verbose:
        print(f"    valid={r.placement.valid}  SWAPs={r.swap_count}  removed={r.gates_removed}")
    return r.final_placement

def get_counts(pub):
    for name in ("c","cr","meas"):
        try: return getattr(pub.data, name).get_counts()
        except AttributeError: pass
    d = pub.data
    return getattr(d, list(vars(d).keys())[0]).get_counts()


# count extraction

def extract_counts(raw, n_bh, t, syn, post_sel=True):
    """
    Returns:
      D_counts  : {d_bits: count}   (t chars, D[0] first)
      R_counts  : {r_bits: count}   (n_bh chars, R[0] first)
      pairs     : {(r,m): {2-bit-str: count}}  matched pairs
    """
    sw      = 2*n_bh if syn else 0
    n_total = sum(raw.values())
    D_c  = {}
    R_c  = {}
    pc   = {(r, r%n_bh): {} for r in range(t)}
    n_a  = 0

    for bs, cnt in raw.items():
        bits = bs[::-1]
        need = sw + n_bh + t
        if len(bits) < need: bits += '0'*(need-len(bits))

        if post_sel and syn:
            if any(bits[k]!='0' for k in range(sw)): continue

        r_b = bits[sw:sw+n_bh]
        d_b = bits[sw+n_bh:sw+n_bh+t]

        D_c[d_b] = D_c.get(d_b,0) + cnt
        R_c[r_b] = R_c.get(r_b,0) + cnt

        for r in range(t):
            m   = r % n_bh
            key = (d_b[r] if r<len(d_b) else '0') + (r_b[m] if m<len(r_b) else '0')
            pc[(r,m)][key] = pc[(r,m)].get(key,0) + cnt

        n_a += cnt

    return D_c, R_c, pc, n_a, n_total


# shadow estimator

def renyi2(shadow_list, n_qubits, n_shadows):
    if n_qubits==0: return 0.0, 1.0
    D = 1<<n_qubits
    F = np.zeros((n_shadows,D))
    for i,c in enumerate(shadow_list):
        tot = sum(c.values())
        if tot==0: continue
        for bs,cnt in c.items():
            F[i, int(bs.zfill(n_qubits)[-n_qubits:],2)] += cnt/tot
    ia = np.arange(D,dtype=np.int32)
    hm = np.zeros((D,D))
    xm = ia[:,None]^ia[None,:]
    for k in range(n_qubits): hm += ((xm>>k)&1).astype(float)
    W  = (-2.0)**(-hm)
    K  = (F@W)@F.T; np.fill_diagonal(K,0)
    np_ = n_shadows*(n_shadows-1)
    if np_==0: return 0.0,1.0
    p = float(np.clip(D*K.sum()/np_,1e-10,1.0))
    return float(-math.log2(p)), p


# matched-pair I2

def compute_I2_pairs(shadow_D_list, shadow_R_list, shadow_pairs_list,
                      t, n_bh, n_shadows):
    """
    For each radiation qubit r, compute I2(D_r : R_{r%n_bh}).
    shadow_pairs_list: list of n_shadows dicts {(r,m): {2-bit: count}}

    Returns list of dicts, one per pair, and mean I2.
    """
    results = []
    for r in range(t):
        m = r % n_bh

        shA = [{k[0]: v for k2, v2 in sh.items()
                for k,v in {k2: v2}.items()
                if isinstance(k2,str) and k2[0]==k[0]}
               for sh in shadow_D_list]
        # simpler: marginalise D_r from D counts
        shA = []
        for c in shadow_D_list:
            mc = {}
            for bs,cnt in c.items():
                bit = bs[r] if r<len(bs) else '0'
                mc[bit] = mc.get(bit,0)+cnt
            shA.append(mc)

        shB = []
        for c in shadow_R_list:
            mc = {}
            for bs,cnt in c.items():
                bit = bs[m] if m<len(bs) else '0'
                mc[bit] = mc.get(bit,0)+cnt
            shB.append(mc)

        shAB = [sh.get((r,m),{}) for sh in shadow_pairs_list]

        s2A,_  = renyi2(shA,  1, n_shadows)
        s2B,_  = renyi2(shB,  1, n_shadows)
        s2AB,_ = renyi2(shAB, 2, n_shadows)
        I2     = s2A + s2B - s2AB

        results.append({"r": r, "m": m, "S2_D": s2A, "S2_R": s2B,
                         "S2_DR": s2AB, "I2": I2})

    mean_I2 = float(np.mean([x["I2"] for x in results])) if results else 0.0
    return results, mean_I2


# print results

def print_results(t_vals, S2_D, S2_R, I2_by_t, mean_I2_by_t, n_bh):
    pt = n_bh//2
    sp = [min(t,n_bh-t)*math.log2(2) for t in t_vals]
    print("\n"+"="*72)
    print(f"  RESULTS  (n_bh={n_bh}, Page time t={pt})")
    print("="*72)
    print(f"  {'t':<4} {'S2(D)':>8} {'page':>8} {'S2(R)':>8} {'mean I2(D:R)':>14}  note")
    print("  "+"-"*60)
    for i,t in enumerate(t_vals):
        m  = ">" if t>pt else " "
        I2 = mean_I2_by_t[i]
        note = "above Page time" if t>pt else ""
        print(f"  {t:<4}{m} {S2_D[i]:>8.4f} {sp[i]:>8.4f} "
              f"{S2_R[i]:>8.4f} {I2:>14.4f}  {note}")
    pre  = [mean_I2_by_t[i] for i,t in enumerate(t_vals) if t<=pt]
    post = [mean_I2_by_t[i] for i,t in enumerate(t_vals) if t> pt]
    print()
    if pre:  print(f"  Mean I2(D:R) t<={pt}: {np.mean(pre):+.4f}")
    if post: print(f"  Mean I2(D:R) t> {pt}: {np.mean(post):+.4f}")
    if pre and post:
        d = np.mean(post)-np.mean(pre)
        print(f"  Delta (post - pre):  {d:+.4f}")
        print()
        if   d > 0.05:  print("  I2 rises after Page time.")
        elif abs(d)<0.05: print("  I2 stays flat.")
        else:             print("  I2 falls after Page time.")
    print("="*72)


# plot

def plot_results(t_vals, S2_D, S2_R, I2_by_t, mean_I2_by_t, n_bh, job_id):
    pt   = n_bh//2
    sp   = [min(t,n_bh-t)*math.log2(2) for t in t_vals]
    ta   = np.array(t_vals,float)
    tf   = np.linspace(0,n_bh,600)
    pf   = np.clip(np.minimum(tf,n_bh-tf)*np.log2(2),0,None)

    BG,PN,GR = '#0d1117','#161b22','#21262d'
    TC,DIM   = '#e6edf3','#8b949e'
    RED,BLUE,AMB,ORG = '#f85149','#58a6ff','#e3b341','#f78166'

    fig,axes = plt.subplots(1,3,figsize=(17,5.5),facecolor=BG)
    for ax in axes:
        ax.set_facecolor(PN)
        ax.tick_params(colors=TC,labelsize=9)
        ax.xaxis.label.set_color(TC); ax.yaxis.label.set_color(TC)
        ax.title.set_color(TC)
        for sp_ in ax.spines.values(): sp_.set_edgecolor(GR)
        ax.grid(True,color=GR,lw=0.6,alpha=0.7)
        ax.set_xlim(-0.3,n_bh+0.3)
        ax.set_xlabel('Evaporation step  t')
        ax.axvline(pt,color=AMB,lw=1.5,ls=':',alpha=0.8)

    # Panel 1: S2(D)
    ax = axes[0]
    ax.plot(tf,pf,'--',color=RED,lw=2,label='Page prediction')
    ax.plot(np.concatenate([[0],ta]),np.concatenate([[0.],S2_D]),
            'o-',color=BLUE,lw=2.5,ms=9,
            markeredgecolor='white',markeredgewidth=1.5,label='Measured S2(D)')
    ax.set_ylabel('S2  (bits)')
    ax.set_title('Radiation entropy  S2(D)')
    ax.legend(fontsize=8,facecolor=GR,labelcolor=TC)

    # Panel 2: S2(R)
    ax = axes[1]
    ax.axhline(n_bh,color='#3fb950',lw=1.5,ls='--',alpha=0.6,
               label=f'Max = {n_bh} bits')
    ax.plot(ta,S2_R,'s-',color='#a5d6ff',lw=2,ms=8,
            markeredgecolor='white',markeredgewidth=1.2,label='Measured S2(R)')
    ax.set_ylabel('S2  (bits)')
    ax.set_title('Reference entropy  S2(R)')
    ax.legend(fontsize=8,facecolor=GR,labelcolor=TC)

    # Panel 3: matched-pair I2
    ax = axes[2]
    ax.axhline(0,color=DIM,lw=1.5,ls='--',label='I2 = 0')
    for t_step,pairs in I2_by_t.items():
        for pr in pairs:
            c = ORG if t_step>pt else BLUE
            ax.scatter(t_step,pr["I2"],color=c,s=70,zorder=5,
                       edgecolors='white',linewidths=0.8,alpha=0.85)
    ax.plot(ta,mean_I2_by_t,'o-',color='white',lw=2,ms=7,
            markeredgecolor='white',label='Mean I2 per step',zorder=6)
    ax.scatter([],[],color=BLUE,s=60,label='t ≤ Page time')
    ax.scatter([],[],color=ORG, s=60,label='t > Page time')
    ax.set_ylabel('I2(D_r : R_matched)  (bits)')
    ax.set_title('Mutual information  I2(D:R)  - matched pairs')
    ax.legend(fontsize=8,facecolor=GR,labelcolor=TC,loc='upper left')

    fig.suptitle(
        f'Hayden-Preskill Experiment - IBM Heron r2 - Job: {job_id}',
        color=TC,fontsize=11,y=1.01)
    plt.tight_layout()
    fname = f"hp_v5_{job_id}.png"
    plt.savefig(fname,dpi=180,bbox_inches='tight',facecolor=BG)
    plt.close()
    return fname


# main

def run():
    if USE_SIMULATOR:
        n_bh,n_sh,shots = SIM_N_BH,SIM_N_SHADOWS,SIM_SHOTS
        mode = "SIMULATOR"
    else:
        n_bh,n_sh,shots = N_BH,N_SHADOWS,SHOTS_PER_CIRC
        mode = f"HARDWARE ({BACKEND_NAME})"

    pt    = n_bh//2
    nc    = n_bh*n_sh

    print("\n"+"="*70)
    if REPLAY_JOB_ID:
        print(f"  HP v5 - REPLAY  job={REPLAY_JOB_ID}")
    else:
        print(f"  HP v5 - {mode}")
    print(f"  N_BH={n_bh}  Page_time={pt}  Circuits={nc}  Shots={nc*shots:,}")
    print("="*70)

    # Parameters and encoding are always needed (deterministic, same seeds)
    print("\n[1/8] Random parameters ...")
    sp = make_scrambler(n_bh, SCRAMBLER_DEPTH, SCRAMBLER_SEED)
    su = make_shadows(n_sh, n_bh, SHADOW_SEED)
    print(f"  Scrambler: {n_bh}q x {SCRAMBLER_DEPTH} layers  seed={SCRAMBLER_SEED}")
    print(f"  Shadows:   shape {su.shape}  {su.nbytes//1024} KB")

    print("\n[2/8] CIQA n=5 encoding ...")
    bh_gates, n_phys_bh, blk, ciqa_obj = encode_bh(n_bh)

    # Build circuit_index - same order as original, no transpilation needed for replay
    cidx = {(t_step, sidx): t_step*n_sh - n_sh + sidx
            for t_step in range(1, n_bh+1)
            for sidx   in range(n_sh)}

    if REPLAY_JOB_ID:
        # Reconnect and fetch completed job - skip steps 3-7
        print(f"\n[3-7/8] Fletching completed job {REPLAY_JOB_ID} ...")
        svc    = QiskitRuntimeService(channel="ibm_quantum_platform",
                                      token=IBM_API_TOKEN,
                                      instance="open-instance")
        job    = svc.job(REPLAY_JOB_ID)
        result = job.result()
        job_id = REPLAY_JOB_ID
        print(f"  Retrieved {len(result)} pub results.")

    else:
        if USE_SIMULATOR:
            print("\n[3/8] Simulation mode.")
            backend, cm_edges = None, []
            n_tot = n_phys_bh + SYN_PER_BLOCK*n_bh + 2*n_bh
            cm,_ = get_sim_cm(n_tot)
            print("\n[4/8] Mock calibration.")
        else:
            print("\n[3/8] Connecting to IBM Quantum ...")
            svc     = QiskitRuntimeService(channel="ibm_quantum_platform",
                                           token=IBM_API_TOKEN)
            backend = svc.backend(BACKEND_NAME)
            print(f"  {backend.name}  ({backend.configuration().n_qubits} qubits)")
            print("\n[4/8] Calibration ...")
            cx_err, ro_err, t1_us = read_cal(backend)
            print(f"  CX edges={len(cx_err)}"
                  f"  avg_err={sum(cx_err.values())/len(cx_err):.5f}"
                  f"  best_T1={max(t1_us.values()):.1f}µs")
            cm, cm_edges = get_ibm_cm(backend)
            print(f"  Coupling map: {len(cm_edges)} edges (actual)")

        print(f"\n[5/8] Layout (t={n_bh}, largest circuit) ...")
        rep = build_circuit(n_bh, 0, bh_gates, n_phys_bh, blk, sp, su, ciqa_obj,
                             n_bh=n_bh, syn=not USE_SIMULATOR)
        pinned = ciqs_layout(rep, cm)

        if not USE_SIMULATOR and cm_edges:
            layout = fix_layout_for_R(pinned, n_bh, cm_edges)
            print(f"  R qubits reassigned adjacent to anchors.")
        else:
            layout = [pinned[i] for i in sorted(pinned.keys())]
        print(f"  Layout (first 12): {layout[:12]}")

        print(f"\n[6/8] Transpiling {nc} circuits ...")
        if not USE_SIMULATOR:
            pm_max = generate_preset_pass_manager(
                backend=backend, optimization_level=1, initial_layout=layout)

        all_tc = []
        t0 = time.time()
        for t_step in range(1, n_bh+1):
            for sidx in range(n_sh):
                qc = build_circuit(t_step, sidx, bh_gates, n_phys_bh, blk, sp, su,
                                    ciqa_obj, n_bh=n_bh, syn=not USE_SIMULATOR)
                if USE_SIMULATOR:
                    tc = qtranspile(qc,
                        basis_gates=['cx','ry','rz','x','h','measure'],
                        optimization_level=0)
                else:
                    l2 = layout[:qc.num_qubits]
                    pm = generate_preset_pass_manager(
                        backend=backend, optimization_level=1, initial_layout=l2)
                    tc = pm.run(qc)
                all_tc.append(tc)
            print(f"  t={t_step}/{n_bh}  ({t_step*n_sh}/{nc})")
        print(f"  Done in {time.time()-t0:.1f}s")

        print(f"\n[7/8] Executing {nc} x {shots} shots ...")
        ts = time.time()
        if USE_SIMULATOR:
            from qiskit_aer.primitives import SamplerV2 as AerSampler
            sampler = AerSampler()
            job     = sampler.run(all_tc, shots=shots)
            result  = job.result()
            job_id  = f"sim_{int(time.time())}"
            print(f"  Done in {time.time()-ts:.1f}s")
        else:
            sampler = IBMSampler(backend)
            job     = sampler.run(all_tc, shots=shots)
            job_id  = job.job_id()
            print(f"  Job: {job_id}  waiting ...")
            result  = job.result()
            print(f"  Received in {time.time()-ts:.1f}s")

    print("\n[8/8] Post-processing ...")
    S2_D, S2_R, mean_I2_by_t = [], [], []
    I2_by_t = {}
    t_vals  = list(range(1, n_bh+1))

    for t_step in t_vals:
        shD, shR, shP = [], [], []

        for sidx in range(n_sh):
            raw = get_counts(result[cidx[(t_step,sidx)]])
            D,R,pc,n_a,_ = extract_counts(raw, n_bh, t_step,
                                            syn=not USE_SIMULATOR,
                                            post_sel=not USE_SIMULATOR)
            if n_a==0:
                D,R,pc,n_a,_ = extract_counts(raw, n_bh, t_step,
                                               syn=False, post_sel=False)
            shD.append(D); shR.append(R); shP.append(pc)

        s2d,_ = renyi2(shD, t_step, n_sh)
        s2r,_ = renyi2(shR, n_bh,   n_sh)
        pairs, mi = compute_I2_pairs(shD, shR, shP, t_step, n_bh, n_sh)

        S2_D.append(s2d); S2_R.append(s2r); mean_I2_by_t.append(mi)
        I2_by_t[t_step] = pairs

        i2s = [f"{x['I2']:+.3f}" for x in pairs]
        print(f"  t={t_step}  S2(D)={s2d:.4f}  S2(R)={s2r:.4f}"
              f"  mean_I2={mi:+.4f}  [{', '.join(i2s)}]")

    print_results(t_vals, S2_D, S2_R, I2_by_t, mean_I2_by_t, n_bh)
    pname = plot_results(t_vals, S2_D, S2_R, I2_by_t, mean_I2_by_t, n_bh, job_id)

    pre  = [mean_I2_by_t[i] for i,t in enumerate(t_vals) if t<=pt]
    post = [mean_I2_by_t[i] for i,t in enumerate(t_vals) if t> pt]

    data = {
        "job_id":            job_id,
        "backend":           BACKEND_NAME,
        "n_bh":              n_bh,
        "blk":               blk,
        "n_phys_bh":         n_phys_bh,
        "n_shadows":         n_sh,
        "shots_per_circ":    shots,
        "total_shots":       nc*shots,
        "scrambler_depth":   SCRAMBLER_DEPTH,
        "scrambler_seed":    SCRAMBLER_SEED,
        "shadow_seed":       SHADOW_SEED,
        "t_values":          t_vals,
        "page_time":         pt,
        "S2_D":              S2_D,
        "S2_R":              S2_R,
        "S2_page":           [min(t,n_bh-t)*math.log2(2) for t in t_vals],
        "mean_I2_by_t":      mean_I2_by_t,
        "I2_matched_pairs":  {
            str(t): [{"D_qubit":p["r"],"R_qubit":p["m"],
                      "S2_D":p["S2_D"],"S2_R":p["S2_R"],
                      "S2_DR":p["S2_DR"],"I2":p["I2"]}
                     for p in pairs]
            for t,pairs in I2_by_t.items()
        },
        "mean_I2_before_page": float(np.mean(pre))  if pre  else None,
        "mean_I2_after_page":  float(np.mean(post)) if post else None,
        "delta_I2":            float(np.mean(post)-np.mean(pre)) if pre and post else None,
        "timestamp":           time.strftime("%Y-%m-%d %H:%M:%S UTC", time.gmtime()),
    }
    jname = f"hp_v5_{job_id}.json"
    with open(jname,"w") as f: json.dump(data,f,indent=2)
    print(f"\n  JSON -> {jname}")
    print(f"  Plot -> {pname}")
    print("\n  Done.\n")

if __name__ == "__main__":
    run()
