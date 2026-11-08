"""Ideal isocratic two-solute model. Lengths in mm, time in seconds."""
from dataclasses import dataclass, asdict
import math
import numpy as np

@dataclass(frozen=True)
class Conditions:
    u: float = 1.0
    A: float = 8.0   # micrometres
    B: float = 6.0   # micrometres * mm / s
    C: float = 6.0   # micrometres * s / mm
    L: float = 150.0
    k: float = 3.0
    alpha: float = 1.12

def calculate(p):
    if min(p.u, p.B, p.C, p.L, p.k) <= 0 or p.A < 0 or p.alpha < 1:
        raise ValueError('Positive velocity, B, C, length and retention are required; alpha >= 1.')
    H = p.A + p.B / p.u + p.C * p.u
    N = p.L * 1000 / H
    t0 = p.L / p.u
    t1, t2 = t0 * (1 + p.k), t0 * (1 + p.alpha * p.k)
    s1, s2 = t1 / math.sqrt(N), t2 / math.sqrt(N)
    return dict(H=H, N=N, t0=t0, t1=t1, t2=t2, s1=s1, s2=s2,
                Rs=(t2-t1)/(2*(s1+s2)), uopt=math.sqrt(p.B/p.C),
                Rs_approx=math.sqrt(N)/4*(p.alpha-1)/p.alpha*(p.alpha*p.k)/(1+p.alpha*p.k))

def chromatogram(p, normalized=False):
    m = calculate(p)
    t = np.linspace(max(0, m['t1']-5*m['s1']), m['t2']+5*m['s2'], 1600)
    y1 = np.exp(-0.5*((t-m['t1'])/m['s1'])**2)/(m['s1']*math.sqrt(2*math.pi))
    y2 = np.exp(-0.5*((t-m['t2'])/m['s2'])**2)/(m['s2']*math.sqrt(2*math.pi))
    return (t/m['t0'], y1*m['t0'], y2*m['t0']) if normalized else (t/60, y1*60, y2*60)
