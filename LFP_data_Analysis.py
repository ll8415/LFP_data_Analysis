#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
LFP GUI Runner (Blackrock) — Compat v4 with Multi-band Decomposition
-------------------------------------------------------------------
在 v3 基础上新增：
- “Multi-band decomposition”按钮：按常见 LFP 展示样式，将原始信号 + 多个频段（delta/theta/alpha/beta/low-gamma/high-gamma）
  逐行堆叠显示（同一时间轴）。
- 兼容无 channel_index 的 Neo 版本；CWT 起点自动夹紧；目录/文件双兼容。

依赖：neo, elephant, quantities, numpy, matplotlib, (pywt 可选, scipy 兜底)
"""

import warnings
import numpy as np
import matplotlib.pyplot as plt
import quantities as pq
from pathlib import Path
import neo
from neo import io
from elephant.signal_processing import butter
from elephant.spectral import welch_psd

# optional imports
try:
    import pywt
    _HAS_PYWT = True
except Exception:
    _HAS_PYWT = False

try:
    from scipy.signal import cwt, morlet2
    _HAS_SCIPY = True
except Exception:
    _HAS_SCIPY = False

import tkinter as tk
from tkinter import ttk, filedialog, messagebox

# ---------------- Open Blackrock (compat) ----------------
_NS_GLOB = ["*.ns1", "*.ns2", "*.ns3", "*.ns4", "*.ns5", "*.ns6", "*.nev", "*.ns*"]

def _pick_a_blackrock_file_in_dir(d: Path):
    for pat in _NS_GLOB:
        files = list(d.glob(pat))
        if files:
            files_sorted = sorted(files, key=lambda p: (p.suffix.lower() != ".ns2", p.name.lower()))
            return files_sorted[0]
    return None

def open_blackrock_compat(path_str: str):
    p = Path(path_str)
    if p.is_file():
        return io.BlackrockIO(filename=str(p))
    try:
        return io.BlackrockIO(dirname=str(p))
    except TypeError:
        f = _pick_a_blackrock_file_in_dir(p)
        if f is None:
            raise FileNotFoundError("No nsx/nev file found in the selected folder")
        return io.BlackrockIO(filename=str(f))

# ---------------- Legacy-compatible helpers ----------------
def time_slice_seconds(a, b):
    def _to_sec(x):
        return x if hasattr(x, 'dimensionality') else x * pq.s
    return (_to_sec(a), _to_sec(b))

def group(r):
    blk = r.read_block(lazy=False)
    return blk.groups[0]

def sigs(r, index):
    seg = r.read_segment(lazy=True)
    return seg.analogsignals[index]

def highpass2Hz(r, index, freq=2.0):
    x = sigs(r, index).load()
    return butter(x, highpass_frequency=freq * pq.Hz)

def channel_select(asig, ch_i):
    """构造单通道 AnalogSignal（兼容无 .channel_index 的 Neo 版本）。"""
    y = asig.magnitude[:, [ch_i]]
    return neo.AnalogSignal(y, units=asig.units, t_start=asig.t_start,
                            sampling_rate=asig.sampling_rate, name=asig.name)

def list_streams(r):
    seg = r.read_segment(lazy=False)
    info = []
    for k, a in enumerate(seg.analogsignals):
        fs = float(a.sampling_rate.rescale('Hz').magnitude)
        n_ch = a.shape[1]
        name = getattr(a, 'name', f'analog[{k}]')
        info.append((k, fs, n_ch, name))
    return info

# ---------------- Plotting (raw vs hp / PSD / CWT) ----------------
def sigs_plot(r, index):
    raw = sigs(r, index).load()
    hp  = highpass2Hz(r, index)
    n_ch = raw.shape[1]
    names = group(r).annotations.get("ch_names", [f"ch{i}" for i in range(n_ch)])
    rows = int(np.ceil(n_ch/4))
    fig, axes = plt.subplots(rows, 4, figsize=(14, 3*rows), tight_layout=True, sharex=True, sharey=True)
    axes = np.ravel(axes)
    for i in range(min(n_ch, len(axes))):
        raw_i = channel_select(raw, i)
        hp_i  = channel_select(hp , i)
        ax = axes[i]
        ax.plot(raw_i.times.rescale('s').magnitude, raw_i.squeeze(), alpha=0.5)
        ax.plot(hp_i .times.rescale('s').magnitude,  hp_i.squeeze(),  alpha=0.5)
        ax.set_title(names[i])
        ax.legend(["raw", "filtered"], fontsize=8)
    for j in range(min(n_ch, len(axes)), len(axes)):
        axes[j].axis('off')
    plt.show()

def psd_plot(r, index):
    raw = sigs(r, index).load()
    hp  = highpass2Hz(r, index)
    n_ch = raw.shape[1]
    rows = int(np.ceil(n_ch/4))
    fig, axes = plt.subplots(rows, 4, figsize=(14, 3*rows), tight_layout=True, sharex=True, sharey=True)
    axes = np.ravel(axes)
    f_raw, p_raw = welch_psd(raw)
    f_hp,  p_hp  = welch_psd(hp)
    names = group(r).annotations.get("ch_names", [f"ch{i}" for i in range(n_ch)])
    for i in range(min(n_ch, len(axes))):
        x1 = channel_select(raw, i)
        x2 = channel_select(hp , i)
        f1, p1 = welch_psd(x1)
        f2, p2 = welch_psd(x2)
        ax = axes[i]
        ax.plot(f1, p1.squeeze(), alpha=.2)
        ax.plot(f2, p2.squeeze(), alpha=.4)
        ax.plot(f_raw, np.mean(p_raw, axis=0), alpha=.4)
        ax.plot(f_hp,  np.mean(p_hp , axis=0), alpha=.4)
        ax.set_title("PSD-" + names[i])
        ax.legend(["raw", "filtered", "raw_mean", "highpass_mean"], fontsize=8)
        ax.set_xlabel("frequency (Hz)")
        ax.set_ylabel("power($\\mu V^2$/Hz)")
        ax.set_xlim(0, 140)
    for j in range(min(n_ch, len(axes)), len(axes)):
        axes[j].axis('off')
    plt.show()

def _cwt_morlet(y, fs, totalscale=2**10):
    if _HAS_PYWT:
        fc = pywt.central_frequency('morl')
        const = fs * fc
        scales = const / np.arange(totalscale - 1, 1, -1)
        coefs, freqs = pywt.cwt(y, scales=scales, wavelet='morl', sampling_period=1.0/fs, method='conv')
        return coefs, freqs
    elif _HAS_SCIPY:
        w = 5.0
        const = fs * (w / (2.0*np.pi))
        scales = const / np.arange(totalscale - 1, 1, -1)
        widths = scales
        coefs = cwt(y, lambda M, s: morlet2(M, s, w=w), widths)
        freqs = const / widths
        return coefs, freqs
    else:
        raise ImportError("Neither pywt nor scipy available for CWT. Please install 'pywavelets' or 'scipy'.")

def wavelet_plot(data_dict, index, t_start_list):
    for i, site in enumerate(data_dict.values()):
        r = open_blackrock_compat(site)
        whole = sigs(r, index).load()
        dur_s = float((whole.t_stop - whole.t_start).rescale('s').magnitude)
        if dur_s <= 0:
            raise ValueError("Signal duration is zero; cannot analyze.")
        t0_req = t_start_list[i][0] if isinstance(t_start_list[i], (list, tuple, np.ndarray)) else t_start_list[i]
        t0 = max(0.0, min(float(t0_req), max(0.0, dur_s - 1e-6)))
        t1 = min(t0 + 410.0, dur_s - 1e-6)
        if t1 <= t0:
            raise ValueError(f"Start t0={t0_req} exceeds duration {dur_s:.3f}s. Please reduce t0.")

        base = whole.t_start
        eps = 1e-9 * pq.s
        t_abs0 = base + t0 * pq.s
        t_abs1 = min(base + t1 * pq.s, whole.t_stop - eps)
        x = whole.time_slice(t_abs0, t_abs1)
        hp = butter(x, highpass_frequency=2.0 * pq.Hz)
        fs = float(hp.sampling_rate.rescale('Hz').magnitude)
        y = np.mean(hp.magnitude, axis=1)
        coefs, freqs = _cwt_morlet(y, fs, totalscale=2**10)
        t = hp.times.rescale('s').magnitude
        plt.figure(figsize=(8, 4.5))
        plt.contourf(t, freqs, np.abs(coefs), levels=60)
        plt.yscale('log'); plt.title(f"Morlet CWT site {i} (t0={t0:.2f}s~{t1:.2f}s)")
        plt.xlabel("Time (s)"); plt.ylabel("Frequency (Hz)")
        plt.tight_layout(); plt.show()

# ---------------- Multi-band Decomposition ----------------
_BANDS = [
    ("Raw (no filter)", None, None),
    ("Delta 0.5–4 Hz", 0.5, 4),
    ("Theta 4–8 Hz", 4, 8),
    ("Alpha 8–12 Hz", 8, 12),
    ("Beta 13–30 Hz", 13, 30),
    ("Low-γ 30–60 Hz", 30, 60),
    ("High-γ 60–120 Hz", 60, 120),
]

def _bandpass(asig, lo, hi):
    kwargs = {}
    if lo is not None:
        kwargs["highpass_frequency"] = lo * pq.Hz
    if hi is not None:
        kwargs["lowpass_frequency"] = hi * pq.Hz
    if not kwargs:
        return asig  # raw
    return butter(asig, **kwargs)

def multiband_plot(r, index, t0=0.0, dur=10.0, ch=0):
    """
    画“原始 + 各频段”的堆叠图（同一时间轴），风格与参考图一致。
    默认单通道（ch=0），以便对比清晰；若要平均多通道可改为 asig.mean(axis=1)。
    """
    whole = sigs(r, index).load()
    dur_total = float((whole.t_stop - whole.t_start).rescale('s').magnitude)
    if dur_total <= 0:
        raise ValueError("Signal duration is zero; cannot analyze.")

    # clamp窗口
    t0 = max(0.0, min(float(t0), max(0.0, dur_total - 1e-6)))
    t1 = min(t0 + float(dur), dur_total - 1e-6)
    if t1 <= t0:
        raise ValueError("窗口长度太短或 t0 越界。")

    base = whole.t_start
    eps = 1e-9 * pq.s
    t_abs0 = base + t0 * pq.s
    t_abs1 = min(base + t1 * pq.s, whole.t_stop - eps)
    asig = whole.time_slice(t_abs0, t_abs1)
    one = channel_select(asig, int(ch))
    t = one.times.rescale('s').magnitude
    ylist, labels = [], []

    for name, lo, hi in _BANDS:
        xf = _bandpass(one, lo, hi)
        ylist.append(xf.squeeze())
        labels.append(name)

    # 绘图：每个频段一行，共享时间轴
    n = len(ylist)
    fig, axes = plt.subplots(n, 1, figsize=(10, max(6, 1.4*n)), sharex=True, tight_layout=True)
    if n == 1:
        axes = [axes]
    for i, (y, lab) in enumerate(zip(ylist, labels)):
        axes[i].plot(t, y, lw=1.0, color="k")
        axes[i].set_ylabel("μV")
        axes[i].set_title(lab, fontsize=10, loc="left")
        axes[i].grid(False)
        axes[i].spines["top"].set_visible(False)
        axes[i].spines["right"].set_visible(False)
    axes[-1].set_xlabel("Time (s)")
    fig.suptitle(f"Multi-band LFP (idx={index}, ch={ch}, {t0:.2f}–{t1:.2f}s)", y=1.02, fontsize=12)
    plt.show()

# ---------------- Tk GUI ----------------
class App(tk.Tk):
    def __init__(self):
        super().__init__()
        self.title('LFP Analysis')
        self.geometry("900x640")
        self.path_var = tk.StringVar()
        self.index_var = tk.IntVar(value=0)
        self.tstart_var = tk.DoubleVar(value=0.0)
        self.tdur_var = tk.DoubleVar(value=10.0)   # 新增：多频段分解窗口长度
        self.chan_var = tk.IntVar(value=0)         # 新增：选择单通道
        self.info = []
        self.r = None

        frm = ttk.Frame(self, padding=10)
        frm.pack(fill="both", expand=True)

        row0 = ttk.Frame(frm); row0.pack(fill="x", pady=4)
        ttk.Label(row0, text="Data folder or file:").pack(side="left")
        ttk.Entry(row0, textvariable=self.path_var, width=70).pack(side="left", padx=6)
        ttk.Button(row0, text="Browse...", command=self.pick_path).pack(side="left")

        row1 = ttk.Frame(frm); row1.pack(fill="both", expand=True, pady=6)
        ttk.Button(row1, text="Scan streams (ns2/ns6)", command=self.scan_streams).pack(side="top", anchor="w")
        self.listbox = tk.Listbox(row1, height=12)
        self.listbox.pack(fill="both", expand=True, padx=0, pady=6)

        row2 = ttk.Frame(frm); row2.pack(fill="x", pady=6)
        ttk.Label(row2, text="Index:").pack(side="left")
        ttk.Entry(row2, textvariable=self.index_var, width=6).pack(side="left", padx=6)
        ttk.Label(row2, text="CWT start t0 (s):").pack(side="left", padx=(10,0))
        ttk.Entry(row2, textvariable=self.tstart_var, width=8).pack(side="left", padx=6)
        ttk.Label(row2, text="Window length (s):").pack(side="left", padx=(10,0))
        ttk.Entry(row2, textvariable=self.tdur_var, width=8).pack(side="left", padx=6)
        ttk.Label(row2, text="Channel:").pack(side="left", padx=(10,0))
        ttk.Entry(row2, textvariable=self.chan_var, width=6).pack(side="left", padx=6)

        row3 = ttk.Frame(frm); row3.pack(fill="x", pady=6)
        ttk.Button(row3, text="Raw vs 2 Hz High-pass", command=self.run_sigs).pack(side="left", padx=4)
        ttk.Button(row3, text="Welch PSD", command=self.run_psd).pack(side="left", padx=4)
        ttk.Button(row3, text="Morlet CWT (410s)", command=self.run_cwt).pack(side="left", padx=4)
        ttk.Button(row3, text="Multi-band decomposition", command=self.run_multiband).pack(side="left", padx=4)

        tip = ("提示：选择目录或文件后，先“扫描流”查看 index。\n"
               "Multi-band decomposition：默认画单通道（通道号可改），窗口长度可调，样式与常见 LFP 展示一致。")
        ttk.Label(frm, text=tip, foreground="#444").pack(side="bottom", anchor="w")

    def pick_path(self):
        p = filedialog.askdirectory(title="Select a folder with Blackrock files")
        if not p:
            p = filedialog.askopenfilename(title="选择 Blackrock 文件(ns2/ns6/nev)",
                                           filetypes=[("Blackrock", "*.ns1 *.ns2 *.ns3 *.ns4 *.ns5 *.ns6 *.nev *.ns*"), ("All", "*.*")])
        if p:
            self.path_var.set(p)
            self.r = None
            self.listbox.delete(0, "end")

    def _open_r(self):
        p = self.path_var.get().strip()
        if not p:
            messagebox.showerror("Error", "Please select a folder or file first"); return None
        try:
            r = open_blackrock_compat(p)
            return r
        except Exception as e:
            messagebox.showerror("Open failed", str(e))
            return None

    def scan_streams(self):
        r = self._open_r()
        if r is None: return
        self.r = r
        try:
            self.info = list_streams(r)
            self.listbox.delete(0, "end")
            self.listbox.insert("end", "index | fs_Hz | n_ch | name")
            for k, fs, n_ch, name in self.info:
                self.listbox.insert("end", f"{k:5d} | {fs:8.1f} | {n_ch:5d} | {name}")
        except Exception as e:
            messagebox.showerror("Read failed", str(e))

    def _selected_index(self):
        try:
            sel = self.listbox.curselection()
            if sel and sel[0] > 0:
                line = self.listbox.get(sel[0])
                idx = int(line.split("|")[0].strip())
                return idx
        except Exception:
            pass
        return int(self.index_var.get())

    def run_sigs(self):
        r = self.r or self._open_r()
        if r is None: return
        idx = self._selected_index()
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            sigs_plot(r, idx)

    def run_psd(self):
        r = self.r or self._open_r()
        if r is None: return
        idx = self._selected_index()
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            psd_plot(r, idx)

    def run_cwt(self):
        r = self.r or self._open_r()
        if r is None: return
        idx = self._selected_index()
        t0 = float(self.tstart_var.get())
        site = self.path_var.get().strip()
        data_dict = {"rec": site}
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            wavelet_plot(data_dict, idx, [t0])

    def run_multiband(self):
        r = self.r or self._open_r()
        if r is None: return
        idx = self._selected_index()
        t0 = float(self.tstart_var.get())
        dur = float(self.tdur_var.get())
        ch  = int(self.chan_var.get())
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            multiband_plot(r, idx, t0=t0, dur=dur, ch=ch)

if __name__ == "__main__":
    app = App()
    app.mainloop()
