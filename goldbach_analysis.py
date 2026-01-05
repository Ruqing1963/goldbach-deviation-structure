#!/usr/bin/env python3
"""
Goldbach Deviation Structure Analysis
=====================================

This script performs fractal-spectral analysis of Goldbach deviations,
computing:
1. Goldbach representation counts G(N)
2. Hardy-Littlewood predictions HL(N)
3. Relative deviations δ(N)
4. Hurst exponent via R/S analysis
5. FFT spectral analysis

Author: Ruqing Chen
Email: ruqing@hotmail.com
License: CC BY 4.0
"""

import numpy as np
from scipy import stats
import matplotlib.pyplot as plt
import pandas as pd
import os

# Twin prime constant
C2 = 0.6601618158

def sieve_of_eratosthenes(limit):
    """Generate prime sieve up to limit."""
    is_prime = [True] * (limit + 1)
    is_prime[0] = is_prime[1] = False
    for i in range(2, int(limit**0.5) + 1):
        if is_prime[i]:
            for j in range(i*i, limit + 1, i):
                is_prime[j] = False
    return is_prime

def compute_goldbach_counts(N_max, step=100):
    """
    Compute Goldbach representation counts and deviations.
    
    Parameters
    ----------
    N_max : int
        Maximum even number to analyze
    step : int
        Step size between consecutive N values
        
    Returns
    -------
    DataFrame with columns: N, G_N, HL_N, delta
    """
    print(f"Computing Goldbach counts up to N = {N_max}...")
    
    is_prime = sieve_of_eratosthenes(N_max)
    primes = [i for i in range(2, N_max) if is_prime[i]]
    prime_set = set(primes)
    
    results = []
    for N in range(4, N_max + 1, step):
        if N % 2 != 0:
            continue
            
        # Count representations
        count = sum(1 for p in primes if p <= N//2 and (N - p) in prime_set)
        
        # Hardy-Littlewood prediction
        ln_N = np.log(N)
        S_N = 1.0
        for p in [3, 5, 7, 11, 13, 17, 19]:
            if N % p == 0:
                S_N *= (p - 1) / (p - 2) if p > 2 else 1
            else:
                S_N *= 1 - 1/((p-1)**2) if p > 2 else 1
        HL = 2 * C2 * S_N * N / (ln_N ** 2)
        
        # Relative deviation
        delta = (count - HL) / HL if HL > 0 else 0
        
        results.append({
            'N': N,
            'ln_N': ln_N,
            'G_N': count,
            'HL_N': HL,
            'delta': delta
        })
    
    return pd.DataFrame(results)

def compute_hurst_exponent(series):
    """
    Compute Hurst exponent using Rescaled Range (R/S) method.
    
    Parameters
    ----------
    series : array-like
        Time series data
        
    Returns
    -------
    H : float
        Hurst exponent
    window_sizes : list
        Window sizes used
    rs_values : list
        R/S values for each window size
    """
    print("Computing Hurst exponent...")
    
    n = len(series)
    window_sizes = []
    rs_values = []
    
    for w in range(10, n // 3):
        num_windows = n // w
        if num_windows < 2:
            continue
            
        rs_list = []
        for i in range(num_windows):
            window = series[i*w:(i+1)*w]
            mean = np.mean(window)
            Y = np.cumsum(window - mean)
            R = np.max(Y) - np.min(Y)
            S = np.std(window, ddof=1)
            if S > 0 and R > 0:
                rs_list.append(R / S)
                
        if len(rs_list) > 0:
            window_sizes.append(w)
            rs_values.append(np.mean(rs_list))
    
    if len(window_sizes) < 5:
        return np.nan, [], []
        
    log_w = np.log10(window_sizes)
    log_rs = np.log10(rs_values)
    H, intercept, r_value, p_value, std_err = stats.linregress(log_w, log_rs)
    
    print(f"  Hurst exponent H = {H:.3f} ± {std_err:.3f}")
    print(f"  R² = {r_value**2:.4f}")
    
    return H, window_sizes, rs_values

def compute_fft_spectrum(delta_series, ln_N_series):
    """
    Compute FFT spectrum of deviation series.
    
    Parameters
    ----------
    delta_series : array-like
        Deviation values
    ln_N_series : array-like
        ln(N) values (for frequency calculation)
        
    Returns
    -------
    frequencies : array
        Frequency values (gamma)
    power : array
        Normalized power spectrum
    """
    print("Computing FFT spectrum...")
    
    # Remove mean
    delta_centered = delta_series - np.mean(delta_series)
    
    # FFT
    fft_result = np.fft.fft(delta_centered)
    d_ln = ln_N_series[1] - ln_N_series[0]
    freqs = np.fft.fftfreq(len(delta_series), d=d_ln)
    
    # Power spectrum (positive frequencies only)
    n_half = len(fft_result) // 2
    power = np.abs(fft_result[:n_half])**2
    power = power / np.max(power)  # Normalize
    
    # Convert to gamma
    gamma = 2 * np.pi * freqs[:n_half]
    
    return gamma, power

def plot_all_figures(df, H, ws, rs, gamma, power, output_dir='../figures'):
    """Generate all analysis figures."""
    
    os.makedirs(output_dir, exist_ok=True)
    
    # Figure 1: Scale Separation
    print("Generating Figure 1: Scale Separation...")
    fig, axes = plt.subplots(1, 2, figsize=(12, 5))
    
    N_range = np.logspace(3, 14, 200)
    signal = N_range / (np.log(N_range)**2)
    noise = np.sqrt(N_range) * np.log(N_range)**2
    
    axes[0].loglog(N_range, signal, 'g-', lw=2, label=r'Signal: $S \sim N/(\ln N)^2$')
    axes[0].loglog(N_range, noise, 'r-', lw=2, label=r'Noise: $\Delta \sim \sqrt{N}(\ln N)^2$')
    axes[0].fill_between(N_range, noise, signal, alpha=0.15, color='green')
    axes[0].set_xlabel('N')
    axes[0].set_ylabel('Magnitude')
    axes[0].set_title('(A) Scale Separation')
    axes[0].legend()
    axes[0].grid(True, alpha=0.3)
    
    relative_error = noise / signal
    axes[1].loglog(N_range, relative_error, 'b-', lw=2)
    axes[1].axhline(y=1, color='red', ls='--', lw=2, label='Failure threshold')
    axes[1].set_xlabel('N')
    axes[1].set_ylabel(r'$|\delta(N)|$')
    axes[1].set_title('(B) Relative Error Decay')
    axes[1].legend()
    axes[1].grid(True, alpha=0.3)
    
    plt.tight_layout()
    plt.savefig(f'{output_dir}/fig1_scale_separation.png', dpi=150)
    plt.close()
    
    # Figure 2: Hurst Analysis
    print("Generating Figure 2: Hurst Analysis...")
    fig, axes = plt.subplots(1, 2, figsize=(12, 5))
    
    axes[0].plot(df['ln_N'], df['delta'], 'b-', lw=0.8, alpha=0.7)
    axes[0].axhline(y=0, color='green', lw=2)
    axes[0].axhline(y=-1, color='red', ls='--', lw=2, label='Failure threshold')
    axes[0].fill_between(df['ln_N'], df['delta'], alpha=0.3)
    axes[0].set_xlabel(r'$\ln N$')
    axes[0].set_ylabel(r'$\delta(N)$')
    axes[0].set_title('(A) Goldbach Residual Series')
    axes[0].legend()
    axes[0].grid(True, alpha=0.3)
    
    if len(ws) > 0:
        axes[1].loglog(ws, rs, 'bo-', ms=4, label='Data')
        log_w = np.log10(ws)
        log_rs = np.log10(rs)
        fit_line = 10**(H * log_w + (log_rs[0] - H * log_w[0]))
        axes[1].loglog(ws, fit_line, 'r-', lw=2, label=f'Fit: H = {H:.3f}')
        random_line = [w**0.5 for w in ws]
        axes[1].loglog(ws, random_line, 'k--', alpha=0.5, label='Random (H=0.5)')
    axes[1].set_xlabel('Window Size')
    axes[1].set_ylabel('R/S Statistic')
    axes[1].set_title(f'(B) R/S Analysis: H = {H:.2f}')
    axes[1].legend()
    axes[1].grid(True, alpha=0.3)
    
    plt.tight_layout()
    plt.savefig(f'{output_dir}/fig2_hurst_analysis.png', dpi=150)
    plt.close()
    
    # Figure 3: FFT Spectrum
    print("Generating Figure 3: FFT Spectrum...")
    fig, ax = plt.subplots(figsize=(10, 6))
    
    ax.semilogy(gamma[:300], power[:300], 'b-', lw=1.5)
    
    L_zeros = {'L(s,χ₃)': 8.04, 'L(s,χ₅)': 6.02, 'ζ(s)': 14.13, 'ζ(s)₂': 21.02}
    colors = ['red', 'green', 'purple', 'orange']
    for (name, zero), color in zip(L_zeros.items(), colors):
        ax.axvline(x=zero, color=color, ls='--', alpha=0.7, label=f'{name}: γ={zero}')
    
    ax.set_xlabel(r'Frequency $\gamma$')
    ax.set_ylabel('Normalized Power')
    ax.set_title('FFT Spectrum of Goldbach Residuals')
    ax.legend()
    ax.grid(True, alpha=0.3)
    ax.set_xlim(0, 35)
    
    plt.tight_layout()
    plt.savefig(f'{output_dir}/fig3_fft_spectrum.png', dpi=150)
    plt.close()
    
    # Figure 4: Safety Gap Diagram
    print("Generating Figure 4: Safety Gap Diagram...")
    fig, ax = plt.subplots(figsize=(10, 7))
    
    ln_N = np.linspace(7, 28, 500)
    envelope = C2 / ln_N
    
    ax.axhline(-1, color='darkred', lw=3, label='Death Line (δ=-1)')
    ax.plot(ln_N, envelope, 'g--', lw=2.5, label=r'Envelope: $+C_2/\ln N$')
    ax.plot(ln_N, -envelope, 'g--', lw=2.5, label=r'Envelope: $-C_2/\ln N$')
    ax.fill_between(ln_N, -envelope, envelope, alpha=0.1, color='green')
    
    # Plot actual data
    ax.scatter(df['ln_N'], df['delta'], c='blue', s=5, alpha=0.5, label='Data')
    
    ax.set_xlabel(r'$\ln N$')
    ax.set_ylabel(r'$\delta(N)$')
    ax.set_title('Safety Gap Diagram')
    ax.legend()
    ax.grid(True, alpha=0.3)
    ax.set_ylim(-1.1, 0.5)
    
    plt.tight_layout()
    plt.savefig(f'{output_dir}/fig4_safety_gap.png', dpi=150)
    plt.close()
    
    print("✅ All figures saved!")

def main():
    """Main analysis pipeline."""
    print("=" * 60)
    print("Goldbach Deviation Structure Analysis")
    print("=" * 60)
    
    # Compute Goldbach data
    df = compute_goldbach_counts(500000, step=100)
    
    # Save data
    os.makedirs('../data', exist_ok=True)
    df.to_csv('../data/goldbach_counts.csv', index=False)
    print(f"✅ Data saved: {len(df)} points")
    
    # Hurst analysis
    H, ws, rs = compute_hurst_exponent(df['delta'].values)
    
    # Save Hurst results
    hurst_df = pd.DataFrame({'window_size': ws, 'rs_value': rs})
    hurst_df.to_csv('../data/hurst_analysis.csv', index=False)
    
    # FFT analysis
    gamma, power = compute_fft_spectrum(df['delta'].values, df['ln_N'].values)
    
    # Generate figures
    plot_all_figures(df, H, ws, rs, gamma, power)
    
    # Summary statistics
    print("\n" + "=" * 60)
    print("SUMMARY")
    print("=" * 60)
    print(f"  N range: {df['N'].min()} to {df['N'].max()}")
    print(f"  Data points: {len(df)}")
    print(f"  Hurst exponent: H = {H:.3f}")
    print(f"  δ range: [{df['delta'].min():.4f}, {df['delta'].max():.4f}]")
    print(f"  Safety margin: {1 + df['delta'].min():.4f}")
    print(f"  Twin prime constant C₂ = {C2}")

if __name__ == '__main__':
    main()
