import numpy as np
import pandas as pd
from pathlib import Path

OUTPUT_DIR = Path(__file__).resolve().parent.parent / "Results"
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)


def generate_synthetic_signals(num_signals=5, length=48, noise_scale=0.06, seed=2026):
    rng = np.random.default_rng(seed)
    time = np.linspace(0.0, 24.0, length)
    signals = []

    for idx in range(num_signals):
        amplitude = 0.75 + 0.10 * rng.standard_normal()
        offset = 0.03 * rng.standard_normal()
        ramp = np.tanh((time - 6.0) * (0.18 + 0.01 * idx))
        shape = amplitude * (0.4 + 0.6 * ramp)
        noise = noise_scale * rng.standard_normal(size=length)
        signal = np.clip(shape + offset + noise, 0.0, 1.0)
        signals.append(signal)

    return time, np.vstack(signals)


def aggregate_signals(signal_matrix):
    return np.mean(signal_matrix, axis=0)


def save_synthetic_dataset(time, signal_matrix, filename="synthetic_signals.csv"):
    columns = [f"signal_{index+1}" for index in range(signal_matrix.shape[0])]
    df = pd.DataFrame(signal_matrix.T, columns=columns)
    df.insert(0, "time", time)
    path = OUTPUT_DIR / filename
    df.to_csv(path, index=False)
    return path, df
