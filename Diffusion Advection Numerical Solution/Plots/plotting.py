import matplotlib.pyplot as plt
import numpy as np
from pathlib import Path

OUTPUT_DIR = Path(__file__).resolve().parent.parent / "Results"
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

try:
    import seaborn as sns
    sns.set_theme(style="whitegrid")
    HAS_SEABORN = True
except ImportError:
    HAS_SEABORN = False


def plot_signal_matrix(time, signal_matrix, title, filename):
    plt.figure(figsize=(10, 6))
    num_signals = signal_matrix.shape[0]
    for idx in range(num_signals):
        plt.plot(time, signal_matrix[idx], marker="o", linestyle="--", alpha=0.55, label=f"Sim {idx + 1}")

    mean_line = np.mean(signal_matrix, axis=0)
    std_line = np.std(signal_matrix, axis=0)
    plt.plot(time, mean_line, color="black", linewidth=2.5, label="Mean signal")
    plt.fill_between(time, mean_line - std_line, mean_line + std_line, color="gray", alpha=0.24, label="±1 std")

    plt.title(title, fontsize=14)
    plt.xlabel("Time step", fontsize=12)
    plt.ylabel("Normalized intensity", fontsize=12)
    plt.ylim(0.0, 1.05)
    plt.grid(True, linestyle="--", linewidth=0.5, alpha=0.7)
    plt.legend(loc="upper left", fontsize=9)
    plt.tight_layout()
    path = OUTPUT_DIR / filename
    plt.savefig(path, dpi=200)
    plt.close()
    return path


def plot_mean_vs_fit(time, mean_signal, fit_signal, params, filename):
    plt.figure(figsize=(10, 6))
    plt.plot(time, mean_signal, marker="o", linestyle="-", color="#1f77b4", label="Synthetic mean")
    plt.plot(time, fit_signal, marker="s", linestyle="--", color="#ff7f0e", label="Model fit")

    plt.title("Mean synthetic signal vs. fitted diffusion model", fontsize=14)
    plt.xlabel("Time step", fontsize=12)
    plt.ylabel("Normalized intensity", fontsize=12)
    plt.ylim(0.0, 1.05)
    plt.grid(True, linestyle="--", linewidth=0.5, alpha=0.7)
    plt.legend(loc="upper left", fontsize=10)

    annotation = (
        f"diffusion={params['diffusion']:.3f}, velocity={params['velocity']:.3f}, "
        f"decay={params['decay']:.3f}, mae={params['mae']:.4f}"
    )
    plt.annotate(annotation, xy=(0.02, 0.02), xycoords="axes fraction", fontsize=10,
                 bbox=dict(boxstyle="round,pad=0.3", fc="white", ec="black", alpha=0.7))

    plt.tight_layout()
    path = OUTPUT_DIR / filename
    plt.savefig(path, dpi=200)
    plt.close()
    return path


def plot_signal_distribution(df, filename):
    plt.figure(figsize=(8, 5))
    if HAS_SEABORN:
        sns.boxplot(data=df.drop(columns=["time"]))
        sns.stripplot(data=df.drop(columns=["time"]), color="black", alpha=0.5, jitter=True)
    else:
        signal_columns = df.columns.drop("time")
        for idx, column in enumerate(signal_columns):
            plt.scatter([idx] * len(df), df[column], alpha=0.6)

    plt.title("Distribution of anonymized synthetic signals", fontsize=14)
    plt.xlabel("Synthetic signal index", fontsize=12)
    plt.ylabel("Normalized intensity", fontsize=12)
    plt.tight_layout()
    path = OUTPUT_DIR / filename
    plt.savefig(path, dpi=200)
    plt.close()
    return path
