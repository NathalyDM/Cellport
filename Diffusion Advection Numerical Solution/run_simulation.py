from DataGeneration.generate_data import generate_synthetic_signals, save_synthetic_dataset, aggregate_signals
from Model.simulation_model import fit_model_to_mean
from Plots.plotting import plot_signal_matrix, plot_mean_vs_fit, plot_signal_distribution


def main():
    time, signal_matrix = generate_synthetic_signals(num_signals=5, length=48, noise_scale=0.06)
    dataset_path, df = save_synthetic_dataset(time, signal_matrix)

    plot_signal_matrix(time, signal_matrix, "Anonymized Synthetic Diffusion-Advection Simulations", "simulation_curves.png")
    plot_signal_distribution(df, "signal_distribution.png")

    mean_signal = aggregate_signals(signal_matrix)
    params, fit_signal = fit_model_to_mean(mean_signal, trials=300)
    plot_mean_vs_fit(time, mean_signal, fit_signal, params, "mean_fit_comparison.png")

    print("Synthetic dataset saved to:", dataset_path)
    print("Generated plots in Results/ folder")


if __name__ == "__main__":
    main()
