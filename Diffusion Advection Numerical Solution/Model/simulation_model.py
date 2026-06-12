import numpy as np
from pathlib import Path

OUTPUT_DIR = Path(__file__).resolve().parent.parent / "Results"
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)


def simulate_diffusion_advection(initial_value, diffusion=0.09, velocity=0.12, decay=0.02, timesteps=48):
    output = np.zeros(timesteps)
    output[0] = float(np.clip(initial_value, 0.0, 1.0))

    for step in range(1, timesteps):
        growth = diffusion * (1.0 - output[step - 1])
        transport = velocity * output[step - 1]
        output[step] = output[step - 1] + growth - transport - decay * output[step - 1]
        output[step] = np.clip(output[step], 0.0, 1.0)

    return output


def mean_absolute_error(target, prediction):
    return np.mean(np.abs(np.array(target) - np.array(prediction)))


def fit_model_to_mean(signal_mean, trials=250, seed=2026):
    rng = np.random.default_rng(seed)
    best_error = None
    best_params = None
    best_fit = None
    timesteps = len(signal_mean)

    for _ in range(trials):
        diffusion = rng.uniform(0.01, 0.18)
        velocity = rng.uniform(0.01, 0.25)
        decay = rng.uniform(0.0, 0.08)
        prediction = simulate_diffusion_advection(signal_mean[0], diffusion=diffusion, velocity=velocity, decay=decay, timesteps=timesteps)
        error = mean_absolute_error(signal_mean, prediction)

        if best_error is None or error < best_error:
            best_error = error
            best_params = {
                "diffusion": float(diffusion),
                "velocity": float(velocity),
                "decay": float(decay),
                "mae": float(error),
            }
            best_fit = prediction

    return best_params, best_fit
