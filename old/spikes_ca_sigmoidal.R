library(ggplot2)
library(patchwork)

# --- Sigmoid forward function (spikes → Ca) ---
sigmoid <- function(x, A=1, k=0.6, x0=8) {
  A / (1 + exp(-k * (x - x0)))
}

# --- Inverse sigmoid (Ca → spikes) ---
inverse_sigmoid <- function(ca, A=1, k=0.6, x0=8) {
  x0 - (1/k) * log((A/ca) - 1)
}

# --- Data for forward curve ---
spike_vals <- seq(0, 20, length.out = 200)
ca_vals <- sigmoid(spike_vals)

df_forward <- data.frame(
  spikes = spike_vals,
  ca = ca_vals
)

# --- Data for inverse curve ---
ca_vals2 <- seq(0.01, 0.99, length.out = 200)
spike_est <- inverse_sigmoid(ca_vals2)

df_inverse <- data.frame(
  ca = ca_vals2,
  spikes = spike_est
)

# --- Panel 1: Forward mapping ---
p1 <- ggplot(df_forward, aes(x = spikes, y = ca)) +
  geom_line(size = 1.5, color = "#2C7BB6") +
  labs(
    x = "Spikes",
    y = "Calcium signal (a.u.)",
    title = "Forward Mapping: Spikes → Ca"
  ) +
  theme_minimal(base_size = 14) +
  theme(plot.title = element_text(face = "bold", size = 16))

# --- Panel 2: Inverse mapping ---
p2 <- ggplot(df_inverse, aes(x = ca, y = spikes)) +
  geom_line(size = 1.5, color = "#D1495B") +
  labs(
    x = "Calcium signal (a.u.)",
    y = "Spike count (equivalent)",
    title = "Inverse Mapping: Ca → Spikes"
  ) +
  theme_minimal(base_size = 14) +
  theme(plot.title = element_text(face = "bold", size = 16))

# --- Combine side by side ---
p1 + p2
