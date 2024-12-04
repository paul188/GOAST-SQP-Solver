import matplotlib.pyplot as plt

# Data points
merit_values = [
    60.4868, 55.0123, 50.6482, 47.1442, 44.3330, 42.3064, 40.7909, 39.5046,
    38.4274, 37.5446, 36.8460, 36.5264, 36.3832, 36.4214, 36.6525, 37.0955,
    37.7800, 38.7488, 40.0639, 41.8119
]

linearization_values = [
    60.4868, 57.4439, 54.4011, 51.3582, 48.3154, 45.2725, 42.2297, 39.1868,
    36.1440, 33.1011, 30.0583, 27.0154, 23.9725, 20.9297, 17.8868, 14.8440,
    11.8011, 8.7583, 5.7154, 2.6726
]

# Generate x-axis values
iterations = list(range(1, len(merit_values) + 1))

# Create the plot
plt.figure(figsize=(10, 6))
plt.plot(iterations, merit_values, label="Merit", marker="o", linestyle="-")
plt.plot(iterations, linearization_values, label="Linearization", marker="s", linestyle="--")

# Label axes and add a title
plt.xlabel("Iteration")
plt.ylabel("Value")
plt.title("Merit vs. Linearization")
plt.legend()
plt.grid(True)

# Show the plot
plt.show()
