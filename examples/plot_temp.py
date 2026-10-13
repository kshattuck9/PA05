import pandas as pd
import matplotlib.pyplot as plt

data = pd.read_csv("build/composite_temperature.csv")
temperature = data.pivot(index="y", columns="x", values="T")
plt.figure(figsize=(9, 5.5))
contour = plt.contourf(
    temperature.columns,
    temperature.index,
    temperature.values,
    levels=30,
    cmap="inferno"
)
plt.colorbar(contour, label="Temperature")
plt.xlabel("x")
plt.ylabel("y")
plt.title("Composite plate temperature (127 × 95)")
plt.tight_layout()
plt.savefig("build/composite_temperature.png", dpi=200)
plt.show()