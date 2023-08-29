import os
import json

from matplotlib import pyplot as plt

J_theory = 54.081

with open("results.json", "r", encoding="UTF-8") as fh:
    results = json.load(fh)

contour_ids_to_set_names = {
    int(set_name.split("_")[-1]): set_name
    for set_name in results["J"].keys()
    if "J_NEAR" in set_name.upper()
}
cf_regions = {
    int(set_name.split("_")[-1]): set_name
    for set_name in results["CF"].keys()
    if "J_NEAR" in set_name.upper()
}
contour_ids = sorted(contour_ids_to_set_names.keys())

J_map = {k: v for k, v in results["J"].items()}
cfx_map = {k: v[0] for k, v in results["CF"].items()}

fig, ax = plt.subplots()  # type: plt.Figure, plt.Axes
fig.set_size_inches(3.15, 2.5)
ax.axvline(21, ls=":", c="k")
ax.text(21 / 2, 55, "region 𝔸", ha="center")
ax.text(21 + (57 - 21) / 2, 55, r"region 𝔹", ha="center")
ax.axhline(J_theory, ls="-", c="g", label=r"$G_{\mathrm{Anderson}}$")
ax.plot(contour_ids, [J_map[contour_ids_to_set_names[idx]] for idx in contour_ids], label=r"$G_{\mathrm{Abaqus}}$")
ax.plot(contour_ids, [-cfx_map[cf_regions[idx]] for idx in contour_ids], label=r"$G_{\mathrm{conforce}}$")
ax.set_xlabel("number of contours")
ax.set_ylabel(r"$G$ [kJ/m²]")
ax.set_xlim(left=0, right=contour_ids[-1])
ax.set_ylim(top=56)
ax.set_xticks([0, 4, 21, 57])
ax.legend()
fig.tight_layout()

fig.savefig("test.png", dpi=600)
plt.show()