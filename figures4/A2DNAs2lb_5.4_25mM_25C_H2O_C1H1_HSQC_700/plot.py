import matplotlib.pyplot as plt
import numpy as np
import nmrglue as ng
import math

contour_start = math.pow(10, 6)
contour_num = 20
contour_factor = 1.4

cl = contour_start * contour_factor ** np.arange(contour_num)

dic, data = ng.pipe.read('test.ft2')
print type(data)
print type(dic)
x_proton = ng.pipe.make_uc(dic, data, dim=1)
x_proton_ppm = x_proton.ppm_scale()
x_proton_ppm_start, x_proton_ppm_end = x_proton.ppm_limits()
y_het = ng.pipe.make_uc(dic, data, dim=0)
y_het_ppm = y_het.ppm_scale()
y_het_ppm_start, y_het_ppm_end = y_het.ppm_limits()

fig, ax = plt.subplots()
ax.contour(data.transpose(), cl, color='k', extent=(y_het_ppm_start, y_het_ppm_end, x_proton_ppm_start, x_proton_ppm_end))
ax.set_xlabel('Proton')
ax.set_ylabel('Carbon')
ax.set_ylim([x_proton_ppm_start, x_proton_ppm_end])
ax.set_xlim([y_het_ppm_start, y_het_ppm_end])
plt.show()
