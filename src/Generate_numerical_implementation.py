# Generation of Numeric Implementations

import subprocess
import os

import Auxiliary_functions as sy_hfkt
import ele_def as ele_def

#########################
# Generate implementation
#########################

# Name of elements for which the implementation should be generated
# The properties of this types have to be defined in ele_def.py
types = (
    'CPE4', 'CPE4R', 'CPE8', 'CPE8R', 'CPE3', 'CPE6',
    # 'C3D6',
    'C3D8', 'C3D8R', 'C3D20', 'C3D20R', 'C3D4', 'C3D10', 'C3D10R'
)

# Define file where the c- implementations should be stored
f_c_functions = open("Conf_Forces.c", "w")
# Define path where the Python wrappers should be stored
f_python_functions = open("Conf_Forces_py.py", "w")

# Write c-headers
f_c_functions.write(sy_hfkt.gen_C_header())
# Write python-headers
f_python_functions.write(sy_hfkt.gen_Python_Wrapper_Header("Conf_Forces"))

for i in range(len(types)):
    poly_power = ele_def.poly_power[types[i]]
    bild_points = ele_def.bild_points[types[i]]
    int_points = ele_def.int_points[types[i]]
    int_weights = ele_def.int_weights[types[i]]

    # Build Implementation
    print("Implementation for " + types[i])
    C_Force_dyn = sy_hfkt.gen_Configurational_Forces_Dynamic(poly_power, bild_points, int_points, int_weights,
                                                             types[i] + "_dynamic")
    C_Force_stat_mbf = sy_hfkt.gen_Configurational_Forces_Static(poly_power, bild_points, int_points, int_weights,
                                                                 types[i] + "_static", method='mbf')
    C_Force_stat_dbf = sy_hfkt.gen_Configurational_Forces_Static(poly_power, bild_points, int_points, int_weights,
                                                                 types[i] + "_static", method='dbf')

    # Write C-Implementation
    f_c_functions.write(C_Force_dyn)
    f_c_functions.write(C_Force_stat_mbf)
    f_c_functions.write(C_Force_stat_dbf)

    # Write Python wrapper
    f_python_functions.write(
        sy_hfkt.gen_Python_Wrapper_static(types[i] + "_static", bild_points.shape[0], int_points.shape[0]))
    f_python_functions.write(
        sy_hfkt.gen_Python_Wrapper_dynamic(types[i] + "_dynamic", bild_points.shape[0], int_points.shape[0]))

f_c_functions.close()
f_python_functions.close()

#############################################
# Compilation of the generated implementation
#############################################
# It is recommended to use clang for compilation.
# On Windows the Build Tools for Visual Studio are also a required.
#
# For a optimized build on a specific cpu use "-march=native".

if os.name == 'nt':
    # Windows
    subprocess.call("clang -shared -O3 -o Conf_Forces.dll Conf_Forces.c", shell=True)
else:
    # Unix based operating systems
    subprocess.call("clang -shared -O3 -fPIC -o  Conf_Forces.so Conf_Forces.c", shell=True)
