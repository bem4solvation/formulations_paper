import bem_electrostatics
import numpy as np
import scipy
import pandas as pd
import pickle


bem_electrostatics.bempp.api.enable_console_logging("info")

folder = sys.argv[1]
molecule_file_path = sys.argv[2]
density = int(sys.argv[3])
formulation = sys.argv[4]
form_type = sys.argv[5]


protein = bem_electrostatics.Solute(molecule_file_path,
                                    mesh_density = density,
                                    mesh_probe_radius = 1.4,
                                    mesh_generator = "nanoshaper",
                                    print_times = False,
                                    save_mesh_build_files = False,
                                    force_field = "parse")

print("Formulation type is: "+formulation)
protein.discrete_form_type = form_type
protein.pb_formulation = formulation
protein.operator_assembler = "dense"
protein.pb_formulation_preconditioning = False
protein.gmres_restart = protein.gmres_max_iterations

protein.calculate_potential(rerun_all=True)
protein.calculate_solvation_energy()

if form_type == "strong":
    A = protein.matrices['A'].strong_form().A
else:
    A = protein.matrices['A'].weak_form().A

cond = np.linalg.cond(A)
print("Condition number is: "+str(cond))
eigenvalues = np.linalg.eigvals(A)
real = eigenvalues.real
imag = eigenvalues.imag

gmres = protein.results["solver_iteration_count"]
energy = protein.results["solvation_energy"]
potential_time = protein.timings["time_compute_potential"]

data = {"id": formulation+"-"+form_type+"-none",
        "condition_number": cond,
        "total_calc_time": protein.timings['time_compute_potential']+protein.timings['time_calc_energy'],
        "gmres_iters": protein.results['solver_iteration_count'],
        "solvation_energy": protein.results['solvation_energy'],
        "eigen_real": pd.Series(real),
        "eigen_imag": pd.Series(imag)}

filename = folder+formulation+"-"+form_type+"-none.pkl"
pickle.dump(data, open(filename,"wb"))
