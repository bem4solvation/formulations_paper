import os.path
import pbj
import resource
import csv
import sys
    
csv_file = sys.argv[1]
solute_file  = sys.argv[2]
density = float(sys.argv[3])

formulation = sys.argv[4]
precon_type = sys.argv[5]

expansion_order = int(sys.argv[6])
ncrit = int(sys.argv[7])

quadrature_points = int(sys.argv[8])


csv_columns = ["molecule",
               "density",
               "mesh_elements",
               "mesh_vertices",
               "mem_usage",
               "total_time",
               "matrix_initialisation",
               "matrix_assembly",
               "rhs_construction",
               "preconditioning_time",
               "gmres_time",
               "potential_time",
               "energy_time",
               "gmres_iter",
               "energy",
               "fmm_order",
               "fmm_ncrit",
               "assembler",
               "quadrature_points",
               "formulation",
               "precon_bool",
               "precon_type"]


pbj.electrostatics.solute.bempp.api.enable_console_logging('info')

pbj.electrostatics.solute.bempp.api.GLOBAL_PARAMETERS.fmm.expansion_order = expansion_order
pbj.electrostatics.solute.bempp.api.GLOBAL_PARAMETERS.fmm.ncrit = ncrit

pbj.electrostatics.solute.bempp.api.GLOBAL_PARAMETERS.quadrature.regular = quadrature_points

protein = pbj.Solute(solute_file,
                     mesh_density = density,
                     mesh_probe_radius = 1.4,
                     mesh_generator = "nanoshaper",
                     print_times = False,
                     force_field = "parse")


protein.ep_in = 4.0 
protein.ep_ex = 80.0
protein.kappa = 0.125

protein.gmres_restart = protein.gmres_max_iterations

protein.pb_formulation = formulation


if precon_type == "none":
    protein.pb_formulation_preconditioning = False
else:
    protein.pb_formulation_preconditioning = True
    protein.pb_formulation_preconditioning_type = precon_type

    
protein.operator_assembler = "fmm"
protein.rhs_constructor = "fmm"

#protein.initialise_matrices()
#protein.assemble_matrices()
#protein.initialise_rhs()
#protein.apply_preconditioning()

protein.calculate_potential()
protein.calculate_solvation_energy()

#protein.calculate_solvation_energy()


result = {"molecule" : protein.solute_name,
          "density" : density,
          "mesh_elements" : protein.mesh.number_of_elements,
          "mesh_vertices" : protein.mesh.number_of_vertices,
          "mem_usage" : resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,
          "total_time" : protein.timings['time_compute_potential']+protein.timings['time_calc_energy'],
          "matrix_initialisation" : protein.timings['time_matrix_initialisation'],
          "matrix_assembly" : protein.timings['time_matrix_assembly'],
          "rhs_construction" : protein.timings['time_rhs_initialisation'],
          "preconditioning_time" : protein.timings['time_preconditioning'],
          "gmres_time" : protein.timings['time_gmres'],
          "potential_time" : protein.timings['time_compute_potential'],
          "energy_time" : protein.timings['time_calc_energy'],
          "gmres_iter" : protein.results['solver_iteration_count'],
          "energy" : protein.results['solvation_energy'],
          "fmm_order" : expansion_order,
          "fmm_ncrit" : ncrit,
          "assembler" : protein.operator_assembler,
          "quadrature_points" : pbj.electrostatics.solute.bempp.api.GLOBAL_PARAMETERS.quadrature.regular,
          "formulation" : formulation,
          "precon_bool" : protein.pb_formulation_preconditioning,
          "precon_type": precon_type}

print(result)

file_exists = os.path.isfile(csv_file)
    
try:
    with open(csv_file, 'a') as csvfile:
        writer = csv.DictWriter(csvfile, fieldnames=csv_columns)
        if not file_exists:
            writer.writeheader()  # file doesn't exist yet, write a header
        writer.writerow(result)
except IOError:
    print("I/O error")
