import os.path
import bem_electrostatics
import resource
import csv
import time
    
csv_file = sys.argv[1]
solute_file  = sys.argv[2]
density = float(sys.argv[3])

formulation = sys.argv[4]
discrete_form = sys.argv[5]
precon_type = sys.argv[6]

expansion_order_main = int(sys.argv[7])
ncrit_main = int(sys.argv[8])
reg_quadrature_points_main = int(sys.argv[9])
sing_quadrature_points_main = int(sys.argv[10])

expansion_order_calderon = int(sys.argv[11])
ncrit_calderon = int(sys.argv[12])
reg_quadrature_points_calderon = int(sys.argv[13])
sing_quadrature_points_calderon = int(sys.argv[14])


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
               "quadrature_points_regular",
               "quadrature_points_singular",
               "discrete_form",
               "formulation",
               "precon_bool",
               "precon_type",
               "fmm_order_calderon",
               "fmm_ncrit_calderon",
               "quadrature_points_calderon_regular",
               "quadrature_points_calderon_singular"]


bem_electrostatics.bempp.api.enable_console_logging('info')

bem_electrostatics.bempp.api.GLOBAL_PARAMETERS.fmm.expansion_order = expansion_order_main
bem_electrostatics.bempp.api.GLOBAL_PARAMETERS.fmm.ncrit = ncrit_main

bem_electrostatics.bempp.api.GLOBAL_PARAMETERS.quadrature.regular = reg_quadrature_points_main
bem_electrostatics.bempp.api.GLOBAL_PARAMETERS.quadrature.singular = sing_quadrature_points_main

protein = bem_electrostatics.Solute(solute_file,
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
protein.discrete_form_type = discrete_form


protein.pb_formulation_preconditioning = True
protein.pb_formulation_preconditioning_type = precon_type

    
protein.operator_assembler = "fmm"
protein.rhs_constructor = "fmm"

protein.initialise_matrices()
protein.assemble_matrices()
protein.initialise_rhs()

precon_start = time.time()

protein.matrices["A"].strong_form()

bem_electrostatics.bempp.api.GLOBAL_PARAMETERS.fmm.expansion_order = expansion_order_calderon
bem_electrostatics.bempp.api.GLOBAL_PARAMETERS.fmm.ncrit = ncrit_calderon

bem_electrostatics.bempp.api.GLOBAL_PARAMETERS.quadrature.regular = reg_quadrature_points_calderon
bem_electrostatics.bempp.api.GLOBAL_PARAMETERS.quadrature.singular = sing_quadrature_points_calderon

protein.matrices["preconditioning_matrix"] = bem_electrostatics.pb_formulation.formulations.lhs.first_kind_external(protein.dirichl_space, protein.neumann_space, protein.ep_in, protein.ep_ex, protein.kappa, protein.operator_assembler)

protein.matrices["preconditioning_matrix"].strong_form()
protein.matrices["A_discrete"] = protein.matrices["preconditioning_matrix"].strong_form() * protein.matrices["A"].strong_form()


bem_electrostatics.bempp.api.GLOBAL_PARAMETERS.fmm.expansion_order = expansion_order_main
bem_electrostatics.bempp.api.GLOBAL_PARAMETERS.fmm.ncrit = ncrit_main

bem_electrostatics.bempp.api.GLOBAL_PARAMETERS.quadrature.regular = reg_quadrature_points_main
bem_electrostatics.bempp.api.GLOBAL_PARAMETERS.quadrature.singular = sing_quadrature_points_main


protein.rhs["rhs_final"] = [protein.rhs["rhs_1"], protein.rhs["rhs_2"]]
protein.rhs["rhs_discrete"] = protein.matrices["preconditioning_matrix"].strong_form() * bem_electrostatics.solute.rhs_to_discrete_form(protein.rhs["rhs_final"], protein.discrete_form_type, protein.matrices["A"])


print(bem_electrostatics.bempp.api.fmm.fmm_assembler._FMM_CACHE)

precon_time = time.time() - precon_start
protein.timings['time_preconditioning'] = precon_time

protein.calculate_potential()
protein.calculate_solvation_energy()


result = {"molecule" : protein.solute_name,
          "density" : density,
          "mesh_elements" : protein.mesh.number_of_elements,
          "mesh_vertices" : protein.mesh.number_of_vertices,
          "mem_usage" : resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,
          "total_time" : protein.timings['time_compute_potential']+protein.timings['time_calc_energy'],
          "matrix_initialisation" : protein.timings['time_matrix_initialisation'],
          "matrix_assembly" : protein.timings['time_matrix_assembly'],
          "rhs_construction" : protein.timings['time_rhs_construction'],
          "preconditioning_time" : protein.timings['time_preconditioning'],
          "gmres_time" : protein.timings['time_gmres'],
          "potential_time" : protein.timings['time_compute_potential'],
          "energy_time" : protein.timings['time_calc_energy'],
          "gmres_iter" : protein.results['solver_iteration_count'],
          "energy" : protein.results['solvation_energy'],
          "fmm_order" : expansion_order_main,
          "fmm_ncrit" : ncrit_main,
          "assembler" : protein.operator_assembler,
          "quadrature_points_regular" : reg_quadrature_points_main,
          "quadrature_points_singular" : sing_quadrature_points_main,
          "discrete_form" : discrete_form,
          "formulation" : formulation,
          "precon_bool" : protein.pb_formulation_preconditioning,
          "precon_type": precon_type,
          "fmm_order_calderon" : expansion_order_calderon,
          "fmm_ncrit_calderon" : ncrit_calderon,
          "quadrature_points_calderon_regular" : reg_quadrature_points_calderon,
          "quadrature_points_calderon_singular" : sing_quadrature_points_calderon}

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
