import highspy
import numpy as np

# Highs h
h = highspy.Highs()

# Load a model from MPS file model.mps
filename = 'model.mps'
status = h.readModel(filename)
status = h.run()
model_status = h.getModelStatus()
print('model_status = ', model_status)
print('Model', filename, 'has return status ', h.modelStatusToString(model_status))
option_value = h.getOptionValue('small_matrix_value')
print("option_value = ", option_value)
small_matrix_value = option_value[1]
print("small_matrix_value = ", small_matrix_value)
[status, output_flag] = h.getOptionValue('output_flag')
[status, solver] = h.getOptionValue('solver')
[status, primal_feasibility_tolerance] = h.getOptionValue('primal_feasibility_tolerance')
[status, simplex_update_limit] = h.getOptionValue('simplex_update_limit')
print('output_flag = ', output_flag, '\nsolver = ', solver, '\nprimal_feasibility_tolerance = ', primal_feasibility_tolerance, '\nsimplex_update_limit = ', simplex_update_limit, '\n')

