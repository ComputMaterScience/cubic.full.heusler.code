import pandas as pd
import pysubgroup as ps

columns = ['X','Y','Z','ne','Tc2','magsum']

target_name = 'Tc2'

# Load your data into a DataFrame
data = pd.read_csv('data_full_stable2.csv',usecols=columns)

# Define your target and search space
target = ps.NumericTarget(target_name)
search_space = ps.create_selectors(data, ignore=[target_name])

# Define the task
task = ps.SubgroupDiscoveryTask(
    data, 
    target=target, 
    search_space=search_space, 
    result_set_size=10,  # the number of resulting subgroups
    depth=3,  # the maximum size of the description length
    qf=ps.StandardQFNumeric(1)  # the quality function
)

# Run the algorithm
result = ps.DFSNumeric().execute(task)

# Print the results
print(result.to_dataframe())
result.to_dataframe().to_csv('pysubgroup.csv')
