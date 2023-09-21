import multiprocessing

# Define the function you want to execute multiple times
def my_function(arg1, arg2):
    # Your function logic here
    result = arg1 + arg2
    return result

# Function to run another function multiple times with given arguments and collect results
def run_function_multiple_times(func, args_list, num_processes):
    # Create a pool of processes
    pool = multiprocessing.Pool(processes=num_processes)

    # Use the pool to map the function to the argument list
    results = pool.starmap(func, args_list)

    # Close the pool to prevent any more tasks from being submitted
    pool.close()

    # Wait for all processes to complete
    pool.join()

    # Print or process the collected results
    for i, result in enumerate(results):
        print(f"Result {i}: {result}")

if __name__ == "__main__":
    num_processes = 4  # Change this to the desired number of processes
    arguments_list = [(1, 2), (3, 4), (5, 6), (7, 8)]  # List of argument tuples

    # Run my_function multiple times with the provided arguments and collect results
    run_function_multiple_times(my_function, arguments_list, num_processes)
