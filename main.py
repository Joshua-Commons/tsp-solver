import csv
import time
import ctypes
import numpy as np


class TSP_Solver:
    '''
        Uses a heuristic to find a valid solution to a metric travelling salesman problem
        (TSP), for network size 10-50, within 30 seconds.

        Usage
        -----
        1. Initialise object with `solver = TSP_Solver()`
        2. Load a distance matrix by either:
            1. Generating a distance matrix by `matrix = solver.generate_matrix(n)` where 
            `n` is the desired size of matrix 10-50 inclusive.
            2. Load a distance matrix from a `.csv` file by `matrix = solver.load_matrix(path)` where
            `path` is the path to your file. Ideally, under the `data` folder.
        3. Validate matrix form (done automatically when `solver.solve(matrix)` is called.)
        4. Find solution to metric TSP with `solver.solve(matrix)` where `matrix` is the
        distance matrix.

        Outputs
        -------
        Prints the found solution to the metric TSP, and its cost.
    '''


    def __init__(self) -> None:
        '''
            Initialise the `TSP_Solver` object. 
            
            Specifically load the C++ library containing the algorithms.
        '''

        # load algorithms
        self.algorithms = ctypes.CDLL('./libs/algorithms.so')

        # define structure to match C++ result struct
        class Result(ctypes.Structure):
            _fields_ = [("cost", ctypes.c_double), ("path", ctypes.POINTER(ctypes.c_int))]

        # define tsp solver types
        self.algorithms.solve_tsp.argtypes = [ctypes.c_int, ctypes.POINTER(ctypes.c_double)]
        self.algorithms.solve_tsp.restype = Result


    def solve(self, matrix: np.ndarray) -> tuple[float, list[int]]:
        '''
            Finds a valid solution to the metric TSP provided.

            Parameters
            ----------
            matrix : np.ndarray
                Distance matrix satisfying the triangle inequality.

            Returns
            -------
            cost : float
                The cost of the found TSP solution.
            path : list[int]
                List of ordered indicies which consist of the TSP tour. 
        '''

        # check for mismatching inputs
        self.validate_distance_matrix(matrix)

        # check that distance matrix satisfies triangle inequality
        if not self.check_triangle_inequality(matrix):
            raise ValueError(f"Error: Distance matrix does not satisfy the triangle inequality.")

        # get size of matrix
        n = matrix.shape[0]

        # flatten the distance matrix for passing to C++
        flat_matrix = matrix.flatten().astype(np.float64)
        flat_matrix_ptr = flat_matrix.ctypes.data_as(ctypes.POINTER(ctypes.c_double))

        # call algorithm
        t1 = time.time()
        result = self.algorithms.solve_tsp(n, flat_matrix_ptr)
        t2 = time.time()

        # extract the cost and path
        cost = result.cost
        path = [result.path[i] for i in range(n)]

        # free the dynamically allocated memory in C++
        self.algorithms.free(result.path)

        # outputting results
        print(f"Cost: {cost:.2f}")
        print(f"Path: {path}")
        print(f"Time Taken: {(t2-t1):.2f}s")

        return cost, path


    def check_triangle_inequality(self, matrix: np.ndarray) -> bool:
        '''
            Ensures that distance matrix satisfies triangle inequality.

            Parameters
            ----------
            matrix : np.ndarray
                Distance matrix.

            Returns
            -------
            cond : Boolean
                Returns `True` if triangle inequality is satisfied, False otherwise.
        '''
        
        # get the size of the matrix
        n = matrix.shape[0]

        # iterate through each combination of 3 nodes
        for i in range(n):
            for j in range(i + 1, n):
                for k in range(n):

                    # skip if there is overlap
                    if k == i or k == j:
                        continue
                    
                    # if triangle inequality is broken return false
                    if matrix[i, j] > matrix[i, k] + matrix[k, j]:
                        return False

        # triangle inequality is satisfied
        return True


    def validate_distance_matrix(self, matrix: np.ndarray) -> bool:
        '''
            Validates that the distance matrix is appropriate to be solved.

            Parameters
            ----------
            matrix : np.ndarray
                Distance matrix.

            Returns : Boolean
                If distance matrix can be safely passed into C++ algorithm and solved.
        '''
        
        # check that it is np.ndarray
        if not isinstance(matrix, np.ndarray):
            raise ValueError(f"Error: Distance matrix is not a np.ndarray. Yours is {type(matrix)}.")
        
        # check that matrix is 2D
        if len(matrix.shape) > 2:
            raise ValueError(f"Error: Distance matrix is not 2D. Yours is {len(matrix.shape)}.")
        
        # check that it is square
        if matrix.shape[0] != matrix.shape[1]:
            raise ValueError(f"Error: Distance matrix is not square. Your dimensions are {matrix.shape}.")

        # check that it is size 10-50
        if matrix.shape[0] < 10 or matrix.shape[0] > 50:
            raise ValueError(f"Error: Distance matrix must be from size 10-50. Your size is {matrix.shape[0]}.")

        # check that it contains only positive floats
        for i in range(matrix.shape[0]):
            for j in range(matrix.shape[1]):
                if not isinstance(matrix[i][j], float):
                    raise ValueError(f"Error: Distance matrix must contain only floats. At ({i}, {j}) the type is {type(matrix[i][j])}.")
                if i != j:
                        if matrix[i][j] <= 0:
                            raise ValueError(f"Error: Distance matrix must contain positive off diagonals. At ({i}, {j}) the value is {matrix[i][j]}.")
                
        # check that it does not contain negatives
        if np.amin(matrix) < 0:
            raise ValueError(f"Error: Distance matrix must not contain negatives. Yours contains the value {np.amin(matrix)}.")

        # check that main diagonal is zeros
        for i in range(matrix.shape[0]):
            if matrix[i][i] != 0:
                raise ValueError(f"Error: Distance matrix must have zero values on the main diagonal. At ({i}, {i}) the value is {matrix[i][i]}.")

        return True

    
    def generate_distance_matrix(self, n: int) -> np.ndarray:
        '''
            Generates triangle inequality compliant matrix suitable for metric TSP solvers.
            Note that it generates points in Euclidean space, a subset of metric space.

            Parameters
            ----------
                n : int
                    Size of matrix to create.

            Returns
            -------
                matrix : np.ndarray
                    Distance matrix solving triangle inequality.
        '''
        
        # error checking
        if not isinstance(n, int):
            raise ValueError(f"Error: Size must be an integer. You put {n}.")
        if n > 50 or n < 10:
            raise ValueError(f"Error: Size must be between 10 and 50. You put {n}.")

        # generate random points in 2D Euclidean space
        points = np.random.rand(n, 2) * 1000

        # calculate pairwise distances
        matrix = np.linalg.norm(points[:, np.newaxis] - points[np.newaxis, :], axis=2)

        # return distance matrix
        return matrix


    def read_distance_matrix(self, path: str) -> np.ndarray:
        '''
            Reads distance matrix stored as a `.csv` file.

            Parameters
            ----------
            path : str
                String containing the path to the csv file relative to the repo.

            Returns
            -------
            matrix np.ndarray
                Distance matrix
        '''

        try:
            # read file line by line
            with open(path, 'r') as file:
                lines = file.readlines()

            # initialise distance matrix as list of lists
            matrix_builder = []
            
            # parse each line into a list of floats
            for line in lines:
                
                try:
                    row = [float(value) for value in line.strip().split(',')]
                    matrix_builder.append(row)
                
                except ValueError:
                    raise ValueError(f"Error: All values must be floats. Line is {line}.")

            # convert the list of lists to a np array
            matrix = np.array(matrix_builder, dtype=float)

            return matrix

        except FileNotFoundError:
            raise FileNotFoundError(f"Error: The file at {path} was not found.")


    def write_distance_matrix(self, matrix: np.ndarray, path: str) -> None:
        '''
            Writes distance matrix to as a `.csv` file.

            Parameters
            ----------
            matrix np.ndarray
                Distance matrix
            path : str
                String containing the path to the csv file relative to the repo.
        '''

        # writing distance matrix to csv
        with open(path, mode='w', newline='') as file:
            writer = csv.writer(file)
            writer.writerows(matrix)
        

# sample workflow
if __name__ == '__main__':
    
    # initialise soler
    solver = TSP_Solver()

    # initialise distance matrix
    read = False
    if read:
        matrix = solver.read_distance_matrix('data/matrix_1.csv')
    else:
        matrix = solver.generate_distance_matrix(50)

    # solve matrix and time function
    cost, path = solver.solve(matrix)