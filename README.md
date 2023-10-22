# DXZ

dxz is a Python package that implements the dxz algorithm that is described in the paper [Dancing with Decision Diagrams: A Combined Approach to Exact Cover](https://ojs.aaai.org/index.php/AAAI/article/view/10662). The algorithms returns all solutions to a given exact cover problem.

## Installation

Use the package manager [pip](https://pip.pypa.io/en/stable/) to install dxz.

```bash
pip install dxz
```

## The Algorithm

dxz is an algorithm for solving an exact cover problem. It build a [zdd](https://en.wikipedia.org/wiki/Zero-suppressed_decision_diagram) that represents all solutions for a given problem. This module converts the zdd to to a list of lists of ints. Each list in the list represents one solution.

For more details about dxz read the paper [Dancing with Decision Diagrams: A Combined Approach to Exact Cover](https://ojs.aaai.org/index.php/AAAI/article/view/10662).

## Usage

```python
import dxz

rows = 8
columns = 5

# A flatten matrix of boolean values. Row 0 is 1,0,0,1,1 row 2 is 0,1,1,0,0 and so on.
matrix = [
    1,0,0,1,1, # 0
    0,1,1,0,0, # 1
    0,0,0,1,1, # 2
    0,1,0,1,0, # 3
    1,0,1,0,0, # 4
    0,0,0,0,1, # 5
    0,0,0,1,0, # 6
    1,0,0,0,0, # 7
]

# Get all solutions to the exact cover problem.
solutions = dxz.dxz_solve(rows, columns, matrix, 0)

print(solutions)
# [
#   [1,0],
#   [1,7,2],
#   [1,7,6,5],
#   [3,4,5]
# ]
# One solutions is row 1 and 2 another solution is row 1, 7 and 2.

# Get only the number of solutions
number_of_solutions = dxz.dxz_count(rows, columns, matrix)


# Get only 3 solutions.
solutions = dxz.dxz_solve(rows, columns, matrix, 3)

```

## Contributing

Pull requests are welcome. For major changes, please open an issue first
to discuss what you would like to change.

## License

[MIT](https://choosealicense.com/licenses/mit/)
