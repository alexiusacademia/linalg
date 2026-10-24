# linalg

Linear algebra for [Tessel](https://tessel-lang.syncster.dev) programs:
a `Matrix` type, vectors, and the factorizations engineering work keeps
needing (LU, Cholesky, QR, eigenvalues of symmetric matrices, the SVD).
When a system has no answer, such as a singular matrix, you get a failure
that says why, never a list of `NaN`s.

```tessel
import linalg

fn main() {
    // Currents (A) in a circuit: 12 V through 2 Ω, then 4 Ω and 6 Ω side
    // by side. Kirchhoff's laws give three equations in I1, I2 and I3.
    let a = linalg.Matrix.fromRows([
        [1, -1, -1],    // I1 = I2 + I3
        [2, 4, 0],      // 2 I1 + 4 I2 = 12
        [0, 4, -6],     // 4 I2 = 6 I3
    ])
    match linalg.solve(a, [0, 12, 0]) {
        .ok(i) -> print("I1 = {i[0]} A, I2 = {i[1]} A, I3 = {i[2]} A")
        .failure(e) -> print("no answer: {e.message}")
    }
}
```

```text
I1 = 2.7272727272727275 A, I2 = 1.6363636363636362 A, I3 = 1.0909090909090908 A
```

## Using it

Add it to a project (Tessel 0.10 or later):

```sh
tessel add linalg
```

then `import linalg` and use its names as `linalg.Matrix`,
`linalg.solve(…)` and so on (or `import linalg as la`, then `la.solve(…)`).

## Matrices

A `Matrix` has `rows`, `columns` and `values`, which holds the values row
after row. Tessel has no operators for your own types, so the arithmetic
is written as methods.

```tessel
let a = linalg.Matrix.fromRows([[1, 2], [3, 4]])
print(a.times(a).describe())
print(a.timesVector([1, 1]))
print(a.transposed().describe(decimals: 1))
```

```text
[  7.0000  10.0000 ]
[ 15.0000  22.0000 ]
[3.0, 7.0]
[ 1.0  3.0 ]
[ 2.0  4.0 ]
```

| Making one | |
|---|---|
| `Matrix.fromRows([[…], […]])` | From its rows, which must be the same length. |
| `Matrix.fromColumns([[…], […]])` | From its columns. |
| `Matrix.zeros(rows:columns:)` | All zeros. |
| `Matrix.identity(n)` | Ones on the diagonal. |
| `Matrix.diagonal([…])` | The given values on the diagonal. |

| Reading it | |
|---|---|
| `at(i, j)`, `set(i, j, x)` | One value (from 0). `set` changes the matrix. |
| `row(i)`, `column(j)`, `diagonalValues()` | As lists. |
| `toRows()` | As a list of rows. |
| `isSquare`, `isSymmetric(tolerance:)` | |
| `trace`, `largest` | The sum of the diagonal; the largest absolute value. |
| `frobeniusNorm`, `norm1`, `normInfinity` | Norms. |
| `describe(decimals:)` | Text with one line per row and the columns lined up (4 decimals unless you say). |
| `isClose(to:tolerance:)` | Whether another matrix has the same size and values within `tolerance`. |

| Arithmetic | |
|---|---|
| `plus(b)`, `minus(b)`, `scaled(by:)` | |
| `times(b)` | The matrix product. |
| `timesVector(x)` | The matrix times a vector. |
| `transposed()` | Rows turned into columns. |

Vectors are plain `[Float]` lists: `linalg.dot`, `norm`, `add`,
`subtract`, `scale(_:by:)`, `normalized`, `cross` (3 values) and `outer`
(a matrix).

## Solving A x = b

| Function | Returns | |
|---|---|---|
| `solve(a, b)` | `Result<[Float]>` | x for a square a. |
| `inverse(a)` | `Result<Matrix>` | To solve a system, `solve` is faster and more accurate. |
| `determinant(a)` | `Float` | 0 when a is singular. |
| `lu(a)` | `Result<LU>` | Factor once, then `solve(b)` many times. It also has `determinant`, `inverse()`, `solveMatrix(B)`, and `lower`, `upper` and `permutation` (P A = L U). |
| `cholesky(a)` | `Result<Cholesky>` | For symmetric positive definite matrices (stiffness, normal equations, covariances): half the work of LU. Has `lower`, `solve(b)`, `determinant` and `inverse()`. |
| `solveTridiagonal(below:diagonal:above:_:)` | `Result<[Float]>` | For the banded systems that splines and finite differences make, in time proportional to n. |

```tessel
let f = try linalg.lu(linalg.Matrix.fromRows([[2, 1], [1, 3]]))
print(f.solve([3, 5]))
print(f.solve([1, 0]))
print(f.determinant)

match linalg.solve(linalg.Matrix.fromRows([[1, 2], [2, 4]]), [1, 2]) {
    .ok(x) -> print(x)
    .failure(e) -> print(e.message)
}
```

```text
[0.8, 1.4]
[0.6, -0.2]
5.0
the matrix is singular: column 2 has no pivot left after the columns before it are eliminated
```

`lu` pivots on the largest value in each column, and counts a matrix as
singular when a column has nothing left larger than rounding error
(n × 2.2e-16 times the matrix's largest value).

## Least squares

With more equations than unknowns, `leastSquares(a, b)` finds the x that
makes the sum of the squares of a x − b smallest. It uses a QR
factorization (`qr(a)`, with `q`, `r` and `solve(b)`), which keeps twice
the digits that solving the normal equations AᵀA x = Aᵀb would.

```tessel
// A spring's extension (mm) under loads of 0 to 4 kN.
let loads = linalg.Matrix.fromRows([[1, 0], [1, 1], [1, 2], [1, 3], [1, 4]])
let fit = try linalg.leastSquares(loads, [0.0, 1.9, 4.1, 6.0, 8.0])
print("extension = {fit[0].formatted(decimals: 3)} + {fit[1].formatted(decimals: 3)} × load")
```

```text
extension = -0.020 + 2.010 × load
```

It fails when the columns aren't independent; `svd(a).solve(b)` then gives
the shortest of the x's that fit equally well.

## Eigenvalues

`symmetricEigen(a)` gives the eigenvalues of a symmetric matrix, smallest
first, and the eigenvectors as the columns of `vectors` (each of length 1,
signed so its largest value is positive). It uses the Jacobi method: slower
than the QR algorithm on big matrices, but very accurate, which suits the
small symmetric matrices of engineering: stress and inertia tensors,
stiffness and mass matrices.

```tessel
let e = try linalg.symmetricEigen(linalg.Matrix.fromRows([[50, 30], [30, -20]]))
print("principal stresses: {e.values[0].formatted(decimals: 2)} and {e.values[1].formatted(decimals: 2)} MPa")
print(e.vectors.describe())
```

```text
principal stresses: -31.10 and 61.10 MPa
[ -0.3469   0.9379 ]
[  0.9379   0.3469 ]
```

## The singular value decomposition

`svd(a)` factors any matrix as A = U Σ Vᵀ: `u`, `values` (largest first)
and `v`. From it come `rank`, `conditionNumber`, `pseudoInverse()` and
`solve(b)`, the shortest least-squares answer. `linalg.rank(a)`,
`conditionNumber(a)` and `pseudoInverse(a)` are shortcuts.

```tessel
print(linalg.rank(linalg.Matrix.fromRows([[1, 2], [2, 4]])))
print(linalg.conditionNumber(linalg.Matrix.fromRows([[1, 1], [1, 1.0001]])))
```

```text
1
40002.000075025
```

A condition number of about 10^k means a system with that matrix can lose
k digits.

## Mistakes and failures

Sizes that don't fit together, such as adding a 2×3 matrix to a 3×2 one
or reading a value outside the matrix, are a mistake in the code, like an
index outside a list. They stop the program with a message saying which
sizes met:

```text
linalg: can't multiply a 2×3 matrix by a 2×3 one: the first needs as many columns as the second has rows
```

Matrices whose values leave no answer (singular, not symmetric, not
positive definite) give a `Result` failure that says why, for the program
to handle.

## Speed

On an Apple M5 Pro, for a 200×200 matrix: LU about 2 ms, a matrix
product about 11 ms, the SVD and symmetric eigenvalues about 0.1 s.

## Tests

```sh
tessel run tests
```

checks every operation against answers known exactly, puts every
factorization back together, runs every example in the documentation, and
exits with 1 if any check fails.

## What's where

| File | |
|---|---|
| `matrix.tsl` | `Matrix`. |
| `vector.tsl` | `dot`, `norm`, `add`, `subtract`, `scale`, `normalized`, `cross`, `outer`. |
| `solve.tsl` | `solve`, `inverse`, `determinant`, `lu` and `LU`, `cholesky` and `Cholesky`, `solveTridiagonal`. |
| `qr.tsl` | `qr` and `QR`, `leastSquares`. |
| `eigen.tsl` | `symmetricEigen` and `Eigen`. |
| `svd.tsl` | `svd` and `SVD`, `rank`, `conditionNumber`, `pseudoInverse`. |
| `tests/main.tsl` | The tests. |
| `tessel.toml` | The library's name, version, description and licence. |

Every function's comment in the code says what it does, its parameters,
what it returns, and gives an example.

## Licence

MIT: see [LICENSE](LICENSE).
