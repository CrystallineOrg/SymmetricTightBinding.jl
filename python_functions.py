import numpy as np
import sympy as sp
from ortools.sat.python import cp_model

fixedGamma2 = (
    [1]
    + list(np.arange(3, 5 + 1, 1))
    + list(np.arange(6, 9 + 1, 1))
    + list(np.arange(25, 46 + 1, 1))
    + list(np.arange(75, 80 + 1, 1))
    + list(np.arange(143, 146 + 1, 1))
    + list(np.arange(156, 161 + 1, 1))
    + list(np.arange(168, 173 + 1, 1))
    + list(np.arange(183, 186 + 1, 1))
)

fixedGamma4 = list(np.arange(99, 110 + 1, 1))


def ir2vec(v: sp.core.add.Add, irreps: list[sp.core.symbol.Symbol]) -> np.ndarray:
    """Convert a sum of the irreps of the space group into a vector of the multiplicities following the order specified in the input `irreps`.

    Args:
        v (sympy.core.add.Add): sum of irreps of the space group.
        irreps (list): list of irreps indicating their order and their names.

    Returns:
        numpy.ndarray: vector of multiplicities of the irreps following the order specified in `irreps`.
    """  # noqa: E501
    y = np.array([v.coeff(irrep) for irrep in irreps])
    return y


def ebr2vec(v: sp.core.add.Add, ebrs: list[sp.core.symbol.Symbol]) -> np.ndarray:
    """Convert a sum of EBRs into a vector of multiplicities of the ebrs of the space group following the order stipulated in `ebrs`.

    Args:
        v (sympy.core.add.Add): sum of EBRs of the space group.
        ebrs (list): list of EBRs indicating their names.
        EBRs (numpy.ndarray): matrix form by columns with the vectors of multiplicities of the irreps for each EBR.

    Returns:
        numpy.ndarray: vector of multiplicities of the ebrs following the order specified in `ebrs`.
    """  # noqa: E501
    vec = np.array([v.coeff(ebr) for ebr in ebrs])
    return vec


def ebr2irvec(
    v: sp.core.add.Add, ebrs: list[sp.core.symbol.Symbol], EBR: np.ndarray
) -> np.ndarray:
    """Convert a sum of EBRs into a vector of multiplicities of the irreps of the space group following the order stipulated in `irreps`.

    Args:
        v (sympy.core.add.Add): sum of EBRs of the space group.
        ebrs (list): list of EBRs indicating their names.
        EBRs (numpy.ndarray): matrix form by columns with the vectors of multiplicities of the irreps for each EBR.

    Returns:
        numpy.ndarray: vector of multiplicities of the irreps following the order specified in `irreps`.
    """  # noqa: E501
    vec = EBR @ np.array([v.coeff(ebr) for ebr in ebrs])
    return vec


def ebr2ir(
    v: sp.core.add.Add,
    irreps: list[sp.core.symbol.Symbol],
    ebrs: list[sp.core.symbol.Symbol],
    EBR: np.ndarray,
) -> sp.core.add.Add:
    """Convert a sum of EBRs into a sum of irreps of the space group.

    Args:
        v (sympy.core.add.Add): sum of EBRs of the space group.
        irreps (list): list of irreps indicating their order and their names.
        ebrs (list): list of EBRs indicating their names.
        EBR (numpy.ndarray): matrix form by columns with the vectors of multiplicities of the irreps for each EBR.

    Returns:
        sympy.core.add.Add: sum of irreps of the space group.
    """  # noqa: E501
    ir = irreps @ ebr2irvec(v, ebrs, EBR)
    return ir


class varArraySolutionObtainer(cp_model.CpSolverSolutionCallback):
    """A class representing a solution callback for the solver.

    Args:
        cp_model (list): list of variables defined in your model.
    """

    def __init__(self, variables: cp_model.IntVar) -> None:
        """Initiate a solution initial state.

        Args:
            variables (ortools.sat.python.cp_model.IntVar): values of the variables.
        """
        cp_model.CpSolverSolutionCallback.__init__(self)
        self.__variables = variables
        self.__solutionCount = 0
        self.__x = []

    def on_solution_callback(self) -> None:
        """Save the current solution in a list and add `1` to the solution count."""
        self.__solutionCount += 1
        self.__x.append([self.Value(v) for v in self.__variables])  # type:ignore

    def solutionCount(self) -> int:
        """Returns the total number of solutions.

        Returns:
            int: Total number of solutions
        """
        return self.__solutionCount

    def solutions(self) -> np.ndarray:
        """Returns all solutions of the model in a matrix form, where each row is a particular solution.

        Returns:
            numpy.ndarray: Numpy array which rows represent a solution to the model.
        """
        return np.array(self.__x)


def searchForAllLongModes(dim: list[int], t: int) -> varArraySolutionObtainer:
    """Construct a model to search for all longitudinal modes up to total dimension `t`.

    Args:
        dim (list): list indicating the dimensions of the EBRs of the space group in the
        order specified by `ebrs`.
        t (int): positive integer indicating up to what dimension of the longitudinal
        modes search.

    Returns:
        varArraySolutionObtainer: returns a container with all solutions in matrix form,
        where each row is a particular solution, and the number of solutions.
    """

    # Creates the model and solver.
    model = cp_model.CpModel()
    solver = cp_model.CpSolver()

    N = len(dim)

    # Creates the variables
    x = [model.NewIntVar(0, 1000, f"x{i}") for i in range(N)]
    solutionObtainer = varArraySolutionObtainer(x)  # type:ignore

    # Create the constraints.
    model.Add(np.dot(x, dim) == t)  # type:ignore

    # Solve.
    solver.SearchForAllSolutions(model, solutionObtainer)

    return solutionObtainer


def vec2ebr(v: np.ndarray, ebrs: list[sp.core.add.Add]) -> list[sp.core.add.Add]:
    """Transform a vector of multiplicities of EBRs into a sum of EBRs.

    Args:
        v (numpy.ndarray): matrix which each row is a vector of multiplicities of EBRs. Each vector must be given in the order
        stipulated by `ebrs`.
        ebrs (list): list of variables defining the posible EBRs of the space group.

    Returns:
        list: list whose elements are the sum of EBRs related to each row of the matrix `v`.
    """
    y = v @ ebrs
    return y.tolist()


def searchForAll_EBRs(
    v: np.ndarray, dim: list[int], EBR: np.ndarray, long_modes
) -> list[np.ndarray]:
    """Looks for all posible linear combinations of EBRs which, after subtracting all
    possible longitudinal modes until dimension `t`, could represent the vector of
    multiplicities `v`.

    Args:
        v (numpy.ndarray): vector of multiplicities of the irreps following the order
        specified in `irreps`.
        dim (list): list indicating the dimensions of the EBRs of the space group in the
        order specified by `ebrs`.
        EBR (numpy.ndarray): matrix form by columns with the vectors of multiplicities
        of the irreps for each EBR.

    Returns:
        list: list of vectors of multiplicities of EBRs which represent a solution to the problem.
    """
    possibleEBRs = []

    N_ebrs = len(dim)
    N_irr = len(v)

    for i in long_modes:
        model = cp_model.CpModel()
        solver = cp_model.CpSolver()
        x = [model.NewIntVar(0, 1000, f"x{i}") for i in range(N_ebrs)]
        solutionObtainer = varArraySolutionObtainer(x)  # type:ignore

        for j in range(N_irr):
            model.Add(EBR[j] @ (x - i) == v[j])

        solver.SearchForAllSolutions(model, solutionObtainer)

        possibleEBRs.append(solutionObtainer.solutions())

    return possibleEBRs


def phys(
    x: sp.core.add.Add,
    v: np.ndarray,
    ebrs: list[sp.core.symbol.Symbol],
    EBR: np.ndarray,
    SG: int,
    irreps: list[sp.core.symbol.Symbol],
) -> int:
    """Determines if a linear combinations of EBRs `x` properly represents the vector of
    multiplicities `v`. We define this situation as *physical*.

    Args:
        x (sp.core.add.Add): linear combinations of EBRs.
        v (numpy.ndarray): vector of multiplicities of the irreps following the order
        specified in `irreps`.
        ebrs (list): list of variables defining the EBRs of the space group.
        EBR (numpy.ndarray): matrix form by columns with the vectors of multiplicities
        of the irreps for each EBR.
        SG (int): number of the Space Group

    Returns:
        int: `0` if it is physical, otherwise in any other case.
    """

    if SG in fixedGamma2:
        phys = (
            ebr2ir(x, ebrs=ebrs, EBR=EBR, irreps=irreps).coeff(irreps[0]) > 0
            or ebr2ir(x, ebrs=ebrs, EBR=EBR, irreps=irreps).coeff(irreps[1]) > 0
        )
    elif SG in fixedGamma4:
        phys = (
            ebr2ir(x, ebrs=ebrs, EBR=EBR, irreps=irreps).coeff(irreps[0]) > 0
            or ebr2ir(x, ebrs=ebrs, EBR=EBR, irreps=irreps).coeff(irreps[3]) > 0
        )
    else:
        phys = (
            np.sum(ebr2irvec(x, ebrs, EBR) - v - np.abs(ebr2irvec(x, ebrs, EBR) - v))
            / 2
            == 0
        )

    return phys


def showAllResults(
    v: np.ndarray,
    t: int,
    dim: list[int],
    EBRs: np.ndarray,
    N_gamma: int,
    ebrs: list[sp.core.symbol.Symbol],
    irreps: list[sp.core.symbol.Symbol],
    SG: int,
) -> list[
    tuple[list[sp.core.add.Add], sp.core.add.Add, list[int], list[sp.core.add.Add]]
]:
    """Obtains all linear combinations of EBRs which represents the symmetry vector `v`.
    In addition, it computes any necessary longitudinal modes to obtain such
    decompositions in EBRs, if they are physical or not, and, the ill-defined content at
    the Gamma point and zero frequency.

    Args:
        v (numpy.ndarray): vector of multiplicities of the irreps following the order
        specified in `irreps`.
        t (int): positive integer indicating up to what dimension of the longitudinal
        modes search.
        dim (list): list indicating the dimensions of the EBRs of the space group in the
        order specified by `ebrs`.
        EBRs (numpy.ndarray): matrix form by columns with the vectors of multiplicities
        of the irreps for each EBR.
        N_gamma (int): number of irreps at Gamma.
        ebrs (list): list of variables defining the EBRs of the space group.
        irreps (list): list of irreps of the space group indicating their order and their names.
        SG (int): number of the Space Group

    Returns:
        list: list of tuples which first component is a list of possible lineal
        combinations of EBRs, second component is the EBRs describing the longitudinal
        modes, third component is a list of numbers indicating if the EBRs
        decompositions are physical or not, and, fourth component is the ill-defined
        symmetry content at Gamma and zero frequency.
    """
    rv = np.delete(v, range(N_gamma))  # type: ignore
    rEBR = np.delete(EBRs, range(N_gamma), axis=0)  # type: ignore

    long_modes = searchForAllLongModes(dim, t).solutions()

    all_EBRs = searchForAll_EBRs(rv, dim, rEBR, long_modes)

    TETB_vs_LM = [
        [(all_EBRs[i] @ ebrs).tolist(), (long_modes[i] @ ebrs)]
        for i in range(len(long_modes))
        if all_EBRs[i].shape[0] > 0
    ]

    physical = [
        [
            phys(TETB_vs_LM[i][0][j], v, ebrs, EBRs, SG, irreps)
            for j in range(len(TETB_vs_LM[i][0]))
        ]
        for i in range(len(TETB_vs_LM))
    ]

    bs = [
        [
            (
                irreps
                @ (
                    ebr2irvec(TETB_vs_LM[i][0][j], ebrs, EBRs)
                    - ebr2irvec(TETB_vs_LM[i][1], ebrs, EBRs)
                    - v
                )
            )
            for j in range(len(TETB_vs_LM[i][0]))
        ]
        for i in range(len(TETB_vs_LM))
    ]

    return [
        (TETB_vs_LM[i][0], TETB_vs_LM[i][1], physical[i], bs[i])
        for i in range(len(TETB_vs_LM))
    ]


def showOnlyPhysical(
    v: np.ndarray,
    t: int,
    dim: list[int],
    EBRs: np.ndarray,
    N_gamma: int,
    ebrs: list[sp.core.symbol.Symbol],
    irreps: list[sp.core.symbol.Symbol],
    SG: int,
) -> list[tuple[list[sp.core.add.Add], sp.core.add.Add, list[sp.core.add.Add]]]:
    """Obtains all physical linear combinations of EBRs which represents the symmetry
    vector `v`. In addition, it computes any necessary longitudinal modes to obtain such
    decompositions in EBRs, and, the ill-defined content at the Gamma point and zero
    frequency.

    Args:
        v (numpy.ndarray): vector of multiplicities of the irreps following the order
        specified in `irreps`.
        t (int): positive integer indicating up to what dimension of the longitudinal
        modes search.
        dim (list): list indicating the dimensions of the EBRs of the space group in the
        order specified by `ebrs`.
        EBRs (numpy.ndarray): matrix form by columns with the vectors of multiplicities
        of the irreps for each EBR.
        N_gamma (int): number of irreps at Gamma.
        ebrs (list): list of variables defining the posible EBRs of the space group.
        irreps (list): list of irreps indicating their order and their names.
        SG (int): number of the Space Group

    Returns:
        list: list of tuples which first component indicate the physical linear
        combination of EBRs, second component is the EBRs representing the longitudinal
        modes used, and, third component the ill-defined symmetry content at Gamma and
        zero frequency.
    """
    rv = np.delete(v, range(N_gamma))  # type: ignore
    rEBR = np.delete(EBRs, range(N_gamma), axis=0)  # type: ignore

    long_modes = searchForAllLongModes(dim, t).solutions()

    all_EBRs = searchForAll_EBRs(rv, dim, rEBR, long_modes)

    TETB_vs_LM = [
        [(all_EBRs[i] @ ebrs).tolist(), (long_modes[i] @ ebrs)]
        for i in range(len(long_modes))
        if all_EBRs[i].shape[0] > 0
    ]

    physical = [
        [
            phys(TETB_vs_LM[i][0][j], v, ebrs, EBRs, SG, irreps)
            for j in range(len(TETB_vs_LM[i][0]))
        ]
        for i in range(len(TETB_vs_LM))
    ]

    bs = [
        [
            (
                irreps
                @ (
                    ebr2irvec(TETB_vs_LM[i][0][j], ebrs, EBRs)
                    - ebr2irvec(TETB_vs_LM[i][1], ebrs, EBRs)
                    - v
                )
            )
            for j in range(len(TETB_vs_LM[i][0]))
        ]
        for i in range(len(TETB_vs_LM))
    ]

    return [
        (TETB_vs_LM[i][0][j], TETB_vs_LM[i][1], bs[i][j])
        for i in range(len(TETB_vs_LM))
        for j in range(len(TETB_vs_LM[i][0]))
        if physical[i][j] == True
    ]
