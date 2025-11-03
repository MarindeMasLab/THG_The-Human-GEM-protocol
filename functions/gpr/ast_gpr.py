from ast import And, BoolOp, Expression, Name, Or
from functools import reduce
from typing import List, Union

# we are only using cobra for the GPR parse
import cobra
from cobra.core.gene import GPR, ast_parse
import logging
import re

_LOGGER = logging.getLogger(__name__)


def divide_gpr_in_ors(ast) -> List:
    """Divide GPR in AST form recursively until a Gene or a AND operation is found."""
    if isinstance(ast, Name) or isinstance(ast, str):
        # base case 1: a gene
        return [ast]
    elif isinstance(ast, BoolOp):
        if isinstance(ast.op, Or):
            return reduce(lambda x, y: x + divide_gpr_in_ors(y), ast.values, [])
        if isinstance(ast.op, And):
            # base case 2: AND is returned as is
            return [ast]


def compare_ast(ast_left, ast_right) -> bool:
    """Compare recursively that 2 AST of genes are the same."""
    # base case 1: the asts are str or Name -> str comparison
    if isinstance(ast_left, str) and isinstance(ast_right, str):
        if ast_left != ast_right:
            return False
    elif isinstance(ast_left, Name) and isinstance(ast_right, Name):
        if ast_left.id != ast_right.id:
            return False
    elif isinstance(ast_left, BoolOp) and isinstance(ast_right, BoolOp):
        for left_leaf in ast_left.values:
            if not any(
                compare_ast(left_leaf, leaf_right) for leaf_right in ast_right.values
            ):
                return False
    else:
        # base case 2: types do not match
        return False
    return True


def deduplicate_gpr(ast: Union[Expression, str]) -> GPR:
    """Merge two GPRs represented as GPR `ast.Expression`."""
    if isinstance(ast, str):
        ast = ast_parse(ast).body[0].value
    or_left = divide_gpr_in_ors(ast)
    to_remove = []
    for i, left_node in enumerate(or_left):
        if i not in to_remove:
            for j, right_node in enumerate(or_left):
                if i != j:
                    if compare_ast(left_node, right_node):
                        to_remove.append(j)

    return GPR(
        Expression(
            BoolOp(Or(), [node for i, node in enumerate(or_left) if i not in to_remove])
        )
    )


def is_or(node: Expression) -> bool:
    if isinstance(node, BoolOp):
        if isinstance(node.op, Or):
            return True
    return False


def naive_reduce_or(ast: Expression):
    """Reduce AST expression by commutative operations at the same level.

    It is not guaranteed to reduce the AST to canonical form, but it
    tackles the usual non-reduced GPRs and it is faster than Quine-McCluskey.
    """
    to_remove = []
    if isinstance(ast, Name) or isinstance(ast, str):
        return
    # merge ORs at the same level
    for i, node in enumerate(ast.values):
        if is_or(node):
            if i != len(ast.values) - 1:
                for j in range(i + 1, len(ast.values)):
                    if is_or(ast.values[j]):
                        to_remove.append(j)
                        node.values += ast.values[j].values
    ast.values = [node for i, node in enumerate(ast.values) if i not in to_remove]
    # walk down the tree
    for node in ast.values:
        if isinstance(node, BoolOp):
            naive_reduce_or(node)


def remove_empty(ast: Expression):
    """Remove nodes with empty identifiers **in-place**."""
    if isinstance(ast, Name) or isinstance(ast, str):
        return ast
    to_remove = []
    for i, node in enumerate(ast.values):
        if isinstance(node, BoolOp):
            remove_empty(node)
        elif isinstance(node, Name):
            if node.id == "ZEMPTYZ":
                to_remove.append(i)
        else:
            if node == "ZEMPTYZ":
                to_remove.append(i)
    if to_remove:
        ast.values = [val for i, val in enumerate(ast.values) if i not in to_remove]


def reduce_gpr(ast_str: str) -> GPR:
    """Sum OR nodes at the same level of a GPR."""
    # first replace empty gene names to allow parsing them
    ast_str = ast_str.replace("''", "ZEMPTYZ")

    # Remove stoichiometric coefficients (*N) as they're not part of standard GPR syntax
    # GPR expressions should only contain gene IDs with boolean operators (and/or)
    import re

    ast_str = re.sub(r"\*\d+", "", ast_str)

    # Fix numeric gene IDs by prefixing them with 'G' to make them valid identifiers
    # This prevents ast.Constant errors when gene IDs are pure numbers
    # Pattern: Match numbers that are NOT preceded or followed by alphanumeric characters
    # This ensures we only replace standalone numbers, not numbers within identifiers
    def prefix_numeric_genes(match):
        num = match.group(1)
        return f"G{num}"

    # Look for numbers that appear after non-word characters or at start,
    # and before non-word characters or at end (but not * which we just removed)
    # (?<![a-zA-Z0-9_]) - not preceded by alphanumeric or underscore
    # (\d+) - one or more digits
    # (?![a-zA-Z0-9_]) - not followed by alphanumeric or underscore
    ast_str = re.sub(
        r"(?<![a-zA-Z0-9_])(\d+)(?![a-zA-Z0-9_])", prefix_numeric_genes, ast_str
    )

    to_reduce = ast_parse(ast_str).body[0].value
    remove_empty(to_reduce)
    naive_reduce_or(to_reduce)
    return cobra.core.gene.GPR(Expression(to_reduce))


def sanitize_gpr(gpr: str) -> str:
    """Remove empty names, reduce SUM-types and deduplicate nodes."""
    original_gpr = gpr

    # If the expression contains only ORs and no ANDs, parentheses are not
    # meaningful for logical grouping (OR is associative). Strip all
    # parentheses/brackets to avoid parse errors from unmatched brackets.
    try:
        lower_gpr = gpr.lower()
        if "and" not in lower_gpr and "or" in lower_gpr:
            _LOGGER.info(
                "GPR appears to contain only ORs (no ANDs). Stripping all brackets for: %r",
                original_gpr,
            )
            gpr = (
                gpr.replace("(", "").replace(")", "").replace("[", "").replace("]", "")
            )
    except Exception:
        _LOGGER.debug("Error while checking for OR-only GPR", exc_info=True)

    # Remove parentheses that wrap a single gene (e.g. `(GENE1)` -> `GENE1`)
    # Do this iteratively to handle nested harmless parentheses.
    try:
        prev = None
        while prev != gpr:
            prev = gpr
            gpr = re.sub(r"\(\s*([A-Za-z0-9_\-]+)\s*\)", r"\1", gpr)
    except Exception:
        _LOGGER.debug("Failed to strip single-gene parentheses", exc_info=True)

    # If parentheses are unmatched, warn and remove them to allow parsing to continue
    if gpr.count("(") != gpr.count(")"):
        _LOGGER.warning(
            "Unmatched parentheses in GPR: %r -- original: %r",
            gpr,
            original_gpr,
        )
        # remove all parentheses (best-effort salvage)
        gpr = gpr.replace("(", "").replace(")", "")

    # Proceed with existing reduction/parsing, but guard against parse failures
    try:
        reduced = reduce_gpr(gpr.replace("]", "").replace("[", "")).to_string()
    except Exception as e:
        _LOGGER.warning(
            "reduce_gpr failed for GPR %r (after pre-clean): %s -- falling back by stripping brackets",
            original_gpr,
            e,
        )
        _LOGGER.debug("traceback:", exc_info=True)
        # fallback: strip all bracket characters and try again
        safe = gpr.replace("(", "").replace(")", "").replace("[", "").replace("]", "")
        reduced = reduce_gpr(safe).to_string()

    return deduplicate_gpr(reduced).to_string()
