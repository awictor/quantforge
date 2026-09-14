"""Arithmetic expression evaluation via the shunting-yard algorithm.

Parses and evaluates infix arithmetic without Python's ``eval`` -- so untrusted input can be
scored safely. `tokenize` splits a string into numbers, operators, and parentheses;
`shunting_yard` converts the token stream to Reverse Polish (postfix) form respecting
precedence and associativity; `eval_rpn` folds a postfix stream to a number; and
`eval_expression` chains all three. Supports ``+ - * / // % ^`` (``^`` = power,
right-associative), unary minus, parentheses, and a small set of named functions and
constants. Pure standard library.
"""

import math

# exponentiation binds tighter than unary minus, so -2^2 == -(2^2) == -4
_PRECEDENCE = {"+": 1, "-": 1, "*": 2, "/": 2, "//": 2, "%": 2, "u-": 3, "^": 4}
_RIGHT_ASSOC = {"^", "u-"}
_FUNCTIONS = {
    "sin": math.sin, "cos": math.cos, "tan": math.tan,
    "exp": math.exp, "log": math.log, "sqrt": math.sqrt, "abs": abs,
}
_CONSTANTS = {"pi": math.pi, "e": math.e}


def tokenize(expr):
    """Split ``expr`` into a token list of numbers (float), operators, parens, and names."""
    tokens = []
    i = 0
    n = len(expr)
    while i < n:
        c = expr[i]
        if c.isspace():
            i += 1
            continue
        if c.isdigit() or c == ".":
            j = i
            while j < n and (expr[j].isdigit() or expr[j] in ".eE"
                             or (expr[j] in "+-" and j > i and expr[j - 1] in "eE")):
                j += 1
            tokens.append(float(expr[i:j]))
            i = j
        elif c.isalpha() or c == "_":
            j = i
            while j < n and (expr[j].isalnum() or expr[j] == "_"):
                j += 1
            tokens.append(expr[i:j])
            i = j
        elif c == "/" and i + 1 < n and expr[i + 1] == "/":
            tokens.append("//")
            i += 2
        elif c in "+-*/%^(),":
            tokens.append(c)
            i += 1
        else:
            raise ValueError("unexpected character %r in expression" % c)
    return tokens


def _is_value_end(tok):
    # a token after which a binary operator (not unary) is expected
    return isinstance(tok, float) or tok == ")" or (isinstance(tok, str) and tok in _CONSTANTS)


def shunting_yard(tokens):
    """Convert an infix token list to a postfix (RPN) list (Dijkstra's shunting-yard)."""
    output = []
    stack = []
    prev = None
    for tok in tokens:
        if isinstance(tok, float):
            output.append(tok)
        elif isinstance(tok, str) and tok in _CONSTANTS:
            output.append(_CONSTANTS[tok])
        elif isinstance(tok, str) and tok in _FUNCTIONS:
            stack.append(tok)
        elif tok == ",":
            while stack and stack[-1] != "(":
                output.append(stack.pop())
            if not stack:
                raise ValueError("misplaced comma or mismatched parentheses")
        elif tok in _PRECEDENCE or tok == "-" and (prev is None or (isinstance(prev, str) and prev in _PRECEDENCE | {"(", ","})):
            op = tok
            if op == "-" and (prev is None or (isinstance(prev, str) and (prev in _PRECEDENCE or prev in "(,"))):
                op = "u-"       # unary minus (prefix): never pops lower-precedence binary ops
            if op != "u-":
                while stack and stack[-1] != "(" and (
                    stack[-1] in _FUNCTIONS
                    or _PRECEDENCE.get(stack[-1], 0) > _PRECEDENCE[op]
                    or (_PRECEDENCE.get(stack[-1], 0) == _PRECEDENCE[op] and op not in _RIGHT_ASSOC)
                ):
                    output.append(stack.pop())
            stack.append(op)
        elif tok == "(":
            stack.append(tok)
        elif tok == ")":
            while stack and stack[-1] != "(":
                output.append(stack.pop())
            if not stack:
                raise ValueError("mismatched parentheses")
            stack.pop()         # discard "("
            if stack and stack[-1] in _FUNCTIONS:
                output.append(stack.pop())
        else:
            raise ValueError("unknown token %r" % (tok,))
        prev = tok
    while stack:
        top = stack.pop()
        if top in "()":
            raise ValueError("mismatched parentheses")
        output.append(top)
    return output


def eval_rpn(rpn):
    """Evaluate a postfix (RPN) token list to a number."""
    stack = []
    for tok in rpn:
        if isinstance(tok, float):
            stack.append(tok)
        elif tok in _FUNCTIONS:
            if not stack:
                raise ValueError("function missing its argument")
            stack.append(float(_FUNCTIONS[tok](stack.pop())))
        elif tok == "u-":
            if not stack:
                raise ValueError("unary minus missing its operand")
            stack.append(-stack.pop())
        else:
            if len(stack) < 2:
                raise ValueError("operator %r missing operands" % tok)
            b = stack.pop()
            a = stack.pop()
            if tok == "+":
                stack.append(a + b)
            elif tok == "-":
                stack.append(a - b)
            elif tok == "*":
                stack.append(a * b)
            elif tok == "/":
                stack.append(a / b)
            elif tok == "//":
                stack.append(float(math.floor(a / b)))
            elif tok == "%":
                stack.append(math.fmod(a, b))
            elif tok == "^":
                stack.append(a ** b)
            else:
                raise ValueError("unknown operator %r" % tok)
    if len(stack) != 1:
        raise ValueError("malformed expression")
    return stack[0]


def eval_expression(expr):
    """Evaluate an infix arithmetic string safely (no Python ``eval``)."""
    return eval_rpn(shunting_yard(tokenize(expr)))
