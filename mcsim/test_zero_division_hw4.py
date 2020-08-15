"""
Tests for Python Standard Library implementation of MC code.
"""

import pytest

import monte_carlo as func1

def test_zero_division_error():
    # set initial parameters
    coordinates = [[0, 0, 0], [0, 0, 0], [0, 0, 1], [1, 1, 1]]
    box_length = 10
    cutoff = 3
    redu_temp = 1.5
    n_steps = 1000

    # run the code
    func1.run_simulation(coordinates, box_length, cutoff, redu_temp, n_steps)