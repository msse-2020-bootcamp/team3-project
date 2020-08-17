"""
Tests for Python Standard Library implementation of MC code.
"""
import math

import pytest

from monte_carlo import calculate_distance as mc_cd
   
def test_calculate_distance():

    point1 = [0, 0, 0]
    point2 = [0, 0, 1]

    expected_value = 1

    observed_value = mc_cd(point1, point2)

    assert expected_value == observed_value

def test_calculate_distance2():

    point3 = [0, 0, 0]
    point4 = [0, 1, 1]

    assert math.sqrt(2) == mc_cd(point3, point4)

def test_calculate_distance_periodic():

    point5 = [0, 0, 0]
    point6 = [0, 0, 8]

    assert 2 == mc_cd(point5, point6, 10)

"""
@pytest.mark.parametrize('variable_name1, variable_name2, ... , variable_nameN', [
    (variable_value1, variable_value2, ... , variable_valueN),
    (variable_value1, variable_value2, ... , variable_valueN),
]
)
def test_function(variable_name1, variable_name2, ... , variable_nameN)
    *** TEST CODE HERE ***
"""
# -@pytest.mark.skip
@pytest.mark.parametrize("point1, point2, expected_distance, box_length", [
    ([0, 0, 0], [1, 0, 0], 1, None),
    ([0, 0, 0], [8, 0, 0], 8, None),
    ([0, 0, 0], [8, 0, 0], 2, 10),
])
def test_calculate_distance_many(point1, point2, expected_distance, box_length):

    observed_distance = mc_cd(point1, point2, box_length)

    assert expected_distance == observed_distance