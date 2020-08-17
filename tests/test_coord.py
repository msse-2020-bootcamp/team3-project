"""
Tests for the coord module
"""
import os
import sys
import pytest

# Add our folder to the system path so python can find our code.
current_location = os.path.dirname(__file__)
sys.path.append(os.path.join(current_location, '..'))

from mcsim.coordinates import calculate_distance

def test_calculate_distance():
    point_1 = [0, 0, 0]
    point_2 = [1, 0, 0]

    expected_distance = 1
    dist1 = calculate_distance(point_1, point_2)
    
    assert dist1 == expected_distance

def test_calculate_distance2():
    point_1 = [0, 0, 0]
    point_2 = [0, 0, 8]
    box_length = 10

    expected_distance = 2
    dist1 = calculate_distance(point_1, point_2, box_length=box_length)
    assert dist1 == expected_distance

@pytest.mark.parametrize("point1, point2, expected_distance, box_length",
    [
        ([0, 0, 0], [1, 0 , 0], 1, None),
        ([0, 0, 0], [8, 0 , 0], 8, None),
        ([0, 0, 0], [8, 0 , 0], 2, 10)
    ],
)
def test_calculate_distance_many(point1, point2, expected_distance, box_length):

    calculated_distance = calculate_distance(point1, point2, box_length=box_length)

    assert calculated_distance == expected_distance