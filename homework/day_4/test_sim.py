"""
Tests for Python Standard Library implementation of Monte Carlo code.
"""
import math, random
import pytest
import os
import mcsim.monte_carlo as mc

@pytest.mark.parametrize("point1, point2, expected_distance, box_length",[
    ([0,0,0], [0,0,8], 8, None),
    ([0,0,0], [1,1,1], math.sqrt(3), None),
    ([0,0,0], [0,0,8], 2, 10),
])
def test_calculate_distance(point1, point2, expected_distance, box_length):
    real_distance = mc.calculate_distance(point1, point2, box_length)
    assert real_distance == expected_distance

@pytest.mark.parametrize("r_ij, expected_LJ",[
    (1, 0),
    (2**(1/6), -1),
])
def test_calculate_LJ(r_ij, expected_LJ):
    real_LJ = mc.calculate_LJ(r_ij)
    assert real_LJ == expected_LJ

def test_calculate_tail_correction():
    num_particles = 1
    test_box_length = 1
    test_cutoff = 1
    assert mc.calculate_tail_correction(num_particles, test_box_length, test_cutoff) == -16/9*math.pi

@pytest.mark.parametrize("seed, delta_e, beta, expected_result",[
    (0, 1, 1, False),
    (1, 1, 1, True),
    (0, -1, 1, True),
    (0, 1, 0.1, True)
])
def test_calculate_pair_energy(seed, delta_e, beta, expected_result):
    random.seed(seed)
    real_result = mc.accept_or_reject(delta_e, beta)
    random.seed()
    assert real_result == expected_result

def test_calculate_pair_energy():
    coordinates = [[0,0,0], [0,0,2**(1/6)], [0,0,2*(2**(1/6))]]
    assert mc.calculate_pair_energy(coordinates, 1, 10, 3) == -2

def test_calculate_total_energy():
    # This is a system test (integration)
    # Test file: lj_sample_config_periodic1.txt
    filepath = os.path.join('..', '..', 'lj_sample_configurations', 'lj_sample_config_periodic1.txt')
    coordinates, box_length = mc.read_xyz(filepath)
    real_etot = mc.calculate_total_energy(coordinates, box_length, 3)
    assert abs(real_etot - (-4351.5)) < 0.1