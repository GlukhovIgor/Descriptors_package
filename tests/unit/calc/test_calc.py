from biodescriptors.calc import *
from pytest import approx


test_structure = 'tests/data/1DB1.pdb'
href = [[127,142], [149,152], [226,246]]


def test_calc_len_of_hel():    
    ans = calc_len_of_hel(test_structure, href)
    correct_ans = [approx(x, abs=0.5) for x in [22.6, 5.6, 30.5]]
    assert ans == correct_ans
    

def test_calc_COM_Calpha_angles():    
    ans = calc_COM_Calpha_angles(test_structure, href)
    correct_ans = [approx(x, abs=0.5) for x in [61.0, 120.5, 14.2]]
    assert ans == correct_ans
    

def test_calc_COM_protein():    
    ans = calc_COM_protein(test_structure)
    correct_ans = [approx(x, abs=0.5) for x in [13.360, 20.776, 42.003]]
    assert ans == correct_ans
    
    
def test_calc_angles_between_hel():    
    ans = [calc_angles_between_hel(test_structure, href)]
    correct_ans = [
        [approx(49.5, abs=0.5), 
         approx(148.37, abs=0.5)
        ], 
        approx(125.10, abs=0.5)
    ]
    assert ans == correct_ans
    
    
def test_calc_COM_helix():    
    ans = calc_COM_helix(test_structure, href)
    correct_ans = [[approx(x, abs=0.5) for x in [19.588, 6.445, 50.056]],
                   [approx(x, abs=0.5) for x in [10.896, 11.219, 26.464]],
                   [approx(x, abs=0.5) for x in [7.358, 15.986, 35.683]]]
    assert ans == correct_ans
    
    
def test_calc_prot_hel_dist():    
    ans = calc_prot_hel_dist(test_structure, href)
    correct_ans = [approx(x, abs=0.5) for x in [17.6, 9.9, 18.4]]
    assert ans == correct_ans