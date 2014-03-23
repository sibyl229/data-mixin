import pytest
from species.base import Species as SP

@pytest.fixture(scope="module")
def species():
    sp = SP('species.yeast','Yeast').load()
    return sp 

def test_sample(species):
    
    assert species.moduleName == 'phl.species.yeast'
    assert 0 in [1,2,0]











