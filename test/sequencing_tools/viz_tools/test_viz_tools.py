from sequencing_tools.viz_tools import mixed_sort, ColorEncoder
from sequencing_tools.utils import SeqUtilsError
import pytest
labels = ['A','A','B','A','B']

def test_mixed_sort():

    a = ["A1", "A10", "A2", "A20"]
    assert mixed_sort(a) == ["A1", "A2", "A10", "A20"]


@pytest.fixture(scope='module')
def color_encoder():
    ce = ColorEncoder()
    ce.fit(labels, ['red', 'blue'])
    return ce


def test_ColorEncoder_encoder(color_encoder):
    assert color_encoder.encoder == {'A':'red', 'B': 'blue'}

def test_ColorEncoder_fit_transform(color_encoder):
    assert all(color_encoder.fit_transform(labels, ['red','blue']) == ['red','red','blue','red','blue'])

def test_ColorEncoder_transform(color_encoder):
    assert all(color_encoder.transform(['A','A','B']) == ['red','red','blue'])


def test_ColorEncoder_transform_error(color_encoder):
    with pytest.raises(SeqUtilsError) as e:
        color_encoder.transform(['A','a','C','B'])
    
    assert "Contain unseen data!!:" in str(e.value)

def test_assert_color_encoder(color_encoder):
    with pytest.raises(SeqUtilsError) as e:
        color_encoder.fit(['A','B','C'], ['red','blue'])
    
    assert "Not enough colors!! 2 colors for 3 categories" in str(e.value)