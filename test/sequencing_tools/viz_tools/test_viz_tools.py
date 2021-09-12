from sequencing_tools.viz_tools import mixed_sort


def test_sort():

    a = ["A1", "A10", "A2", "A20"]
    assert mixed_sort(a) == ["A1", "A2", "A10", "A20"]
