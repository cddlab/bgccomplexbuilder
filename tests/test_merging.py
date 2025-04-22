import pytest

from complexbuilder.merging import make_list_of_directories


def test_make_list_of_directories():
    fastafile = "/data2/moriwaki/BGCcomplex/newfastas/newBGC0000001.fasta"
    make_list_of_directories(fastafile)
    assert make_list_of_directories(fastafile)[0] == "aek75490.1_aek75490.1"
    assert len(make_list_of_directories(fastafile)) == 378
