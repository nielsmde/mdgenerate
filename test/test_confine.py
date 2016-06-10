
from mdgenerate.mdgenerate import run_gmx
import mdgenerate.confine as conf

import mdevaluate as md
import os


def test_generate_spherical_water():
    """
    Check if a correct topology is created, by runing grompp with it.
    """
    simdir = '/data/robin/sim/nvt/12kwater/200_8_0_NVT'
    testdir = '/tmp/jdahfjncckjke/'
    if not os.path.exists(simdir):
        print('Skipping test of generate_spherical_water.\nFile not found:{}'.format(simdir))
        return
    tr = md.open(simdir)
    conf.generate_spherical_water(testdir, tr, 99, 3.0)

    os.chdir(testdir)
    with open('grompp.mdp', 'a'):
        pass
    stdout, stderr = run_gmx('gmx grompp -c water.gro -p water.top')
    assert 'error' not in stdout.lower()
    assert 'error' not in stderr.lower()
