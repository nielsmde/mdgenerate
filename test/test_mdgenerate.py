from mdgenerate import mdgenerate

import os
import pytest


@pytest.fixture
def meta_grompp():
    meta = mdgenerate.MetaDict({
        'name': 'test',
        'time': '1ns',
        'timestep': '1fs',
        'topology-files': {'topol.top', 'conf.gro'}
    })


def test_time_to_ps():
    assert mdgenerate.time_to_ps('1s') == 1e12
    assert mdgenerate.time_to_ps('1ms') == 1e9
    assert mdgenerate.time_to_ps('1us') == 1e6
    assert mdgenerate.time_to_ps('1ns') == 1e3
    assert mdgenerate.time_to_ps('1ps') == 1e0
    assert mdgenerate.time_to_ps('1fs') == 1e-3


def test_locate_template():
    mdgenerate.TEMPLATE_DIR = template_dir = '/tmp/.mdgenerate/templates'

    assert mdgenerate.locate_template('base') == os.path.join(template_dir, 'base.yaml')
    assert mdgenerate.locate_template('base.ext') == os.path.join(template_dir, 'base.ext.yaml')
    assert mdgenerate.locate_template('base/pore') == os.path.join(template_dir, 'base/pore.yaml')
    assert mdgenerate.locate_template('/full/path/template.yaml') == '/full/path/template.yaml'
