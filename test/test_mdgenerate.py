from mdgenerate import mdgenerate

import os
import pytest


@pytest.fixture(scope='session')
def yaml(request):
    with open('/tmp/base.yaml', 'w') as f:
        f.write(
            """
            test-key: t
            queue: short
            timestep: 1fs
            mdp-options:
              mdp1: t
            """
        )
    with open('/tmp/meta.yaml', 'w') as f:
        f.write(
            """
            template: /tmp/base.yaml
            queue: nodes
            mdp-options:
              mdp2: t
            """
        )

    def finalize():
        os.remove('/tmp/base.yaml')
        os.remove('/tmp/meta.yaml')
    request.addfinalizer(finalize)


def test_time_to_ps():
    assert mdgenerate.time_to_ps('1s')  == 1e12
    assert mdgenerate.time_to_ps('1ms') == 1e9
    assert mdgenerate.time_to_ps('1us') == 1e6
    assert mdgenerate.time_to_ps('1ns') == 1e3
    assert mdgenerate.time_to_ps('1ps') == 1e0
    assert mdgenerate.time_to_ps('1fs') == 1e-3


def test_locate_template():
    home = os.environ['HOME']
    assert mdgenerate.locate_template('base') == '{}/.mdgenerate/templates/base.yaml'.format(home)
    assert mdgenerate.locate_template('base.ext') == '{}/.mdgenerate/templates/base.ext.yaml'.format(home)
    assert mdgenerate.locate_template('base/pore') == '{}/.mdgenerate/templates/base/pore.yaml'.format(home)
    assert mdgenerate.locate_template('/full/path/template') == '/full/path/template.yaml'
