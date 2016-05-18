
from mdgenerate import mdgenerate

import os
import pytest

tmp_templates = []


@pytest.fixture(scope='session')
def yaml(request):
    with open('/tmp/base.yaml', 'w') as f:
        f.write(
            """
            test-key: t
            queue: short
            mdp:
              mdp1: t
            """
        )
    with open('/tmp/meta.yaml', 'w') as f:
        f.write(
            """
            template: /tmp/base.yaml
            time: 1ns
            timestep: 1fs
            T: 300
            queue: nodes
            mdp:
              mdp2: t
            """
        )

    def finalize():
        for t in tmp_templates:
            os.remove(t)
    request.addfinalizer(finalize)


@pytest.fixture(scope='session')
def yaml_multi():
    with open('/tmp/base2.yaml', 'w') as f:
        f.write(
            """
            test-key-2: t
            mdp:
              mdp1_2: t
            """
        )
    with open('/tmp/meta2.yaml', 'w') as f:
        f.write(
            """
            template:
              - /tmp/base.yaml
              - /tmp/base2.yaml
            queue: nodes
            mdp:
              mdp2: t
            """
        )
    global tmp_templates
    tmp_templates.append('/tmp/base2.yaml')
    tmp_templates.append('/tmp/meta2.yaml')


@pytest.fixture(scope='session')
def yaml_dicts(yaml):
    with open('/tmp/base3.yaml', 'w') as f:
        f.write(
            """
            slurm:
              skey: t
            """
        )
    with open('/tmp/meta3.yaml', 'w') as f:
        f.write(
            """
            template:
              - /tmp/base.yaml
              - /tmp/base3.yaml
            mdrun:
              mkey: t
            """
        )


def test_inherit(yaml):
    base = mdgenerate.MetaDict('/tmp/base.yaml')
    meta = mdgenerate.MetaDict('/tmp/meta.yaml')
    assert base.directory == '/tmp'
    assert base.filename == 'base.yaml'
    assert meta.filename == 'meta.yaml'

    assert base['test-key']
    assert base['queue'] == 'short'
    assert meta['test-key']
    assert meta['queue'] == 'nodes'

    mdp_base = base['mdp']
    mdp_meta = meta['mdp']
    assert mdp_base['mdp1']
    assert mdp_meta['mdp1']
    assert mdp_meta['mdp2']


def test_multinherit(yaml_multi):
    meta2 = mdgenerate.MetaDict('/tmp/meta2.yaml')

    assert meta2['test-key']
    assert meta2['test-key-2']
    assert meta2['queue'] == 'nodes'
    mdp = meta2['mdp']
    assert mdp['mdp1']
    assert mdp['mdp1_2']
    assert mdp['mdp2']


def test_multidicts(yaml_dicts):
    meta3 = mdgenerate.MetaDict('/tmp/meta3.yaml')

    assert meta3['mdp']
    assert meta3['slurm']
    assert meta3['mdrun']


def test_mdp(yaml):
    meta = mdgenerate.MetaDict('/tmp/meta.yaml')
    mdp = meta.as_mdp
    assert isinstance(mdp, str)
    assert 'mdp1' in mdp
    assert 'mdp2' in mdp
