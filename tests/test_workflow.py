import pytest

from mbworkbench.workflow.workflow import Workflow
from mbworkbench.lib.block import Block




def test_need_inputfile(caplog):
    
    workflow = Workflow(None)
    assert "Workflow was not given an input file" in caplog.text
    

# TODO: test that we get the expected output as well
@pytest.mark.parametrize('expected', (
    ('pass.yaml', 'pass'),
    ('fail.yaml', 'fail'),
))
def test_workflow(expected, monkeypatch, caplog):
    '''
    test if workflow will run / fail properly
    '''

    class Recorder(object):
        called = False

    def fake_write_chk(arg):
        Recorder.called = True

    infile = expected[0]
    
    workflow = Workflow(infile)
    assert workflow.current_block == 0
    
    workflow.build()
    assert len(workflow.blocks) > 0
    
    monkeypatch.setattr('mbworkbench.workflow.workflow.Workflow.write_chk', fake_write_chk)
    workflow.run()
    assert Recorder.called

