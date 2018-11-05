import unittest
from nose.tools import raises
from kqml import KQMLList
from bioagents.tests.util import ekb_from_text
from bioagents.biosense.biosense_module import BioSense_Module, _get_agent
from bioagents.biosense.biosense import BioSense, InvalidAgentError, \
    InvalidCollectionError, UnknownCategoryError, \
    CollectionNotFamilyOrComplexError


# example ekb terms
mek1_ekb = ekb_from_text('MAP2K1')  # agent
dusp_ekb = ekb_from_text('DUSP6')  # agent
mek_ekb = ekb_from_text('MEK')  # family
foo_ekb = ekb_from_text('foo')  # invalid
mek1 = _get_agent(mek1_ekb)
mek = _get_agent(mek_ekb)
dusp = _get_agent(dusp_ekb)


# BioSense python API unit tests
def test_choose_sense():
    bs = BioSense()
    cases = [(mek1_ekb, 'MAP2K1', 'ONT::GENE'),
             (dusp_ekb, 'DUSP6', 'ONT::GENE'),
             (mek_ekb, 'MEK', 'ONT::PROTEIN-FAMILY')]
    for case in cases:
        agents, _ = bs.choose_sense(case[0])
        agent, ont_type, _ = list(agents.values())[0]
        assert agent.name == case[1]
        assert ont_type == case[2]


@raises(InvalidAgentError)
def test_choose_sense_invalid_agent():
    """should raise InvalidAgentError if the agent is not recognized"""
    bs = BioSense()
    invalid_case = foo_ekb
    bs.choose_sense(invalid_case)


def test_choose_nonsense():
    """ekb terms that aren't biological agents should have ont-type None

    BAGEL is from ONT::BAGELS-BISCUITS
    """
    bs = BioSense()
    case = ekb_from_text('bagel')
    agents, _ = bs.choose_sense(case)
    _, ont_type, _ = list(agents.values())[0]
    assert ont_type is None


def test_choose_sense_category():
    bs = BioSense()
    cases = [(mek1, [('kinase activity', 'TRUE'),
                     ('enzyme', 'TRUE'),
                     ('kinase', 'TRUE'),
                     ('transcription-factor', 'FALSE'),
                     ('W::KINASE', 'TRUE'),
                     ('phosphatase', 'FALSE')]),
             (dusp, [('phosphatase', 'TRUE'), ('enzyme', 'TRUE')]),
             (_get_agent(ekb_from_text('BRAF')), [('kinase', 'TRUE')])]
    for ag, result_tuples in cases:
        for cat, result in result_tuples:
            print('Testing: %s. Expect result %s.' % (cat, result))
            in_category = bs.choose_sense_category(ag, cat)
            assert in_category == (result == 'TRUE')


def test_choose_sense_category_unknown_category():
    """should raise UnknownCategoryError if the category is not recognized"""
    bs = BioSense()
    assert bs.choose_sense_category(mek, 'foo') is None


def test_choose_sense_is_member():
    bs = BioSense()
    cases = [(mek1, mek, True),
             (dusp, mek, False)]
    for (agent, collection, result) in cases:
        assert bs.choose_sense_is_member(agent, collection) == result


def test_choose_sense_what_member():
    bs = BioSense()
    members = bs.choose_sense_what_member(mek)
    member_names = [agent.name for agent in members]
    result = ['MAP2K1', 'MAP2K2']
    assert member_names == result


def test_get_synonyms():
    bs = BioSense()
    example_synonyms = {'PRKMK1', 'MEK1', 'MAP2K1', 'ERK activator kinase 1',
                        'MKK1', 'MEK 1'}
    synonyms = set(bs.get_synonyms(mek1))
    assert example_synonyms.issubset(synonyms)


# BioSense module unit tests
def test_respond_choose_sense():
    bs = BioSense_Module(testing=True)
    msg_content = KQMLList('CHOOSE-SENSE')
    msg_content.sets('ekb-term', mek1_ekb)
    res = bs.respond_choose_sense(msg_content)
    print(res)
    agents = res.get('agents')
    assert agents and agents.data
    agent = agents[0]
    name = agent.gets('name')
    assert name == 'MAP2K1'
    ont_type = agent.get('ont-type')
    assert ont_type == 'ONT::GENE'


def test_respond_choose_nonsense():
    bs = BioSense_Module(testing=True)
    msg_content = KQMLList('CHOOSE-SENSE')
    msg_content.sets('ekb-term', ekb_from_text('bagel'))
    res = bs.respond_choose_sense(msg_content)
    print(res)
    assert res.head() == 'SUCCESS'
    assert res.get('agents')[0].gets('ont-type') is None


@unittest.skip('No ambiguity reported here yet')
def test_respond_choose_sense_ambiguity():
    bs = BioSense_Module(testing=True)
    msg_content = KQMLList('CHOOSE-SENSE')
    pdk1_ekb = ekb_from_text('PDK1')
    msg_content.sets('ekb-term', pdk1_ekb)
    res = bs.respond_choose_sense(msg_content)
    print(res)
    agents = res.get('agents')
    assert agents and agents.data
    agent = agents[0]
    name = agent.gets('name')
    assert name == 'PDK1'
    ont_type = agent.get('ont-type')
    assert ont_type == 'ONT::GENE'


def test_respond_choose_sense_category():
    bs = BioSense_Module(testing=True)
    cases = [(mek1_ekb, [('kinase activity', 'TRUE'),
                         ('enzyme', 'TRUE'),
                         ('kinase', 'TRUE'),
                         ('transcription-factor', 'FALSE'),
                         ('W::KINASE', 'TRUE'),
                         ('phosphatase', 'FALSE')]),
             (dusp_ekb, [('phosphatase', 'TRUE'), ('enzyme', 'TRUE')]),
             (ekb_from_text('BRAF'), [('kinase', 'TRUE')])]
    for ekb, result_tuples in cases:
        msg_content = KQMLList('CHOOSE-SENSE-CATEGORY')
        msg_content.sets('ekb-term', ekb)
        for cat, result in result_tuples:
            print('Testing: %s. Excpet result %s.' % (cat, result))
            msg_content.sets('category', cat)
            res = bs.respond_choose_sense_category(msg_content)
            print(res)
            print(res.head())
            assert(res.head() == 'SUCCESS')
            assert(res.get('in-category') == result)


def test_respond_choose_sense_is_member():
    bs = BioSense_Module(testing=True)
    msg_content = KQMLList('CHOOSE-SENSE-IS-MEMBER')
    msg_content.sets('ekb-term', mek1_ekb)
    msg_content.sets('collection', mek_ekb)
    print(msg_content)
    res = bs.respond_choose_sense_is_member(msg_content)
    print(res)
    print(res.head())
    assert(res.head() == 'SUCCESS')
    assert(res.get('is-member') == 'TRUE')


def test_respond_choose_sense_what_member():
    bs = BioSense_Module(testing=True)
    msg_content = KQMLList('CHOOSE-SENSE-WHAT-MEMBER')
    msg_content.sets('collection', mek_ekb)
    print(msg_content)
    res = bs.respond_choose_sense_what_member(msg_content)
    print(res)
    print(res.head())
    assert(res.head() == 'SUCCESS')
    assert(len(res.get('members')) == 2)
    m1 = res.get('members')[0]
    m2 = res.get('members')[1]
    assert m1.gets('name') == 'MAP2K1', m1.gets('name')
    assert m2.gets('name') == 'MAP2K2', m2.gets('name')


def test_respond_get_synonyms():
    bs = BioSense_Module(testing=True)
    msg_content = KQMLList('GET-SYNONYMS')
    msg_content.sets('entity', mek1_ekb)
    res = bs.respond_get_synonyms(msg_content)
    assert res.head() == 'SUCCESS'
    syns = res.get('synonyms')
    syn_strs = [s.gets(':name') for s in syns]
    assert 'MAP2K1' in syn_strs
    assert 'MEK1' in syn_strs
    assert 'MKK1' in syn_strs
