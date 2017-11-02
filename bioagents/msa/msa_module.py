import os
import sys
import logging
import re
from bioagents import Bioagent
from indra.sources.trips.processor import TripsProcessor
from kqml import KQMLPerformative, KQMLList
import pickle


logging.basicConfig(format='%(levelname)s: %(name)s - %(message)s',
                    level=logging.INFO)
logger = logging.getLogger('MSA')


def _read_signor_afs():
    path = os.path.dirname(os.path.abspath(__file__)) + \
            '/../resources/signor_active_forms.pkl'
    with open(path, 'rb') as pkl_file:
        stmts = pickle.load(pkl_file)
    if isinstance(stmts, dict):
        signor_afs = []
        for _, stmt_list in stmts.iteritems():
            signor_afs += stmt_list
    else:
        signor_afs = stmts
    return signor_afs


class MSA_Module(Bioagent):
    name = 'MSA'
    tasks = ['PHOSPHORYLATION-ACTIVATING']
    signor_afs = _read_signor_afs()

    def receive_tell(self, msg, content):
        tell_content = content[0].to_string().upper()
        if tell_content == 'START-CONVERSATION':
            logger.info('MSA resetting')

    def respond_phosphorylation_activating(self, content):
        """Return response content to phosphorylation_activating request."""
        heading = content.head()
        m = re.match('(\w+)-(\w+)', heading)
        if m is None:
            return self.make_failure('UNKNOWN_ACTION')
        action, polarity = [s.lower() for s in m.groups()]
        target_ekb = content.gets('target')
        if target_ekb is None or target_ekb == '':
            return self.make_failure('MISSING_TARGET')
        agent = self._get_agent(target_ekb)
        logger.debug('Found agent (target): %s.' % agent.name)
        residue = content.gets('residue')
        position = content.gets('position')
        related_results = [
            s for s in self.signor_afs
            if self._matching(s, agent, residue, position, action, polarity)
            ]
        if not len(related_results):
            return self.make_failure(
                'MISSING_MECHANISM',
                "Could not find statement matching phosphorylation activating "
                "%s, %s, %s, %s." % (agent.name, residue, position, 'phosphorylation')
                )
        else:
            self._add_provenance(related_results)
            msg = KQMLPerformative('SUCCESS')
            msg.set('is-activating', 'TRUE')
            return msg

    def tell(self, content):
        """Send a tell message."""
        msg = KQMLPerformative('tell')
        msg.set('content', content)
        return self.send(msg)

    def _add_provenance(self, related_results):
        """Creates the content for an add-provenance tell message.

        The message is used to provide evidence supporting the conclusion.
        """
        url_base = 'https://www.ncbi.nlm.nih.gov/pubmed/?term'
        stmt_evidence_fmt = ('Found at pmid <a href={url}={pmid} '
                             'target="_blank">{pmid}</a>:\n{evidence}\n\n')
        content_fmt = "<text>Supporting evidence:\n%s</text><hr>"
        content = KQMLList('add-provenance')
        evidence_list = [stmt.evidence[0] for stmt in related_results]
        pmid_text_dict = {
            pmid: [e.text for e in evidence_list if e.pmid == pmid]
            for pmid in set([e.pmid for e in evidence_list])
            }
        content.sets(
            'html',
            content_fmt % '\n'.join([
                stmt_evidence_fmt.format(
                    url=url_base,
                    pmid=pmid,
                    evidence='\n'.join(['<i>%s</i>' % e for e in elist])
                    )
                for pmid, elist in pmid_text_dict.items()
                ])
            )
        return self.tell(content)

    @staticmethod
    def _get_agent(agent_ekb):
        tp = TripsProcessor(agent_ekb)
        terms = tp.tree.findall('TERM')
        term_id = terms[0].attrib['id']
        agent = tp._get_agent_by_id(term_id, None)
        return agent

    def _matching(self, stmt, agent, residue, position, action, polarity):
        if stmt.agent.name != agent.name:
            return False
        if stmt.is_active is not (polarity == 'activating'):
            return False
        matching_residues = any([
            m.residue == residue
            and m.position == position
            and m.mod_type == action
            for m in stmt.agent.mods])
        return matching_residues


if __name__ == "__main__":
    MSA_Module(argv=sys.argv[1:])
