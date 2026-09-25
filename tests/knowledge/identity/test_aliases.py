"""Exact lookup keys and attribution remain scientific input data."""
from simforge.knowledge.identity import ExternalIdentifierKey
from registry_contract import command, alias, eid


def test_lookup_nfc_but_no_provider_specific_normalization(registry):
    supplied = alias(accession=' e\u0301Ab.1 ')
    registry.write(command(aliases=(supplied,)))
    key = ExternalIdentifierKey('uniprot', ' éAb.1 ')
    assert registry.lookup_external(key)[0].alias == supplied
    for accession in ('éAb.1', ' éab.1 ', ' éAb ', ' éAb.2 '):
        assert registry.lookup_external(ExternalIdentifierKey('uniprot', accession)) == ()
    assert registry.lookup_external(ExternalIdentifierKey('other', key.accession)) == ()


def test_source_supersession_does_not_redirect_or_merge(registry):
    registry.write(command(aliases=(alias(lifecycle='DEPRECATED', release='R2'),)))
    registry.write(command(n=2, op=2, aliases=(alias(2, release='R2'),)))
    matches = registry.lookup_external(ExternalIdentifierKey('uniprot', 'P00533'))
    assert [m.entity.sf_id for m in matches] == [eid(), eid(2)]
    assert [m.alias.lifecycle for m in matches] == ['DEPRECATED', 'ACTIVE']
    assert registry.get(eid()).entity.lifecycle == 'ACTIVE'
