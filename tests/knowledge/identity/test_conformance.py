"""Run the exact reusable contract against the backend selected by fixtures."""
from registry_contract import RegistryBehavior, AliasBehavior, ConcurrencyBehavior


class TestRegistryConformance(RegistryBehavior, AliasBehavior, ConcurrencyBehavior):
    pass
