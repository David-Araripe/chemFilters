import unittest
from pathlib import Path

from chemFilters.filters import RDKitFilters


class TestRDKitFilters(unittest.TestCase):
    def setUp(self) -> None:
        self.testroot = Path(__file__).parent
        self.filterFunc = RDKitFilters()
        self.availFilters = self.filterFunc.availableFilters
        self.testSmiles = (
            (self.testroot / "resources/testSmiles.smi").read_text().splitlines()
        )

    def test_loadting_filters(self):
        for _filter in self.availFilters:
            try:
                newFilter = RDKitFilters(_filter)
            except ValueError:
                self.fail(f"Filter type {_filter} not available.")
            self.assertIsInstance(newFilter, RDKitFilters)

    def tearDown(self) -> None:
        return super().tearDown()
