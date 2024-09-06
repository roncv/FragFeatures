import pytest
# from unittest.mock import MagicMock
from pandas import DataFrame
from FragFeatures.fragalysis.target_parser import TargetParser  # Adjust the import path as necessary
import logging

@pytest.fixture
def mock_csv(mocker):
    # Mock pd.read_csv to return a DataFrame
    data = {
        "Code": ["comp1", "comp2", "comp3"],
        "Property1": [100, 200, 300],
        "Property2": ["A", "B", "C"]
    }
    df = DataFrame(data)
    mocker.patch('pandas.read_csv', return_value=df)
    return df

@pytest.fixture
def target_parser(mocker, mock_csv):
    # Mock logging to avoid side effects on the file system or console
    mocker.patch.object(logging.Logger, 'info')
    mocker.patch.object(logging.Logger, 'debug')

    # Instantiate TargetParser with mocked directory and verbosity turned off
    return TargetParser("/fake/directory", verbose=False)

def test_initialization(target_parser):
    # Test correct initialization
    assert target_parser.target_dir == "/fake/directory"
    assert target_parser.metadata is not None

def test_get_all_compounds(target_parser):
    # Test get_all_compounds method
    expected_compounds = ["comp1", "comp2", "comp3"]
    assert target_parser.get_all_compounds() == expected_compounds

def test_metadata_columns(target_parser, mock_csv):
    # Test if metadata.columns contains the expected columns
    expected_columns = list(mock_csv.columns)
    assert list(target_parser.metadata.columns) == expected_columns