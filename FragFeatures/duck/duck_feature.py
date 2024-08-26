"""
Process the duck output for a feature of a compound.
"""





class DUckFeature():
    """
    Process the FragFeature metadata for a compound's feature.
    """
    def __init__(self, feature_dir):
        self.feature_dir = feature_dir

    def get_feature(self):
        """
        Return the feature of the compound.
        """
        pass

    def get_feature_dict(self):
        """
        Return the feature of the compound as a dictionary.
        """
        pass

    def get_feature_json(self):
        """
        Return the feature of the compound as a json file.
        """
        pass


# Test
if __name__ == '__main__':
    duck_feature = DUckFeature('/Users/nfo24278/Documents/dphil/diamond/DuCK/structures/CHIKV_Mac_simulations/Experiment/cx0270a/A_ARG_144_NH1')
    # duck_feature.get_feature()
    # duck_feature.get_feature_dict()
    # duck_feature.get_feature_json()