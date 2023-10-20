from gammapy.datasets import Datasets
from gammapy.modeling.models import Models


def save_datasets_with_models(datasets, path_data, path_models):
    datasets.write(path_data, overwrite=True)
    datasets.models.write(path_models, overwrite=True)


def load_datasets_with_models(path_data, path_model):
    data = Datasets.read(path_data)
    data.models = Models.read(path_model)
    return data
