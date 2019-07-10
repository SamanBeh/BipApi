from joblib import load
import json


model = load("SVMModel.joblib")


jmodel = json.dumps(str(model))


with open("SVMModel.json", 'w') as data_file:
    data_file.write(jmodel)
